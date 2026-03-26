#!/usr/bin/env python3
"""
Convert GIZMO from AoS (Array of Structures) to SoA (Structure of Arrays) layout.

This script:
1. Parses particle_data.h and cell_data.h to extract field definitions
2. Generates SoA struct headers (particle_soa.h, cell_soa.h)
3. Transforms P[i].Field -> P.Field[i] and CellP[i].Field -> CellP.Field[i] across all source files
4. Generates allocation and helper functions
"""

import re
import os
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent

def parse_struct_fields(filepath):
    """Parse a struct definition file and extract field declarations with their preprocessor context."""
    with open(filepath) as f:
        lines = f.readlines()

    fields = []  # list of (preprocessor_context, type, name, array_suffix, comment)
    ifdef_stack = []  # stack of #ifdef/#if conditions
    in_nested_struct = None  # name of nested struct if inside one
    nested_fields = []

    for line in lines:
        stripped = line.strip()

        # Track preprocessor conditionals
        if stripped.startswith('#if'):
            cond = stripped
            ifdef_stack.append(cond)
            if in_nested_struct:
                nested_fields.append(('IFDEF', cond, None, None, None))
            else:
                fields.append(('IFDEF', cond, None, None, None))
            continue
        if stripped.startswith('#else'):
            if in_nested_struct:
                nested_fields.append(('ELSE', None, None, None, None))
            else:
                fields.append(('ELSE', None, None, None, None))
            continue
        if stripped.startswith('#endif'):
            if ifdef_stack:
                ifdef_stack.pop()
            if in_nested_struct:
                nested_fields.append(('ENDIF', None, None, None, None))
            else:
                fields.append(('ENDIF', None, None, None, None))
            continue
        if stripped.startswith('#define'):
            if in_nested_struct:
                nested_fields.append(('DEFINE', stripped, None, None, None))
            else:
                fields.append(('DEFINE', stripped, None, None, None))
            continue

        # Detect nested struct start
        if re.match(r'^\s*struct\s*$', stripped) or re.match(r'^\s*struct\s*\{', stripped):
            in_nested_struct = '__pending__'
            nested_fields = []
            continue
        if stripped == '{' and in_nested_struct == '__pending__':
            continue

        # Detect nested struct end with name
        if in_nested_struct and stripped.startswith('}'):
            m = re.match(r'\}\s*(\w+)\s*;', stripped)
            if m:
                struct_name = m.group(1)
                fields.append(('NESTED_STRUCT', struct_name, nested_fields, None, None))
                in_nested_struct = None
                nested_fields = []
                continue

        # Detect union start
        if re.match(r'^\s*union\s*$', stripped) or re.match(r'^\s*union\s*\{', stripped):
            in_nested_struct = '__union_pending__'
            nested_fields = []
            continue
        if stripped == '{' and in_nested_struct == '__union_pending__':
            continue
        if in_nested_struct and in_nested_struct.startswith('__union') and stripped.startswith('}'):
            m = re.match(r'\}\s*(\w+)\s*;', stripped)
            if m:
                struct_name = m.group(1)
                fields.append(('UNION', struct_name, nested_fields, None, None))
                in_nested_struct = None
                nested_fields = []
                continue

        # Parse field declarations
        # Match patterns like: Type Name; or Type Name[N]; or Type Name[N][M];
        # Also handles: Vec3<MyDouble> Name;
        field_match = re.match(
            r'^\s*'
            r'((?:ALIGN\(\d+\)\s+)?'  # optional ALIGN
            r'(?:(?:const|volatile|unsigned|signed|short|long|static|extern)\s+)*'  # qualifiers
            r'(?:\w+(?:<[^>]+>)?(?:\s*\*)*)'  # base type (possibly templated, possibly pointer)
            r'(?:\s+(?:const|volatile|unsigned|signed|short|long)\s*)*'  # more qualifiers
            r')'
            r'\s+'
            r'(\w+)'  # field name
            r'((?:\[[^\]]*\])*)'  # array dimensions
            r'\s*;'  # semicolon
            r'(.*)',  # rest (comment)
            stripped
        )

        if field_match:
            field_type = field_match.group(1).strip()
            field_name = field_match.group(2)
            array_dims = field_match.group(3)
            comment = field_match.group(4).strip()

            entry = ('FIELD', field_type, field_name, array_dims, comment)
            if in_nested_struct:
                nested_fields.append(entry)
            else:
                fields.append(entry)
            continue

        # Handle multi-declaration lines: "Type name1, name2, name3;"
        multi_match = re.match(
            r'^\s*'
            r'((?:(?:const|volatile|unsigned|signed|short|long|static|extern)\s+)*'
            r'(?:\w+(?:<[^>]+>)?(?:\s*\*)*)'
            r'(?:\s+(?:const|volatile|unsigned|signed|short|long)\s*)*'
            r')'
            r'\s+'
            r'(\w+(?:\s*,\s*\w+)+)'  # comma-separated names
            r'\s*;'
            r'(.*)',
            stripped
        )
        if multi_match:
            field_type = multi_match.group(1).strip()
            names_str = multi_match.group(2)
            comment = multi_match.group(3).strip()
            for name in re.split(r'\s*,\s*', names_str):
                entry = ('FIELD', field_type, name.strip(), '', comment)
                comment = ''  # only attach comment to first
                if in_nested_struct:
                    nested_fields.append(entry)
                else:
                    fields.append(entry)
            continue

    return fields


def field_to_soa_pointer(field_type, field_name, array_dims):
    """Convert a struct field declaration to SoA pointer declaration.

    Type Name;          -> Type *Name;
    Type Name[N];       -> Type (*Name)[N];  (pointer to array)
    Type Name[N][M];    -> Type (*Name)[N][M];
    ALIGN(32) Type Name; -> Type *Name;  (drop alignment, will align arrays separately)
    """
    # Strip ALIGN
    clean_type = re.sub(r'ALIGN\(\d+\)\s+', '', field_type)

    if array_dims:
        return f"{clean_type} (*{field_name}){array_dims};"
    else:
        return f"{clean_type} *{field_name};"


def generate_soa_header(fields, struct_name, var_name, buf_name, orig_header):
    """Generate SoA struct header from parsed fields."""
    lines = []
    lines.append(f"/* SoA (Structure of Arrays) version of {orig_header}")
    lines.append(f" * Auto-generated by convert_to_soa.py")
    lines.append(f" * Original struct fields become per-field arrays.")
    lines.append(f" * Access pattern: {var_name}[i].Field -> {var_name}.Field[i]")
    lines.append(f" */")
    lines.append(f"")
    lines.append(f"#pragma once")
    lines.append(f"")

    # Forward declare the helper
    lines.append(f"struct {struct_name};")
    lines.append(f"void allocate_{var_name}({struct_name}& s, int count, const char* label);")
    lines.append(f"void copy_particle_{var_name}({struct_name}& s, int dest, int src);")
    lines.append(f"void move_particles_{var_name}({struct_name}& s, int dest, int src, int count);")
    lines.append(f"")

    lines.append(f"struct {struct_name}")
    lines.append(f"{{")

    def emit_fields(field_list, indent="    "):
        for entry in field_list:
            kind = entry[0]
            if kind == 'IFDEF':
                lines.append(f"{entry[1]}")
            elif kind == 'ELSE':
                lines.append(f"#else")
            elif kind == 'ENDIF':
                lines.append(f"#endif")
            elif kind == 'DEFINE':
                lines.append(f"{entry[1]}")
            elif kind == 'FIELD':
                _, ftype, fname, fdims, fcomment = entry
                decl = field_to_soa_pointer(ftype, fname, fdims)
                if fcomment:
                    lines.append(f"{indent}{decl} {fcomment}")
                else:
                    lines.append(f"{indent}{decl}")
            elif kind == 'NESTED_STRUCT':
                _, sname, nested, _, _ = entry
                lines.append(f"{indent}struct {{")
                emit_fields(nested, indent + "    ")
                lines.append(f"{indent}}} {sname};")
            elif kind == 'UNION':
                _, uname, nested, _, _ = entry
                # For SoA, union fields become union of pointers
                # so they alias to same memory (same pointer)
                lines.append(f"{indent}union {{")
                emit_fields(nested, indent + "    ")
                lines.append(f"{indent}}} {uname};")

    emit_fields(fields)

    lines.append(f"}};")
    lines.append(f"")
    lines.append(f"extern struct {struct_name} {var_name}; /* SoA data on local processor */")
    lines.append(f"extern struct {struct_name} {buf_name}; /* buffer for domain decomposition */")
    lines.append(f"")

    return "\n".join(lines)


def generate_alloc_function(fields, struct_name, var_name):
    """Generate allocation function for SoA struct."""
    lines = []
    lines.append(f"void allocate_{var_name}({struct_name}& s, int count, const char* label)")
    lines.append(f"{{")
    lines.append(f"    if(count <= 0) return;")

    def emit_allocs(field_list, prefix="s."):
        for entry in field_list:
            kind = entry[0]
            if kind == 'IFDEF':
                lines.append(f"{entry[1]}")
            elif kind == 'ELSE':
                lines.append(f"#else")
            elif kind == 'ENDIF':
                lines.append(f"#endif")
            elif kind == 'DEFINE':
                lines.append(f"{entry[1]}")
            elif kind == 'FIELD':
                _, ftype, fname, fdims, _ = entry
                clean_type = re.sub(r'ALIGN\(\d+\)\s+', '', ftype)
                if fdims:
                    # Array field: allocate as array of arrays
                    lines.append(f"    {prefix}{fname} = ({clean_type} (*){fdims}) mymalloc(\"{var_name}.{fname}\", count * sizeof({clean_type}{fdims}));")
                else:
                    lines.append(f"    {prefix}{fname} = ({clean_type} *) mymalloc(\"{var_name}.{fname}\", count * sizeof({clean_type}));")
            elif kind == 'NESTED_STRUCT':
                _, sname, nested, _, _ = entry
                emit_allocs(nested, prefix=f"{prefix}{sname}.")
            elif kind == 'UNION':
                _, uname, nested, _, _ = entry
                # For unions, allocate the first field and alias the rest
                first = True
                first_name = None
                first_type = None
                for sub in nested:
                    if sub[0] == 'FIELD':
                        if first:
                            _, ft, fn, fd, _ = sub
                            ct = re.sub(r'ALIGN\(\d+\)\s+', '', ft)
                            if fd:
                                lines.append(f"    {prefix}{uname}.{fn} = ({ct} (*){fd}) mymalloc(\"{var_name}.{uname}.{fn}\", count * sizeof({ct}{fd}));")
                            else:
                                lines.append(f"    {prefix}{uname}.{fn} = ({ct} *) mymalloc(\"{var_name}.{uname}.{fn}\", count * sizeof({ct}));")
                            first_name = fn
                            first_type = ct
                            first = False
                        else:
                            _, ft, fn, fd, _ = sub
                            ct = re.sub(r'ALIGN\(\d+\)\s+', '', ft)
                            # Alias to first field
                            if fd:
                                lines.append(f"    {prefix}{uname}.{fn} = ({ct} (*){fd}) {prefix}{uname}.{first_name};")
                            else:
                                lines.append(f"    {prefix}{uname}.{fn} = ({ct} *) {prefix}{uname}.{first_name};")
                    elif sub[0] in ('IFDEF', 'ELSE', 'ENDIF', 'DEFINE'):
                        if sub[0] == 'IFDEF':
                            lines.append(f"{sub[1]}")
                        elif sub[0] == 'ELSE':
                            lines.append(f"#else")
                        elif sub[0] == 'ENDIF':
                            lines.append(f"#endif")
                        elif sub[0] == 'DEFINE':
                            lines.append(f"{sub[1]}")

    emit_allocs(fields)
    lines.append(f"}}")
    return "\n".join(lines)


def generate_copy_function(fields, struct_name, var_name):
    """Generate per-field copy function."""
    lines = []
    lines.append(f"void copy_particle_{var_name}({struct_name}& s, int dest, int src)")
    lines.append(f"{{")

    def emit_copies(field_list, prefix="s."):
        for entry in field_list:
            kind = entry[0]
            if kind == 'IFDEF':
                lines.append(f"{entry[1]}")
            elif kind == 'ELSE':
                lines.append(f"#else")
            elif kind == 'ENDIF':
                lines.append(f"#endif")
            elif kind == 'DEFINE':
                pass  # skip defines in copy
            elif kind == 'FIELD':
                _, ftype, fname, fdims, _ = entry
                if fdims:
                    clean_type = re.sub(r'ALIGN\(\d+\)\s+', '', ftype)
                    lines.append(f"    memcpy(&{prefix}{fname}[dest], &{prefix}{fname}[src], sizeof({clean_type}{fdims}));")
                else:
                    lines.append(f"    {prefix}{fname}[dest] = {prefix}{fname}[src];")
            elif kind == 'NESTED_STRUCT':
                _, sname, nested, _, _ = entry
                emit_copies(nested, prefix=f"{prefix}{sname}.")
            elif kind == 'UNION':
                _, uname, nested, _, _ = entry
                # Only copy the first field (they alias)
                first = True
                for sub in nested:
                    if sub[0] == 'FIELD':
                        if first:
                            _, ft, fn, fd, _ = sub
                            if fd:
                                ct = re.sub(r'ALIGN\(\d+\)\s+', '', ft)
                                lines.append(f"    memcpy(&{prefix}{uname}.{fn}[dest], &{prefix}{uname}.{fn}[src], sizeof({ct}{fd}));")
                            else:
                                lines.append(f"    {prefix}{uname}.{fn}[dest] = {prefix}{uname}.{fn}[src];")
                            first = False
                    elif sub[0] in ('IFDEF', 'ELSE', 'ENDIF'):
                        if sub[0] == 'IFDEF': lines.append(f"{sub[1]}")
                        elif sub[0] == 'ELSE': lines.append(f"#else")
                        elif sub[0] == 'ENDIF': lines.append(f"#endif")

    emit_copies(fields)
    lines.append(f"}}")
    return "\n".join(lines)


def generate_move_function(fields, struct_name, var_name):
    """Generate per-field memmove function for moving ranges of particles."""
    lines = []
    lines.append(f"void move_particles_{var_name}({struct_name}& s, int dest, int src, int count)")
    lines.append(f"{{")
    lines.append(f"    if(count <= 0 || dest == src) return;")

    def emit_moves(field_list, prefix="s."):
        for entry in field_list:
            kind = entry[0]
            if kind == 'IFDEF':
                lines.append(f"{entry[1]}")
            elif kind == 'ELSE':
                lines.append(f"#else")
            elif kind == 'ENDIF':
                lines.append(f"#endif")
            elif kind == 'DEFINE':
                pass
            elif kind == 'FIELD':
                _, ftype, fname, fdims, _ = entry
                clean_type = re.sub(r'ALIGN\(\d+\)\s+', '', ftype)
                if fdims:
                    lines.append(f"    memmove(&{prefix}{fname}[dest], &{prefix}{fname}[src], count * sizeof({clean_type}{fdims}));")
                else:
                    lines.append(f"    memmove(&{prefix}{fname}[dest], &{prefix}{fname}[src], count * sizeof({clean_type}));")
            elif kind == 'NESTED_STRUCT':
                _, sname, nested, _, _ = entry
                emit_moves(nested, prefix=f"{prefix}{sname}.")
            elif kind == 'UNION':
                _, uname, nested, _, _ = entry
                first = True
                for sub in nested:
                    if sub[0] == 'FIELD' and first:
                        _, ft, fn, fd, _ = sub
                        ct = re.sub(r'ALIGN\(\d+\)\s+', '', ft)
                        if fd:
                            lines.append(f"    memmove(&{prefix}{uname}.{fn}[dest], &{prefix}{uname}.{fn}[src], count * sizeof({ct}{fd}));")
                        else:
                            lines.append(f"    memmove(&{prefix}{uname}.{fn}[dest], &{prefix}{uname}.{fn}[src], count * sizeof({ct}));")
                        first = False
                    elif sub[0] in ('IFDEF', 'ELSE', 'ENDIF'):
                        if sub[0] == 'IFDEF': lines.append(f"{sub[1]}")
                        elif sub[0] == 'ELSE': lines.append(f"#else")
                        elif sub[0] == 'ENDIF': lines.append(f"#endif")

    emit_moves(fields)
    lines.append(f"}}")
    return "\n".join(lines)


def transform_access_patterns(filepath):
    """Transform P[expr].Field -> P.Field[expr] and CellP[expr].Field -> CellP.Field[expr]
    in a source file. Handles nested brackets in the index expression."""

    with open(filepath) as f:
        content = f.read()

    original = content

    def find_balanced_bracket(text, start):
        """Find the matching ] for a [ at position start. Returns index of ]."""
        depth = 0
        i = start
        while i < len(text):
            if text[i] == '[':
                depth += 1
            elif text[i] == ']':
                depth -= 1
                if depth == 0:
                    return i
            i += 1
        return -1

    def transform_var(content, var_name):
        """Transform var[expr].member_chain to var.member_chain[expr]

        The member chain stops when we encounter an identifier followed by '('
        (i.e., a method call on the element type, like Vec3::norm()).
        This ensures P[i].GravAccel.norm() -> P.GravAccel[i].norm() (correct)
        rather than P.GravAccel.norm[i]() (incorrect).
        """
        result = []
        i = 0
        var_len = len(var_name)

        while i < len(content):
            # Check if we're at a potential match
            if content[i:i+var_len] == var_name and i + var_len < len(content) and content[i+var_len] == '[':
                # Check it's a word boundary before var_name
                if i > 0 and (content[i-1].isalnum() or content[i-1] == '_'):
                    result.append(content[i])
                    i += 1
                    continue

                # Find the matching ]
                bracket_start = i + var_len
                bracket_end = find_balanced_bracket(content, bracket_start)

                if bracket_end == -1:
                    result.append(content[i])
                    i += 1
                    continue

                # Extract the index expression
                idx_expr = content[bracket_start+1:bracket_end]

                # Check for . after the ]
                pos_after = bracket_end + 1
                if pos_after < len(content) and content[pos_after] == '.':
                    # Extract member chain: identifiers separated by dots,
                    # but STOP if the next identifier is followed by '(' (method call)
                    member_parts = []
                    pos = pos_after + 1

                    while pos < len(content):
                        # Extract an identifier
                        id_start = pos
                        while pos < len(content) and (content[pos].isalnum() or content[pos] == '_'):
                            pos += 1
                        if pos == id_start:
                            break  # no identifier found

                        identifier = content[id_start:pos]

                        # Check what follows the identifier
                        if pos < len(content) and content[pos] == '(':
                            # It's a method call - don't include this identifier in the chain
                            # Rewind pos to before this identifier
                            pos = id_start - 1  # back to the '.' before this identifier
                            break

                        member_parts.append(identifier)

                        # Check for more dots
                        if pos < len(content) and content[pos] == '.':
                            next_pos = pos + 1
                            if next_pos < len(content) and (content[next_pos].isalnum() or content[next_pos] == '_'):
                                pos = next_pos
                                continue
                        break

                    member_chain = '.'.join(member_parts)

                    if member_chain:
                        # Transform: var[idx].member -> var.member[idx]
                        result.append(f"{var_name}.{member_chain}[{idx_expr}]")
                        i = pos
                        continue

                # Not followed by . member access - leave as is
                # (This handles cases like pointer arithmetic: P + offset)
                result.append(content[i:bracket_end+1])
                i = bracket_end + 1
                continue

            result.append(content[i])
            i += 1

        return ''.join(result)

    content = transform_var(content, 'P')
    content = transform_var(content, 'CellP')

    if content != original:
        with open(filepath, 'w') as f:
            f.write(content)
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: convert_to_soa.py [parse|transform|all]")
        sys.exit(1)

    action = sys.argv[1]

    if action in ('parse', 'all'):
        print("Parsing particle_data.h...")
        p_fields = parse_struct_fields(ROOT / 'declarations' / 'particle_data.h')
        print(f"  Found {sum(1 for f in p_fields if f[0] == 'FIELD')} direct fields")

        print("Parsing cell_data.h...")
        c_fields = parse_struct_fields(ROOT / 'declarations' / 'cell_data.h')
        print(f"  Found {sum(1 for f in c_fields if f[0] == 'FIELD')} direct fields")

        # Generate SoA headers
        p_header = generate_soa_header(p_fields, 'particle_soa', 'P', 'DomainPartBuf', 'particle_data.h')
        with open(ROOT / 'declarations' / 'particle_soa.h', 'w') as f:
            f.write(p_header + '\n')
        print("  Generated declarations/particle_soa.h")

        c_header = generate_soa_header(c_fields, 'cell_soa', 'CellP', 'DomainGasBuf', 'cell_data.h')
        with open(ROOT / 'declarations' / 'cell_soa.h', 'w') as f:
            f.write(c_header + '\n')
        print("  Generated declarations/cell_soa.h")

        # Generate helper functions
        helpers = []
        helpers.append("/* SoA helper functions - auto-generated by convert_to_soa.py */")
        helpers.append('#include "../declarations/allvars.h"')
        helpers.append('#include "../core/proto.h"')
        helpers.append('#include <string.h>')
        helpers.append('')
        helpers.append(generate_alloc_function(p_fields, 'particle_soa', 'P'))
        helpers.append('')
        helpers.append(generate_alloc_function(c_fields, 'cell_soa', 'CellP'))
        helpers.append('')
        helpers.append(generate_copy_function(p_fields, 'particle_soa', 'P'))
        helpers.append('')
        helpers.append(generate_copy_function(c_fields, 'cell_soa', 'CellP'))
        helpers.append('')
        helpers.append(generate_move_function(p_fields, 'particle_soa', 'P'))
        helpers.append('')
        helpers.append(generate_move_function(c_fields, 'cell_soa', 'CellP'))
        helpers.append('')

        with open(ROOT / 'system' / 'soa_helpers.cc', 'w') as f:
            f.write('\n'.join(helpers))
        print("  Generated system/soa_helpers.cc")

    if action in ('transform', 'all'):
        print("\nTransforming access patterns...")
        # Find all .cc and .h files
        count = 0
        for ext in ('*.cc', '*.h'):
            for filepath in ROOT.rglob(ext):
                # Skip generated files and test files we don't need
                rel = filepath.relative_to(ROOT)
                if str(rel).startswith('.claude') or str(rel).startswith('scripts/convert'):
                    continue
                if transform_access_patterns(filepath):
                    count += 1
                    print(f"  Transformed: {rel}")
        print(f"\nTransformed {count} files total.")


if __name__ == '__main__':
    main()
