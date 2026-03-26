#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../system/vector.h"

/* Tree walk timing instrumentation */
static double _treewalk_time_total = 0;
static long _treewalk_call_count = 0;
#ifdef _OPENMP
#pragma omp threadprivate(_treewalk_time_total, _treewalk_call_count)
#endif
void treewalk_timing_reset(void) { _treewalk_time_total = 0; _treewalk_call_count = 0; }
void treewalk_timing_report(void) {
    double total = _treewalk_time_total;
    long count = _treewalk_call_count;
    double total_all = 0; long count_all = 0;
    MPI_Reduce(&total, &total_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&count, &count_all, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0) {
        printf("TREEWALK TIMING: %.3f sec avg/rank, %ld calls avg/rank\n",
               total_all / NTask, count_all / NTask);
    }
}

/*!
 * This file contains routines for neighbour finding.  We use the gravity-tree and a range-searching technique to find neighbours.
 */
/*!
 * This file was originally part of the GADGET3 code developed by Volker Springel. The code has been heavily modified
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO (adding/consolidating some of the search routines as needed for different fluids).
 * the modules now are more modular and primarily run on a generic structure which was built entirely for GIZMO, as opposed to the
 * GADGET3-style neighbor finding, to allow for greater flexibility, more stable memory use, and efficient multi-threading.
 * More importantly general routines for neighbor-searching and walking trees have been completely rewritten into a new
 * set of modular files and structured portions of routines designed to be called by a high-level evaluation routine, so that all code can use
 * consistent definitions and these only need to be modified in a single place.
 */


/*! This function constructs the neighbour tree. To this end, we actually need to construct the gravitational tree, 
 *  because we use it now for the neighbour search.
 */
void ngb_treebuild(void)
{
    if(ThisTask == 0) {printf("Begin Ngb-tree construction.\n");}
    CPU_Step[CPU_MISC] += measure_time();
    force_treebuild(NumPart, NULL);
    /* Note: NgbCache (from tree_optimization) is NOT needed with SoA layout.
       Individual fields like P.Type[p], P.Mass[p], P.Pos[p] are already in
       contiguous arrays — no cache pollution from loading unrelated fields. */
    CPU_Step[CPU_TREEBUILD] += measure_time();
    if(ThisTask == 0) {printf("Ngb-Tree contruction finished \n");}
}


/*! This routine finds all neighbours `j' that can interact with the particle `i' in the communication buffer.
 *  Note that an interaction can take place if: \f$ r_{ij} < h_i \f$  OR if  \f$ r_{ij} < h_j \f$.
 *
 *  In the range-search this is taken into account, i.e. it is guaranteed that all particles are found that fulfill this condition, 
 *  including the (more difficult) second part of it. For this purpose, each node knows the maximum h occuring among the particles it represents.
 */
int ngb_treefind_pairs_threads(MyDouble searchcenter[3], MyFloat rkern, int target, int *startnode,
                               int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
#include "../system/ngb_codeblock_before_condition.h"
    if(P.Type[p] > 0) continue; // skip particles with non-gas types
    if(P.Mass[p] <= 0) continue; // skip zero-mass particles
#define NGB_ONLY_OPEN_NODES_CONTAINING_GAS // only want gas
#define SEARCHBOTHWAYS 1 // need neighbors that can -mutually- see one another, not just single-directional searching here
#include "../system/ngb_codeblock_after_condition_threaded.h"
#undef SEARCHBOTHWAYS // must be undefined after code block inserted, or compiler will crash
#undef NGB_ONLY_OPEN_NODES_CONTAINING_GAS
}


/*! This function returns neighbours with distance <= rkern and returns them in Ngblist. Actually, particles in a box of half side length rkern are
 *  returned, i.e. the reduction to a sphere still needs to be done in the calling routine.
 */
int ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat rkern, int target, int *startnode,
				  int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
#include "../system/ngb_codeblock_before_condition.h"
    if(P.Type[p] > 0) continue; // skip particles with non-gas types
    if(P.Mass[p] <= 0) continue; // skip zero-mass particles
#define NGB_ONLY_OPEN_NODES_CONTAINING_GAS // only want gas
#define SEARCHBOTHWAYS 0 // only need neighbors inside of search radius, not particles 'looking at' primary
#include "../system/ngb_codeblock_after_condition_threaded.h"
#undef SEARCHBOTHWAYS
#undef NGB_ONLY_OPEN_NODES_CONTAINING_GAS
}


/*! Thread-local pre-found neighbor list for batched evaluation.
 *  When set, ngb_treefind_optimized returns this list instead of doing a tree walk. */
static int *_prefound_ngblist = NULL;
static int  _prefound_numngb = -1;
#ifdef _OPENMP
#pragma omp threadprivate(_prefound_ngblist, _prefound_numngb)
#endif

void ngb_set_prefound_list(int *list, int count) { _prefound_ngblist = list; _prefound_numngb = count; }
void ngb_clear_prefound_list(void) { _prefound_ngblist = NULL; _prefound_numngb = -1; }


/*! Optimized single-query neighbor search with prefound-list fast path.
 *  If a pre-found neighbor list has been set, returns it immediately without tree walk.
 *  Otherwise does a full tree walk reading P[] directly (not NgbCache). */
int ngb_treefind_optimized(MyDouble searchcenter[3], MyFloat rkern, int target, int *startnode,
                           int mode, int *exportflag, int *exportnodecount, int *exportindex,
                           int *ngblist, int search_both_ways)
{
    /* fast path: return pre-found neighbor list from batched tree walk */
    if(_prefound_numngb >= 0) {
        int n = _prefound_numngb;
        if(_prefound_ngblist != ngblist) {for(int i = 0; i < n; i++) ngblist[i] = _prefound_ngblist[i];}
        _prefound_numngb = -1;
        *startnode = -1;
        return n;
    }
    /* fall through to regular tree walk */
    if(search_both_ways) {
        return ngb_treefind_pairs_threads(searchcenter, rkern, target, startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
    } else {
        return ngb_treefind_variable_threads(searchcenter, rkern, target, startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
    }
}


/* this is the same as above, but the simpler un-threaded version, useful for historical reasons and because some sub-routines use 
    this model without threading because there is no real performance gain. note the slightly different construction of the subroutine below.
    TARGET_BITMASK should be set as a bitmask, i.e. SUM(2^n), where n are all the particle types desired for neighbor finding,
    so e.g. if you want particle types 0 and 4, set TARGET_BITMASK = 17 = 1 + 16 = 2^0 + 2^4
*/
int ngb_treefind_variable_targeted(MyDouble searchcenter[3], MyFloat rkern, int target, int *startnode, int mode, int *nexport, int *nsend_local, int TARGET_BITMASK)
{
    long nexport_save = *nexport; /* this line must be here in the un-threaded versions */
#include "../system/ngb_codeblock_before_condition.h" // call the same variable/initialization block
    if(!((1 << P.Type[p]) & (TARGET_BITMASK))) continue; // skip anything not of the desired type
    if(P.Mass[p] <= 0) continue; // skip zero-mass particles
#define SEARCHBOTHWAYS 0 // only need neighbors inside of search radius, not particles 'looking at' primary
#include "../system/ngb_codeblock_after_condition_unthreaded.h" // call the main loop block as above, but this time the -unthreaded- version
#undef SEARCHBOTHWAYS
}
/* identical to above but includes 'both ways' search for interacting neighbors */
int ngb_treefind_pairs_targeted(MyDouble searchcenter[3], MyFloat rkern, int target, int *startnode, int mode, int *nexport, int *nsend_local, int TARGET_BITMASK)
{
    long nexport_save = *nexport; /* this line must be here in the un-threaded versions */
#include "../system/ngb_codeblock_before_condition.h" // call the same variable/initialization block
    if(!((1 << P.Type[p]) & (TARGET_BITMASK))) continue; // skip anything not of the desired type
    if(P.Mass[p] <= 0) continue; // skip zero-mass particles
#define SEARCHBOTHWAYS 1 // only need neighbors inside of search radius, not particles 'looking at' primary
#include "../system/ngb_codeblock_after_condition_unthreaded.h" // call the main loop block as above, but this time the -unthreaded- version
#undef SEARCHBOTHWAYS
}


/*  slightly modified version of treefind that searches for one or more types of particles: 
        TARGET_BITMASK should be set as a bitmask, i.e. SUM(2^n), where n are all the particle types desired for neighbor finding,
        so e.g. if you want particle types 0 and 4, set TARGET_BITMASK = 17 = 1 + 16 = 2^0 + 2^4
 */
int ngb_treefind_variable_threads_targeted(MyDouble searchcenter[3], MyFloat rkern, int target, int *startnode,
                                           int mode, int *exportflag, int *exportnodecount, int *exportindex,
                                           int *ngblist, int TARGET_BITMASK)
{
#include "../system/ngb_codeblock_before_condition.h"
    if(!((1 << P.Type[p]) & (TARGET_BITMASK))) continue; // skip anything not of the desired type
    if(P.Mass[p] <= 0) continue; // skip zero-mass particles
#define SEARCHBOTHWAYS 0 // only need neighbors inside of search radius, not particles 'looking at' primary
#include "../system/ngb_codeblock_after_condition_threaded.h"
#undef SEARCHBOTHWAYS
}
/* identical to above but includes 'both ways' search for interacting neighbors */
int ngb_treefind_pairs_threads_targeted(MyDouble searchcenter[3], MyFloat rkern, int target, int *startnode,
                                           int mode, int *exportflag, int *exportnodecount, int *exportindex,
                                           int *ngblist, int TARGET_BITMASK)
{
#include "../system/ngb_codeblock_before_condition.h"
    if(!((1 << P.Type[p]) & (TARGET_BITMASK))) continue; // skip anything not of the desired type
    if(P.Mass[p] <= 0) continue; // skip zero-mass particles
#define SEARCHBOTHWAYS 1 // also want particles 'looking at' primary
#include "../system/ngb_codeblock_after_condition_threaded.h"
#undef SEARCHBOTHWAYS
}






/*! Batched neighbor search: walks the tree ONCE for up to NGB_BATCH_SIZE nearby queries simultaneously.
 *  For each candidate particle, checks all B queries — this inner loop auto-vectorizes.
 *  Supports MPI export via per-query export tracking. Returns 0 on success, -1 on buffer overflow.
 */
int ngb_treefind_pairs_threads_batched(
    int batch_size,
    MyDouble searchcenters[][3],
    MyFloat rkerns[],
    int targets[],
    int *exportflag,
    int *exportnodecount,
    int *exportindex,
    int *ngblists[],
    int numngb_out[],
    int search_both_ways)
{
    MyDouble xtmp; xtmp=0; /* needed by NGB_PERIODIC_BOX_LONG macros */
    int maxPart = All.MaxPart;
    int maxNodes = MaxNodes;
    integertime ti_Current = All.Ti_Current;
    int b;

    int numngb[NGB_BATCH_SIZE];
    for(b = 0; b < batch_size; b++) {numngb[b] = 0;}

    int no = maxPart; /* root node */

    while(no >= 0)
    {
        if(no < maxPart) /* single particle */
        {
            int p = no;
            no = Nextnode[no];

            if(P.Type[p] > 0) continue;
            if(P.Mass[p] <= 0) continue;

            if(__builtin_expect(P.Ti_current[p] != ti_Current, 0))
            {
#ifdef _OPENMP
#pragma omp critical(_partdriftngb_)
#endif
                {drift_particle(p, ti_Current);}
            }

            MyDouble cp0 = P.Pos[p][0], cp1 = P.Pos[p][1], cp2 = P.Pos[p][2];
            MyFloat cp_hr = P.KernelRadius[p];

            if(search_both_ways) {
                for(b = 0; b < batch_size; b++)
                {
                    MyFloat dist_b = (cp_hr > rkerns[b]) ? cp_hr : rkerns[b];
                    MyDouble ddx = NGB_PERIODIC_BOX_LONG_X(cp0 - searchcenters[b][0], cp1 - searchcenters[b][1], cp2 - searchcenters[b][2], -1);
                    if(ddx > dist_b) continue;
                    MyDouble ddy = NGB_PERIODIC_BOX_LONG_Y(cp0 - searchcenters[b][0], cp1 - searchcenters[b][1], cp2 - searchcenters[b][2], -1);
                    if(ddy > dist_b) continue;
                    MyDouble ddz = NGB_PERIODIC_BOX_LONG_Z(cp0 - searchcenters[b][0], cp1 - searchcenters[b][1], cp2 - searchcenters[b][2], -1);
                    if(ddz > dist_b) continue;
                    if(ddx*ddx + ddy*ddy + ddz*ddz > dist_b*dist_b) continue;
                    ngblists[b][numngb[b]++] = p;
                }
            } else {
                for(b = 0; b < batch_size; b++) {
                    MyDouble ddx = NGB_PERIODIC_BOX_LONG_X(cp0 - searchcenters[b][0], cp1 - searchcenters[b][1], cp2 - searchcenters[b][2], -1);
                    if(ddx > rkerns[b]) continue;
                    MyDouble ddy = NGB_PERIODIC_BOX_LONG_Y(cp0 - searchcenters[b][0], cp1 - searchcenters[b][1], cp2 - searchcenters[b][2], -1);
                    if(ddy > rkerns[b]) continue;
                    MyDouble ddz = NGB_PERIODIC_BOX_LONG_Z(cp0 - searchcenters[b][0], cp1 - searchcenters[b][1], cp2 - searchcenters[b][2], -1);
                    if(ddz > rkerns[b]) continue;
                    if(ddx*ddx + ddy*ddy + ddz*ddz > rkerns[b]*rkerns[b]) continue;
                    ngblists[b][numngb[b]++] = p;
                }
            }
        }
        else if(no >= maxPart + maxNodes) /* pseudo particle */
        {
            int task = DomainTask[no - (maxPart + maxNodes)];
            for(b = 0; b < batch_size; b++)
            {
                int tgt = targets[b];
                if(tgt < 0) continue;
                if(exportflag[task] != tgt)
                {
                    exportflag[task] = tgt;
                    exportnodecount[task] = NODELISTLENGTH;
                }
                if(exportnodecount[task] == NODELISTLENGTH)
                {
                    int exitFlag = 0, nexp;
#ifdef _OPENMP
#pragma omp critical(_nexportngb_)
#endif
                    {
                        if(Nexport >= All.BunchSize) {BufferFullFlag = 1; exitFlag = 1;}
                        else {nexp = Nexport; Nexport++;}
                    }
                    if(exitFlag) {for(int bb=0;bb<batch_size;bb++) numngb_out[bb]=numngb[bb]; return -1;}
                    exportnodecount[task] = 0;
                    exportindex[task] = nexp;
                    DataIndexTable[nexp].Task = task;
                    DataIndexTable[nexp].Index = tgt;
                    DataIndexTable[nexp].IndexGet = nexp;
                }
                DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] = DomainNodeIndex[no - (maxPart + maxNodes)];
                if(exportnodecount[task] < NODELISTLENGTH)
                    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
            }
            no = Nextnode[no - maxNodes];
            continue;
        }
        else /* internal node */
        {
            struct NODE *current = &Nodes[no];

            if(__builtin_expect(current->Ti_current != ti_Current, 0))
            {
#ifdef _OPENMP
#pragma omp critical(_nodedriftngb_)
#endif
                {force_drift_node(no, ti_Current);}
            }

            if(__builtin_expect(current->N_part <= 1, 0))
            {
                if(current->u.d.mass) {no = current->u.d.nextnode; continue;}
            }

            double node_hmax = Extnodes[no].hmax;
            if(All.Time > All.TimeBegin) {if(node_hmax <= 0) {no = current->u.d.sibling; continue;}}

            MyFloat nc0 = current->center[0], nc1 = current->center[1], nc2 = current->center[2];
            MyFloat half_len = 0.5 * current->len;
            MyFloat sphere_extra = CUBE_EDGEFACTOR_1 * current->len;

            int open_node = 0;
            for(b = 0; b < batch_size; b++)
            {
                MyFloat dist_b = (search_both_ways ? ((node_hmax > rkerns[b]) ? node_hmax : rkerns[b]) : rkerns[b]) + half_len;
                MyDouble ddx = NGB_PERIODIC_BOX_LONG_X(nc0 - searchcenters[b][0], nc1 - searchcenters[b][1], nc2 - searchcenters[b][2], -1);
                if(ddx > dist_b) continue;
                MyDouble ddy = NGB_PERIODIC_BOX_LONG_Y(nc0 - searchcenters[b][0], nc1 - searchcenters[b][1], nc2 - searchcenters[b][2], -1);
                if(ddy > dist_b) continue;
                MyDouble ddz = NGB_PERIODIC_BOX_LONG_Z(nc0 - searchcenters[b][0], nc1 - searchcenters[b][1], nc2 - searchcenters[b][2], -1);
                if(ddz > dist_b) continue;
                MyFloat dist_sph = dist_b + sphere_extra;
                if(ddx*ddx + ddy*ddy + ddz*ddz > dist_sph*dist_sph) continue;
                open_node = 1;
                break;
            }

            if(open_node) {no = current->u.d.nextnode;} else {no = current->u.d.sibling;}
            continue;
        }
    }

    for(b = 0; b < batch_size; b++) {numngb_out[b] = numngb[b];}
    return 0;
}


/*
    custom code for FOF finder -- needs to be able to deal with complications like pure node-linkages and hard-codes a
    local requirement, so we can't use our simple routines above. this is a customized version of the "ngb_treefind_variable" routine above.
    as a result, updates to the core neighbor search routine will not alter this subroutine
 */
int ngb_treefind_fof_primary(MyDouble searchcenter[3], MyFloat rkern, int target, int *startnode, int mode, int *nexport, int *nsend_local, int MyFOF_PRIMARY_LINK_TYPES)
{
    int numngb, no, p, task, nexport_save;
    struct NODE *current;
    MyDouble dx, dy, dz, dist, r2;
    // cache some global vars locally for improved compiler alias analysis
    int maxPart = All.MaxPart;
    int maxNodes = MaxNodes;
    long bunchSize = All.BunchSize;
    MyDouble xtmp; xtmp=0;

    nexport_save = *nexport;
    
    numngb = 0;
    no = *startnode;
    
    while(no >= 0)
    {
        if(no < maxPart)		/* single particle */
        {
            p = no;
            no = Nextnode[no];
            
            if(!((1 << P.Type[p]) & (MyFOF_PRIMARY_LINK_TYPES)))
                continue;
            
            if(mode == 0)
                continue;
            
            dist = rkern;
            dx = NGB_PERIODIC_BOX_LONG_X(P.Pos[p][0] - searchcenter[0], P.Pos[p][1] - searchcenter[1], P.Pos[p][2] - searchcenter[2],-1);
            if(dx > dist) continue;
            dy = NGB_PERIODIC_BOX_LONG_Y(P.Pos[p][0] - searchcenter[0], P.Pos[p][1] - searchcenter[1], P.Pos[p][2] - searchcenter[2],-1);
            if(dy > dist) continue;
            dz = NGB_PERIODIC_BOX_LONG_Z(P.Pos[p][0] - searchcenter[0], P.Pos[p][1] - searchcenter[1], P.Pos[p][2] - searchcenter[2],-1);
            if(dz > dist) continue;
            if(dx * dx + dy * dy + dz * dz > dist * dist) continue;

            Ngblist[numngb++] = p;
        }
        else
        {
            if(no >= maxPart + maxNodes)	/* pseudo particle */
            {
                if(mode == 1)
                    endrun(123125);
                
                if(mode == 0)
                {
                    if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
                    {
                        Exportflag[task] = target;
                        Exportnodecount[task] = NODELISTLENGTH;
                    }
                    
                    if(Exportnodecount[task] == NODELISTLENGTH)
                    {
                        if(*nexport >= bunchSize)
                        {
                            *nexport = nexport_save;
                            if(nexport_save == 0) {endrun(13005);}	/* in this case, the buffer is too small to process even a single particle */
                            for(task = 0; task < NTask; task++) {nsend_local[task] = 0;}
                            for(no = 0; no < nexport_save; no++) {nsend_local[DataIndexTable[no].Task]++;}
                            return -1; /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */
                        }
                        Exportnodecount[task] = 0;
                        Exportindex[task] = *nexport;
                        DataIndexTable[*nexport].Task = task;
                        DataIndexTable[*nexport].Index = target;
                        DataIndexTable[*nexport].IndexGet = *nexport;
                        *nexport = *nexport + 1;
                        nsend_local[task]++;
                    }
                    
                    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
                    DomainNodeIndex[no - (maxPart + maxNodes)];
                    
                    if(Exportnodecount[task] < NODELISTLENGTH)
                        DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                }
                
                if(mode == -1)
                {
                    *nexport = 1;
                }
                
                no = Nextnode[no - maxNodes];
                continue;
                
            }
            
            current = &Nodes[no];
            
            if(mode == 1)
            {
                if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                {
                    *startnode = -1;
                    return numngb;
                }
            }
            
            if(mode == 0)
            {
                if(!(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL)))	/* we have a node with only local particles, can skip branch */
                {
                    no = current->u.d.sibling;
                    continue;
                }
            }
            
            no = current->u.d.sibling;	/* in case the node can be discarded */
            
            dist = rkern + 0.5 * current->len;
            dx = NGB_PERIODIC_BOX_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
            if(dx > dist) continue;
            dy = NGB_PERIODIC_BOX_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
            if(dy > dist) continue;
            dz = NGB_PERIODIC_BOX_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
            if(dz > dist) continue;
            /* now test against the minimal sphere enclosing everything */
            dist += CUBE_EDGEFACTOR_1 * current->len;
            if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist) continue;
            
            if((current->u.d.bitflags & ((1 << BITFLAG_TOPLEVEL) + (1 << BITFLAG_DEPENDS_ON_LOCAL_ELEMENT))) == 0)	/* only use fully local nodes */
            {
                /* test whether the node is contained within the sphere */
                dist = rkern - CUBE_EDGEFACTOR_2 * current->len;
                if(dist > 0)
                    if(r2 < dist * dist)
                    {
                        if(current->u.d.bitflags & (1 << BITFLAG_INSIDE_LINKINGLENGTH))	/* already flagged */
                        {
                            /* sufficient to return only one particle inside this cell */
                            
                            p = current->u.d.nextnode;
                            while(p >= 0)
                            {
                                if(p < maxPart)
                                {
                                    if(((1 << P.Type[p]) & (MyFOF_PRIMARY_LINK_TYPES)))
                                    {
                                        dx = NGB_PERIODIC_BOX_LONG_X(P.Pos[p][0] - searchcenter[0], P.Pos[p][1] - searchcenter[1], P.Pos[p][2] - searchcenter[2],-1);
                                        dy = NGB_PERIODIC_BOX_LONG_Y(P.Pos[p][0] - searchcenter[0], P.Pos[p][1] - searchcenter[1], P.Pos[p][2] - searchcenter[2],-1);
                                        dz = NGB_PERIODIC_BOX_LONG_Z(P.Pos[p][0] - searchcenter[0], P.Pos[p][1] - searchcenter[1], P.Pos[p][2] - searchcenter[2],-1);
                                        if(dx * dx + dy * dy + dz * dz > rkern * rkern) break;

                                        Ngblist[numngb++] = p;
                                        break;
                                    }
                                    p = Nextnode[p];
                                }
                                else if(p >= maxPart + maxNodes)
                                    p = Nextnode[p - maxNodes];
                                else
                                    p = Nodes[p].u.d.nextnode;
                            }
                            continue;
                        }
                        else
                        {
                            /* flag it now */
                            current->u.d.bitflags |= (1 << BITFLAG_INSIDE_LINKINGLENGTH);
                        }
                    }
            }
            
            no = current->u.d.nextnode;	/* ok, we need to open the node */
        }
    }
    
    *startnode = -1;
    return numngb;
}





