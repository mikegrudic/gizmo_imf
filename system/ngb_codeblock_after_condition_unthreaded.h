if(P.Ti_current[p] != ti_Current)
{
    drift_particle(p, ti_Current);
}

#if (SEARCHBOTHWAYS==1)
dist = DMAX(P.KernelRadius[p], rkern);
#if defined(BOX_PERIODIC) && !(defined(BOX_REFLECT_X) || defined(BOX_OUTFLOW_X))
dx = boxHalf_X - fabs(fabs(P.Pos[p][0] - searchcenter[0]) - boxHalf_X);
#else
dx = fabs(P.Pos[p][0] - searchcenter[0]);
#endif
if(dx > dist) continue;
#if defined(BOX_PERIODIC) && !(defined(BOX_REFLECT_Y) || defined(BOX_OUTFLOW_Y))
dy = boxHalf_Y - fabs(fabs(P.Pos[p][1] - searchcenter[1]) - boxHalf_Y);
#else
dy = fabs(P.Pos[p][1] - searchcenter[1]);
#endif
if(dy > dist) continue;
#if defined(BOX_PERIODIC) && !(defined(BOX_REFLECT_Z) || defined(BOX_OUTFLOW_Z))
dz = boxHalf_Z - fabs(fabs(P.Pos[p][2] - searchcenter[2]) - boxHalf_Z);
#else
dz = fabs(P.Pos[p][2] - searchcenter[2]);
#endif
if(dz > dist) continue;
if(dx * dx + dy * dy + dz * dz > dist * dist) continue;
#endif

Ngblist[numngb++] = p;
}
else
{
    if(no >= maxPart + maxNodes)	/* pseudo particle */
    {
        if(mode == 1) {endrun(123129);}
        if(target >= 0)	/* if no target is given, export will not occur */
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
                    if(nexport_save == 0) {endrun(13004);} /* in this case, the buffer is too small to process even a single particle */
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

            DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] = DomainNodeIndex[no - (maxPart + maxNodes)];
            if(Exportnodecount[task] < NODELISTLENGTH) {DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;}
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

    if(current->Ti_current != ti_Current) {force_drift_node(no, ti_Current);}

    if(current->N_part <= 1)
    {
        if(current->u.d.mass)	/* open cell */
        {
            no = current->u.d.nextnode;
            continue;
        }
    }

    double hmax = Extnodes[no].hmax;
#if (SEARCHBOTHWAYS==1)
    dist = DMAX(hmax, rkern) + 0.5 * current->len;
#else
    dist = rkern + 0.5 * current->len;
#endif
    no = current->u.d.sibling;	// in case the node can be discarded //
#include "ngb_codeblock_checknode.h"
    no = current->u.d.nextnode;	// ok, we need to open the node //
}
}

*startnode = -1;
return numngb;
