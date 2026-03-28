/*
 this defines a code-block to be inserted in the neighbor search routines after the conditions for neighbor-validity are applied
 (valid particle types checked)
 */
if(P.Ti_current[p] != ti_Current)
{
#ifdef _OPENMP
#pragma omp critical(_partdriftngb_)
#endif
    {
        drift_particle(p, ti_Current);
    }
}

#if (SEARCHBOTHWAYS==1)
dist = DMAX(P.KernelRadius[p], rkern);
/* branchless periodic wrap where applicable — avoids shared xtmp variable */
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
ngblist[numngb++] = p;
}
else
{
    if(no >= maxPart + maxNodes)	/* pseudo particle */
    {
        if(mode == 1) {endrun(123128);}

        if(target >= 0)
        {
            if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
            {
                exportflag[task] = target;
                exportnodecount[task] = NODELISTLENGTH;
            }

            if(exportnodecount[task] == NODELISTLENGTH)
            {
                int exitFlag = 0, nexp;
#ifdef _OPENMP
#pragma omp critical(_nexportngb_)
#endif
                {
                    if(Nexport >= bunchSize)
                    {
                        BufferFullFlag = 1;
                        exitFlag = 1;
                    }
                    else
                    {
                        nexp = Nexport;
                        Nexport++;
                    }
                }
                if(exitFlag) {return -1;}

                exportnodecount[task] = 0;
                exportindex[task] = nexp;
                DataIndexTable[nexp].Task = task;
                DataIndexTable[nexp].Index = target;
                DataIndexTable[nexp].IndexGet = nexp;
            }
            DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] = DomainNodeIndex[no - (maxPart + maxNodes)];
            if(exportnodecount[task] < NODELISTLENGTH)
                DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
                }

        no = Nextnode[no - maxNodes];
        continue;
    }

    current = &Nodes[no];

    if(mode == 1)
    {
        if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))
        {
            *startnode = -1;
            return numngb;
        }
    }

    if(current->Ti_current != ti_Current)
    {
#ifdef _OPENMP
#pragma omp critical(_nodedriftngb_)
#endif
        {
            force_drift_node(no, ti_Current);
        }
    }

    if(current->N_part <= 1)
    {
        if(current->u.d.mass)
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
    no = current->u.d.sibling;
#include "ngb_codeblock_checknode.h"
    no = current->u.d.nextnode;
}
}

*startnode = -1;
return numngb;
