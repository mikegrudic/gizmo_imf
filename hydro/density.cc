#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"
#include "../mesh/mesh_motion.h"
#include "density_tile.h"

/*! \file density.c
 *  \brief hydro kernel size and neighbor determination, volumetric quantities calculated
 *
 *  This file contains the "first hydro loop", where the gas densities and some
 *  auxiliary quantities are computed.  There is also functionality that corrects the kernel length if needed.
 */
/*!
 * This file was originally part of the GADGET3 code developed by Volker Springel.
 * The code has been modified substantially (condensed, different criteria for kernel lengths, optimizatins,
 * rewritten parallelism, new physics included, new variable/memory conventions added, fundamentally different
 * criteria and conditioning and calcuilations actually being done for the modular hydro solvers)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

struct kernel_density /*! defines a number of useful variables we will use below */
{
  Vec3<double> dp; Vec3<double> dv; double r, wk, dwk, hinv, hinv3, hinv4, mj_wk, mj_dwk_r;
};


/*! routine to determine if a given element is actually going to be active in the density subroutines below */
int density_isactive(int n)
{
    /* first check our 'marker' for particles which have finished iterating to an KernelRadius solution (if they have, dont do them again) */
    if(P.TimeBin[n] < 0) {return 0;}
    if(P.Type[n] == 0) {if(CellP.recent_refinement_flag[n] == 1) return 1;}
    
#if defined(GRAIN_FLUID)
    if((1 << P.Type[n]) & (GRAIN_PTYPES)) {return 1;} /* any of the particle types flagged as a valid grain-type is active here */
#endif

#if defined(SINK_INTERACT_ON_GAS_TIMESTEP)
    if(P.Type[n] == 5){if(!P.do_gas_search_this_timestep[n] && All.Ti_Current > 0) return 0;} /* not enough time has elapsed since the last gas interaction */
#endif
#if defined(RT_SOURCE_INJECTION)
    if((1 << P.Type[n]) & (RT_SOURCES))
    {
#if defined(GALSF)
       if(((P.Type[n] == 4)||((All.ComovingIntegrationOn==0)&&((P.Type[n] == 2)||(P.Type[n]==3))))&&(P.Mass[n]>0))
        {
            double star_age = evaluate_stellar_age_Gyr(n);
            if((star_age < 0.1)&&(star_age > 0)&&(!isnan(star_age))) return 1;
        }
#else
        if(Flag_FullStep) {return 1;} // only do on full timesteps
#endif
    }
#endif

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    if(((P.Type[n] == 4)||((All.ComovingIntegrationOn==0)&&((P.Type[n] == 2)||(P.Type[n]==3))))&&(P.Mass[n]>0))
    {
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
        /* check if there is going to be a SNe this timestep, in which case, we want the density info! */
        if(P.SNe_ThisTimeStep[n]>0) return 1;
#if defined(GALSF_FB_FIRE_STELLAREVOLUTION)
        if(P.MassReturn_ThisTimeStep[n]>0) return 1;
#ifdef GALSF_FB_FIRE_RPROCESS
        if(P.RProcessEvent_ThisTimeStep[n]>0) return 1;
#endif
#if defined(GALSF_FB_FIRE_AGE_TRACERS)
        if(P.AgeDeposition_ThisTimeStep[n]>0) return 1;
#endif
#endif
#endif
        
#if defined(GALSF)
        if(P.DensityAroundParticle[n] <= 0) return 1;
        if(All.ComovingIntegrationOn == 0) // only do stellar age evaluation if we have to //
        {
            double star_age = evaluate_stellar_age_Gyr(n);
            if(star_age < 0.035) return 1;
        }
#endif
#if (defined(GRAIN_FLUID) || defined(RADTRANSFER)) && (!defined(GALSF) && !(defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)))
        return 1;
#endif
    }
#endif

#ifdef SINK_PARTICLES
    if(P.Type[n] == 5) return 1;
#endif

    if(P.Type[n] == 0 && P.Mass[n] > 0) return 1;
    return 0; /* default to 0 if no check passed */
}




//#define NGB_SEARCH_BOTH_WAYS 0 /* opt in to batched neighbor search (SEARCHBOTHWAYS=0 for density) */
#define CORE_FUNCTION_NAME density_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME hydrokerneldensity_particle2in    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME hydrokerneldensity_out2particle  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(density_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P.Mass[i]>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/*! this structure defines the variables that need to be sent -from- the 'searching' element */
static struct INPUT_STRUCT_NAME
{
  Vec3<MyDouble> Pos;
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
  Vec3<MyFloat> Accel;
#endif
  Vec3<MyFloat> Vel;
  MyFloat KernelRadius;
#ifdef GALSF_SUBGRID_WINDS
  MyFloat DelayTime;
#endif
  int NodeList[NODELISTLENGTH];
  int Type;
}
 *DATAIN_NAME, *DATAGET_NAME;

/*! this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
void hydrokerneldensity_particle2in(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k;
    in->Type = P.Type[i];
    in->KernelRadius = P.KernelRadius[i];
    in->Pos = P.Pos[i];
    if(P.Type[i]==0) {in->Vel=CellP.VelPred[i];} else {in->Vel=P.Vel[i];}
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
    if(P.Type[i] == 0) {in->Accel = All.cf_a2inv * P.GravAccel[i] + CellP.HydroAccel[i];} // PHYSICAL units //
#endif
#ifdef GALSF_SUBGRID_WINDS
    if(P.Type[i]==0) {in->DelayTime = CellP.DelayTime[i];} else {in->DelayTime=0;}
#endif
}

/*! this structure defines the variables that need to be sent -back to- the 'searching' element */
static struct OUTPUT_STRUCT_NAME
{
    MyDouble Ngb;
    MyDouble Rho;
    MyDouble DrkernNgb;
    MyDouble Particle_DivVel;
    SymmetricTensor2<MyDouble> NV_T;
    Vec3<MyDouble> NV_T_face_weights; /*!< weighted first moments sum(wk*dp[k]); used for face area estimation */
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
    Vec3<MyDouble> ParticleVel;
#endif
#ifdef HYDRO_SPH
    MyDouble DrkernHydroSumFactor;
#endif
#ifdef RT_SOURCE_INJECTION
    MyDouble KernelSum_Around_RT_Source;
#endif
#ifdef HYDRO_PRESSURE_SPH
    MyDouble EgyRho;
#endif
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
    MyFloat NV_D[3][3];
    MyFloat NV_A[3][3];
#endif
#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    Vec3<MyFloat> GradRho;
#endif
#if defined(SINK_PARTICLES)
    int Sink_TimeBinGasNeighbor;
#if defined(BH_ACCRETE_NEARESTFIRST) || defined(SINGLE_STAR_TIMESTEPPING)
    MyDouble Sink_dr_to_NearestGasNeighbor;
#endif
#endif
#if defined(TURB_DRIVING) || defined(GRAIN_FLUID)
    Vec3<MyDouble> GasVel;
#endif
#if defined(GRAIN_FLUID)
    MyDouble Gas_InternalEnergy;
#if defined(GRAIN_LORENTZFORCE)
    Vec3<MyDouble> Gas_B;
#endif
#endif
}
 *DATARESULT_NAME, *DATAOUT_NAME;

/*! this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
void hydrokerneldensity_out2particle(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    int j,k;
    ASSIGN_ADD(P.NumNgb[i], out->Ngb, mode);
    ASSIGN_ADD(P.DrkernNgbFactor[i], out->DrkernNgb, mode);
    ASSIGN_ADD(P.Particle_DivVel[i], out->Particle_DivVel,   mode);

    if(P.Type[i] == 0)
    {
        ASSIGN_ADD(CellP.Density[i], out->Rho, mode);
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
        ASSIGN_ADD(CellP.ParticleVel[i], out->ParticleVel, mode);
#endif
        for(k=0;k<6;k++) {ASSIGN_ADD(CellP.NV_T[i].data[k], out->NV_T.data[k], mode);}
        ASSIGN_ADD(CellP.NV_T_face_weights[i], out->NV_T_face_weights, mode);

#ifdef HYDRO_SPH
        ASSIGN_ADD(CellP.DrkernHydroSumFactor[i], out->DrkernHydroSumFactor, mode);
#endif
#ifdef HYDRO_PRESSURE_SPH
        ASSIGN_ADD(CellP.EgyWtDensity[i],   out->EgyRho,   mode);
#endif
#if defined(TURB_DRIVING)
        ASSIGN_ADD(CellP.SmoothedVel[i], out->GasVel, mode);
#endif
#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        for(k = 0; k < 3; k++)
            for(j = 0; j < 3; j++)
            {
                ASSIGN_ADD(CellP.NV_D[i][k][j], out->NV_D[k][j], mode);
                ASSIGN_ADD(CellP.NV_A[i][k][j], out->NV_A[k][j], mode);
            }
#endif
    } // P.Type[i] == 0 //

#if defined(GRAIN_FLUID)
    if((1 << P.Type[i]) & (GRAIN_PTYPES))
    {
        ASSIGN_ADD(P.Gas_Density[i], out->Rho, mode);
        ASSIGN_ADD(P.Gas_InternalEnergy[i], out->Gas_InternalEnergy, mode);
        ASSIGN_ADD(P.Gas_Velocity[i], out->GasVel, mode);
#if defined(GRAIN_LORENTZFORCE)
        ASSIGN_ADD(P.Gas_B[i], out->Gas_B, mode);
#endif
    }
#endif

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    ASSIGN_ADD(P.DensityAroundParticle[i], out->Rho, mode);
    ASSIGN_ADD(P.GradRho[i], out->GradRho, mode);
#endif

#if defined(RT_SOURCE_INJECTION)
#if defined(RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION)
    if(All.TimeStep == 0) // we only do this on the 0'th timestep, since we haven't done a sink loop yet to get the angle weights we'll use normally
#endif
    if((1 << P.Type[i]) & (RT_SOURCES)) {ASSIGN_ADD(P.KernelSum_Around_RT_Source[i], out->KernelSum_Around_RT_Source, mode);}
#endif

#ifdef SINK_PARTICLES
    if(P.Type[i] == 5)
    {
        if(mode == 0) {P.Sink_TimeBinGasNeighbor[i] = out->Sink_TimeBinGasNeighbor;} else {if(P.Sink_TimeBinGasNeighbor[i] > out->Sink_TimeBinGasNeighbor) {P.Sink_TimeBinGasNeighbor[i] = out->Sink_TimeBinGasNeighbor;}}
#if defined(SINGLE_STAR_TIMESTEPPING)
        if(mode == 0) {P.Sink_dr_to_NearestGasNeighbor[i] = out->Sink_dr_to_NearestGasNeighbor;} else {if(P.Sink_dr_to_NearestGasNeighbor[i] > out->Sink_dr_to_NearestGasNeighbor) {P.Sink_dr_to_NearestGasNeighbor[i] = out->Sink_dr_to_NearestGasNeighbor;}}
#endif
    } /* if(P.Type[i] == 5) */
#endif
}

/*! declare this utility function here now that the relevant structures it uses have been defined */
void density_evaluate_extra_physics_gas(struct INPUT_STRUCT_NAME *local, struct OUTPUT_STRUCT_NAME *out, struct kernel_density *kernel, int j);


/*! This function represents the core of the initial hydro kernel-identification and volume computation. The target particle may either be local, or reside in the communication buffer. */
/*!   -- this subroutine should in general contain no writes to shared memory. for optimization reasons, a couple of such writes have been included here in the sub-code for some sink routines -- those need to be handled with special care, both for thread safety and because of iteration. in general writes to shared memory in density.c are strongly discouraged -- */
int density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int j, n, startnode, numngb_inbox, listindex = 0; double r2, h2, u, mass_j, wk;
    struct kernel_density kernel; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    if(mode == 0) {hydrokerneldensity_particle2in(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    h2 = local.KernelRadius * local.KernelRadius; kernel_hinv(local.KernelRadius, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
#if defined(SINK_PARTICLES)
    out.Sink_TimeBinGasNeighbor = TIMEBINS;
#endif
    DensityNeighborTile dtile;
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads(local.Pos.data_ptr(), local.KernelRadius, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}
            density_tile_gather(dtile, ngblist, numngb_inbox);

            /* ---- Loop fission: Phase 1 ---- */
            /* Pre-compute distances and kernel values for all tile neighbors in a
               tight, branch-free loop that the compiler can auto-vectorize.  When
               the tile path is not active we fall back to the original fused loop. */
            double   r2_pre[DENSITY_TILE_NGB_MAX];
            double   r_pre[DENSITY_TILE_NGB_MAX];
            double   u_pre[DENSITY_TILE_NGB_MAX];
            double   wk_pre[DENSITY_TILE_NGB_MAX];
            double   dwk_pre[DENSITY_TILE_NGB_MAX];
            Vec3<MyDouble> dp_pre[DENSITY_TILE_NGB_MAX];
            int      valid_pre[DENSITY_TILE_NGB_MAX];

            if(dtile.use_tile && numngb_inbox > 0)
            {
#if defined(__clang__)
#pragma clang loop vectorize(enable) interleave(enable)
#elif defined(__GNUC__)
#pragma GCC ivdep
#endif
                for(n = 0; n < numngb_inbox; n++)
                {
                    dp_pre[n] = local.Pos - dtile.pos[n];
                    nearest_xyz(dp_pre[n]);
                    r2_pre[n] = dp_pre[n].norm_sq();
                    int v = (dtile.mass[n] > 0) & (r2_pre[n] < h2);
#ifdef GALSF_SUBGRID_WINDS
                    if(dtile.delay_time[n] > 0 && local.DelayTime <= 0) {v = 0;}
#endif
                    valid_pre[n] = v;
                    double r_n = v ? sqrt(r2_pre[n]) : 0;
                    r_pre[n] = r_n;
                    u_pre[n] = r_n * kernel.hinv;
                    kernel_main_branchless(u_pre[n], kernel.hinv3, kernel.hinv4, &wk_pre[n], &dwk_pre[n]);
                    /* zero out kernel values for invalid neighbors so they contribute nothing */
                    if(!v) { wk_pre[n] = 0; dwk_pre[n] = 0; }
                }
            }

            /* ---- Loop fission: Phase 2 ---- */
            /* Accumulation pass: uses pre-computed distances and kernel values
               from Phase 1 when the tile path is active, otherwise computes
               everything inline (original code path). */
            for(n = 0; n < numngb_inbox; n++)
            {
                if(dtile.use_tile)
                {
                    /* fast path: use pre-computed values */
                    if(!valid_pre[n]) continue;
                    j = dtile.ngb_idx[n];
                    kernel.dp = dp_pre[n];
                    r2 = r2_pre[n];
                    kernel.r = r_pre[n];
                    u = u_pre[n];
                    kernel.wk = wk_pre[n];
                    kernel.dwk = dwk_pre[n];
                    mass_j = dtile.mass[n];
                }
                else
                {
                    /* fallback path: compute inline (no tile) */
                    j = ngblist[n];
#ifdef GALSF_SUBGRID_WINDS
                    if(CellP.DelayTime[j] > 0) {if(local.DelayTime <= 0) {continue;}}
#endif
                    if(P.Mass[j] <= 0) continue;
                    kernel.dp = local.Pos - P.Pos[j];
                    nearest_xyz(kernel.dp);
                    r2 = kernel.dp.norm_sq();
                    if(r2 >= h2) continue;
                    kernel.r = sqrt(r2);
                    u = kernel.r * kernel.hinv;
                    kernel_main_branchless(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk);
                    mass_j = P.Mass[j];
                }

                    kernel.mj_wk = (mass_j * kernel.wk);

                    out.Ngb += kernel.wk;
                    out.Rho += kernel.mj_wk;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                    if(local.Type == 0 && kernel.r==0) {out.ParticleVel += kernel.mj_wk * (dtile.use_tile ? dtile.velpred[n] : CellP.VelPred[j]);} // just the self-contribution //
#endif
#if defined(RT_SOURCE_INJECTION)
#if defined(RT_SINK_ANGLEWEIGHT_PHOTON_INJECTION)
                    if(All.TimeStep == 0) // we only do this on the 0'th timestep, since we haven't done a sink loop yet to get the angle weights we'll use normally
#endif
                    if((1 << local.Type) & (RT_SOURCES)) {out.KernelSum_Around_RT_Source += 1.-u*u;}
#endif
                    out.DrkernNgb += -(NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
#ifdef HYDRO_SPH
                    double mass_eff = mass_j;
#ifdef HYDRO_PRESSURE_SPH
                    mass_eff *= dtile.use_tile ? dtile.ie_pred[n] : CellP.InternalEnergyPred[j];
                    out.EgyRho += kernel.wk * mass_eff;
#endif
                    out.DrkernHydroSumFactor += -mass_eff * (NUMDIMS * kernel.hinv * kernel.wk + u * kernel.dwk);
#endif
                    /* for everything below, we do NOT include the particle self-contribution! */
                    if(kernel.r > 0)
                    {
                        if(local.Type == 0)
                        {
                            wk = kernel.wk; /* MAKE SURE THIS MATCHES CHOICE IN GRADIENTS.c! */
                            /* the weights for the MLS tensor used for gradient estimation */
                            out.NV_T += wk * outer_product(kernel.dp);
                            /* weighted first moments, used for face area estimation */
                            out.NV_T_face_weights += wk * kernel.dp;
                        }
                        Vec3<MyDouble> velpred_j = dtile.use_tile ? dtile.velpred[n] : CellP.VelPred[j];
                        Vec3<MyDouble> pos_j = dtile.use_tile ? dtile.pos[n] : P.Pos[j];
                        kernel.dv = local.Vel - velpred_j;
                        NGB_SHEARBOX_BOUNDARY_VELCORR_(local.Pos,pos_j,kernel.dv,1); /* wrap velocities for shearing boxes if needed */
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                        // do neighbor contribution to smoothed particle velocity here, after wrap, so can account for shearing boxes correctly //
                        out.ParticleVel += kernel.mj_wk * (local.Vel - kernel.dv);
#endif
                        out.Particle_DivVel -= kernel.dwk * dot(kernel.dp, kernel.dv) / kernel.r;
                        /* this is the -particle- divv estimator, which determines how KernelRadius will evolve (particle drift) */

                        density_evaluate_extra_physics_gas(&local, &out, &kernel, j);
                    } // kernel.r > 0
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {hydrokerneldensity_out2particle(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}



/*! this is an extra function to simplify additional computations within the kernel that need to be done as part of the evaluation above */
void density_evaluate_extra_physics_gas(struct INPUT_STRUCT_NAME *local, struct OUTPUT_STRUCT_NAME *out, struct kernel_density *kernel, int j)
{
    kernel->mj_dwk_r = P.Mass[j] * kernel->dwk / kernel->r;

    if(local->Type != 0)
    {
#if defined(GRAIN_FLUID)
        if((1 << local->Type) & (GRAIN_PTYPES))
        {
            out->Gas_InternalEnergy += kernel->mj_wk * CellP.InternalEnergyPred[j];
            out->GasVel += kernel->mj_wk * (local->Vel - kernel->dv);
#if defined(GRAIN_LORENTZFORCE)
            out->Gas_B += kernel->wk * CellP.BPred[j];
#endif
        }
#endif

#if defined(SINK_PARTICLES)
        /* note, we will have some writes to shared memory here for some initializations of 'j' quantities. fortunately these do not depend on previous values of those quantities, so can be done thread-safely with minor edits using the constructs below */
        if(local->Type == 5)
        {
            #pragma omp atomic write
            P.SwallowID[j] = 0;  // this way we don't have to do a global loop over local particles in sink_accretion() to reset these quantities...
#ifdef SINGLE_STAR_SINK_DYNAMICS
            #pragma omp atomic write
            P.SwallowTime[j] = MAX_REAL_NUMBER; // initialize as a large number before looking
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
            if(kernel->r < DMAX(P.KernelRadius[j], SinkParticle_GravityKernelRadius)) {
                #pragma omp atomic write
                P.Sink_Ngb_Flag[j] = 1; // note that this particle is inside of a sink's kernel function
            }
#endif
            short int TimeBin_j = P.TimeBin[j]; if(TimeBin_j < 0) {TimeBin_j = -TimeBin_j - 1;} // need to make sure we correct for the fact that TimeBin is used as a 'switch' here to determine if a particle is active for iteration, otherwise this gives nonsense!
            if(out->Sink_TimeBinGasNeighbor > TimeBin_j) {out->Sink_TimeBinGasNeighbor = TimeBin_j;}
#if defined(SINGLE_STAR_TIMESTEPPING)
            double dr_eff_wtd = Get_Particle_Size(j);
            dr_eff_wtd=sqrt(dr_eff_wtd*dr_eff_wtd + (kernel->r)*(kernel->r)); /* effective distance for Gaussian-type kernel, weighted by density */
            if((dr_eff_wtd < out->Sink_dr_to_NearestGasNeighbor) && (P.Mass[j] > 0)) {out->Sink_dr_to_NearestGasNeighbor = dr_eff_wtd;}
#endif
        }
#endif // SINK_PARTICLES
        
    } else { /* local.Type == 0 */

#if defined(TURB_DRIVING)
        out->GasVel += kernel->mj_wk * (local->Vel - kernel->dv);
#endif

#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
        double wk = kernel->wk;
        out->NV_A[0][0] += (local->Accel[0] - All.cf_a2inv*P.GravAccel[j][0] - CellP.HydroAccel[j][0]) * kernel->dp[0] * wk;
        out->NV_A[0][1] += (local->Accel[0] - All.cf_a2inv*P.GravAccel[j][0] - CellP.HydroAccel[j][0]) * kernel->dp[1] * wk;
        out->NV_A[0][2] += (local->Accel[0] - All.cf_a2inv*P.GravAccel[j][0] - CellP.HydroAccel[j][0]) * kernel->dp[2] * wk;
        out->NV_A[1][0] += (local->Accel[1] - All.cf_a2inv*P.GravAccel[j][1] - CellP.HydroAccel[j][1]) * kernel->dp[0] * wk;
        out->NV_A[1][1] += (local->Accel[1] - All.cf_a2inv*P.GravAccel[j][1] - CellP.HydroAccel[j][1]) * kernel->dp[1] * wk;
        out->NV_A[1][2] += (local->Accel[1] - All.cf_a2inv*P.GravAccel[j][1] - CellP.HydroAccel[j][1]) * kernel->dp[2] * wk;
        out->NV_A[2][0] += (local->Accel[2] - All.cf_a2inv*P.GravAccel[j][2] - CellP.HydroAccel[j][2]) * kernel->dp[0] * wk;
        out->NV_A[2][1] += (local->Accel[2] - All.cf_a2inv*P.GravAccel[j][2] - CellP.HydroAccel[j][2]) * kernel->dp[1] * wk;
        out->NV_A[2][2] += (local->Accel[2] - All.cf_a2inv*P.GravAccel[j][2] - CellP.HydroAccel[j][2]) * kernel->dp[2] * wk;

        out->NV_D[0][0] += kernel->dv[0] * kernel->dp[0] * wk;
        out->NV_D[0][1] += kernel->dv[0] * kernel->dp[1] * wk;
        out->NV_D[0][2] += kernel->dv[0] * kernel->dp[2] * wk;
        out->NV_D[1][0] += kernel->dv[1] * kernel->dp[0] * wk;
        out->NV_D[1][1] += kernel->dv[1] * kernel->dp[1] * wk;
        out->NV_D[1][2] += kernel->dv[1] * kernel->dp[2] * wk;
        out->NV_D[2][0] += kernel->dv[2] * kernel->dp[0] * wk;
        out->NV_D[2][1] += kernel->dv[2] * kernel->dp[1] * wk;
        out->NV_D[2][2] += kernel->dv[2] * kernel->dp[2] * wk;
#endif

    } // Type = 0 check

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
    /* this is here because for the models of BH growth and self-shielding of stars, we
     just need a quick-and-dirty, single-pass approximation for the gradients (the error from
     using this as opposed to the higher-order gradient estimators is small compared to the
     Sobolev approximation): use only for -non-gas- particles */
    out->GradRho += kernel->mj_dwk_r * kernel->dp;
#endif

}




/*! This function computes the local neighbor kernel for each active hydro element, the number of neighbours in the current kernel radius, and the divergence
 * and rotation of the velocity field.  This is used then to compute the effective volume of the element in MFM/MFV/SPH-type methods, which is then used to
 * update volumetric quantities like density and pressure. The routine iterates to attempt to find a target kernel size set adaptively -- see code user guide for details
 */
void density(void)
{
    /* initialize variables used below, in particlar the structures we need to call throughout the iteration */
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second(); MyFloat *Left, *Right; double fac, fac_lim, desnumngb, desnumngbdev; long long ntot;
    int i, npleft, iter=0, redo_particle, particle_set_to_minrkern_flag = 0, particle_set_to_maxrkern_flag = 0;
    Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
    
#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES /* define a variable for below to know which stellar types qualify here */
    int valid_stellar_types = 2+4+8+16, invalid_stellar_types = 1+32; // allow types 1,2,3,4 here //
#if (defined(GRAIN_FLUID) || defined(RADTRANSFER)) && (!defined(GALSF) && !(defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)))
    valid_stellar_types = 16; invalid_stellar_types = 1+2+4+8+32; // -only- type-4 sources in these special problems
#ifdef RADTRANSFER
    invalid_stellar_types = 64; valid_stellar_types = RT_SOURCES; // any valid 'injection' source is allowed
#endif
#ifdef GRAIN_FLUID
    invalid_stellar_types = GRAIN_PTYPES;
#endif
#endif
#endif
    
    /* initialize anything we need to about the active particles before their loop */
    for (int i : ActiveParticleList) {
        if(density_isactive(i)) {
            Left[i] = Right[i] = 0;
#ifdef SINK_PARTICLES
            P.SwallowID[i] = 0;
#ifdef SINGLE_STAR_SINK_DYNAMICS
            P.SwallowTime[i] = MAX_REAL_NUMBER;
#endif
#if (SINGLE_STAR_SINK_FORMATION & 8)
            P.Sink_Ngb_Flag[i] = 0;
#endif
#endif
            double maxsoft = All.MaxKernelRadius; /* before the first pass, need to ensure the particles do not exceed the maximum KernelRadius allowed */
#if defined(DO_DENSITY_AROUND_NONGAS_PARTICLES) && defined(GALSF)
            if( ((1 << P.Type[i]) & (valid_stellar_types)) && !((1 << P.Type[i]) & (invalid_stellar_types)) ) {maxsoft = 2.0 / (UNIT_LENGTH_IN_KPC*All.cf_atime);}
#endif
#ifdef SINK_PARTICLES
            if(P.Type[i] == 5) {maxsoft = All.SinkMaxAccretionRadius / All.cf_atime;}  // MaxAccretionRadius is now defined in params.txt in PHYSICAL units
#endif
            if((P.KernelRadius[i] < 0) || !isfinite(P.KernelRadius[i]) || (P.KernelRadius[i] > 0.99*maxsoft)) {P.KernelRadius[i] = 0.99*maxsoft;} /* don't set to exactly maxsoft because our looping below won't treat this correctly */
            
        }} /* done with intial zero-out loop */

    /* allocate buffers to arrange communication */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do
    {
        #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */

        /* do check on whether we have enough neighbors, and iterate for density-rkern solution */
        double tstart = my_second(), tend;
        npleft = 0; for (int i : ActiveParticleList)
        {
            desnumngb = All.DesNumNgb; desnumngbdev = All.MaxNumNgbDeviation;
            /* in the initial timestep and iteration, use a much more strict tolerance for the neighbor number */
            if(All.Time==All.TimeBegin) {if(All.MaxNumNgbDeviation > 0.05) desnumngbdev=0.05;}
            MyDouble desnumngbdev_0 = desnumngbdev, Tinv[3][3], ConditionNumber=0; int k,k1,k2; k=0;
            if(density_isactive(i))
            {
                if(P.NumNgb[i] > 0)
                {
                    P.DrkernNgbFactor[i] *= P.KernelRadius[i] / (NUMDIMS * P.NumNgb[i]);
                    P.Particle_DivVel[i] /= P.NumNgb[i];
                    /* spherical volume of the Kernel (use this to normalize 'effective neighbor number') */
                    P.NumNgb[i] *= VOLUME_NORM_COEFF_FOR_NDIMS * pow(P.KernelRadius[i],NUMDIMS);
                } else {
                    P.NumNgb[i] = P.DrkernNgbFactor[i] = P.Particle_DivVel[i] = 0;
                }
#if defined(ADAPTIVE_GRAVSOFT_FORALL) /* if particle is AGS-active and non-gas, set DivVel to zero because it will be reset in ags_rkern routine */
                if(ags_density_isactive(i) && (P.Type[i] > 0)) {P.Particle_DivVel[i] = 0;}
#endif

                // inverse of fluid volume element (to satisfy constraint implicit in Lagrange multipliers)
                if(P.DrkernNgbFactor[i] > -0.9) {P.DrkernNgbFactor[i] = 1 / (1 + P.DrkernNgbFactor[i]);} else {P.DrkernNgbFactor[i] = 1;} /* note: this would be -1 if only a single particle at zero lag is found */
                P.Particle_DivVel[i] *= P.DrkernNgbFactor[i];

                double dimless_face_leak=0; MyDouble NV_T_prev[6]; NV_T_prev[0]=CellP.NV_T[i][0][0]; NV_T_prev[1]=CellP.NV_T[i][1][1]; NV_T_prev[2]=CellP.NV_T[i][2][2]; NV_T_prev[3]=CellP.NV_T[i][0][1]; NV_T_prev[4]=CellP.NV_T[i][0][2]; NV_T_prev[5]=CellP.NV_T[i][1][2];
                if(P.Type[i] == 0) /* invert the NV_T matrix we just measured */
                {
                    /* use the single-moment terms of NV_T to construct the faces one would have if the system were perfectly symmetric in reconstruction 'from both sides' */
                    double V_i = VOLUME_NORM_COEFF_FOR_NDIMS * pow(P.KernelRadius[i],NUMDIMS) / P.NumNgb[i], dx_i = pow(V_i , 1./NUMDIMS); // this is the effective volume which will be used below
                    dx_i = sqrt(V_i * CellP.NV_T[i].trace()); // this is the sqrt of the weighted sum of (w*r^2)
                    double Face_Area_OneSided_Estimator_in[3]={0}, Face_Area_OneSided_Estimator_out[3]={0}; Face_Area_OneSided_Estimator_in[0]=CellP.NV_T_face_weights[i][0]; Face_Area_OneSided_Estimator_in[1]=CellP.NV_T_face_weights[i][1]; Face_Area_OneSided_Estimator_in[2]=CellP.NV_T_face_weights[i][2];
                    double dimensional_NV_T_normalizer = pow( P.KernelRadius[i] , 2-NUMDIMS ); /* this has the same dimensions as NV_T here */
                    double NV_T_local[3][3]; /* local working copy for inversion: avoids passing SymmetricTensor2 to matrix_invert_ndims and avoids double-applying normalizer to off-diagonal elements */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {NV_T_local[k1][k2] = CellP.NV_T[i][k1][k2] / dimensional_NV_T_normalizer;}} /* dimensionless copy */
                    /* Also, we want to be able to calculate the condition number of the matrix to be inverted, since
                        this will tell us how robust our procedure is (and let us know if we need to expand the neighbor number */
                    double ConditionNumber_threshold = 10. * CONDITION_NUMBER_DANGER; /* set a threshold condition number - above this we will 'pre-condition' the matrix for better behavior */
                    double trace_initial = NV_T_local[0][0] + NV_T_local[1][1] + NV_T_local[2][2]; /* initial trace of this symmetric, positive-definite matrix; used below as a characteristic value for adding the identity */
                    double conditioning_term_to_add = 1.05 * (trace_initial / NUMDIMS) / ConditionNumber_threshold; /* this will be added as a test value if the code does not reach the desired condition number */
                    while(1)
                    {
                        ConditionNumber = matrix_invert_ndims(NV_T_local, Tinv);
                        if(ConditionNumber < ConditionNumber_threshold) {break;}
                        for(k1=0;k1<NUMDIMS;k1++) {NV_T_local[k1][k1] += conditioning_term_to_add;} /* add the conditioning term which should make the matrix better-conditioned for subsequent use */
                        conditioning_term_to_add *= 1.2; /* multiply the conditioning term so it will grow and eventually satisfy our criteria */
                    }
                    for(k1=0;k1<3;k1++) {for(k2=k1;k2<3;k2++) {CellP.NV_T[i][k1][k2] = Tinv[k1][k2] / dimensional_NV_T_normalizer;}} /* re-insert normalization correctly */
                    /* now NV_T holds the inverted matrix elements, for use in hydro */
                    for(k1=0;k1<3;k1++) {for(k2=0;k2<3;k2++) {Face_Area_OneSided_Estimator_out[k1] += 2.*V_i*CellP.NV_T[i][k1][k2]*Face_Area_OneSided_Estimator_in[k2];}} /* calculate mfm/mfv areas that we would have by default, if both sides of reconstruction were symmetric */
                    for(k1=0;k1<3;k1++) {dimless_face_leak += fabs(Face_Area_OneSided_Estimator_out[k1]) / NUMDIMS;} // average of absolute values
#ifdef HYDRO_KERNEL_SURFACE_VOLCORR
                    double closure_asymm=0; for(k1=0;k1<3;k1++) {closure_asymm += Face_Area_OneSided_Estimator_in[k1]*Face_Area_OneSided_Estimator_in[k1];}
                    double particle_inverse_volume = P.NumNgb[i] / ( VOLUME_NORM_COEFF_FOR_NDIMS * pow(P.KernelRadius[i],NUMDIMS) );
                    closure_asymm = sqrt(closure_asymm) / (P.KernelRadius[i] * particle_inverse_volume); // dimensionnless measure of asymmetry in kernel
                    CellP.FaceClosureError[i] = DMIN(DMAX(1.0259-2.52444*closure_asymm,0.344301),1.); // correction factor for 'missing' volume assuming a wendland C2 kernel and a sharp surface from Reinhardt & Stadel 2017 (arXiv:1701.08296)
#else
                    CellP.FaceClosureError[i] = dimless_face_leak / (2.*NUMDIMS*pow(dx_i,NUMDIMS-1));
#endif
                } // P.Type[i] == 0 //

                /* now check whether we had enough neighbours */
                double ncorr_ngb = 1.0;
                double cn=1;
                double c0 = 0.1 * (double)CONDITION_NUMBER_DANGER;
                if(P.Type[i]==0)
                {
                    /* use the previous timestep condition number to correct how many neighbors we should use for stability */
                    if((iter==0)&&(ConditionNumber>CellP.ConditionNumber[i])&&(CellP.ConditionNumber[i]>0))
                    {
                        /* if we find ourselves with a sudden increase in condition number - check if we have a reasonable
                            neighbor number for the previous iteration, and if so, use the new (larger) correction */
                        ncorr_ngb=1; cn=CellP.ConditionNumber[i]; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
                        double dn_ngb = fabs(P.NumNgb[i]-All.DesNumNgb*ncorr_ngb)/(desnumngbdev_0*ncorr_ngb);
                        ncorr_ngb=1; cn=ConditionNumber; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
                        double dn_ngb_alt = fabs(P.NumNgb[i]-All.DesNumNgb*ncorr_ngb)/(desnumngbdev_0*ncorr_ngb);
                        dn_ngb = DMIN(dn_ngb,dn_ngb_alt);
                        if(dn_ngb < 10.0) CellP.ConditionNumber[i] = ConditionNumber;
                    }
                    ncorr_ngb=1; cn=CellP.ConditionNumber[i]; if(cn>c0) {ncorr_ngb=sqrt(1.0+(cn-c0)/((double)CONDITION_NUMBER_DANGER));} if(ncorr_ngb>2) ncorr_ngb=2;
#if !defined(HYDRO_KERNEL_SURFACE_VOLCORR)
                    double d00=0.35; if(CellP.FaceClosureError[i] > d00) {ncorr_ngb = DMAX(ncorr_ngb , DMIN(CellP.FaceClosureError[i]/d00 , 2.));}
#endif
                }
                desnumngb = All.DesNumNgb * ncorr_ngb;
                desnumngbdev = desnumngbdev_0 * ncorr_ngb;
                /* allow the neighbor tolerance to gradually grow as we iterate, so that we don't spend forever trapped in a narrow iteration */
#if !defined(EOS_ELASTIC)
                if(iter > 1) {desnumngbdev = DMIN( 0.25*desnumngb , desnumngbdev * exp(0.1*log(desnumngb/(16.*desnumngbdev))*(double)iter) );}
#endif

#ifdef SINK_PARTICLES
                if(P.Type[i] == 5)
                {
                    desnumngb = All.DesNumNgb * All.SinkNgbFactor;
#ifdef SINGLE_STAR_SINK_DYNAMICS
                    desnumngbdev = (All.SinkNgbFactor+1);
#else
                    desnumngbdev = 4 * (All.SinkNgbFactor+1);
#endif
                }
#endif

#ifdef GRAIN_FLUID /* for the grains, we only need to estimate neighboring gas properties, we don't need to worry about condition numbers or conserving an exact neighbor number */
                if((1 << P.Type[i]) & (GRAIN_PTYPES))
                {
                    desnumngb = All.DesNumNgb; desnumngbdev = All.DesNumNgb / 4;
#if defined(GRAIN_BACKREACTION)
                    desnumngbdev = desnumngbdev_0;
#endif
                }
#endif

                double minsoft = All.MinKernelRadius;
                double maxsoft = All.MaxKernelRadius;

#ifdef DO_DENSITY_AROUND_NONGAS_PARTICLES
                /* use a much looser check for N_neighbors when the central point is a star particle,
                 since the accuracy is limited anyways to the coupling efficiency -- the routines use their
                 own estimators+neighbor loops, anyways, so this is just to get some nearby particles */
                if( ((1 << P.Type[i]) & (valid_stellar_types)) && !((1 << P.Type[i]) & (invalid_stellar_types)) )
                {
                    desnumngb = All.DesNumNgb;
#if defined(RT_SOURCE_INJECTION)
                    if(desnumngb < 64.0) {desnumngb = 64.0;} // we do want a decent number to ensure the area around the particle is 'covered'
#endif
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                    desnumngb = DMAX(desnumngb , 128); // we do want a decent number to ensure the area around the particle is 'covered'
                    if(KERNEL_FUNCTION > 3) {desnumngb = DMAX(desnumngb, 256);}
#endif
#ifdef GALSF
                    if(desnumngb < 64.0) {desnumngb = 64.0;} // we do want a decent number to ensure the area around the particle is 'covered'
                    // if we're finding this for feedback routines, there isn't any good reason to search beyond a modest physical radius //
                    double unitlength_in_kpc=UNIT_LENGTH_IN_KPC*All.cf_atime;
                    maxsoft = 2.0 / unitlength_in_kpc;
#if defined(GALSF_FB_FIRE_STELLAREVOLUTION) && defined(SINK_PARTICLES) && (defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL))
                    if((P.SNe_ThisTimeStep[i]>0) || (P.MassReturn_ThisTimeStep[i]>0) || (All.Time==All.TimeBegin)) {maxsoft=2.0/unitlength_in_kpc;} else {maxsoft=0.1/unitlength_in_kpc;};
#endif
#endif
                    desnumngbdev = desnumngb / 2; // enforcing exact number not important
                }
#endif

#ifdef SINK_PARTICLES
                if(P.Type[i] == 5) {maxsoft = All.SinkMaxAccretionRadius / All.cf_atime;}  // MaxAccretionRadius is now defined in params.txt in PHYSICAL units
#ifdef SINGLE_STAR_SINK_DYNAMICS
		        if(P.Type[i] == 5) {minsoft = SinkParticle_GravityKernelRadius;} // we should always find all neighbours within the softening kernel/accretion radius, which is a lower bound on the accretion radius
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
                if(P.Type[i] == 5) {minsoft = DMAX(minsoft, DMIN(P.SinkRadius[i] , 0.1*SinkParticle_GravityKernelRadius));}
#endif
#endif
#endif

                redo_particle = 0;

                /* check if we are in the 'normal' range between the max/min allowed values */
                if((P.NumNgb[i] < (desnumngb - desnumngbdev) && P.KernelRadius[i] < 0.999*maxsoft) ||
                   (P.NumNgb[i] > (desnumngb + desnumngbdev) && P.KernelRadius[i] > 1.001*minsoft))
                    {redo_particle = 1;}

                /* check maximum kernel size allowed */
                particle_set_to_maxrkern_flag = 0;
                if((P.KernelRadius[i] >= 0.999*maxsoft) && (P.NumNgb[i] < (desnumngb - desnumngbdev)))
                {
                    redo_particle = 0;
                    if(P.KernelRadius[i] == maxsoft)
                    {
                        /* iteration at the maximum value is already complete */
                        particle_set_to_maxrkern_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the maximum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        P.KernelRadius[i] = maxsoft;
                        particle_set_to_maxrkern_flag = 1;
                    }
                }

                /* check minimum kernel size allowed */
                particle_set_to_minrkern_flag = 0;
                if((P.KernelRadius[i] <= 1.001*minsoft) && (P.NumNgb[i] > (desnumngb + desnumngbdev)))
                {
                    redo_particle = 0;
                    if(P.KernelRadius[i] == minsoft)
                    {
                        /* this means we've already done an iteration with the MinKernelRadius value, so the
                         neighbor weights, etc, are not going to be wrong; thus we simply stop iterating */
                        particle_set_to_minrkern_flag = 0;
                    } else {
                        /* ok, the particle needs to be set to the minimum, and (if gas) iterated one more time */
                        redo_particle = 1;
                        P.KernelRadius[i] = minsoft;
                        particle_set_to_minrkern_flag = 1;
                    }
                }

#ifdef GALSF
                if((All.ComovingIntegrationOn)&&(All.Time>All.TimeBegin))
                {
                    if((P.Type[i]==4)&&(iter>1)&&(P.NumNgb[i]>4)&&(P.NumNgb[i]<100)&&(redo_particle==1)) {redo_particle=0;}
                }
#endif

                if((redo_particle==0)&&(P.Type[i] == 0))
                {
                    /* ok we have reached the desired number of neighbors: save the condition number for next timestep */
                    if(ConditionNumber > 1e6 * (double)CONDITION_NUMBER_DANGER) {
                        PRINT_WARNING("Condition number=%g CNum_prevtimestep=%g CNum_danger=%g iter=%d Num_Ngb=%g desnumngb=%g KernelRadius=%g KernelRadius_min=%g KernelRadius_max=%g \n i=%d task=%d ID=%llu Type=%d KernelRadius=%g Drkern=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)  \n NVT=%.17g/%.17g/%.17g %.17g/%.17g/%.17g %.17g/%.17g/%.17g NVT_inv=%.17g/%.17g/%.17g %.17g/%.17g/%.17g %.17g/%.17g/%.17g ",
                               ConditionNumber,CellP.ConditionNumber[i],CONDITION_NUMBER_DANGER,iter,P.NumNgb[i],desnumngb,P.KernelRadius[i],All.MinKernelRadius,All.MaxKernelRadius, i, ThisTask,
                               (unsigned long long) P.ID[i], P.Type[i], P.KernelRadius[i], P.DrkernNgbFactor[i], Left[i], Right[i],
                               (float) P.NumNgb[i], Right[i] - Left[i], particle_set_to_maxrkern_flag, particle_set_to_minrkern_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P.Pos[i][0], P.Pos[i][1], P.Pos[i][2],
                               CellP.NV_T[i][0][0],CellP.NV_T[i][0][1],CellP.NV_T[i][0][2],CellP.NV_T[i][1][0],CellP.NV_T[i][1][1],CellP.NV_T[i][1][2],CellP.NV_T[i][2][0],CellP.NV_T[i][2][1],CellP.NV_T[i][2][2],
                               NV_T_prev[0],NV_T_prev[3],NV_T_prev[4],NV_T_prev[3],NV_T_prev[1],NV_T_prev[5],NV_T_prev[4],NV_T_prev[5],NV_T_prev[2]);}
                    CellP.ConditionNumber[i] = ConditionNumber;
                }

                if(redo_particle)
                {
                    if(iter >= MAXITER - 10)
                    {
                        PRINT_WARNING("i=%d task=%d ID=%llu iter=%d Type=%d KernelRadius=%g Drkern=%g Left=%g Right=%g Ngbs=%g Right-Left=%g maxh_flag=%d minh_flag=%d  minsoft=%g maxsoft=%g desnum=%g desnumtol=%g redo=%d pos=(%g|%g|%g)",
                               i, ThisTask, (unsigned long long) P.ID[i], iter, P.Type[i], P.KernelRadius[i], P.DrkernNgbFactor[i], Left[i], Right[i],
                               (float) P.NumNgb[i], Right[i] - Left[i], particle_set_to_maxrkern_flag, particle_set_to_minrkern_flag, minsoft,
                               maxsoft, desnumngb, desnumngbdev, redo_particle, P.Pos[i][0], P.Pos[i][1], P.Pos[i][2]);
                    }

                    /* need to redo this particle */
                    npleft++;

                    if(Left[i] > 0 && Right[i] > 0)
                        if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                        {
                            /* this one should be ok */
                            npleft--;
                            P.TimeBin[i] = -P.TimeBin[i] - 1;	/* Mark as inactive */
                            CellP.ConditionNumber[i] = ConditionNumber;
                            continue;
                        }

                    if((particle_set_to_maxrkern_flag==0)&&(particle_set_to_minrkern_flag==0))
                    {
                        if(P.NumNgb[i] < (desnumngb - desnumngbdev)) {Left[i] = DMAX(P.KernelRadius[i], Left[i]);}
                        else
                        {
                            if(Right[i] != 0) {if(P.KernelRadius[i] < Right[i]) {Right[i] = P.KernelRadius[i];}} else {Right[i] = P.KernelRadius[i];}
                        }

                        // right/left define upper/lower bounds from previous iterations
                        if(Right[i] > 0 && Left[i] > 0)
                        {
                            // geometric interpolation between right/left //
                            double maxjump=0;
                            if(iter>1) {maxjump = 0.2*log(Right[i]/Left[i]);}
                            if(P.NumNgb[i] > 1)
                            {
                                double jumpvar = P.DrkernNgbFactor[i] * log( desnumngb / P.NumNgb[i] ) / NUMDIMS;
                                if(iter>1) {if(fabs(jumpvar) < maxjump) {if(jumpvar<0) {jumpvar=-maxjump;} else {jumpvar=maxjump;}}}
                                P.KernelRadius[i] *= exp(jumpvar);
                            } else {
                                P.KernelRadius[i] *= 2.0;
                            }
                            if((P.KernelRadius[i]<Right[i])&&(P.KernelRadius[i]>Left[i]))
                            {
                                if(iter > 1)
                                {
                                    double hfac = exp(maxjump);
                                    if(P.KernelRadius[i] > Right[i] / hfac) {P.KernelRadius[i] = Right[i] / hfac;}
                                    if(P.KernelRadius[i] < Left[i] * hfac) {P.KernelRadius[i] = Left[i] * hfac;}
                                }
                            } else {
                                if(P.KernelRadius[i]>Right[i]) P.KernelRadius[i]=Right[i];
                                if(P.KernelRadius[i]<Left[i]) P.KernelRadius[i]=Left[i];
                                P.KernelRadius[i] = pow(P.KernelRadius[i] * Left[i] * Right[i] , 1.0/3.0);
                            }
                        }
                        else
                        {
                            if(Right[i] == 0 && Left[i] == 0)
                            {
                                char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
                                snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "Right[i] == 0 && Left[i] == 0 && P.KernelRadius[i]=%g\n", P.KernelRadius[i]); terminate(buf);
                            }

                            if(Right[i] == 0 && Left[i] > 0)
                            {
                                if (P.NumNgb[i] > 1)
                                    {fac_lim = log( desnumngb / P.NumNgb[i] ) / NUMDIMS;} // this would give desnumgb if constant density (+0.231=2x desnumngb)
                                else
                                    {fac_lim = 1.4;} // factor ~66 increase in N_NGB in constant-density medium

                                if((P.NumNgb[i] < 2*desnumngb)&&(P.NumNgb[i] > 0.1*desnumngb))
                                {
                                    double slope = P.DrkernNgbFactor[i];
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4) {if(P.DrkernNgbFactor[i]==1) {fac *= 10;}} // tries to help with being trapped in small steps

                                    if(fac < fac_lim+0.231)
                                    {
                                        P.KernelRadius[i] *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                    {
                                        P.KernelRadius[i] *= exp(fac_lim+0.231);
                                        // fac~0.26 leads to expected doubling of number if density is constant,
                                        //   insert this limiter here b/c we don't want to get *too* far from the answer (which we're close to)
                                    }
                                }
                                else
                                    {P.KernelRadius[i] *= exp(fac_lim);} // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }

                            if(Right[i] > 0 && Left[i] == 0)
                            {
                                if(P.NumNgb[i] > 1)
                                    {fac_lim = log( desnumngb / P.NumNgb[i] ) / NUMDIMS;} // this would give desnumgb if constant density (-0.231=0.5x desnumngb)
                                else
                                    {fac_lim = 1.4;} // factor ~66 increase in N_NGB in constant-density medium

                                if(fac_lim < -1.535) {fac_lim = -1.535;} // decreasing N_ngb by factor ~100

                                if((P.NumNgb[i] < 2*desnumngb)&&(P.NumNgb[i] > 0.1*desnumngb))
                                {
                                    double slope = P.DrkernNgbFactor[i];
                                    if(iter>2 && slope<1) slope = 0.5*(slope+1);
                                    fac = fac_lim * slope; // account for derivative in making the 'corrected' guess
                                    if(iter>=4) {if(P.DrkernNgbFactor[i]==1) {fac *= 10;}} // tries to help with being trapped in small steps

                                    if(fac > fac_lim-0.231)
                                    {
                                        P.KernelRadius[i] *= exp(fac); // more expensive function, but faster convergence
                                    }
                                    else
                                        {P.KernelRadius[i] *= exp(fac_lim-0.231);} // limiter to prevent --too-- far a jump in a single iteration
                                }
                                else
                                    {P.KernelRadius[i] *= exp(fac_lim);} // here we're not very close to the 'right' answer, so don't trust the (local) derivatives
                            }
                        } // closes if[particle_set_to_max/minrkern_flag]
                    } // closes redo_particle
                    /* resets for max/min values */
                    if(P.KernelRadius[i] < minsoft) {P.KernelRadius[i] = minsoft;}
                    if(particle_set_to_minrkern_flag==1) {P.KernelRadius[i] = minsoft;}
                    if(P.KernelRadius[i] > maxsoft) {P.KernelRadius[i] = maxsoft;}
                    if(particle_set_to_maxrkern_flag==1) {P.KernelRadius[i] = maxsoft;}
                }
                else {P.TimeBin[i] = -P.TimeBin[i] - 1;}	/* Mark as inactive */
            } //  if(density_isactive(i))
        } // npleft = 0; for (int i : ActiveParticleList)

        tend = my_second();
        timecomp += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 10) {PRINT_STATUS("ngb iteration %d: need to repeat for %d%09d particles", iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
            if(iter > MAXITER) {printf("failed to converge in neighbour iteration in density()\n"); fflush(stdout); endrun(1155);}
        }
    }
    while(ntot > 0);

    /* iteration is done - de-malloc everything now */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    myfree(Right); myfree(Left);

    /* mark as active again */
    for (int i : ActiveParticleList)
    {
        if(P.TimeBin[i] < 0) {P.TimeBin[i] = -P.TimeBin[i] - 1;}
    }


    /* now that we are DONE iterating to find rkern, we can do the REAL final operations on the results
     ( any quantities that only need to be evaluated once, on the final iteration --
     won't save much b/c the real cost is in the neighbor loop for each particle, but it's something )
     -- also, some results (for example, viscosity suppression below) should not be calculated unless
     the quantities are 'stabilized' at their final values -- */
    for (int i : ActiveParticleList)
    {
        if(density_isactive(i))
        {
            if(P.Type[i] == 0 && P.Mass[i] > 0)
            {
                if(CellP.Density[i] > 0)
                {
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
                    /* set motion of the mesh-generating points */
#if (HYDRO_FIX_MESH_MOTION==4)
                    set_mesh_motion(i); // use user-specified analytic function to define mesh motions //
#elif ((HYDRO_FIX_MESH_MOTION==5)||(HYDRO_FIX_MESH_MOTION==6))
                    double eps_pvel = 0.3; // normalization for how much 'weight' to give to neighbors (unstable if >=0.5)
                    CellP.ParticleVel[i] = CellP.VelPred[i] * (1.-eps_pvel) + CellP.ParticleVel[i] * (eps_pvel/CellP.Density[i]); // assign mixture velocity
#elif (HYDRO_FIX_MESH_MOTION==7)
                    CellP.ParticleVel[i] = CellP.VelPred[i]; // move with fluid
#endif
#endif

#ifdef HYDRO_SPH
#ifdef HYDRO_PRESSURE_SPH
                    if(CellP.InternalEnergyPred[i] > 0)
                    {
                        CellP.EgyWtDensity[i] /= CellP.InternalEnergyPred[i];
                    } else {
                        CellP.EgyWtDensity[i] = 0;
                    }
#endif
                    /* need to divide by the sum of x_tilde=1, i.e. numden_ngb */
                    if((P.KernelRadius[i] > 0)&&(P.NumNgb[i] > 0))
                    {
                        double numden_ngb = P.NumNgb[i] / ( VOLUME_NORM_COEFF_FOR_NDIMS * pow(P.KernelRadius[i],NUMDIMS) );
                        CellP.DrkernHydroSumFactor[i] *= P.KernelRadius[i] / (NUMDIMS * numden_ngb);
                        CellP.DrkernHydroSumFactor[i] *= -P.DrkernNgbFactor[i]; /* now this is ready to be called in hydro routine */
                    } else {
                        CellP.DrkernHydroSumFactor[i] = 0;
                    }
#endif


#if defined(SPHAV_CD10_VISCOSITY_SWITCH)
                    int k1, k2;
                    for(k1 = 0; k1 < 3; k1++)
                        for(k2 = 0; k2 < 3; k2++)
                        {
                            CellP.NV_D[i][k2][k1] *= All.cf_a2inv; // converts to physical velocity/length
                            CellP.NV_A[i][k2][k1] /= All.cf_atime; // converts to physical accel/length
                        }
                    // all quantities below in this block should now be in proper PHYSICAL units, for subsequent operations //
                    double dtDV[3][3], A[3][3], V[3][3], S[3][3];
                    for(k1=0;k1<3;k1++)
                        for(k2=0;k2<3;k2++)
                        {
                            V[k1][k2] = CellP.NV_D[i][k1][0]*CellP.NV_T[i][0][k2] + CellP.NV_D[i][k1][1]*CellP.NV_T[i][1][k2] + CellP.NV_D[i][k1][2]*CellP.NV_T[i][2][k2];
                            A[k1][k2] = CellP.NV_A[i][k1][0]*CellP.NV_T[i][0][k2] + CellP.NV_A[i][k1][1]*CellP.NV_T[i][1][k2] + CellP.NV_A[i][k1][2]*CellP.NV_T[i][2][k2];
                        }
                    CellP.NV_DivVel[i] = V[0][0] + V[1][1] + V[2][2];
                    CellP.NV_trSSt[i] = 0;
                    for(k1=0;k1<3;k1++)
                        for(k2=0;k2<3;k2++)
                        {
                            dtDV[k1][k2] = A[k1][k2] - (V[k1][0]*V[0][k2] + V[k1][1]*V[1][k2] + V[k1][2]*V[2][k2]);
                            /* S = 0.5*(V+V_transpose) - delta_ij*div_v/3 */
                            S[k1][k2] = 0.5 * (V[k1][k2] + V[k2][k1]);
                            if(k2==k1) S[k1][k2] -= CellP.NV_DivVel[i] / NUMDIMS;
                            /* Trace[S*S_transpose] = SSt[0][0]+SSt[1][1]+SSt[2][2] = |S|^2 = sum(Sij^2) */
                            CellP.NV_trSSt[i] += S[k1][k2]*S[k1][k2];
                        }
                    CellP.NV_dt_DivVel[i] = dtDV[0][0] + dtDV[1][1] + dtDV[2][2];
#endif


#if defined(TURB_DRIVING)
                    if(CellP.Density[i] > 0)
                    {
                        CellP.SmoothedVel[i] /= CellP.Density[i];
                    } else {
                        CellP.SmoothedVel[i] = {};
                    }
#endif
                }

#ifndef HYDRO_SPH
                if((P.KernelRadius[i] > 0)&&(P.NumNgb[i] > 0))
                {
                    CellP.Density[i] = P.Mass[i] * P.NumNgb[i] / ( VOLUME_NORM_COEFF_FOR_NDIMS * pow(P.KernelRadius[i],NUMDIMS) ); // divide mass by volume
                } else {
                    if(P.KernelRadius[i] <= 0)
                    {
                        CellP.Density[i] = 0; // in this case, give up, no meaningful volume
                    } else {
                        CellP.Density[i] = P.Mass[i] / ( VOLUME_NORM_COEFF_FOR_NDIMS * pow(P.KernelRadius[i],NUMDIMS) ); // divide mass (lone particle) by volume
                    }
                }
#endif
                double Volume_0; Volume_0 = P.Mass[i] / CellP.Density[i]; // save for potential later use
#if defined(HYDRO_KERNEL_SURFACE_VOLCORR)
                CellP.Density[i] /= CellP.FaceClosureError[i]; // correct volume of the cell based on the free surface correction above
                CellP.FaceClosureError[i] = Volume_0;
#endif
#ifdef HYDRO_EXPLICITLY_INTEGRATE_VOLUME
                Volume_0 = P.Mass[i] / CellP.Density[i];
                if(All.Time == All.TimeBegin) {CellP.Density_ExplicitInt[i] = CellP.Density[i];} // set initial value to density calculated above
                    else {CellP.Density[i] = CellP.Density_ExplicitInt[i];} // set to explicitly-evolved density field
                CellP.FaceClosureError[i] = Volume_0;
#endif
#ifdef HYDRO_VOLUME_CORRECTIONS
                CellP.Volume_1[i] = CellP.Volume_0[i] = Volume_0; // initialize this value for use in the correction loop, and in case this is not set in the subsequent loop because of inactivity, set this first to the zeroth-order estimator
#endif
                set_eos_pressure(i);		// should account for density independent pressure

            } // P.Type[i] == 0


#if defined(GRAIN_FLUID)
            if((1 << P.Type[i]) & (GRAIN_PTYPES))
            {
                int k;
                if(P.Gas_Density[i] > 0)
                {
                    P.Gas_InternalEnergy[i] /= P.Gas_Density[i];
                    P.Gas_Velocity[i] /= P.Gas_Density[i];
                } else {
                    P.Gas_InternalEnergy[i] = 0;
                    P.Gas_Velocity[i] = {};
#if defined(GRAIN_LORENTZFORCE)
                    P.Gas_B[i] = {};
#endif
                }
            }
#endif

         /* finally, convert NGB to the more useful format, NumNgb^(1/NDIMS),
            which we can use to obtain the corrected particle sizes. Because of how this number is used above, we --must-- make
            sure that this operation is the last in the loop here */
            if(P.NumNgb[i] > 0) {P.NumNgb[i]=pow(P.NumNgb[i],1./NUMDIMS);} else {P.NumNgb[i]=0;}

#if defined(MAGNETIC)
            if(P.Type[i] == 0) {
                if(CellP.recent_refinement_flag[i] == 1) {
                    CellP.BPred[i] = CellP.B[i] = CellP.BField_prerefinement[i] * (P.Mass[i] / CellP.Density[i]); // reset B-fields to desired values given the conserved variable is VB, after refinement or de-refinement step
                    CellP.BField_prerefinement[i] = {}; // reset this variable to null
                    }}
#endif
            if(P.Type[i] == 0) {CellP.recent_refinement_flag[i] = 0;} // reset this flag after density re-computation
            
        } // density_isactive(i)
        
#if defined(SINK_WIND_SPAWN_SET_BFIELD_POLTOR) /* re-assign magnetic fields after getting the correct density for newly-spawned cells when these options are enabled */
        if(P.Type[i]==0) {if(P.ID[i]==All.SpawnedWindCellID && CellP.IniDen[i]<0) {CellP.IniDen[i]=CellP.Density[i]; CellP.BPred[i]=CellP.B[i]=CellP.IniB[i]*((All.UnitMagneticField_in_gauss/UNIT_B_IN_GAUSS)*(P.Mass[i]/(All.cf_a2inv*CellP.Density[i])));}}
#endif
        
    } // for (int i : ActiveParticleList)

    /* collect some timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_DENSCOMPUTE] += timecomp; CPU_Step[CPU_DENSWAIT] += timewait;
    CPU_Step[CPU_DENSCOMM] += timecomm; CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */



/* Routines for a loop after the iterative density loop needed to find neighbors, etc, once all have converged, to apply additional correction terms to the cell volumes and faces (for those needed -before- the gradients loop because they alter primitive quantities needed for gradients, such as particle densities, pressures, etc.)
    This was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO. */
#ifdef HYDRO_VOLUME_CORRECTIONS

#define CORE_FUNCTION_NAME cellcorrections_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME particle2in_cellcorrections    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME out2particle_cellcorrections  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(GasGrad_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. for current implementation, only cells eligible for gradients and hydro should be called */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/* define structures to use below */
struct INPUT_STRUCT_NAME {Vec3<MyDouble> Pos; MyDouble KernelRadius, Volume_0; int NodeList[NODELISTLENGTH];} *DATAIN_NAME, *DATAGET_NAME;

/* define properties to be sent to nodes */
void particle2in_cellcorrections(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{in->Volume_0=CellP.Volume_0[i]; in->KernelRadius=P.KernelRadius[i]; in->Pos=P.Pos[i];}

/* define output structure to use below */
struct OUTPUT_STRUCT_NAME {MyFloat Volume_1;} *DATARESULT_NAME, *DATAOUT_NAME;

/* define properties to be collected from nodes */
void out2particle_cellcorrections(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{ASSIGN_ADD(CellP.Volume_1[i], out->Volume_1, mode);}

/* core subroutine. this does not write to shared memory. */
int cellcorrections_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int j, n, k, startnode, numngb_inbox, listindex = 0; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* generic variables we always use, and set initial memory */
    if(mode == 0) {particle2in_cellcorrections(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode; /* open it */} /* start usual neighbor tree search */
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.KernelRadius, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}
            for(n=0; n<numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                Vec3<double> dp = local.Pos - P.Pos[j];
                nearest_xyz(dp); // find the closest image in the given box size  //
                double r2 = dp.norm_sq(); // distance
                if(r2 >= P.KernelRadius[j]*P.KernelRadius[j]) {continue;} // need to be inside of 'j's kernel search
                double u,hinv,hinv3,hinv4,wk,dwk; kernel_hinv(P.KernelRadius[j], &hinv, &hinv3, &hinv4); u=sqrt(r2)*hinv; wk=0; dwk=0; // define kernel-needed variables
                kernel_main(u, hinv3, hinv4, &wk, &dwk, -1); // calculate the normal kernel weight 'wk'
                out.Volume_1 += CellP.Volume_0[j] * CellP.Volume_0[j] * wk; // this is the next-order correction to the cell volume quadrature
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode=DATAGET_NAME[target].NodeList[listindex]; if(startnode>=0) {startnode=Nodes[startnode].u.d.nextnode;}}} // handle opening nodes
    } // closes while(startnode >= 0)
    if(mode == 0) {out2particle_cellcorrections(&out, target, 0, 0);} else {DATARESULT_NAME[target] = out;} /* collect the result at the right place */
    return 0; /* done */
}

/* final operations for after the updates are computed */
void cellcorrections_final_operations_and_cleanup(void)
{
    int i; for (int i : ActiveParticleList) { /* check all active elements */
        CONDITIONFUNCTION_FOR_EVALUATION /* ensures only the ones which met our criteria above are actually treated here */
        {
            if(CellP.Volume_1[i] > 0) {CellP.Density[i] = P.Mass[i] / CellP.Volume_1[i];} else {CellP.Volume_1[i] = CellP.Volume_0[i];} // set the updated density. other variables that need volumes will all scale off this, so we can rely on it to inform everything else [if bad value here, revert to the 0th-order volume quadrature]
            set_eos_pressure(i);
        }}
}

/* parent routine which calls the work loop above */
void cellcorrections_calc(void)
{
    CPU_Step[CPU_DENSMISC] += measure_time(); double t00_truestart = my_second();
    PRINT_STATUS(" ..calculating first-order corrections to cell sizes/faces");
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    cellcorrections_final_operations_and_cleanup(); /* do final operations on results */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_DENSCOMPUTE] += timecomp; CPU_Step[CPU_DENSWAIT] += timewait; CPU_Step[CPU_DENSCOMM] += timecomm;
    CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm); /* collect timings and reset clock for next timing */
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */

#endif // parent if statement for all code in the HYDRO_VOLUME_CORRECTIONS block
