#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"

/*! Routines for the drift/predict step */

/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially in detail (although the highest-level algorithm
 * structure remains essentially the same)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO, and many new
 * options and subroutines added for flexibility with different
 * hydro solvers, timestepping schemes, boundary conditions, and
 * different mesh-motion options added.
 */

void reconstruct_timebins(void)
{
    int i, bin;
    long long glob_sum;
    
    for(bin = 0; bin < TIMEBINS; bin++)
    {
        TimeBinCount[bin] = 0;
        TimeBinCountGas[bin] = 0;
        FirstInTimeBin[bin] = -1;
        LastInTimeBin[bin] = -1;
#ifdef GALSF
        TimeBinSfr[bin] = 0;
#endif
#ifdef SINK_PARTICLES
        TimeBin_Sink_mass[bin] = 0;
        TimeBin_Sink_dynamicalmass[bin] = 0;
        TimeBin_Sink_Mdot[bin] = 0;
        TimeBin_Sink_Medd[bin] = 0;
#endif
    }
    
    for(i = 0; i < NumPart; i++)
    {
        bin = P.TimeBin[i];
        
        if(TimeBinCount[bin] > 0)
        {
            PrevInTimeBin[i] = LastInTimeBin[bin];
            NextInTimeBin[i] = -1;
            NextInTimeBin[LastInTimeBin[bin]] = i;
            LastInTimeBin[bin] = i;
        }
        else
        {
            FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
            PrevInTimeBin[i] = NextInTimeBin[i] = -1;
        }
        TimeBinCount[bin]++;
        if(P.Type[i] == 0)
            TimeBinCountGas[bin]++;
        
#ifdef GALSF
        if(P.Type[i] == 0)
            TimeBinSfr[bin] += CellP.Sfr[i];
#endif
#ifdef SINK_PARTICLES
        if(P.Type[i] == 5)
        {
            TimeBin_Sink_mass[bin] += P.Sink_Mass[i];
            TimeBin_Sink_dynamicalmass[bin] += P.Mass[i];
            TimeBin_Sink_Mdot[bin] += P.Sink_Mdot[i];
            TimeBin_Sink_Medd[bin] += P.Sink_Mdot[i] / P.Sink_Mass[i];
        }
#endif
    }
    
    make_list_of_active_particles();
    
    NumForceUpdate = 0;
    for (int i : ActiveParticleList)
    {
        NumForceUpdate++;
        if(i >= NumPart)
        {
            printf("Bummer i=%d\n", i);
            terminate("inconsistent list");
        }
    }
    
    sumup_large_ints(1, &NumForceUpdate, &glob_sum);
    GlobNumForceUpdate = glob_sum;
}





void drift_particle(int i, integertime time1)
{
    int j __attribute__((unused)); double dt_drift; integertime time0 = P.Ti_current[i];
    if(time1 < time0)
    {
        printf("i=%d time0=%lld time1=%lld\n", i, (long long)time0, (long long)time1);
        terminate("no prediction into past allowed");
    }
    if(time1 == time0) {return;}
    
    dt_drift = get_drift_factor(time0, time1, i, 0);
        
#if !defined(FREEZE_HYDRO)
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
    if(P.Type[i]==0) {advect_mesh_point(i,dt_drift);} else {P.Pos[i] += P.Vel[i] * dt_drift;}
#elif (SINGLE_STAR_TIMESTEPPING > 0)
    Vec3<double> fewbody_drift_dx, fewbody_kick_dv; // if super-timestepping, the updates above account for COM motion of the binary; now we account for the internal motion
    if( (P.Type[i] == 5) && (P.SuperTimestepFlag[i]>=2) )
    {
        Vec3<double> COM_Vel = P.Vel[i] + P.comp_dv[i] * (P.comp_Mass[i]/(P.Mass[i]+P.comp_Mass[i])); //center of mass velocity
        P.Pos[i] += COM_Vel * dt_drift; //center of mass drift
        odeint_super_timestep(i, dt_drift, fewbody_kick_dv, fewbody_drift_dx, 1); // do_fewbody_drift
        P.GravAccel[i] = P.COM_GravAccel[i]; //Overwrite the acceleration with center of mass value
        P.Pos[i] += fewbody_drift_dx; //Keplerian evolution
        P.Vel[i] += fewbody_kick_dv; //move on binary.orbit
    } else {
       P.Pos[i] += P.Vel[i] * dt_drift;
    }
#else
    P.Pos[i] += P.Vel[i] * dt_drift;
#endif
#endif // FREEZE_HYDRO clause
#if (NUMDIMS==1)
    P.Pos[i][1]=P.Pos[i][2]=0; // force zero-ing
#endif
#if (NUMDIMS==2)
    P.Pos[i][2]=0; // force zero-ing
#endif

#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
    double a_fac = return_timestep_dilation_factor(i,0);
    if(a_fac > 1.) {
        double cfac = dt_drift * (1. - 1./a_fac);
        P.Pos[i] += P.vel_of_nearest_special[i] * cfac; /* add back in the mean drift of the surrounding stuff to dilate only the local dynamics */
    }
#endif

    double divv_fac = P.Particle_DivVel[i] * dt_drift;
    double divv_fac_max = 0.3; //1.5; // don't allow KernelRadius to change too much in predict-step //
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(ags_density_isactive(i) && P.Type[i]>0) {divv_fac_max=4;} // can [should] allow larger changes when using adapting soft for all
#endif
    if(divv_fac > +divv_fac_max) divv_fac = +divv_fac_max;
    if(divv_fac < -divv_fac_max) divv_fac = -divv_fac_max;
    
#ifdef GRAIN_FLUID
    if((1 << P.Type[i]) & (GRAIN_PTYPES))
    {
        P.KernelRadius[i] *= exp((double)divv_fac / ((double)NUMDIMS));
        if(P.KernelRadius[i] < All.MinKernelRadius) {P.KernelRadius[i] = All.MinKernelRadius;}
        if(P.KernelRadius[i] > All.MaxKernelRadius) {P.KernelRadius[i] = All.MaxKernelRadius;}
    }
#endif

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(ags_density_isactive(i) && (dt_drift>0)) /* particle is AGS-active */
    {
        double minsoft = ags_return_minsoft(i), maxsoft = ags_return_maxsoft(i);
        P.AGS_KernelRadius[i] *= exp((double)divv_fac / ((double)NUMDIMS));
        if(P.AGS_KernelRadius[i] < minsoft) {P.AGS_KernelRadius[i] = minsoft;}
        if(P.AGS_KernelRadius[i] > maxsoft) {P.AGS_KernelRadius[i] = maxsoft;}
    } else {P.AGS_KernelRadius[i] = ForceSoftening_KernelRadius(i);} /* non-AGS-active particles use fixed softening */
#endif
    
#ifdef DM_FUZZY
    do_dm_fuzzy_drift_kick(i, dt_drift, 1);
#endif
    
    if((P.Type[i] == 0) && (P.Mass[i] > 0))
        {
            double dt_gravkick, dt_gravkick_pm, dt_hydrokick, dt_entr;
            dt_entr = dt_hydrokick = (time1 - time0) * UNIT_INTEGERTIME_IN_PHYSICAL(i);
            dt_gravkick = get_gravkick_factor(time0, time1, i, 0);
            
#ifdef PMGRID
            dt_gravkick_pm = get_gravkick_factor(time0, time1, -1, 0);
            CellP.VelPred[i] += P.GravAccel[i]*dt_gravkick + P.GravPM[i]*dt_gravkick_pm + CellP.HydroAccel[i]*(dt_hydrokick*All.cf_atime); /* make sure v is in code units */
#else
            CellP.VelPred[i] += P.GravAccel[i] * dt_gravkick + CellP.HydroAccel[i] * (dt_hydrokick*All.cf_atime); /* make sure v is in code units */
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
	        if((P.Type[i] == 5) && (P.SuperTimestepFlag[i]>=2)) {CellP.VelPred[i] += Vec3<double>{fewbody_kick_dv[0], fewbody_kick_dv[1], fewbody_kick_dv[2]};}
#endif	    
            
#if defined(TURB_DRIVING)
            CellP.VelPred[i] += CellP.TurbAccel[i] * dt_gravkick;
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            CellP.VelPred[i] += CellP.Rad_Accel[i] * (All.cf_atime * dt_hydrokick);
#endif
            
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            P.Mass[i] = DMAX(P.Mass[i] + CellP.DtMass[i] * dt_entr, 0.5 * CellP.MassTrue[i]);
#endif
            
            CellP.Density[i] *= exp(-divv_fac);
            double etmp = CellP.InternalEnergyPred[i] + CellP.DtInternalEnergy[i] * dt_entr;
#if defined(RADTRANSFER) && defined(RT_EVOLVE_ENERGY) /* block here to deal with tricky cases where radiation energy density is -much- larger than thermal */ 
            int kfreq; double erad_tot=0,tot_e_min=0,enew=0,int_e_min=0,dErad=0,rsol_fac=C_LIGHT_CODE_REDUCED(i)/C_LIGHT_CODE; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad_tot+=CellP.Rad_E_gamma_Pred[i][kfreq];}
            if(erad_tot > 0)
            {
                int_e_min=0.025*CellP.InternalEnergyPred[i]; tot_e_min=0.025*(erad_tot/rsol_fac+CellP.InternalEnergyPred[i]*P.Mass[i]);
                enew=DMAX(erad_tot/rsol_fac+etmp*P.Mass[i],tot_e_min); etmp=(enew-erad_tot/rsol_fac)/P.Mass[i]; if(etmp<int_e_min) {dErad=rsol_fac*(etmp-int_e_min); etmp=int_e_min;}
                if(dErad<-0.975*erad_tot) {dErad=-0.975*erad_tot;} CellP.InternalEnergyPred[i] = etmp; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {CellP.Rad_E_gamma_Pred[i][kfreq] *= 1 + dErad/erad_tot;}
            } else {
                if(etmp<0.5*CellP.InternalEnergyPred[i]) {CellP.InternalEnergyPred[i] *= 0.5;} else {CellP.InternalEnergyPred[i]=etmp;}
            }
#else
            if(etmp<0.5*CellP.InternalEnergyPred[i]) {CellP.InternalEnergyPred[i] *= 0.5;} else {CellP.InternalEnergyPred[i]=etmp;}
#endif
            if(CellP.InternalEnergyPred[i]<All.MinEgySpec) CellP.InternalEnergyPred[i]=All.MinEgySpec;
            
#ifdef HYDRO_PRESSURE_SPH
            CellP.EgyWtDensity[i] *= exp(-divv_fac);
#endif
            
#if (HYDRO_FIX_MESH_MOTION > 0)
            P.KernelRadius[i] *= exp((double)divv_fac / ((double)NUMDIMS));
            if(P.KernelRadius[i] < All.MinKernelRadius) {P.KernelRadius[i] = All.MinKernelRadius;}
            if(P.KernelRadius[i] > All.MaxKernelRadius) {P.KernelRadius[i] = All.MaxKernelRadius;}
#ifdef ADAPTIVE_GRAVSOFT_FORALL
            if(1 & ADAPTIVE_GRAVSOFT_FORALL) {P.AGS_KernelRadius[i] = P.KernelRadius[i];} /* gas is AGS-active, so needs to be set here to match updated KernelRadius */
#endif
#endif
            drift_extra_physics(i, time0, time1, dt_entr);

            set_eos_pressure(i);
        }
    
    /* check for reflecting or outflow or otherwise special boundaries: if so, do the reflection/boundary! */
    apply_special_boundary_conditions(i,P.Mass[i],0);

    P.Ti_current[i] = time1;
}





void move_particles(integertime time1)
{
    int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(i=0; i<NumPart; i++) {drift_particle(i, time1);}
}





void drift_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
#ifdef MAGNETIC
    double BphysVolphys_to_BcodeVolCode = 1 / All.cf_atime;
    CellP.BPred[i] += CellP.DtB[i] * (dt_entr * BphysVolphys_to_BcodeVolCode); // fluxes are always physical, convert to code units //
#ifdef DIVBCLEANING_DEDNER
    double PhiphysVolphys_to_PhicodeVolCode = 1 / All.cf_a3inv; // for mass-based phi fluxes (otherwise coefficient is 1)
    double dtphi_code = (PhiphysVolphys_to_PhicodeVolCode) * CellP.DtPhi[i];
    CellP.PhiPred[i] += dtphi_code  * dt_entr;
    double t_damp = Get_Gas_PhiField_DampingTimeInv(i);
    if((t_damp>0) && (!isnan(t_damp)))
    {
        CellP.PhiPred[i] *= exp( -dt_entr * t_damp );
    }
#endif
#ifdef MHD_ALTERNATIVE_LEAPFROG_SCHEME
    CellP.B[i] = CellP.BPred[i];
#ifdef DIVBCLEANING_DEDNER
    CellP.Phi[i]=CellP.PhiPred[i];
#endif
#endif
#endif
#ifdef COSMIC_RAY_FLUID
    CosmicRay_Update_DriftKick(i,dt_entr,1);
#endif
#ifdef RADTRANSFER
    rt_update_driftkick(i,dt_entr,1);
#endif
#ifdef EOS_ELASTIC
    elastic_body_update_driftkick(i,dt_entr,1);
#endif
}





/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef BOX_PERIODIC
void do_box_wrapping(void)
{
    int i, j;
    double boxsize[3];
    boxsize[0] = boxSize_X;
    boxsize[1] = boxSize_Y;
    boxsize[2] = boxSize_Z;
    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(j)
#endif
    for(i = 0; i < NumPart; i++)
    {
        for(j = 0; j < 3; j++)
        {
            while(P.Pos[i][j] < 0)
            {
                P.Pos[i][j] += boxsize[j];
#ifdef BOX_SHEARING
                if(j==0)
                {
                    P.Vel[i][BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
                    P.dp[i][BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset * P.Mass[i];
                    if(P.Type[i]==0)
                    {
                        CellP.VelPred[i][BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) // if have moving cells need to wrap them, too (if cells aren't moving, should never reach this wrap) //
                        CellP.ParticleVel[i][BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Vel_Offset;
#endif
                    }
#if (BOX_SHEARING > 1)
                    /* if we're not assuming axisymmetry, we need to shift the coordinates for the shear flow at the boundary */
                    P.Pos[i][BOX_SHEARING_PHI_COORDINATE] -= Shearing_Box_Pos_Offset;
#endif
                }
#endif
            }
            
            while(P.Pos[i][j] >= boxsize[j])
            {
                P.Pos[i][j] -= boxsize[j];
#ifdef BOX_SHEARING
                if(j==0)
                {
                    P.Vel[i][BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
                    P.dp[i][BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset * P.Mass[i];
                    if(P.Type[i]==0)
                    {
                        CellP.VelPred[i][BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) // if have moving cells need to wrap them, too (if cells aren't moving, should never reach this wrap) //
                        CellP.ParticleVel[i][BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Vel_Offset;
#endif
                    }
#if (BOX_SHEARING > 1)
                    /* if we're not assuming axisymmetry, we need to shift the coordinates for the shear flow at the boundary */
                    P.Pos[i][BOX_SHEARING_PHI_COORDINATE] += Shearing_Box_Pos_Offset;
#endif
                }
#endif
            }
        }
    }
}
#endif




/* ====================================================================== */
/* ================== Functions for physical information ================ */
/* ====================================================================== */


/* this function returns the effective (grid-equivalent) particle 'size'; useful for things like 
    time-stepping and limiter functions */
double INLINE_FUNC Get_Particle_Size(int i)
{
    /* in previous versions of the code, we took NumNgb^(1/NDIMS) here; however, now we 
        take that when NumNgb is computed (at the end of the density routine), so we 
        don't have to re-compute it each time. That makes this function fast enough to 
        call -inside- of loops (e.g. hydro computations) */
#if (NUMDIMS == 1)
    return 2.00000 * P.KernelRadius[i] / P.NumNgb[i]; // (2)^(1/1)
#endif
#if (NUMDIMS == 2)
    return 1.77245 * P.KernelRadius[i] / P.NumNgb[i]; // (pi)^(1/2)
#endif
#if (NUMDIMS == 3)
    return 1.61199 * P.KernelRadius[i] / P.NumNgb[i]; // (4pi/3)^(1/3)
#endif
}



double INLINE_FUNC Get_Particle_Expected_Area(double h)
{
#if (NUMDIMS == 1)
    return 2;
#endif
#if (NUMDIMS == 2)
    return 2 * M_PI * h;
#endif
#if (NUMDIMS == 3)
    return 4 * M_PI * h * h;
#endif
}


/* return the estimated local column (physical units) from a local Sobolev approximation, or using the 'treecol' approximation from the gravity tree if the relevant config flag options are enabled */
double evaluate_NH_from_GradRho(MyFloat gradrho[3], double rkern, double rho, double numngb_ndim, double include_h, int target)
{
    double gradrho_mag=0;
    if(rho>0)
    {
#ifdef RT_USE_TREECOL_FOR_NH        
        gradrho_mag = include_h * rho * rkern / numngb_ndim; if(target>=0) {gradrho_mag += P.SigmaEff[target];}
#else             
        gradrho_mag = sqrt(gradrho[0]*gradrho[0]+gradrho[1]*gradrho[1]+gradrho[2]*gradrho[2]);
        if(gradrho_mag > 0) {gradrho_mag = rho*rho/gradrho_mag;} else {gradrho_mag=0;}
        if(include_h > 0) if(numngb_ndim > 0) gradrho_mag += include_h * rho * rkern / numngb_ndim; // quick-and-dirty approximation to the effective neighbor number needed here
#endif        
    }
    return gradrho_mag * All.cf_a2inv; // (physical units) // *(Z/Zsolar) add metallicity dependence
}







#ifdef MAGNETIC
/* this function is needed to control volume fluxes of the normal components of B and phi in the 
    -bad- situation where the meshless method 'faces' do not properly close (usually means you are 
    using boundary conditions that you should not) */
double Get_DtB_FaceArea_Limiter(int i)
{
#ifdef HYDRO_SPH
    return 1;
#else
    /* define some variables */
    double dt_entr = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
    /* check the magnitude of the predicted change in B-fields, vs. B-magnitude */
    Vec3<double> dB = CellP.DtB[i] * (dt_entr / All.cf_atime); /* converts to code units of Vol_code*B_code = Vol_phys*B_phys/a */
    double dBmag = dB.norm(), Bmag = CellP.BPred[i].norm();
    /* also make sure to check the actual pressure, since if P>>B, we will need to allow larger changes in B per timestep */
    double P_BV_units = sqrt(2.*CellP.Pressure[i]*All.cf_a3inv)*P.Mass[i]/CellP.Density[i] / All.cf_a2inv;
    /* the above should be in CODE Bcode*Vol_code units! */
    double Bmag_max = DMAX(Bmag, DMIN( P_BV_units, 10.*Bmag ));
    /* now check how accurately the cell is 'closed': the face areas are ideally zero */
    double area_sum = fabs(CellP.Face_Area[i][0])+fabs(CellP.Face_Area[i][1])+fabs(CellP.Face_Area[i][2]);
    /* but this needs to be normalized to the 'expected' area given KernelRadius */
    double area_norm = Get_Particle_Expected_Area(P.KernelRadius[i] * All.cf_atime);
    /* ok, with that in hand, define an error tolerance based on this */
    if(area_norm>0)
    {
        double area_norm_min_threshold = 0.001;
        double area_norm_weight = 200.0;
        if(area_sum/area_norm > area_norm_min_threshold)
        {
            double tol = (All.CourantFac/0.2) * DMAX( 0.01, area_norm/(area_norm_weight * area_sum) );
            tol *= Bmag_max; /* give the limiter dimensions */
            if(dBmag > tol) {return tol/dBmag;} /* now actually check if we exceed this */
        }
    }
    return 1;
#endif
}


#ifdef DIVBCLEANING_DEDNER
double INLINE_FUNC Get_Gas_PhiField(int i_particle_id)
{
    //return CellP.PhiPred[i_particle_id] * CellP.Density[i_particle_id] / P.Mass[i_particle_id]; // volumetric phy-flux (requires extra term compared to mass-based flux)
    return CellP.PhiPred[i_particle_id] / P.Mass[i_particle_id]; // mass-based phi-flux
}

double INLINE_FUNC Get_Gas_PhiField_DampingTimeInv(int i_particle_id)
{
    /* this timescale should always be returned as a -physical- time */
#ifdef HYDRO_SPH
    /* PFH: add simple damping (-phi/tau) term */
    double damping_tinv = 0.5 * All.DivBcleanParabolicSigma * (CellP.MaxSignalVel[i_particle_id] / (All.cf_atime*Get_Particle_Size(i_particle_id)));
#else
    double damping_tinv;
#ifdef SELFGRAVITY_OFF
    damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveSpeed / (All.cf_atime*Get_Particle_Size(i_particle_id)); // fastest wavespeed has units of [vphys]
    //double damping_tinv = All.DivBcleanParabolicSigma * All.FastestWaveDecay * All.cf_a2inv; // no improvement over fastestwavespeed; decay has units [vphys/rphys]
#else
    // only see a small performance drop from fastestwavespeed above to maxsignalvel below, despite the fact that below is purely local (so allows more flexible adapting to high dynamic range)
    damping_tinv = 0.0;
    
    if(P.KernelRadius[i_particle_id] > 0)
    {
        double h_eff = Get_Particle_Size(i_particle_id);
        double vsig2 = 0.5 * fabs(CellP.MaxSignalVel[i_particle_id]);
        double phi_B_eff = 0.0;
        if(vsig2 > 0) {phi_B_eff = Get_Gas_PhiField(i_particle_id) / (All.cf_atime * vsig2);}
        double vsig1 = 0.0;
        if(CellP.Density[i_particle_id] > 0)
        {
            vsig1 = sqrt( Get_Gas_effective_soundspeed_i(i_particle_id)*Get_Gas_effective_soundspeed_i(i_particle_id) +
                 (1. / All.cf_atime) *
                 (Get_Gas_BField(i_particle_id,0)*Get_Gas_BField(i_particle_id,0) +
                  Get_Gas_BField(i_particle_id,1)*Get_Gas_BField(i_particle_id,1) +
                  Get_Gas_BField(i_particle_id,2)*Get_Gas_BField(i_particle_id,2) +
                  phi_B_eff*phi_B_eff) / CellP.Density[i_particle_id] );
        }
        vsig1 = DMAX(vsig1, vsig2);
        vsig2 = 0.0;
        vsig2 = CellP.Gradients.Velocity[i_particle_id].frobenius_norm();
        vsig2 = 3.0 * h_eff * DMAX( vsig2, fabs(P.Particle_DivVel[i_particle_id])) / All.cf_atime;
        double prefac_fastest = 0.1;
        double prefac_tinv = 0.5;
        double area_0 = 0.1;
#ifdef MHD_CONSTRAINED_GRADIENT
        prefac_fastest = 1.0;
        prefac_tinv = 2.0;
        area_0 = 0.05;
        vsig2 *= 5.0;
        if(CellP.FlagForConstrainedGradients[i_particle_id] <= 0) prefac_tinv *= 30;
#endif
        prefac_tinv *= sqrt(1. + CellP.ConditionNumber[i_particle_id]/100.);
        double area = fabs(CellP.Face_Area[i_particle_id][0]) + fabs(CellP.Face_Area[i_particle_id][1]) + fabs(CellP.Face_Area[i_particle_id][2]);
        area /= Get_Particle_Expected_Area(P.KernelRadius[i_particle_id]);
        prefac_tinv *= (1. + area/area_0)*(1. + area/area_0);
        
        double vsig_max = DMAX( DMAX(vsig1,vsig2) , prefac_fastest * All.FastestWaveSpeed );
        damping_tinv = prefac_tinv * All.DivBcleanParabolicSigma * (vsig_max / (All.cf_atime * h_eff));
    }
#endif
#endif
    return damping_tinv;
}

#endif // dedner
#endif // magnetic





/* -------------------------------------------------------------------------------------------------------------------------------------
 ------------------- the following routines are not setting the velocity, but instead are useful routines for computation of
 -------------------  various quantities needed in the mesh motion for different coordinate systems or assumed mesh shapes
 ------------------------------------------------------------------------------------------------------------------------------------- */

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
/* time-step the positions of the mesh points. this is trivial except if we are evolving the mesh points in non-cartesian coordinates
    (cylindrical or spherical) based on assumed fixed initial velocities (if HYDRO_FIX_MESH_MOTION=2 or 3),
    in which case we have to convert back and forth. */
void advect_mesh_point(int i, double dt)
{
#if (HYDRO_FIX_MESH_MOTION == 2) || (HYDRO_FIX_MESH_MOTION == 3) // cylindrical or spherical coordinates
    // define the location relative to the origin (needed in these coordinate systems)
    Vec3<double> dp = P.Pos[i]; Vec3<double> dp_offset = {}; // assume center is at coordinate origin
#if defined(GRAVITY_ANALYTIC_ANCHOR_TO_PARTICLE) // unless we use a special anchor, to define the center
    dp_offset = P.Pos[i] - P.Min_xyz_to_Sink[i];
#elif defined(BOX_PERIODIC) // or if periodic, the box mid-point is instead the center
#if (NUMDIMS==1)
    dp_offset[0] = -boxHalf_X;
#elif (NUMDIMS==2)
    dp_offset[0] = -boxHalf_X; dp_offset[1] = -boxHalf_Y;
#else
    dp_offset = Vec3<double>{-boxHalf_X, -boxHalf_Y, -boxHalf_Z};
#endif
#endif
    dp += dp_offset;
#if (HYDRO_FIX_MESH_MOTION == 2) // cylindrical
    double r2=dp[0]*dp[0]+dp[1]*dp[1], r=sqrt(r2), c0=dp[0]/r, s0=dp[1]/r, z=dp[2]; // get r, sin/cos theta, z
    double vr=c0*CellP.ParticleVel[i][0] + s0*CellP.ParticleVel[i][1], vt=s0*CellP.ParticleVel[i][0] - c0*CellP.ParticleVel[i][1], vz=CellP.ParticleVel[i][2]; // velocities in these directions
    double r_n=r+vr*dt, z_n=z+vz*dt, c_n=c0-s0*(vt/r)*dt, s_n=s0+c0*(vt/r)*dt; // updated cylindrical values
    dp[0] = c_n*r_n; dp[1] = s_n*r_n; dp[2] = z_n; // back to coordinates
    CellP.ParticleVel[i][0] = c_n*vr + s_n*vt; // re-set velocities in these coordinates //
    CellP.ParticleVel[i][1] = s_n*vr - c_n*vt;
    CellP.ParticleVel[i][2] = vz;
    return;
#elif (HYDRO_FIX_MESH_MOTION == 3) // spherical
    Vec3<double> v = CellP.ParticleVel[i]; double r2=dp.norm_sq(); // assume center is at coordinate origin
    double r=sqrt(r2), rxy=sqrt(dp[0]*dp[0]+dp[1]*dp[1]), vr=dot(dp,v)/r; // updated r is easy
    double ct = 1./sqrt(1.+dp[1]*dp[1]/(dp[0]*dp[0])), st = (dp[1]/dp[0])*ct; // cos and sin theta
    double cp = sqrt(1.-dp[2]*dp[2]/(r*r)), sp = dp[2]/r; // cos and sin phi
    double t_dot = (v[0]*dp[1]-v[1]*dp[0])/(rxy*rxy), p_dot = (dp[2]*(dp[0]*v[0]+dp[1]*v[1])-rxy*rxy*v[2])/(r*r*rxy); // theta, phi derivatives
    double r_n=r+vr*dt, ct_n=ct-st*t_dot, st_n=st+ct*t_dot, cp_n=cp-sp*t_dot, sp_n=sp+cp*t_dot; // updated angles and positions in spherical
    dp[0] = r_n * ct_n * cp_n; dp[1] = r_n * st_n * cp_n; dp[2] = r_n * sp_n; // back to coordinates
    rxy = sqrt(dp[0]*dp[0] + dp[1]*dp[1]); // updated rxy
    CellP.ParticleVel[i][0] = (dp[0]/r_n) * vr + dp[1] * t_dot + dp[0]*dp[2]/rxy * p_dot; // back to cartesian velocities
    CellP.ParticleVel[i][1] = (dp[1]/r_n) * vr - dp[0] * t_dot + dp[1]*dp[2]/rxy * p_dot; // back to cartesian velocities
    CellP.ParticleVel[i][2] = (dp[2]/r_n) * vr - rxy * p_dot; // back to cartesian velocities
    return;
#endif
    // ok now have the updated x/y/z positions relative to the origin, convert these back to the simulation coordinate frame
    P.Pos[i] = dp - dp_offset;
#endif // ok done with cylindrical/spherical coordinates
    
    
    // ok anything else ('normal' coordinates), does down here
    P.Pos[i] += CellP.ParticleVel[i] * dt; // for standard grid velocities, this is trivial //
    return;
}




/* routine to calculate the overlapping face area of two cuboids in NDIMS dimensions based on their relative positions */
double calculate_face_area_for_cartesian_mesh(double *dp, double rinv, double l_side, double *Face_Area_Vec)
{
    Face_Area_Vec[0]=Face_Area_Vec[1]=Face_Area_Vec[2]=0; double Face_Area_Norm;
#if (NUMDIMS==1)
    Face_Area_Norm = 1; Face_Area_Vec[0] = Face_Area_Norm * dp[0]/fabs(dp[0]);
#elif (NUMDIMS==2)
    if(fabs(dp[0]) > fabs(dp[1])) {Face_Area_Vec[0] = Face_Area_Norm = DMAX(0,l_side-fabs(dp[1])) * dp[0]/fabs(dp[0]) * All.cf_atime;} else {Face_Area_Vec[1] = Face_Area_Norm = DMAX(0,l_side-fabs(dp[0])) * dp[1]/fabs(dp[1]) * All.cf_atime;}
#else
    double dp_abs[3]; int k,kdir; for(k=0;k<3;k++) {dp_abs[k] = fabs(dp[k]);}
    if((dp_abs[0]>=dp_abs[1])&&(dp_abs[0]>=dp_abs[2])) {kdir=0;} else if ((dp_abs[1]>=dp_abs[0])&&(dp_abs[1]>=dp_abs[2])) {kdir=1;} else {kdir=2;}
    Face_Area_Norm=1; for(k=0;k<3;k++) {if(k!=kdir) {Face_Area_Norm *= DMAX(0,l_side-dp_abs[k]) * All.cf_atime*All.cf_atime;}}
    Face_Area_Vec[kdir] = Face_Area_Norm * dp[kdir]/fabs(dp[kdir]);
#endif
    return fabs(Face_Area_Norm);
}

#endif




#ifdef SPECIAL_POINT_WEIGHTED_MOTION
double weight_function_for_weighted_motion_smoothing(double r, int mode)
{
    double wt = 1, amax = 1.e2, rmax_pc = 100., slope_index = 1, r_pc = r * All.cf_atime * UNIT_LENGTH_IN_PC;
    if(r_pc < rmax_pc) {wt = DMAX(pow(r_pc / rmax_pc , slope_index) , 1./amax);}
    if(mode) {return 1 - sqrt(wt);} else {return wt;}
}
#endif
