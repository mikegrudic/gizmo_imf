#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"

/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO 
 * (added energy/entropy switch, terms for explicit mass conservation in mass fluxes, 
 *  and updates to additional fluid variables, options for different hydro solvers,
 *  libraries and architectures for different grid/cell structures, rewritten
 *  to allow boundary conditions other than periodic and open, etc.)
 */

void apply_long_range_kick(integertime tstart, integertime tend);

void do_first_halfstep_kick(void)
{
    int i; integertime ti_step, tstart=0, tend=0;
    
#ifdef TURB_DRIVING
    do_turb_driving_step_first_half();
#endif
    
#ifdef PMGRID
    if(All.PM_Ti_begstep == All.Ti_Current)	/* need to do long-range kick */
    {
        ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
        tstart = All.PM_Ti_begstep;
        tend = tstart + ti_step / 2;
        apply_long_range_kick(tstart, tend);
    }
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME    
    /* as currently written with some revisions to MFV methods, should only update on active timesteps */
    for(i = 0; i < NumPart; i++)
    {
        if((TimeBinActive[P.TimeBin[i]]) || (P.Type[i]==0)) /* active OR gas, need to check each timestep to ensure manifest conservation */
#else
    for (int i : ActiveParticleList) /* 'full' kick for active particles */
    {	    
#endif
        {
            if(P.Mass[i] > 0)
            {
                ti_step = GET_PARTICLE_INTEGERTIME(i);
                tstart = P.Ti_begstep[i];	/* beginning of step */
                tend = P.Ti_begstep[i] + ti_step / 2;	/* midpoint of step */
                do_the_kick(i, tstart, tend, P.Ti_current[i], 0);
            }
        }
    } // for(i = 0; i < NumPart; i++) //
}

void do_second_halfstep_kick(void)
{
    int i; integertime ti_step, tstart=0, tend=0;
    
#ifdef PMGRID
    if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
        ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
        tstart = All.PM_Ti_begstep + ti_step / 2;
        tend = tstart + ti_step / 2;
        apply_long_range_kick(tstart, tend);
    }
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    for(i = 0; i < NumPart; i++)
    {
        if((TimeBinActive[P.TimeBin[i]]) || (P.Type[i]==0)) /* active OR gas, need to check each timestep to ensure manifest conservation */
#else
    for (int i : ActiveParticleList) /* 'full' kick for active particles */
    {
#endif
        {
            if(P.Mass[i] > 0)
            {
                ti_step = GET_PARTICLE_INTEGERTIME(i);
                tstart = P.Ti_begstep[i] + ti_step / 2;	/* midpoint of step */
                tend = P.Ti_begstep[i] + ti_step;	/* end of step */
                do_the_kick(i, tstart, tend, P.Ti_current[i], 1);
                set_predicted_quantities_for_extra_physics(i);
            }
        }
    } // for(i = 0; i < NumPart; i++) //
    
#ifdef TURB_DRIVING
    do_turb_driving_step_second_half();
#endif
}

#ifdef HERMITE_INTEGRATION
int eligible_for_hermite(int i)
{
    if(!(HERMITE_INTEGRATION & (1<<P.Type[i]))) {return 0;} // hermite flag said to not include these types
#if defined(DM_FUZZY) || defined(CBE_INTEGRATOR)
    if(P.Type[i]==1) {return 0;} // not compatible with these flags for these types
#endif
#if defined(GRAIN_FLUID)
    if((1 << P.Type[i]) & (GRAIN_PTYPES)) {return 0;} // not compatible with these flags for these types
#endif
#if defined(SINK_PARTICLES) || defined(GALSF)
    if(P.StellarAge[i] >= DMAX(All.Time - 2*(GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i)*All.cf_hubble_a), 0)) {return 0;} // if we were literally born yesterday then let things settle down a bit with the less-accurate, but more-robust regular integration
    if(P.AccretedThisTimestep[i]) {return 0;}
#endif
#if (SINGLE_STAR_TIMESTEPPING > 0)
    if(P.SuperTimestepFlag[i] >= 2) {return 0;}
#endif   
    return 1;
}

// Initial "prediction" step of Hermite integration, performed after the initial force evaluation 
// Note: the below routines only account for gravitational acceleration - only appropriate for stars or collisionless particles
void do_hermite_prediction(void)
{
    int i; integertime ti_step, tstart=0, tend=0;
    for (int i : ActiveParticleList) {
	if(eligible_for_hermite(i)) { /* check if we're actually eligible */
	    if(P.Mass[i] > 0) { /* skip massless particles scheduled for deletion */
		ti_step = GET_PARTICLE_INTEGERTIME(i);
		tstart = P.Ti_begstep[i];    /* beginning of step */
		tend = P.Ti_begstep[i] + ti_step;    /* end of step */
            double dt_grav = get_gravkick_factor(tstart, tend, i, 0);
#ifdef PMGRID
            //Add the long-range kick from the first half-step, if necessary (since we are overwriting the previous kick operations with the Hermite scheme)
            if(All.PM_Ti_begstep == All.Ti_Current)	/* need to do long-range kick */
            {
                double dt_grav_pm = get_gravkick_factor(All.PM_Ti_begstep, All.PM_Ti_begstep + (All.PM_Ti_endstep - All.PM_Ti_begstep)/2, i, 0);
                P.OldVel[i] += P.GravPM[i] * dt_grav_pm;
            }
#endif
            P.Pos[i] = P.OldPos[i] + (P.OldVel[i] + (P.Hermite_OldAcc[i] + P.OldJerk[i] * (dt_grav/3)) * (dt_grav/2)) * dt_grav;
            P.Vel[i] = P.OldVel[i] + (P.Hermite_OldAcc[i] + P.OldJerk[i] * (dt_grav/2)) * dt_grav;
		}}} // for (int i : ActiveParticleList) 
}

void do_hermite_correction(void) // corrector step
{
    int i; integertime ti_step, tstart=0, tend=0;
    for (int i : ActiveParticleList) {
	if(eligible_for_hermite(i)){
                if(P.Mass[i] > 0) {
                    ti_step = GET_PARTICLE_INTEGERTIME(i);
                    tstart = P.Ti_begstep[i];    /* beginning of step */
                    tend = P.Ti_begstep[i] + ti_step;    /* end of step */
                    double dt_grav = get_gravkick_factor(tstart, tend, i, 0);
                    P.Vel[i] = P.OldVel[i] + (P.Hermite_OldAcc[i] + P.GravAccel[i]) * (dt_grav * 0.5) + (P.OldJerk[i] - P.GravJerk[i]) * (dt_grav * dt_grav / 12);
                    P.Pos[i] = P.OldPos[i] + (P.Vel[i] + P.OldVel[i]) * (dt_grav * 0.5) + (P.Hermite_OldAcc[i] - P.GravAccel[i]) * (dt_grav * dt_grav / 12);
#ifdef PMGRID
                    //Add the long-range kick from the second half-step, if necessary (since we are overwriting the previous kick operations with the Hermite scheme)
                    if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
                    {
                        double dt_grav_pm = get_gravkick_factor(All.PM_Ti_begstep + (All.PM_Ti_endstep - All.PM_Ti_begstep)/2, All.PM_Ti_endstep, i, 0);
                        P.OldVel[i] += P.GravPM[i] * dt_grav_pm;
                    }
#endif
		}}} //     for (int i : ActiveParticleList)
}
#endif // HERMITE_INTEGRATION


#ifdef PMGRID
void apply_long_range_kick(integertime tstart, integertime tend)
{
    int i;
    double dt_gravkick = get_gravkick_factor(tstart, tend, -1, 0);
    for(i = 0; i < NumPart; i++)
    {
        if(P.Mass[i] > 0)
        {
            Vec3<double> dvel = P.GravPM[i] * dt_gravkick; /* do the kick, only collisionless particles */
            P.Vel[i] += dvel;
            P.dp[i] += dvel * P.Mass[i];
        }
    }
}
#endif


void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent, int mode)
{
    Vec3<double> dp; double dt_entr, dt_gravkick, dt_hydrokick;
    double mass_old, mass_pred, mass_new;
    mass_old = mass_pred = mass_new = P.Mass[i];    
    
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    /* need to do the slightly more complicated update scheme to maintain exact mass conservation */
    if(P.Type[i]==0)
    {
        if(CellP.dMass[i] != 0) //ent_old = CellP.InternalEnergy[i]; for(j=0;j<3;j++) v_old[j] = P.Vel[i][j];
        {
            double dMass=0; // fraction of delta_conserved to couple per kick step (each 'kick' is 1/2-timestep) // double dv[3], v_old[3], dMass, ent_old=0, d_inc = 0.5;
            if(mode != 0) // update the --conserved-- variables of each particle //
            {
                dMass = ((tend - tstart) * UNIT_INTEGERTIME_IN_PHYSICAL(i)) * CellP.DtMass[i]; if(dMass * CellP.dMass[i] < 0) {dMass = 0;} // slope-limit: no opposite reconstruction! //
                if((fabs(dMass) > fabs(CellP.dMass[i]))) {dMass = CellP.dMass[i];} // try to get close to what the time-integration scheme would give //
                CellP.dMass[i] -= dMass;
            } else {dMass = CellP.dMass[i];}
            if(dMass < -0.99*CellP.MassTrue[i]) {dMass = -0.99*CellP.MassTrue[i];} // limiter to prevent madness //

            /* load and update the particle masses : particle mass update here, from hydro fluxes */
            mass_old = CellP.MassTrue[i]; mass_pred = P.Mass[i]; mass_new = mass_old + dMass; CellP.MassTrue[i] = mass_new; // UNITS: remember all time derivatives (DtX, dX) are in -physical- units; as are mass, entropy/internal energy, but -not- velocity //
            /* double e_old = mass_old * CellP.InternalEnergy[i]; for(j = 0; j< 3; j++) e_old += 0.5*mass_old * (P.Vel[i][j]/All.cf_atime)*(P.Vel[i][j]/All.cf_atime); // physical //
            for(j = 0; j < 3; j++) // momentum-space-kick
            {
                dp[j] = d_inc * CellP.dMomentum[i][j]; // now update the velocity based on the total momentum change
                P.Vel[i][j] = (mass_old*P.Vel[i][j] + dp[j]*All.cf_atime) / mass_new; // call after tabulating dP[j] //
            } // kick for gas internal energy/entropy
            e_old += d_inc * CellP.dInternalEnergy[i]; // for(j = 0; j< 3; j++) e_old -= 0.5*mass_new * (P.Vel[i][j]/All.cf_atime)*(P.Vel[i][j]/All.cf_atime); // increment of total (thermal+kinetic) energy; subtract off the new kinetic energy //
            CellP.InternalEnergy[i] = e_old / mass_new; check_particle_for_temperature_minimum(i); // obtain the new internal energy per unit mass, check floor // */
             
            // at the end of this kick, need to re-zero the dInternalEnergy, and other conserved-variable gas/fluid quantities set in the hydro loop, to avoid double-counting them
            if(mode==0) {CellP.dMass[i]=0;} /* CellP.dInternalEnergy[i]=0; CellP.dMomentum[i][0]=CellP.dMomentum[i][1]=CellP.dMomentum[i][2]=0; */
        }
    } // if(P.Type[i]==0) //
#endif
    
    /* only enter the 'normal' kick loop below for genuinely active particles */
    if(TimeBinActive[P.TimeBin[i]])
    {
        /* get the timestep (physical units for dt_entr and dt_hydrokick) */
        dt_entr = dt_hydrokick = (tend - tstart) * UNIT_INTEGERTIME_IN_PHYSICAL(i);
        dt_gravkick = get_gravkick_factor(tstart, tend, i, 0);
        
        if(P.Type[i]==0)
        {
            Vec3<double> grav_acc; double dEnt_Gravity = 0;
            grav_acc = P.GravAccel[i] * All.cf_a2inv;
#ifdef PMGRID
            grav_acc += P.GravPM[i] * All.cf_a2inv;
#endif

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
            dEnt_Gravity += -dot(CellP.GravWorkTerm[i], grav_acc) * All.cf_atime * dt_hydrokick;
#endif
            double du_tot = CellP.DtInternalEnergy[i] * dt_hydrokick + dEnt_Gravity;
#if defined(COOLING) && !defined(COOLING_OPERATOR_SPLIT)
            if((mode == 1) && (du_tot != 0) && (dt_hydrokick > 0)) { /* if about to consider second-halfstep kick (just after hydro), decide if we need to split this particular cell on this particular timestep, since this un-split solver can lead to energy conservation problems if the mechanical heating is much larger than cooling */
                CellP.CoolingIsOperatorSplitThisTimestep[i]=1; /* default to assume split */
                double DtInternalEnergyEffCGS = (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (PROTONMASS_CGS/HYDROGEN_MASSFRAC) * (du_tot/dt_hydrokick), DtInternalEnergyReference = 1.e-23*CellP.Density[i]*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS; /* define the effective work term in cgs and a reference typical cooling time */
                if(DtInternalEnergyEffCGS < DtInternalEnergyReference) {CellP.CoolingIsOperatorSplitThisTimestep[i]=0;} /* cooling is fast compared to the hydro work term, or the hydro term is negative [cooling], so un-split the operation */
            }
            if(CellP.CoolingIsOperatorSplitThisTimestep[i]==0) {du_tot=0;} /* cooling in unsplit, so zero contribution here */
#endif
            double dEnt = CellP.InternalEnergy[i] + du_tot;
            
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            /* if we're using a Riemann solver, we include an energy/entropy-type switch to ensure
                that we don't corrupt the temperature evolution of extremely cold, adiabatic flows */
            /* MHD tests suggest that this switch often does more harm than good: we will
             pay the price of noisier temperature fields (corrupting them when c_s << v_A << v_bulk)
             and they are dynamically irrelevant, in exchange for avoiding potentially much more
             serious errors if this tripped when the B-fields were important */
            double e_thermal,e_kinetic,e_potential;
            e_potential = grav_acc.norm_sq();
            e_potential = P.Mass[i] * sqrt(e_potential) * (Get_Particle_Size(i)*All.cf_atime); // = M*|a_grav|*h (physical)
            e_kinetic = 0.5 * P.Mass[i] * All.cf_a2inv * CellP.MaxKineticEnergyNgb[i];
            e_thermal = DMAX(0.5*CellP.InternalEnergy[i], dEnt) * P.Mass[i];
#ifdef MAGNETIC
            e_thermal += 0.5*CellP.B[i].norm_sq()*CellP.Density[i]/(All.cf_atime*P.Mass[i]);
#endif
            int do_entropy = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
            if(0.01*(e_thermal+e_kinetic) > e_thermal) {do_entropy=1;}
#else
            if(0.005*(e_thermal+e_kinetic) > e_thermal) {do_entropy=1;}
#endif
            do_entropy = 0;
            if(0.01*e_potential > e_thermal) {do_entropy=1;}
            // note that for the Zeldovich problem, either the gravity or kinetic energy switch is sufficient for good resolution;
            //  both are not needed. we find slightly cleaner results on that test keeping the gravity and removing the KE switch
            
            // also check for flows which are totally dominated by the adiabatic component of their temperature evolution //
            // double mach = fabs(CellP.MaxSignalVel[i]/Get_Gas_effective_soundspeed_i(i) - 2.0); //
            // if(mach < 1.1) {do_entropy=1;} // (actually, this switch tends to do more harm than good!) //
            //do_entropy = 0; // seems unstable in tests like interacting blastwaves... //
            if(do_entropy)
            {
                /* use the pure-SPH entropy equation, which is exact up to the mass flux, for adiabatic flows */
                CellP.DtInternalEnergy[i] = -(CellP.Pressure[i]/CellP.Density[i]) * P.Particle_DivVel[i]*All.cf_a2inv;
#ifdef MAGNETIC
                CellP.DtB[i] = CellP.B[i] * ((1./3.) * All.cf_atime * P.Particle_DivVel[i]*All.cf_a2inv);
#ifdef DIVBCLEANING_DEDNER
                CellP.DtPhi[i] = (1./3.) * (CellP.Phi[i]*All.cf_a3inv) * P.Particle_DivVel[i]*All.cf_a2inv; // cf_a3inv from mass-based phi-fluxes
#endif
#endif
                if(All.ComovingIntegrationOn) {CellP.DtInternalEnergy[i] -= 3*(GAMMA(i)-1) * CellP.InternalEnergyPred[i] * All.cf_hubble_a;}
                dEnt = CellP.InternalEnergy[i] + CellP.DtInternalEnergy[i] * dt_hydrokick; /* gravity term not included here, as it makes this unstable */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                CellP.dMass[i] = CellP.DtMass[i] = 0;
#endif
            }
#endif // closes ENERGY_ENTROPY_SWITCH_IS_ACTIVE
            
#ifdef HYDRO_EXPLICITLY_INTEGRATE_VOLUME
            CellP.Density_ExplicitInt[i] *= exp(-DMIN(1.5,DMAX(-1.5,P.Particle_DivVel[i]*All.cf_a2inv * dt_hydrokick))); /*!< explicitly integrated volume/density variable to be used if integrating the SPH-like form of the continuity directly */
            if(CellP.FaceClosureError[i] > 0) {double drho2 = CellP.Gradients.Density[i].norm_sq(); /* the evolved density evolves back to the explicit density on a relaxation time of order the sound-crossing or tension wave-crossing time across the density gradient length */
                if(drho2>0 && CellP.Density_ExplicitInt[i]>0 && CellP.Density[i]>0) {
                    double Lgrad = CellP.Density[i] / sqrt(drho2); Lgrad=DMAX(Lgrad,P.KernelRadius[i]); double cs_eff_forrestoringforce=Get_Gas_effective_soundspeed_i(i); /* gradient scale length and sound speed */
#if defined(EOS_TILLOTSON)
                    cs_eff_forrestoringforce=DMIN(cs_eff_forrestoringforce , sqrt(All.Tillotson_EOS_params[CellP.CompositionType[i]][10] / CellP.Density[i])); /* speed of deviatoric waves, which is most relevant, if defined */
#endif
                    double delta = 0.1 * dt_hydrokick * cs_eff_forrestoringforce / Lgrad, q0=log(CellP.Density_ExplicitInt[i]), q1=log(P.Mass[i]/CellP.FaceClosureError[i]), qn=0; if(delta > 0.005) {qn=q0*exp(-delta) + q1*(1.-exp(-delta));} else {qn=q0 + (q1-q0)*delta*(1.-0.5*delta);} /* evolves in log-space across this span */
                    CellP.Density_ExplicitInt[i] = exp(q0); /* set final density */
                }}
#endif

#ifdef RADTRANSFER /* block here to deal with tricky cases where radiation energy density is -much- larger than thermal, re-distribute the energy that would have taken us negative in gas back into radiation */
            int kfreq; double erad_tot=0,emin=0,enew=0,demin=0,dErad=0,rsol_fac=C_LIGHT_CODE_REDUCED(i)/C_LIGHT_CODE;  for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {erad_tot+=CellP.Rad_E_gamma[i][kfreq];}
            if(erad_tot > 0) // do some checks if this helps or hurts (identical setup in predict) - seems relatively ok for now, in new form
            {
                demin=0.025*CellP.InternalEnergy[i]; emin=0.025*(erad_tot/rsol_fac + CellP.InternalEnergy[i]*P.Mass[i]); enew=DMAX(erad_tot/rsol_fac + dEnt*P.Mass[i], emin);
                dEnt=(enew - erad_tot/rsol_fac) / P.Mass[i]; if(dEnt < demin) {dErad=rsol_fac*(dEnt-demin); dEnt=demin;}
                if(dErad<-0.975*erad_tot) {dErad=-0.975*erad_tot;} CellP.InternalEnergy[i] = dEnt; for(kfreq=0;kfreq<N_RT_FREQ_BINS;kfreq++) {CellP.Rad_E_gamma[i][kfreq] *= 1 + dErad/erad_tot;}
            } else {
                if(dEnt < 0.5*CellP.InternalEnergy[i]) {CellP.InternalEnergy[i] *= 0.5;} else {CellP.InternalEnergy[i] = dEnt;}
            }
#else
            if(dEnt < 0.5*CellP.InternalEnergy[i]) {CellP.InternalEnergy[i] *= 0.5;} else {CellP.InternalEnergy[i] = dEnt;}
#endif
            check_particle_for_temperature_minimum(i); /* if we've fallen below the minimum temperature, force the 'floor' */
        }
        
        /* now, kick for non-gas/fluid quantities (accounting for momentum conservation if masses are changing) */
        dp = {};
        if(P.Type[i]==0)
        {
            dp += P.GravAccel[i] * (mass_pred * dt_gravkick);
            dp += CellP.HydroAccel[i] * (mass_pred * All.cf_atime * dt_hydrokick); // convert to code units
#ifdef TURB_DRIVING
            dp += CellP.TurbAccel[i] * (mass_pred * dt_gravkick);
#endif
#ifdef RT_RAD_PRESSURE_OUTPUT
            dp += CellP.Rad_Accel[i] * (mass_pred * All.cf_atime * dt_hydrokick);
#endif
        } else {
            dp += P.GravAccel[i] * (mass_pred * dt_gravkick);
        }
#if (SINGLE_STAR_TIMESTEPPING > 0)  //if we're super-timestepping, the above accounts for the change in COM velocity. Now we do the internal binary velocity change
        if((P.Type[i] == 5) && (P.SuperTimestepFlag[i]>=2)) {dp += (P.COM_GravAccel[i] - P.GravAccel[i]) * (mass_pred * dt_gravkick);}
#endif
#ifdef HERMITE_INTEGRATION
        // we augment this to a whole-step kick for the initial Hermite prediction step, which is done alongside the first half-step kick.
        if((1<<P.Type[i]) & HERMITE_INTEGRATION)
        {
            if(mode == 0)
            {
                P.OldVel[i] = P.Vel[i];
                P.OldPos[i] = P.Pos[i];
                P.OldJerk[i] = P.GravJerk[i];
                P.Hermite_OldAcc[i] = P.GravAccel[i]; // this is the value from the first Hermite tree pass for this timestep
            }
        }
#endif
        P.Vel[i] += dp * (1.0 / mass_new); /* correctly accounts for mass change if its allowed */

#ifdef DILATION_FOR_STELLAR_KINEMATICS_ONLY
        double a_fac = return_timestep_dilation_factor(i,0);
        if(a_fac > 1.) {
            double cfac = dt_gravkick * (1. - 1./a_fac);
            P.Vel[i] += P.acc_of_nearest_special[i] * cfac; /* add back in the mean kick of the surrounding stuff to dilate only the local dynamics */
        }
#endif

 
        /* check for reflecting or outflow or otherwise special boundaries: if so, do the reflection/boundary! */
        apply_special_boundary_conditions(i,mass_new,1);
        if(P.Mass[i] <= 0 || !isfinite(P.Mass[i])) {return;} /* exit if we have zero'd the particle mass, to avoid errors with dividing by zero */

        /* any other gas-specific kicks (e.g. B-fields, radiation) go here */
        if(P.Type[i]==0)
        {
            do_kick_for_extra_physics(i, tstart, tend, dt_entr);

            /* after completion of a full step, set the predicted values of gas/fluid quantities
             * to the current values. They will then predicted further along in drift operations */
            if(mode==1)
            {
#ifdef HYDRO_GENERATE_TARGET_MESH // it is often desirable to damp transient velocities when setting up a stable mesh: do so here by un-commenting the line below //
                //for(j=0;j<3;j++) {P.Vel[i][j] *= exp(-0.15);} // coefficient is constant per-timestep: adjust to make as aggressive or weak as desired //
#endif

                CellP.VelPred[i] = P.Vel[i]; //(mass_old*v_old[j] + dp[j]) / mass_new;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                P.Mass[i] = CellP.MassTrue[i]; //mass_old + CellP.DtMass[i] * dt_hydrokick;
#endif
                CellP.InternalEnergyPred[i] = CellP.InternalEnergy[i]; //ent_old + CellP.DtInternalEnergy[i] * dt_entr;
#ifdef HYDRO_EXPLICITLY_INTEGRATE_VOLUME
                CellP.Density[i] = CellP.Density_ExplicitInt[i]; /*!< explicitly integrated volume/density variable to be used if integrating the SPH-like form of the continuity directly */
#endif
            }
        }
        
        /* set the momentum shift so we know how to move the tree! */
        P.dp[i] += dp;
#ifdef DM_FUZZY
        do_dm_fuzzy_drift_kick(i, dt_entr, 0); /* kicks for fuzzy-dm integration */
#endif
#ifdef CBE_INTEGRATOR
        do_cbe_drift_kick(i, dt_entr); /* kicks for cbe integration of phase-space distribution function */
#endif
        
    } // if(TimeBinActive[P.TimeBin[i]]) //
}


void set_predicted_quantities_for_extra_physics(int i)
{
    if(P.Type[i] == 0 && P.Mass[i] > 0)
    {
        int k, kf; k=0, kf=0;
#if defined(MAGNETIC)
#ifndef MHD_ALTERNATIVE_LEAPFROG_SCHEME
        CellP.BPred[i] = CellP.B[i];
#if defined(DIVBCLEANING_DEDNER)
        CellP.PhiPred[i] = CellP.Phi[i];
#endif
#endif
#endif
#ifdef COSMIC_RAY_FLUID
        for(kf=0;kf<N_CR_PARTICLE_BINS;kf++)
        {
            CellP.CosmicRayEnergyPred[i][kf] = CellP.CosmicRayEnergy[i][kf];
            CellP.CosmicRayFluxPred[i][kf] = CellP.CosmicRayFlux[i][kf];
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(k=0;k<2;k++) {CellP.CosmicRayAlfvenEnergyPred[i][kf][k] = CellP.CosmicRayAlfvenEnergy[i][kf][k];}
#endif
        }
#endif
        
#if defined(RT_EVOLVE_ENERGY)
        for(kf=0;kf<N_RT_FREQ_BINS;kf++)
        {
            CellP.Rad_E_gamma_Pred[i][kf] = CellP.Rad_E_gamma[i][kf];
#if defined(RT_EVOLVE_FLUX)
            CellP.Rad_Flux_Pred[i][kf] = CellP.Rad_Flux[i][kf];
#endif
        }
        rt_eddington_update_calculation(i);
#endif
#ifdef RT_EVOLVE_INTENSITIES
        for(kf=0;kf<N_RT_FREQ_BINS;kf++) {for(k=0;k<N_RT_INTENSITY_BINS;k++) {CellP.Rad_Intensity_Pred[i][kf][k] = CellP.Rad_Intensity[i][kf][k];}}
#endif

#ifdef EOS_ELASTIC
        for(k=0;k<3;k++) {for(kf=0;kf<3;kf++) {CellP.Elastic_Stress_Tensor_Pred[i][k][kf]=CellP.Elastic_Stress_Tensor[i][k][kf];}}
#endif
        
        set_eos_pressure(i);
    }
}



void do_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr)
{
    int j; j=0;
#ifdef MAGNETIC
#ifndef MHD_ALTERNATIVE_LEAPFROG_SCHEME
    double BphysVolphys_to_BcodeVolCode = 1 / All.cf_atime;
    CellP.B[i] += CellP.DtB[i] * (dt_entr * BphysVolphys_to_BcodeVolCode); // fluxes are always physical, convert to code units //
#ifdef DIVBCLEANING_DEDNER
    double PhiphysVolphys_to_PhicodeVolCode = 1 / All.cf_a3inv; // for mass-based phi-fluxes (otherwise is just "1")
    /* phi units are [vcode][Bcode]=a^3 * vphys*Bphys */
    if(CellP.Density[i] > 0)
    {
        /* now we're going to check for physically reasonable phi values */
        double cs_phys = Get_Gas_effective_soundspeed_i(i);
        double b_phys = 0.0;
        for(j = 0; j < 3; j++) {b_phys += Get_Gas_BField(i,j)*Get_Gas_BField(i,j);}
        b_phys = sqrt(b_phys)*All.cf_a2inv;
        double vsig1 = sqrt(cs_phys*cs_phys + b_phys*b_phys/(CellP.Density[i]*All.cf_a3inv));
        double vsig2 = 0.5 * fabs(CellP.MaxSignalVel[i]);
        double vsig_max = DMAX( DMAX(vsig1,vsig2) , All.FastestWaveSpeed );
        double phi_phys_abs = fabs(Get_Gas_PhiField(i)) * All.cf_a3inv;
        double vb_phy_abs = vsig_max * b_phys;

        if((!isnan(CellP.DtPhi[i]))&&(phi_phys_abs>0)&&(vb_phy_abs>0)&&(!isnan(phi_phys_abs))&&(!isnan(vb_phy_abs)))
        {
            double phi_max_tolerance = 10.0;
            if(phi_phys_abs > 1000. * phi_max_tolerance * vb_phy_abs)
            {
                /* this can indicate a problem! issue a warning and zero phi */
                if(phi_phys_abs > 1.0e6 * phi_max_tolerance * vb_phy_abs) {
                    PRINT_WARNING("significant growth detected in phi-field: phi_phys_abs=%g vb_phy_abs=%g vsig_max=%g b_phys=%g particle_id_i=%d dtphi_code=%g Pressure=%g rho=%g x/y/z=%g/%g/%g vx/vy/vz=%g/%g/%g Bx/By/Bz=%g/%g/%g h=%g u=%g m=%g phi=%g bin=%d SigVel=%g a=%g \n",
                       phi_phys_abs,vb_phy_abs,vsig_max,b_phys,i,CellP.DtPhi[i],CellP.Pressure[i],CellP.Density[i],P.Pos[i][0],P.Pos[i][1],P.Pos[i][2],
                       P.Vel[i][0],P.Vel[i][1],P.Vel[i][2],CellP.B[i][0],CellP.B[i][1],CellP.B[i][2],
                       P.KernelRadius[i],CellP.InternalEnergy[i],P.Mass[i],CellP.Phi[i],P.TimeBin[i],CellP.MaxSignalVel[i],All.cf_atime);}
                CellP.PhiPred[i] = CellP.Phi[i] = CellP.DtPhi[i] = 0;
            } else {
                if(phi_phys_abs > phi_max_tolerance * vb_phy_abs)
                {
                    /* in this limit, only allow for decay of phi: to avoid over-shooting, we apply the force as damping */
                    if(CellP.Phi[i] > 0) {CellP.DtPhi[i]=DMIN(CellP.DtPhi[i],0);} else {CellP.DtPhi[i]=DMAX(CellP.DtPhi[i],0);}
                    double dtphi_code = dt_entr * PhiphysVolphys_to_PhicodeVolCode * CellP.DtPhi[i];
                    if(CellP.Phi[i] != 0) {CellP.Phi[i] *= exp( - fabs(dtphi_code) / fabs(CellP.Phi[i]) );}
                } else {
                    /* ok, in this regime, we're safe to apply the 'normal' time evolution */
                    double dtphi_code = dt_entr * PhiphysVolphys_to_PhicodeVolCode * CellP.DtPhi[i];
                    CellP.Phi[i] += dtphi_code;
                }
            }
        }
    } else {
        CellP.Phi[i] = CellP.PhiPred[i] = CellP.DtPhi[i] = 0;
    }
    /* now apply the usual damping term */
    double t_damp = Get_Gas_PhiField_DampingTimeInv(i);
    if((t_damp>0) && (!isnan(t_damp)) && (dt_entr>0))
    {
        CellP.Phi[i] *= exp( -dt_entr * t_damp );
    }
    if(isnan(CellP.DtPhi[i])) {CellP.DtPhi[i]=0;}
    if(isnan(CellP.Phi[i])) {CellP.Phi[i]=0;}
    if(isnan(CellP.PhiPred[i])) {CellP.PhiPred[i]=CellP.Phi[i];}
#endif
#endif
#endif
    
#ifdef COSMIC_RAY_FLUID
    CosmicRay_Update_DriftKick(i,dt_entr,0);
#endif
    
#ifdef RADTRANSFER
    rt_update_driftkick(i,dt_entr,0);
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
    if(P.Pos[i][2] > DMIN(19., DMAX(1.1*All.Time*C_LIGHT_CODE_REDUCED(i), DMIN(18.*boxSize_X + (All.Vertical_Grain_Accel*All.Dust_to_Gas_Mass_Ratio - All.Vertical_Gravity_Strength)*All.Time*All.Time/2., 19.)))) {for(j=0;j<N_RT_FREQ_BINS;j++) {CellP.Rad_E_gamma[i][j]*=0.5; CellP.Rad_E_gamma_Pred[i][j]*=0.5;
#ifdef RT_EVOLVE_FLUX
        if(CellP.Rad_Flux[i][j][2] < 0) {CellP.Rad_Flux[i][j][2]=-CellP.Rad_Flux[i][j][2]; CellP.Rad_Flux_Pred[i][j][2]=CellP.Rad_Flux[i][j][2];}
#endif
    }}
#endif
#endif

#ifdef EOS_ELASTIC
    elastic_body_update_driftkick(i,dt_entr,0);
#endif
}

    
    
void apply_special_boundary_conditions(int i, double mass_for_dp, int mode)
{
#if BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE
    double box_upper[3]; int j;
    box_upper[0]=boxSize_X; box_upper[1]=boxSize_Y; box_upper[2]=boxSize_Z;
    for(j=0; j<3; j++)
    {
        if(P.Pos[i][j] <= 0)
        {
            if(special_boundary_condition_xyz_def_reflect[j] == 0 || special_boundary_condition_xyz_def_reflect[j] == -1)
            {
                if(P.Vel[i][j]<0) {P.Vel[i][j]=-P.Vel[i][j]; if(P.Type[i]==0) {CellP.VelPred[i][j]=P.Vel[i][j]; CellP.HydroAccel[i][j]=0;} if(mode==1) {P.dp[i][j]+=2*P.Vel[i][j]*mass_for_dp;}}
                P.Pos[i][j]=DMAX((0.+((double)P.ID[i])*2.e-8)*box_upper[j], 0.1*P.Pos[i][j]); // old  was 1e-9, safer on some problems, but can artificially lead to 'trapping' in some low-res tests
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                P.Pos[i][j]+=3.e-3*boxSize_X; P.Vel[i][j] += 0.1; /* special because of our wierd boundary condition for this problem, sorry to have so many hacks for this! */
#endif
#ifdef RT_EVOLVE_FLUX
                if(P.Type[i]==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {if(CellP.Rad_Flux[i][kf][j]<0) {CellP.Rad_Flux[i][kf][j]=-CellP.Rad_Flux[i][kf][j]; CellP.Rad_Flux_Pred[i][kf][j]=CellP.Rad_Flux[i][kf][j];}}}
#endif

#ifdef COSMIC_RAY_FLUID
                if(P.Type[i]==0) {int kf; for(kf=0;kf<N_CR_PARTICLE_BINS;kf++) {if(CellP.CosmicRayFlux[i][kf][j]<0) {CellP.CosmicRayFlux[i][kf][j]=-CellP.CosmicRayFlux[i][kf][j]; CellP.CosmicRayFluxPred[i][kf][j]=CellP.CosmicRayFlux[i][kf][j];}}}
#endif
            }
            if(special_boundary_condition_xyz_def_outflow[j] == 0 || special_boundary_condition_xyz_def_outflow[j] == -1) {P.Mass[i]=0; if(mode==1) {P.dp[i][0]=P.dp[i][1]=P.dp[i][2]=0;}}
        }
        else if (P.Pos[i][j] >= box_upper[j])
        {
            if(special_boundary_condition_xyz_def_reflect[j] == 0 || special_boundary_condition_xyz_def_reflect[j] == 1)
            {
                if(P.Vel[i][j]>0) {P.Vel[i][j]=-P.Vel[i][j]; if(P.Type[i]==0) {CellP.VelPred[i][j]=P.Vel[i][j]; CellP.HydroAccel[i][j]=0;} if(mode==1) {P.dp[i][j]+=2*P.Vel[i][j]*mass_for_dp;}}
                P.Pos[i][j]=box_upper[j]*(1.-((double)P.ID[i])*2.e-8);
#ifdef RT_EVOLVE_FLUX
                if(P.Type[i]==0) {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {if(CellP.Rad_Flux[i][kf][j]>0) {CellP.Rad_Flux[i][kf][j]=-CellP.Rad_Flux[i][kf][j]; CellP.Rad_Flux_Pred[i][kf][j]=CellP.Rad_Flux[i][kf][j];}}}
#endif
#ifdef COSMIC_RAY_FLUID
                if(P.Type[i]==0) {int kf; for(kf=0;kf<N_CR_PARTICLE_BINS;kf++) {if(CellP.CosmicRayFlux[i][kf][j]>0) {CellP.CosmicRayFlux[i][kf][j]=-CellP.CosmicRayFlux[i][kf][j]; CellP.CosmicRayFluxPred[i][kf][j]=CellP.CosmicRayFlux[i][kf][j];}}}
#endif
            }
            if(special_boundary_condition_xyz_def_outflow[j] == 0 || special_boundary_condition_xyz_def_outflow[j] == 1) {P.Mass[i]=0; if(mode==1) {P.dp[i][0]=P.dp[i][1]=P.dp[i][2]=0;}}
        }
    }
#endif
    return;
}
