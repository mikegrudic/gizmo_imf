#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "../declarations/allvars.h"
#include "../core/proto.h"
#include "../mesh/kernel.h"


/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */
/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified heavily
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO: initializing
 * new/modified variables, rewriting order-of-operations, standardizing
 * where some calls occur, and adding new subroutines that need to be called
 * on startup now.
 */

/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If gas cells are present,
 *  the initial gas kernel lengths are determined.
 */
void init(void)
{
    int i, j; double a3, atime, a2_fac;

#ifdef MAGNETIC
    double gauss2gizmo = All.UnitMagneticField_in_gauss / UNIT_B_IN_GAUSS;
    /* NOTE: we will always work -internally- in code units where MU_0 = 1; hence the 4pi here; [much simpler, but be sure of your conversions!] */
#endif

#ifdef SINK_PARTICLES
    int count_holes = 0;
#endif

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();

    if(RestartFlag != 1) {All.MinMassForParticleMerger = 0; All.MaxMassForParticleSplit = 0;}

    if(RestartFlag == 3 && RestartSnapNum < 0)
    {
        if(ThisTask == 0) {printf("Need to give the snapshot number if FOF/SUBFIND is selected for output\n");}
        endrun(0);
    }

    if(RestartFlag == 4 && RestartSnapNum < 0)
    {
        if(ThisTask == 0) {printf("Need to give the snapshot number if snapshot should be converted\n");}
        endrun(0);
    }

    if(RestartFlag == 5 && RestartSnapNum < 0)
    {
        if(ThisTask == 0) {printf("Need to give the snapshot number if power spectrum and two-point correlation function should be calculated\n");}
        endrun(0);
    }

    if(RestartFlag == 6 && RestartSnapNum < 0)
    {
        if(ThisTask == 0) {printf("Need to give the snapshot number if velocity power spectrum for the gas cells should be calculated\n");}
        endrun(0);
    }


    switch (All.ICFormat)
    {
        case 1:
        case 2:
        case 3:
        case 4:
            if(RestartFlag >= 2 && RestartSnapNum >= 0)
            {
                char fname[MAX_PATH_BUFFERSIZE_TOUSE];
                if(All.NumFilesPerSnapshot > 1) {snprintf(fname, MAX_PATH_BUFFERSIZE_TOUSE, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase, RestartSnapNum);}
                    else {snprintf(fname, MAX_PATH_BUFFERSIZE_TOUSE, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);}
                read_ic(fname);

            }
            else {read_ic(All.InitCondFile);}
            break;

        default:
            if(ThisTask == 0) {printf("ICFormat=%d not supported.\n", All.ICFormat);}
            endrun(0);
    }

#ifdef CHIMES_INITIALISE_IN_EQM
    for (i = 0; i < N_gas; i++) {allocate_gas_abundances_memory(&(ChimesGasVars[i]), &ChimesGlobalVars);}
#endif

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();


#if defined(COOLING) && !defined(CHIMES)
    IonizeParams();
#endif

    All.Ti_Current = 0;
    if(All.ComovingIntegrationOn)
    {
        All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
        a3 = All.Time * All.Time * All.Time; atime = All.Time; a2_fac = (All.Time * All.Time);
    }
    else
    {
        All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
        a3 = atime = a2_fac = 1;
    }

    set_softenings();

    All.NumCurrentTiStep = 0;	/* setup some counters */
    All.SnapshotFileCount = 0;
    if(RestartFlag == 2)
    {
        if(RestartSnapNum < 0)
        {
            char *underscore = strrchr(All.InitCondFile, '_');
            if(!underscore)
            {
                char buf[DEFAULT_PATH_BUFFERSIZE_TOUSE];
                snprintf(buf, DEFAULT_PATH_BUFFERSIZE_TOUSE, "Your input file '%s' lacks an underscore. Cannot infer next snapshot number.\n", All.InitCondFile);
                terminate(buf);
            }
            else {All.SnapshotFileCount = atoi(underscore + 1) + 1;}
        }
        else {All.SnapshotFileCount = RestartSnapNum + 1;}
    }

#ifdef OUTPUT_LINEOFSIGHT
    All.Ti_nextlineofsight = (int) (log(All.TimeFirstLineOfSight / All.TimeBegin) / All.Timebase_interval);
    if(RestartFlag == 2) {endrun(78787);}
#endif

    All.TotNumOfForces = 0;
    All.TopNodeAllocFactor = 0.008; /* this will start from a low value and be iteratively increased until it is well-behaved */
#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    All.TopNodeAllocFactor = 0.1; /* for optimization on startup this needs to be increased for these extreme dynamic range runs */
#endif
    All.TreeAllocFactor = 0.45; /* this will also iteratively increase to fit the particle distribution */
    /* To construct the BH-tree for N particles, somewhat less than N
     internal tree-nodes are necessary for ‘normal’ particle distributions.
     TreeAllocFactor sets the number of internal tree-nodes allocated in units of the particle number.
     By experience, space for ~0.65N internal nodes is usually fully sufficient for typical clustered
     particle distributions, so a value of 0.7 should put you on the safe side. If the employed particle
     number per processor is very small (less than a thousand or so), or if there are many particle pairs
     with identical or nearly identical coordinates, a higher value may be required. Since the number of
     particles on a given processor may be higher by a factor PartAllocFactor than the average particle
     number, the total amount of memory requested for the Barnes-Hut tree on a single processor scales proportional
     to PartAllocFactor*TreeAllocFactor. */

#ifdef SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM
    for(i = 0; i < SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM; i++)
    {
        All.SpecialParticle_Position_ForRefinement[i][0]=All.SpecialParticle_Position_ForRefinement[i][1]=All.SpecialParticle_Position_ForRefinement[i][2]=0;
        All.Mass_Accreted_By_SpecialParticle[i]=0; All.Mass_of_SpecialParticle[i]=0;
    }
#endif
    
#ifdef GALSF_LIMIT_FBTIMESTEPS_FROM_BELOW
    if(RestartFlag != 1) {All.Dt_Since_LastFBCalc_Gyr=0; All.Dt_Min_Between_FBCalc_Gyr=((double)(GALSF_LIMIT_FBTIMESTEPS_FROM_BELOW))/1.e9;}
#endif

#ifdef BOX_PERIODIC
    if(All.ComovingIntegrationOn) {check_omega();}
#endif
    All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#if (defined(SINK_PARTICLES) || defined(GALSF_SUBGRID_WINDS)) && defined(FOF)
    All.TimeNextOnTheFlyFoF = All.TimeBegin;
#endif

    for(i = 0; i < GRAVCOSTLEVELS; i++) {All.LevelToTimeBin[i] = 0;}

    for(i = 0; i < NumPart; i++) {for(j = 0; j < GRAVCOSTLEVELS; j++) {P.GravCost[i][j] = 0;}}

    if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
        {for(i=0;i<NumPart;i++) {P.Vel[i] *= sqrt(All.Time)*All.Time;}}

#ifdef DM_SIDM
    init_self_interactions();
#endif

#ifdef METALS
    for(j=0;j<NUM_METAL_SPECIES;j++) {All.SolarAbundances[j]=0;} // initialize all to zero
    All.SolarAbundances[0]=0.02;        // all metals (by mass); present photospheric abundances from Asplund et al. 2009 (Z=0.0134, proto-solar=0.0142) in notes;
    //   also Anders+Grevesse 1989 (older, but hugely-cited compilation; their Z=0.0201, proto-solar=0.0213)
#ifdef COOL_METAL_LINES_BY_SPECIES
    All.SolarAbundances[1]=0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314), with proto-solar Y=0.27
    All.SolarAbundances[2]=3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3); proto-solar from Asplund=8.47 -> 2.53e-3
    All.SolarAbundances[3]=1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3); PS=7.87->7.41e-4
    All.SolarAbundances[4]=8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3); PS=8.73->6.13e-3
    All.SolarAbundances[5]=2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3); PS=7.97->1.34e-3
    All.SolarAbundances[6]=9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4); PS=7.64->7.57e-4
    All.SolarAbundances[7]=1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4); PS=7.55->7.12e-4
    All.SolarAbundances[8]=6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4); PS=7.16->3.31e-4
    All.SolarAbundances[9]=1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4); PS=6.38->6.87e-5
    All.SolarAbundances[10]=1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3); PS=7.54->1.38e-3
#if (GALSF_FB_FIRE_STELLAREVOLUTION > 2) // new default abundances; using Asplund et al. 2009 proto-solar abundances
    All.SolarAbundances[0]=0.0142; if(NUM_METAL_SPECIES>=10) {
        All.SolarAbundances[1]=0.27030; All.SolarAbundances[2]=2.53e-3; All.SolarAbundances[3]=7.41e-4; All.SolarAbundances[4]=6.13e-3; All.SolarAbundances[5]=1.34e-3;
        All.SolarAbundances[6]=7.57e-4; All.SolarAbundances[7]=7.12e-4; All.SolarAbundances[8]=3.31e-4; All.SolarAbundances[9]=6.87e-5; All.SolarAbundances[10]=1.38e-3;}
#endif
#endif
#endif


    for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
        P.GravAccel[i] = {};

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE /* init tidal tensor for first output (not used for calculation) */
        for(j=0;j<3;j++) {int kt; for(kt=0;kt<3;kt++) {P.tidal_tensorps[i][j][kt]=0;}}
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        P.tidal_tensor_mag_prev[i] = 0; P.tidal_zeta[i]=0; for(j=0;j<3;j++) {int kt; for(kt=0;kt<3;kt++) {P.tidal_tensorps_prevstep[i][j][kt]=0;}}
#endif
#ifdef PMGRID
        for(j=0;j<3;j++) {int kt; for(kt=0;kt<3;kt++) {P.tidal_tensorpsPM[i][j][kt]=0;}}
#endif
#endif

#ifdef ADAPTIVE_TREEFORCE_UPDATE
        P.time_since_last_treeforce[i] = 0;
        P.tdyn_step_for_treeforce[i] = 0;
#endif        
        

#ifdef PMGRID
        P.GravPM[i] = {};
#endif
        P.Ti_begstep[i] = 0;
        P.Ti_current[i] = (integertime)0;
        P.TimeBin[i] = 0;
        if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS) {P.OldAcc[i] = 0;}	/* Do not zero in 2lpt case as masses are stored here */

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)
        P.Potential[i] = 0;
#endif
#ifdef GALSF
        if(RestartFlag == 0) {P.StellarAge[i] = 0;}
#ifdef GALSF_SFR_IMF_VARIATION
        if(RestartFlag == 0) {P.IMF_Mturnover[i] = 2.0;} /* gives a solar-type IMF for our calculations in current code */
#endif
#ifdef GALSF_SFR_IMF_SAMPLING
        if(RestartFlag == 0) {P.IMF_NumMassiveStars[i] = 0;}
#if defined(SINGLE_STAR_AND_SSP_HYBRID_MODEL) && defined(SINGLE_STAR_RESTART_FROM_FIRESIM)
        if(RestartFlag == 2) {P.IMF_NumMassiveStars[i] = 0;}
#endif
#endif
#endif
        
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION
        if(RestartFlag == 0) {P.ProtoStellarStage[i] = 0;}
#endif

        if(RestartFlag != 1)
        {
#if defined(DO_DENSITY_AROUND_NONGAS_PARTICLES)
            P.DensityAroundParticle[i] = 0;
            P.GradRho[i][0]=0;
            P.GradRho[i][1]=0;
            P.GradRho[i][2]=1;
#endif
#if defined(SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION)
#if defined(SINGLE_STAR_FB_SNE)
            P.Mass_final[i] = P.Mass[i]; // best guess, only matters if we restart in the middle of spawning an SN
#endif
#if defined(SINGLE_STAR_FB_WINDS)
            P.wind_mode[i] = 0; // this will make single_star_wind_mdot reset it
            Vec3<double> nx,ny,nz; int kw; get_random_orthonormal_basis(P.ID[i],nx,ny,nz); for(kw=0;kw<3;kw++) {P.Wind_direction[i][kw] = nx[kw]; P.Wind_direction[i][kw+3] = ny[kw];}
#endif
#endif
#if defined(GALSF_FB_FIRE_RT_LOCALRP)
            P.NewStar_Momentum_For_JetFeedback[i] = 0;
#endif
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
            P.SNe_ThisTimeStep[i] = 0;
#endif
#ifdef GALSF_FB_MECHANICAL
            int k; for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {P.Area_weighted_sum[i][k] = 0;}
#ifdef GALSF_FB_FIRE_STELLAREVOLUTION
            P.MassReturn_ThisTimeStep[i] = 0;
#ifdef GALSF_FB_FIRE_RPROCESS
            P.RProcessEvent_ThisTimeStep[i] = 0;
#endif
#ifdef GALSF_FB_FIRE_AGE_TRACERS
            P.AgeDeposition_ThisTimeStep[i] = 0;
#endif
#endif
#endif
#if defined(FIRE_SUPERLAGRANGIAN_JEANS_REFINEMENT) || defined(SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM)
            P.Time_Of_Last_MergeSplit[i] = All.TimeBegin;
#endif
#ifdef SPECIAL_POINT_WEIGHTED_MOTION
            P.Time_Of_Last_SmoothedVelUpdate[i] = All.TimeBegin;
#endif
        }

#if defined(INIT_STELLAR_METALS_AGES_DEFINED) && defined(GALSF)
        if(RestartFlag == 0) {P.StellarAge[i] = -2.0 * All.InitStellarAgeinGyr / (UNIT_TIME_IN_GYR) * get_random_number(P.ID[i] + 3);}
#endif
        
#ifdef GRAIN_FLUID
        if((RestartFlag == 0) && ((1 << P.Type[i]) & (GRAIN_PTYPES)))
        {
            int grain_subtype = 1; P.Grain_Size[i] = 0; /* default assumption about particulate sub-type for operations below */
#if defined(PIC_MHD)
            grain_subtype = P.MHD_PIC_SubType[i]; /* check if the 'grains' are really PIC elements */
#endif
            /* Change grain mass to change the distribution of sizes.  Grain_Size_Spectrum_Powerlaw parameter sets d\mu/dln(R_d) ~ R_d^Grain_Size_Spectrum_Powerlaw */
            if(grain_subtype <= 2)
            {
                P.Grain_Size[i] = All.Grain_Size_Min * exp( gsl_rng_uniform(random_generator) * log(All.Grain_Size_Max/All.Grain_Size_Min) );
                if(All.Grain_Size_Max > All.Grain_Size_Min*1.0001 && fabs(All.Grain_Size_Spectrum_Powerlaw) != 0) {
                    P.Mass[i] *= (All.Grain_Size_Spectrum_Powerlaw/(pow(All.Grain_Size_Max/All.Grain_Size_Min,All.Grain_Size_Spectrum_Powerlaw)-1.)) *
                    pow(P.Grain_Size[i]/All.Grain_Size_Min,All.Grain_Size_Spectrum_Powerlaw) * log(All.Grain_Size_Max/All.Grain_Size_Min);}
#ifdef GRAIN_RDI_TESTPROBLEM /* initialize various quantities for test problems from parameters set in the ICs */
                P.Mass[i] *= All.Dust_to_Gas_Mass_Ratio;
                int k, non_gdir=1; Vec3<double> A={}, B={}, A_cross_B={}; double amag, rho_gas_expected, cs_gas_expected, acc_ang=All.Vertical_Grain_Accel_Angle * M_PI/180., tS0, a0, ct=1, tau2=0, ct2=0, w0, agamma=9.*M_PI/128.; B[2]=1; if(GRAV_DIRECTION_RDI==1) {non_gdir=2;}
                rho_gas_expected = 1.*UNIT_DENSITY_IN_CGS, cs_gas_expected = 1.*UNIT_VEL_IN_CGS; /* guess for the gas density here [set custom for e.g. stratified problems */
                tS0 = ((0.626657 * sqrt(GAMMA_DEFAULT) * P.Grain_Size[i] * All.Grain_Internal_Density) / (rho_gas_expected * cs_gas_expected)) / UNIT_TIME_IN_CGS; /* stopping time [Epstein] for driftvel->0 [cgs->code units] */
                A[GRAV_DIRECTION_RDI]=cos(acc_ang)*All.Vertical_Grain_Accel - All.Vertical_Gravity_Strength; A[0]=sin(acc_ang)*All.Vertical_Grain_Accel; /* define angles/direction of external acceleration */
                amag=sqrt(A.norm_sq()+MIN_REAL_NUMBER); A /= amag;
                a0 = (tS0 * amag / (1.+All.Dust_to_Gas_Mass_Ratio)) / (cs_gas_expected/UNIT_VEL_IN_CGS) ; /* acc * tS0 / (1+mu) -- we're assuming that the code unit velocity equals the sound speed, for simplicity here */
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
                a0 *= All.Grain_Size_Max / P.Grain_Size[i];
#endif
#ifdef GRAIN_RDI_TESTPROBLEM_LIVE_RADIATION_INJECTION
                double q_a = (0.75*All.Grain_Q_at_MaxGrainSize) / (All.Grain_Internal_Density*All.Grain_Size_Max), kappa_0 = All.Grain_Absorbed_Fraction_vs_Total_Extinction * q_a * All.Dust_to_Gas_Mass_Ratio; // this will be in cgs here
                double rho_base_setup = 1.*UNIT_DENSITY_IN_CGS, H_scale_setup = 1.*boxSize_X*UNIT_LENGTH_IN_CGS; // define in code units the -assumed- initial scaling of the base gas density and vertical scale-length (PROBLEM SPECIFIC HERE!)
#ifdef GRAIN_RDI_TESTPROBLEM_ACCEL_DEPENDS_ON_SIZE
                kappa_0 *= sqrt(All.Grain_Size_Max / All.Grain_Size_Min); // opacity must be corrected for dependence of Q on grainsize or lack thereof
#endif
                a0 *= exp(-kappa_0*rho_base_setup*H_scale_setup*(1.-exp(-P.Pos[i][2]/H_scale_setup))); // attenuate incident flux (and reduce acceleration) according to equilibrium expectation, if we're using single-scattering radiation pressure [otherwise comment this line out] //
#endif
                w0=sqrt((sqrt(1.+4.*agamma*a0*a0)-1.)/(2.*agamma)); // exact solution if no Lorentz forces and Epstein drag //
#ifdef GRAIN_LORENTZFORCE
                double Bmag, tL_i=0, tau2_0=0, f_tau_guess2=0; B = {All.BiniX, All.BiniY, All.BiniZ}; Bmag=B.norm(); B /= Bmag;
                tL_i = (All.Grain_Charge_Parameter*sqrt(rho_gas_expected/UNIT_DENSITY_IN_CGS)/((All.Grain_Internal_Density/UNIT_DENSITY_IN_CGS)*(All.Grain_Size_Max/UNIT_LENGTH_IN_CGS))) * pow(All.Grain_Size_Max/P.Grain_Size[i],2) * Bmag; // 1/t_Lorentz in code units
                ct=dot(A,B); ct2=ct*ct; tau2_0=pow(tS0*tL_i,2); // variables for below //
                for(k=0;k<20;k++)
                {
                   tau2 = tau2_0 / (1. + agamma*w0*w0); // guess tau [including velocity dependence] //
                   f_tau_guess2 = (1.+tau2*ct2) / (1.+tau2); // what the projection factor (reduction in w from projection) would be //
                   w0 = sqrt((sqrt(1.+4.*agamma*a0*a0*f_tau_guess2)-1.)/(2.*agamma)); // re-calculate w0 with this //
                }
#endif
                w0 /= sqrt((1.+tau2)*(1.+tau2*ct2)); // ensures normalization to unity with convention below //
                A_cross_B = cross(A, B);
                P.Vel[i] = (A + A_cross_B*sqrt(tau2) + B*(tau2*ct)) * (w0 * (cs_gas_expected/UNIT_VEL_IN_CGS));
#ifdef BOX_SHEARING
                // now add linearly the NHS drift solution for our shearing box setup
                double v00 = -All.Pressure_Gradient_Accel / (2. * BOX_SHEARING_OMEGA_BOX_CENTER);
                double v_K = -(P.Pos[i][0]-boxHalf_X) * BOX_SHEARING_Q*BOX_SHEARING_OMEGA_BOX_CENTER;
                double tau_s = tS0 * BOX_SHEARING_OMEGA_BOX_CENTER;
                v00 /= (1. + tau_s*tau_s); // appears in both terms here //
                P.Vel[i][0] += v00 * 2.*tau_s; // radial drift
                P.Vel[i][BOX_SHEARING_PHI_COORDINATE] = v_K + v00; // azimuthal drift relative to keplerian frame
#endif
#endif // closes rdi_testproblem
            }
            P.Gas_Density[i] = P.Gas_InternalEnergy[i] = 0; P.Gas_Velocity[i] = {}; P.Grain_AccelTimeMin[i] = MAX_REAL_NUMBER;
#if defined(GRAIN_BACKREACTION)
            P.Grain_DeltaMomentum[i] = {};
#endif
#if defined(GRAIN_LORENTZFORCE)
            P.Gas_B[i] = {};
#endif
        } // closes check on restartflag and particle type
#endif // closes grain_fluid


#ifdef METALS
        if(RestartFlag == 0) {
#if defined(INIT_STELLAR_METALS_AGES_DEFINED)
            P.Metallicity[i][0] = All.InitMetallicityinSolar*All.SolarAbundances[0];
#else
            P.Metallicity[i][0] = 0;
#endif
            /* initialize abundance ratios. for now, assume solar */
            for(j=0;j<NUM_METAL_SPECIES;j++) {P.Metallicity[i][j]=All.SolarAbundances[j]*(P.Metallicity[i][0]/All.SolarAbundances[0]);}
            /* need to allow for a primordial He abundance */
            if(NUM_LIVE_SPECIES_FOR_COOLTABLES>=10) P.Metallicity[i][1]=(1.-HYDROGEN_MASSFRAC)+(All.SolarAbundances[1]-(1.-HYDROGEN_MASSFRAC))*P.Metallicity[i][0]/All.SolarAbundances[0];
        } // if(RestartFlag == 0)

#if defined(GALSF_ISMDUSTCHEM_MODEL)
        Initialize_ISMDustChem_Variables(i);
#endif

#ifdef CHIMES
#ifdef COOL_METAL_LINES_BY_SPECIES
	if (P.Type[i] == 0)
	  {
	    double H_mass_fraction = 1.0 - (P.Metallicity[i][0] + P.Metallicity[i][1]);
	    ChimesGasVars[i].element_abundances[0] = (ChimesFloat) (P.Metallicity[i][1] / (4.0 * H_mass_fraction));   // He
	    ChimesGasVars[i].element_abundances[1] = (ChimesFloat) (P.Metallicity[i][2] / (12.0 * H_mass_fraction));  // C
	    ChimesGasVars[i].element_abundances[2] = (ChimesFloat) (P.Metallicity[i][3] / (14.0 * H_mass_fraction));  // N
	    ChimesGasVars[i].element_abundances[3] = (ChimesFloat) (P.Metallicity[i][4] / (16.0 * H_mass_fraction));  // O
	    ChimesGasVars[i].element_abundances[4] = (ChimesFloat) (P.Metallicity[i][5] / (20.0 * H_mass_fraction));  // Ne
	    ChimesGasVars[i].element_abundances[5] = (ChimesFloat) (P.Metallicity[i][6] / (24.0 * H_mass_fraction));  // Mg
	    ChimesGasVars[i].element_abundances[6] = (ChimesFloat) (P.Metallicity[i][7] / (28.0 * H_mass_fraction));  // Si
	    ChimesGasVars[i].element_abundances[7] = (ChimesFloat) (P.Metallicity[i][8] / (32.0 * H_mass_fraction));  // S
	    ChimesGasVars[i].element_abundances[8] = (ChimesFloat) (P.Metallicity[i][9] / (40.0 * H_mass_fraction));  // Ca
	    ChimesGasVars[i].element_abundances[9] = (ChimesFloat) (P.Metallicity[i][10] / (56.0 * H_mass_fraction)); // Fe

	    ChimesGasVars[i].metallicity = (ChimesFloat) (P.Metallicity[i][0] / 0.0129);  // In Zsol. CHIMES uses Zsol = 0.0129.
	    ChimesGasVars[i].dust_ratio = ChimesGasVars[i].metallicity;
	  }
#else
	if (ThisTask == 0)
	  {
	    printf("ERROR: Config flags CHIMES and METALS are switched on, but COOL_METAL_LINES_BY_SPECIES is switched off. \n");
	    printf("If you want to include metals with CHIMES, you will also need to switch on COOL_METAL_LINES_BY_SPECIES. Aborting. \n");
	    endrun(202);
	  }
#endif // COOL_METAL_LINES_BY_SPECIES
#endif // CHIMES
#else
#ifdef CHIMES
	if (P.Type[i] == 0)
	  {
	    double H_mass_fraction = HYDROGEN_MASSFRAC;
	    ChimesGasVars[i].element_abundances[0] = (ChimesFloat) ((1.0 - H_mass_fraction) / (4.0 * H_mass_fraction));  // He
	    for (j = 1; j < 10; j++) {ChimesGasVars[i].element_abundances[j] = 0.0;}
	    ChimesGasVars[i].metallicity = 0.0;
	    ChimesGasVars[i].dust_ratio = 0.0;
	  }
#endif // CHIMES
#endif // METALS



#ifdef SINK_PARTICLES
#if (SINGLE_STAR_SINK_FORMATION & 8)
        P.Sink_Ngb_Flag[i] = 0;
#endif
#ifdef SINGLE_STAR_FB_TIMESTEP_LIMIT
 // start with a large value (> plausible values v_ejecta or v_wind) as a conservative choice when starting up a simulation with an active feedback-emmiting star - this will get updated to a more reasonable value once the particle walks the gravity tree, but need this to ensure the first timestep is stable.
        P.MaxFeedbackVel[i] = 1e4 / UNIT_VEL_IN_KMS;
#endif
#ifdef SINGLE_STAR_TIMESTEPPING
	    P.Min_Sink_Approach_Time[i] = P.Min_Sink_Freefall_time[i] = MAX_REAL_NUMBER;
#if (SINGLE_STAR_TIMESTEPPING > 0)
	    P.SuperTimestepFlag[i] = 0;
#endif
#endif
        if(P.Type[i] == 5)
        {
            count_holes++;
            if(RestartFlag == 0)
            {
                P.Sink_Mass[i] = All.SeedSinkMass;
                P.Sink_Formation_Mass[i] = P.Mass[i];
#ifdef SINK_RIAF_SUBEDDINGTON_MODEL
                P.Sink_Mdot_ROI[i] = 0;
                P.Sink_ROI[i] = 0;
#endif
#ifdef SINGLE_STAR_SINK_DYNAMICS
                P.Sink_Mass[i] = P.Mass[i];
#endif
#ifdef SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION // properly initialize luminosity
                singlestar_subgrid_protostellar_evolution_update_track(i,0,0);
#if (SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION == 2)
                P.ZAMS_Mass[i] = P.Sink_Mass[i];
                calculate_individual_stellar_luminosity(P.Sink_Mdot[i], P.Sink_Mass[i], i);
#endif                
#endif
#ifdef GRAIN_FLUID
                P.Sink_Dust_Mass[i] = 0;
#endif
#ifdef SINK_GRAVCAPTURE_FIXEDSINKRADIUS
                P.SinkRadius[i] = KERNEL_FAC_FROM_FORCESOFT_TO_PLUMMER * SinkParticle_GravityKernelRadius;
#endif
#ifdef SINK_ALPHADISK_ACCRETION
                P.Sink_Mass_Reservoir[i] = All.SeedReservoirMass;
#endif
#ifdef SINK_FOLLOW_ACCRETED_ANGMOM
                double sink_mu=2*get_random_number(P.ID[i]+3)-1, sink_phi=2*M_PI*get_random_number(P.ID[i]+4), sink_sin=sqrt(1-sink_mu*sink_mu);
                double spin_prefac = All.G * P.Sink_Mass[i] / C_LIGHT_CODE; // assume initially maximally-spinning black hole with random orientation
                P.Sink_Specific_AngMom[i][0]=spin_prefac*sink_sin*cos(sink_phi); P.Sink_Specific_AngMom[i][1]=spin_prefac*sink_sin*sin(sink_phi); P.Sink_Specific_AngMom[i][2]=spin_prefac*sink_mu;
#endif
#ifdef SINK_COUNTPROGS
                P.Sink_CountProgs[i] = 1;
#endif
            }
#ifdef SINK_INTERACT_ON_GAS_TIMESTEP
            P.dt_since_last_gas_search[i] = 0;
            P.do_gas_search_this_timestep[i] = 1;
#endif 
#if defined(SINK_SWALLOWGAS) && !defined(SINK_GRAVCAPTURE_GAS)
            if(RestartFlag != 1) {P.Sink_AccretionDeficit[i] = 0;}
#endif
        }
#endif
    }

#ifdef SINK_PARTICLES
    MPI_Allreduce(&count_holes, &All.TotSinks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

    for(i = 0; i < TIMEBINS; i++) {TimeBinActive[i] = 1;}

    reconstruct_timebins();

#ifdef PMGRID
    All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

    for(i = 0; i < N_gas; i++)	/* initialize gas/fluid cell properties */
    {
        CellP.InternalEnergyPred[i] = CellP.InternalEnergy[i];

        CellP.VelPred[i] = P.Vel[i];
        CellP.HydroAccel[i] = {};

        P.Particle_DivVel[i] = 0;
        CellP.ConditionNumber[i] = 1;
        CellP.DtInternalEnergy[i] = 0;
        CellP.FaceClosureError[i] = 0;
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        CellP.MaxKineticEnergyNgb[i] = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        CellP.dMass[i] = 0;
        CellP.DtMass[i] = 0;
        CellP.MassTrue[i] = P.Mass[i];
        CellP.GravWorkTerm[i] = {};
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
        P.AGS_zeta[i] = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        if(1 & ADAPTIVE_GRAVSOFT_FORALL) {P.AGS_KernelRadius[i] = P.KernelRadius[i];} else {P.AGS_KernelRadius[i] = All.ForceSoftening[P.Type[i]];}
#endif
#endif

#ifdef CONDUCTION
        CellP.Kappa_Conduction[i] = 0;
#endif
#ifdef MHD_NON_IDEAL
        CellP.Eta_MHD_OhmicResistivity_Coeff[i] = 0;
        CellP.Eta_MHD_HallEffect_Coeff[i] = 0;
        CellP.Eta_MHD_AmbiPolarDiffusion_Coeff[i] = 0;
#endif
#ifdef VISCOSITY
        CellP.Eta_ShearViscosity[i] = 0;
        CellP.Zeta_BulkViscosity[i] = 0;
#endif


#ifdef TURB_DIFFUSION
        CellP.TD_DiffCoeff[i] = 0;

#ifdef TURB_DIFF_DYNAMIC
        int u, v; /* start with the standard Smagorinsky-Lilly constant from Kolmogorov theory */
        CellP.TD_DynDiffCoeff[i] = 0.01;
        CellP.h_turb[i] = 0;
        CellP.FilterWidth_bar[i] = 0;
        CellP.MagShear_bar[i] = 0;
        CellP.MagShear[i] = 0;
        CellP.Norm_hat[i] = 0;
        CellP.Dynamic_numerator[i] = 0;
        CellP.Dynamic_denominator[i] = 0;
#ifdef OUTPUT_TURB_DIFF_DYNAMIC_ERROR
        CellP.TD_DynDiffCoeff_error[i] = 0;
#endif
        if (RestartFlag != 7) {
            CellP.Velocity_bar[i] = {};
            CellP.Velocity_hat[i] = {};
        }
        for (u = 0; u < 3; u++) {
            for (v = 0; v < 3; v++) {
                CellP.VelShear_bar[i][u][v] = 0;
            }
        }
#endif
#endif

        if(RestartFlag == 0)
        {
#ifndef INPUT_READ_KERNELRADIUS
            P.KernelRadius[i] = 0;
#endif
            CellP.Density[i] = -1;
#ifdef COOLING
#ifndef CHIMES
            CellP.Ne[i] = 1.0;
#endif
#if defined(COOL_MOLECFRAC_NONEQM)
            CellP.MolecularMassFraction[i] = 0.0; CellP.MolecularMassFraction_perNeutralH[i] = 0.0; // start atomic
#endif
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
            CellP.Rad_Flux_UV[i] = 0;
            CellP.Rad_Flux_EUV[i] = 0;
#endif
#ifdef CHIMES_STELLAR_FLUXES
	    int kc; for (kc = 0; kc < CHIMES_LOCAL_UV_NBINS; kc++) {CellP.Chimes_fluxPhotIon[i][kc] = 0; CellP.Chimes_G0[i][kc] = 0;}
#endif
#ifdef SINK_COMPTON_HEATING
            CellP.Rad_Flux_AGN[i] = 0;
#endif
        }
#ifdef GALSF_SUBGRID_WINDS
        if(RestartFlag == 0) {CellP.DelayTime[i] = 0;}
#if (GALSF_SUBGRID_WIND_SCALING==1)
        CellP.HostHaloMass[i] = 0;
#endif
#endif
#if defined(GALSF_FB_FIRE_RT_HIIHEATING)
        CellP.DelayTimeHII[i] = 0;
#endif
#ifdef GALSF_FB_TURNOFF_COOLING
        CellP.DelayTimeCoolingSNe[i] = 0;
#endif
#ifdef GALSF
        CellP.Sfr[i] = 0;
#if defined(GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED)
        CellP.AlphaVirial_SF_TimeSmoothed[i] = 0;
#endif
#endif
#ifdef COSMIC_RAY_FLUID
        if(RestartFlag == 0) {for(j=0;j<N_CR_PARTICLE_BINS;j++) {CellP.CosmicRayEnergy[i][j] = 0;}}
#if defined(CRFLUID_INJECTION_AT_SHOCKS)
        if(RestartFlag != 1) {CellP.DtCREgyNewInjectionFromShocks[i] = 0;}
#endif
#if defined(SINK_CR_INJECTION_AT_TERMINATION)
        if(RestartFlag != 1) {CellP.Sink_CR_Energy_Available_For_Injection[i] = 0;}
#endif
#if defined(CRFLUID_EVOLVE_SPECTRUM)
        //if(RestartFlag == 0) {for(j=0;j<N_CR_PARTICLE_BINS;j++) {CellP.CosmicRay_PwrLaw_Slopes_in_Bin[i][j] = -2.5; CellP.CosmicRay_Number_in_Bin[i][j] = 0; CellP.DtCosmicRay_Number_in_Bin[i][j] = 0;}} // initialize a flat spectrum in each bin
        if(RestartFlag == 0) {for(j=0;j<N_CR_PARTICLE_BINS;j++) {CellP.CosmicRay_Number_in_Bin[i][j] = 0; CellP.DtCosmicRay_Number_in_Bin[i][j] = 0;}} // initialize the number in each bin
        if(RestartFlag == 2) { // if we don't directly evolve the slopes, we do write them out and read them in, so need to re-construct the correct number in bin from what we actually read, which was the -slope- information
#if defined(CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART)
            double e0 = 0.5 * P.Mass[i] * CellP.InternalEnergy[i]; // snapshot from which we read does not have the full CR info, so we need to initialize it from something //
#if (CRFLUID_ALT_SPECTRUM_SPECIALSNAPRESTART==1)
            e0 = CellP.CosmicRayEnergy[i][0]; // we had one value of energy (total) available to read in //
#endif
            // now define the desired spectrum //
            for(j=0;j<N_CR_PARTICLE_BINS;j++) {
                int species = CR_species_ID_in_bin[j];
                double f_norm = 1.e-20;
                double f_elec = 0.05; // fraction of the energy to put into e- as opposed to p+ at injection [early experiments with 'observed'  fraction ~ 1% give lower e-/p+ actually observed in the end, so tentative favoring closer to equal at injection? but not run to z=0, so U_rad high from CMB; still experimenting here]
                if(species == -1) {f_norm = f_elec;} // e-
                if(species == +1) {f_norm = 1.-f_elec;} // p
                if(species == -2) {f_norm = 1.e-10 * f_elec;} // e+ (initialize to negligible since want to start with primary)
                if(species > 1)
                {
                    double Zfac = P.Metallicity[i][0]/All.SolarAbundances[0]; // scale heavier elements to the metallicity of the gas into which CRs are being accelerated
                    Zfac *= pow(return_CRbin_CRmass_in_mp(-1,j) / fabs(return_CRbin_CR_charge_in_e(-1,j)) , 1.5); // approximate injection factor for a constant-beta distribution at a given R_GV needed below
                    if(species == 2) {f_norm = 3.7e-9 * Zfac;} // B (for standard elements initialize to solar ratios assuming similar energy/nucleon)
                    if(species == 3) {f_norm = 2.4e-3 * Zfac;} // C
                    if(species == 4) {f_norm = 1.4e-10 * Zfac;} // Be7+9 (stable)
                    if(species == 5) {f_norm = 1.4e-20 * Zfac;} // Be10 (radioactive)
                    if(species == 6) {f_norm = 0.0094 * Zfac;} // CNO (combined bin)
                }
                double e_tmp = 1.0 * f_norm * e0, x_RGV = CR_global_rigidity_at_bin_center[j], a_0=0.715197, x_0=1.7, gamma_0=-2.3, gamma_1=0.8, fac, slope; // hadrons
                if(species < 0) {a_0=0.506309; x_0=1.0; gamma_0=-0.6; gamma_1=1.3;} // leptons
                double dlnR = log(CR_global_max_rigidity_in_bin[j]/CR_global_min_rigidity_in_bin[j]); // bin-width needed below for assignment
                double qx0=pow(x_RGV/x_0,gamma_0), qx1=pow(x_RGV/x_0,gamma_1); fac=a_0*dlnR/(qx0+qx1); slope=-(gamma_0*qx0 + gamma_1*qx1)/(qx0+qx1)-2.; // adopt an extremely simple two-power law spectrum, identical in E space for everything, except normalization, to initialize
                //double fac = 2.3 / (pow(x_RGV,-0.6) + pow(x_RGV,0.8)), slope = (3. - 4.*pow(x_RGV,1.4)) / (5. + 5.*pow(x_RGV,1.4)) - 2.;
                if(CR_check_if_bin_is_nonrelativistic(j)) {slope -= 1.;} // correct for NR terms
                CellP.CosmicRayEnergy[i][j] = e_tmp * fac; CellP.CosmicRay_Number_in_Bin[i][j] = slope; // actually assign the energy and power-law slope
            }
#endif
            /* now initialize the number in each bin from the slopes that we have either read in or assumed */
            for(j=0;j<N_CR_PARTICLE_BINS;j++) {
                double slope_from_snapshot = CellP.CosmicRay_Number_in_Bin[i][j];
                CellP.CosmicRay_Number_in_Bin[i][j] = CR_get_number_in_bin_from_slope(i,j,CellP.CosmicRayEnergy[i][j],slope_from_snapshot);
            }
        }
#endif
#endif
#ifdef MAGNETIC
#if defined MHD_B_SET_IN_PARAMS
        if(RestartFlag == 0)
        {			/* Set only when starting from ICs */
            CellP.B[i][0]=CellP.BPred[i][0] = All.BiniX;
            CellP.B[i][1]=CellP.BPred[i][1] = All.BiniY;
            CellP.B[i][2]=CellP.BPred[i][2] = All.BiniZ;
        }
#endif /*MHD_B_SET_IN_PARAMS*/
        CellP.BPred[i] *= a2_fac * gauss2gizmo;
        CellP.B[i] = CellP.BPred[i];
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
        CellP.Balpha[i] = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
        CellP.Phi[i] = CellP.PhiPred[i] = CellP.DtPhi[i] = 0;
#endif
#ifdef SINK_RETURN_BFLUX
        P.B[i][0] = P.B[i][1] = P.B[i][2] = 0;
#endif
#endif
#ifdef SPHAV_CD10_VISCOSITY_SWITCH
        CellP.alpha[i] = 0.0;
#endif
#if defined(SINK_THERMALFEEDBACK)
        CellP.Injected_Sink_Energy[i] = 0;
#endif
    }

#ifndef BOX_SHEARING
#if (NUMDIMS==2)
    for(i = 0; i < NumPart; i++)
    {
        P.Pos[i][2] = 0;
        //P.Vel[i][2] = 0; // this should be set in the ICs, not here //

        P.GravAccel[i][2] = 0;

        if(P.Type[i] == 0)
        {
            CellP.VelPred[i][2] = 0;
            CellP.HydroAccel[i][2] = 0;
        }
    }
#endif
#endif

#if (NUMDIMS==1)
    for(i = 0; i < NumPart; i++)
    {
        P.Pos[i][1] = P.Pos[i][2] = 0;
        //P.Vel[i][1] = P.Vel[i][2] = 0; // this should be set in the ICs, not here //

        P.GravAccel[i][1] = P.GravAccel[i][2] = 0;

        if(P.Type[i] == 0)
        {
            CellP.VelPred[i][1] = CellP.VelPred[i][2] = 0;
            CellP.HydroAccel[i][1] = CellP.HydroAccel[i][2] = 0;
        }
    }
#endif

#ifdef ASSIGN_NEW_IDS
    assign_unique_ids();
#endif
    /* assign other ID parameters needed */

    if(RestartFlag==0) {for(i = 0; i < NumPart; i++) {P.ID_child_number[i] = 0; P.ID_generation[i] = 0;}}
#ifdef NO_CHILD_IDS_IN_ICS
    if(RestartFlag != 1) {for(i = 0; i < NumPart; i++) {P.ID_child_number[i] = 0; P.ID_generation[i] = 0;}}
#endif

#ifdef TEST_FOR_IDUNIQUENESS
    test_id_uniqueness();
#endif

    Flag_FullStep = 1;		/* to ensure that Peano-Hilbert order is done */
    TreeReconstructFlag = 1;

#ifdef SINK_WIND_SPAWN
    Max_Unspawned_MassUnits_fromSink = 0;
#endif

#ifdef SHIFT_BY_HALF_BOX
    for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
        {
            double boxtmp = 0;
            if(j==0) {boxtmp = boxSize_X;}
            if(j==1) {boxtmp = boxSize_Y;}
            if(j==2) {boxtmp = boxSize_Z;}
            P.Pos[i][j] += 0.5 * boxtmp;
        }
#endif


    Gas_split = 0;
#ifdef GALSF
    Stars_converted = 0;
#endif
    domain_Decomposition(0, 0, 0);	/* do initial domain decomposition (gives equal numbers of particles) */

    set_softenings();

    /* will build tree */
    ngb_treebuild();

    All.Ti_Current = 0;

    if(RestartFlag != 3 && RestartFlag != 5) {setup_smoothinglengths();}

#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
    if(RestartFlag != 3 && RestartFlag != 5) {ags_setup_smoothinglengths();}
#endif
#ifdef CBE_INTEGRATOR
    do_cbe_initialization();
#endif

#ifdef GALSF_SUBGRID_WINDS
#if (GALSF_SUBGRID_WIND_SCALING==2)
    if(RestartFlag != 3 && RestartFlag != 5) {disp_setup_smoothinglengths();}
#endif
#endif

#if defined GALSF_SFR_IMF_VARIATION
    for(i = 0; i < NumPart; i++) {P.IMF_Mturnover[i] = 2.0;} // reset to normal IMF
#endif

#if defined(WAKEUP) && defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
    for(i=0;i<NumPart;i++) {P.wakeup[i]=0;}
#endif


    /* HELLO! This here is where you should insert custom code for hard-wiring the ICs of various test problems */



    density();
    for(i = 0; i < N_gas; i++)	/* initialize gas/fluid cell properties */
    {
        int k __attribute__((unused)) = 0;
        CellP.InternalEnergyPred[i] = CellP.InternalEnergy[i];
        CellP.recent_refinement_flag[i] = 0; // always initialize to zero, no recent refinements

        // re-match the predicted and initial velocities and B-field values, just to be sure //
        CellP.VelPred[i] = P.Vel[i];
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) && (HYDRO_FIX_MESH_MOTION==0)
        CellP.ParticleVel[i] = {}; // set these to zero and forget them, for the rest of the run //
#endif

#ifdef MAGNETIC
        CellP.B[i] = CellP.BPred[i] * (P.Mass[i] / CellP.Density[i]); // convert to the conserved unit V*B //
        CellP.BPred[i] = CellP.B[i]; CellP.DtB[i] = {}; CellP.BField_prerefinement[i] = {};
#endif
#ifdef COSMIC_RAY_FLUID
        for(k=0;k<N_CR_PARTICLE_BINS;k++)
        {
            CellP.CosmicRayEnergyPred[i][k]=CellP.CosmicRayEnergy[i][k]; CellP.CosmicRayDiffusionCoeff[i][k]=0; CellP.DtCosmicRayEnergy[i][k]=0;
            CellP.CosmicRayFlux[i][k] = {}; CellP.CosmicRayFluxPred[i][k] = {};
#ifdef CRFLUID_EVOLVE_SCATTERINGWAVES
            for(j=0;j<2;j++) {CellP.CosmicRayAlfvenEnergy[i][k][j]=0; CellP.CosmicRayAlfvenEnergyPred[i][k][j]=0; CellP.DtCosmicRayAlfvenEnergy[i][k][j]=0;}
#endif
        }
#endif
#if defined(EOS_ELASTIC)
        if(RestartFlag != 1)
        {
            for(k=0;k<3;k++) {for(j=0;j<3;j++) {CellP.Dt_Elastic_Stress_Tensor[i][j][k] = CellP.Elastic_Stress_Tensor_Pred[i][j][k] = CellP.Elastic_Stress_Tensor[i][j][k] = 0;}}
        } else {
            for(k=0;k<3;k++) {for(j=0;j<3;j++) {CellP.Elastic_Stress_Tensor_Pred[i][j][k] = CellP.Elastic_Stress_Tensor[i][j][k]; CellP.Dt_Elastic_Stress_Tensor[i][j][k] = 0;}}
        }
#endif
        CellP.DtInternalEnergy[i] = 0;
#if defined(COOLING) && !defined(COOLING_OPERATOR_SPLIT)
        CellP.CoolingIsOperatorSplitThisTimestep[i] = 1; /* default to more conservative split */
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        CellP.dMass[i] = 0;
        CellP.DtMass[i] = 0;
        CellP.MassTrue[i] = P.Mass[i];
        CellP.GravWorkTerm[i] = {};
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE)
        P.AGS_zeta[i] = 0;
#endif
#ifdef WAKEUP
        if(RestartFlag!=0) {P.wakeup[i]=0;}
        NeedToWakeupParticles = 0;
        NeedToWakeupParticles_local = 0;
#endif
#ifdef SUPER_TIMESTEP_DIFFUSION
        CellP.Super_Timestep_Dt_Explicit[i] = 0;
        CellP.Super_Timestep_j[i] = 0;
#endif
#ifdef GALSF_FB_FIRE_RT_LONGRANGE
        CellP.Rad_Flux_UV[i] = 0;
        CellP.Rad_Flux_EUV[i] = 0;
#endif
#ifdef SINK_COMPTON_HEATING
        CellP.Rad_Flux_AGN[i] = 0;
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP.Rad_E_gamma[i][kf]=0;}}
#endif
#if defined(RT_USE_GRAVTREE_SAVE_RAD_FLUX)
        {int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {CellP.Rad_Flux[i][kf] = {};}}
#endif
#if defined(COSMIC_RAY_SUBGRID_LEBRON)
        CellP.SubGrid_CosmicRayEnergyDensity[i] = 0;
#endif
        
#if (SINGLE_STAR_AND_SSP_NUCLEAR_ZOOM_SPECIALBOUNDARIES>=2)
        if(RestartFlag != 1) {if(P.ID[i] == All.SpawnedWindCellID) {P.ID[i] += 1;}} // reset any of these so can obey desired merge-split rules
#endif

#ifdef COOL_GRACKLE
        if(RestartFlag == 0)
        {
#if (COOL_GRACKLE_CHEMISTRY >= 1)
            CellP.grHI[i]    = HYDROGEN_MASSFRAC;
            CellP.grHII[i]   = 1.0e-20;
            CellP.grHM[i]    = 1.0e-20;
            CellP.grHeI[i]   = 1.0 - HYDROGEN_MASSFRAC;
            CellP.grHeII[i]  = 1.0e-20;
            CellP.grHeIII[i] = 1.0e-20;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 2)
            CellP.grH2I[i]   = 1.0e-20;
            CellP.grH2II[i]  = 1.0e-20;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 3)
            CellP.grDI[i]    = 2.0 * 3.4e-5;
            CellP.grDII[i]   = 1.0e-20;
            CellP.grHDI[i]   = 1.0e-20;
#endif
        }
#endif

    }


    /* we should define the maximum and minimum particle masses
        below/above which particles are merged/split */
    if(RestartFlag != 1)
    {
        double mass_min = MAX_REAL_NUMBER;
        double mass_max = -MAX_REAL_NUMBER;
        double mass_tot = 0;
        for(i = 0; i < N_gas; i++)	/* initialize gas/fluid cell properties */
        {
            mass_tot += P.Mass[i];
            if(P.Mass[i] > mass_max) mass_max = P.Mass[i];
            if(P.Mass[i] < mass_min) mass_min = P.Mass[i];
        }
        /* broadcast this and get the min and max values over all processors */
        double mpi_mass_min, mpi_mass_max;
        MPI_Allreduce(&mass_min, &mpi_mass_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&mass_max, &mpi_mass_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double mpi_splitmerge_readmin, mpi_splitmerge_readmax; /* check if this has been initialized by broadcasting to all processors */
        MPI_Allreduce(&All.MinMassForParticleMerger, &mpi_splitmerge_readmin, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&All.MaxMassForParticleSplit, &mpi_splitmerge_readmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if(mpi_splitmerge_readmin <= 0) { /* initialize if this isn't saved in the ICs */
            All.MinMassForParticleMerger = 0.49 * mpi_mass_min;
#ifdef GALSF_GENERATIONS
            All.MinMassForParticleMerger /= (float)GALSF_GENERATIONS;
#endif
        } else {All.MinMassForParticleMerger = mpi_splitmerge_readmin;} /* use the version from the ICs */
        if(mpi_splitmerge_readmax <= 0) {All.MaxMassForParticleSplit  = 3.01 * mpi_mass_max;} else {All.MaxMassForParticleSplit = mpi_splitmerge_readmax;}
#ifdef MERGESPLIT_HARDCODE_MAX_MASS
        All.MaxMassForParticleSplit = MERGESPLIT_HARDCODE_MAX_MASS;
#endif
#ifdef MERGESPLIT_HARDCODE_MIN_MASS
        All.MinMassForParticleMerger = MERGESPLIT_HARDCODE_MIN_MASS;
#endif

#ifdef SINGLE_STAR_SINK_DYNAMICS /* Get mean gas mass, used in various subroutines */
        double mpi_mass_tot; long mpi_Ngas; long Ngas_l = (long) N_gas;
        MPI_Allreduce(&mass_tot, &mpi_mass_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Ngas_l, &mpi_Ngas, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        All.MeanGasParticleMass = mpi_mass_tot/( (double)mpi_Ngas );
        if(RestartFlag==0){
            for(i=0; i<NumPart; i++){
                if(P.Type[i]==5){P.Sink_Formation_Mass[i] = All.MeanGasParticleMass;} // will behave as if this sink formed from a gas cell with the average mass
            }
        }
#endif
    }

    if(RestartFlag == 3)
    {
#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
        if(ThisTask == 0) {printf("*AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE* Computation of softening lengths... \n");}
        ags_setup_smoothinglengths();
        if(ThisTask == 0) {printf("*AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE* Computation of softening lengths done. \n");}
#endif

#ifdef FOF
        fof_fof(RestartSnapNum);
#endif
        endrun(0);
    }

#ifdef OUTPUT_TWOPOINT_ENABLED
    if(RestartFlag == 5)
    {
        /* calculating powerspec and twopoint function */
#ifdef PMGRID
        long_range_init_regionsize();
#ifdef BOX_PERIODIC
        /* determine global and local particle numbers */
        int n, n_type[6]; long long ntot_type_all[6];
        for(n = 0; n < 6; n++) {n_type[n] = 0;}
        for(n = 0; n < NumPart; n++) {n_type[P.Type[n]]++;}
        sumup_large_ints(6, n_type, ntot_type_all);
        calculate_power_spectra(RestartSnapNum, ntot_type_all);
#endif
#endif
        force_treebuild(NumPart, NULL);
        twopoint();
        endrun(0);
    }
#endif


    if(RestartFlag == 4)
    {
        All.Time = All.TimeBegin = header.time;
        snprintf(All.SnapshotFileBase, 100, "%s_converted", All.SnapshotFileBase);
        if(ThisTask == 0) {printf("Start writing file %s\n", All.SnapshotFileBase);}
        printf("RestartSnapNum %d\n", RestartSnapNum);

        All.TopNodeAllocFactor = 0.008;

        savepositions(RestartSnapNum);
        endrun(0);
    }
    

#if defined(COOL_MOLECFRAC_NONEQM)
    if(RestartFlag == 2) // should have read in CellP.MolecularMassFraction_perNeutralH[i]
    {
        CellP.MolecularMassFraction_perNeutralH[i] = DMIN(1,DMAX(0,CellP.MolecularMassFraction_perNeutralH[i]));
        CellP.MolecularMassFraction[i] = DMIN(1,DMAX(0, 1.-CellP.Ne[i]/1.25)) * CellP.MolecularMassFraction_perNeutralH[i];
    }
#endif
    

#ifdef CHIMES_INITIALISE_IN_EQM
    if (RestartFlag != 1)
      {
	/* Note that stellar fluxes computed through the
	 * gravity tree are all zero at this stage,
	 * because the gravitational forces have not yet
	 * been computed. So the equilibrium abundances
	 * computed here include only the extragalactic UVB. */
	if (ThisTask == 0)
	  printf("Computing equilibrium CHIMES abundances. \n");

	int iter_number;

#ifdef _OPENMP
	int ThisThread;

#pragma omp parallel private(i, iter_number, ThisThread)
	{
	  ThisThread = omp_get_thread_num();

#pragma omp for schedule(dynamic)
#endif
	  for(i = 0; i < N_gas; i++)
	    {
	      initialise_gas_abundances(&(ChimesGasVars[i]), &ChimesGlobalVars);

#ifdef CHIMES_TURB_DIFF_IONS
	      chimes_update_turbulent_abundances(i, 1);
#endif

	      chimes_update_gas_vars(i);

	      // Evolve the chemistry for (1 / nH) Myr (limited to 1 Gyr) ten times at fixed temperature.
	      ChimesGasVars[i].hydro_timestep = (ChimesFloat) DMIN(3.16e13 / ChimesGasVars[i].nH_tot, 3.16e16);
	      ChimesGasVars[i].ThermEvolOn = 0;

	      for (iter_number = 0; iter_number < 10; iter_number++) chimes_network(&(ChimesGasVars[i]), &ChimesGlobalVars);


#ifdef CHIMES_TURB_DIFF_IONS
	      chimes_update_turbulent_abundances(i, 1);
#endif
	    }
#ifdef _OPENMP
	} // End of parallel block
#endif
      } // RestartFlag != 1
#endif // CHIMES_INITIALISE_IN_EQM
}



/*! This routine computes the mass content of the box and compares it to the specified value of Omega-matter.  If discrepant, the run is terminated. */
#ifdef BOX_PERIODIC
void check_omega(void)
{
    double mass = 0, masstot, omega; int i;
    for(i = 0; i < NumPart; i++) {mass += P.Mass[i];}
    MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    omega = masstot / (boxSize_X*boxSize_Y*boxSize_Z) / (3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G));
#ifdef GR_TABULATED_COSMOLOGY_G
    omega *= All.Gini / All.G;
#endif
    if(fabs(omega - All.OmegaMatter) > 1.0e-2) // look for a 1% tolerance of omega-matter
        {PRINT_WARNING("\n\nMass content in the ICs accounts only for Omega_M=%g,\nbut you specified Omega_M=%g in the parameterfile.\nRun will stop.\n",omega, All.OmegaMatter); endrun(1);}
}
#endif


/*! This function is used to find an initial kernel length (what used to be called the
 *  'smoothing length' for SPH, but is just the kernel size for the mesh-free methods) for each gas
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the kernel length is provided to the function density(), which will
 *  then iterate if needed to find the right kernel length.
 */
void setup_smoothinglengths(void)
{
    int i, no, p;
    if((RestartFlag == 0)||(RestartFlag==2)) // best for stability if we re-calc KernelRadius for snapshot restarts //
    {
#if defined(DO_DENSITY_AROUND_NONGAS_PARTICLES) || defined(GRAIN_FLUID)
        for(i = 0; i < NumPart; i++)
#else
        for(i = 0; i < N_gas; i++)
#endif
        {
                no = Father[i];
                while(2 * All.DesNumNgb * P.Mass[i] > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0) {break;}
                    no = p;
                }

                if((RestartFlag == 0)||(P.Type[i] != 0)) // if Restartflag==2, use the saved KernelRadius of the gas as initial guess //
                {
#ifndef INPUT_READ_KERNELRADIUS
#if NUMDIMS == 3
                    P.KernelRadius[i] = pow(3.0 / (4 * M_PI) * All.DesNumNgb * P.Mass[i] / Nodes[no].u.d.mass, 0.333333) * Nodes[no].len;
#endif
#if NUMDIMS == 2
                    P.KernelRadius[i] = pow(1.0 / (M_PI) * All.DesNumNgb * P.Mass[i] / Nodes[no].u.d.mass, 0.5) * Nodes[no].len;
#endif
#if NUMDIMS == 1
                    P.KernelRadius[i] = All.DesNumNgb * (P.Mass[i] / Nodes[no].u.d.mass) * Nodes[no].len;
#endif
#ifndef SELFGRAVITY_OFF
                    double soft = All.ForceSoftening[P.Type[i]];
                    if(soft != 0) {if((P.KernelRadius[i]>100.*soft)||(P.KernelRadius[i]<=0.01*soft)||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {P.KernelRadius[i] = soft;}}
#else
                    if((Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {
#if (defined(BOX_PERIODIC) || defined(BOX_SHEARING) || defined(BOX_DEFINED_SPECIAL_XYZ_BOUNDARY_CONDITIONS_ARE_ACTIVE) || defined(BOX_LONG_X) || defined(BOX_LONG_Y) || defined(BOX_LONG_Z))
                        P.KernelRadius[i] = 0.05 * All.BoxSize;
#else
                        P.KernelRadius[i] = 1;
#endif
                    }
#endif
#endif // INPUT_READ_KERNELRADIUS
                } // closes if((RestartFlag == 0)||(P.Type[i] != 0))
            }
    }
    if((RestartFlag==0 || RestartFlag==2) && All.ComovingIntegrationOn) {for(i=0;i<N_gas;i++) {P.KernelRadius[i] *= pow(All.OmegaMatter/All.OmegaBaryon,1./NUMDIMS);}} /* correct (crudely) for baryon fraction, used in the estimate above for KernelRadius */

#ifdef SINK_PARTICLES
    if(RestartFlag==0 || RestartFlag==2) {for(i=0;i<NumPart;i++) {if(P.Type[i] == 5) {P.KernelRadius[i] = All.ForceSoftening[P.Type[i]];}}}
#endif

#ifdef GRAIN_FLUID
    if(RestartFlag==0 || RestartFlag==2) {for(i=0;i<NumPart;i++) {P.KernelRadius[i] *= pow(2.,1./NUMDIMS);}} /* very rough correction assuming comparable numbers of dust and gas elements */
#endif

    density();
}


void assign_unique_ids(void)
{
    int i, *numpartlist;
    MyIDType idfirst;

    numpartlist = (int *) mymalloc("numpartlist", NTask * sizeof(int));

    MPI_Allgather(&NumPart, 1, MPI_INT, numpartlist, 1, MPI_INT, MPI_COMM_WORLD);

    idfirst = 1;

    for(i = 0; i < ThisTask; i++)
        idfirst += numpartlist[i];

    for(i = 0; i < NumPart; i++)
    {
        P.ID[i] = idfirst;
        idfirst++;
    }

    myfree(numpartlist);
}


#ifdef AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE
void ags_setup_smoothinglengths(void)
{
    int i, no, p;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            P.Particle_DivVel[i] = 0;
            P.AGS_zeta[i] = 0;
            if(ags_density_isactive(i) || P.Type[i]==0) // type is AGS-active //
            {
                if(P.Type[i] > 0)
                {
                    no = Father[i];
                    while(10 * All.AGS_DesNumNgb * P.Mass[i] > Nodes[no].u.d.mass)
                    {
                        p = Nodes[no].u.d.father;
                        if(p < 0) break;
                        no = p;
                    }
                    P.AGS_KernelRadius[i] = 2. * pow(1.0/VOLUME_NORM_COEFF_FOR_NDIMS * All.AGS_DesNumNgb * P.Mass[i] / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                    double soft = All.ForceSoftening[P.Type[i]];
                    if(soft != 0)
                    {
                        if((P.AGS_KernelRadius[i]>1e6*soft)||(P.AGS_KernelRadius[i]<=1e-3*soft)||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {P.AGS_KernelRadius[i] = 1.e2*soft;} /* random guess to get things started here, thats all */
                    }
                } else {
                    P.AGS_KernelRadius[i] = P.KernelRadius[i];
                }
            } else {
                P.AGS_KernelRadius[i] = All.ForceSoftening[P.Type[i]]; /* not AGS-active, use fixed softening */
            }
        }
    }
    ags_density();
#ifdef DM_FUZZY
    do_dm_fuzzy_initialization();
#endif
}
#endif // AGS_KERNELRADIUS_CALCULATION_IS_ACTIVE


#if defined(GALSF_SUBGRID_WINDS)
#if (GALSF_SUBGRID_WIND_SCALING==2)
void disp_setup_smoothinglengths(void)
{
    int i, no, p;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            if(P.Type[i] == 0)
            {
                no = Father[i];
                while(10 * 2.0 * 64 * P.Mass[i] > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0) {break;}
                    no = p;
                }
                CellP.KernelRadiusDM[i] = pow(1.0/VOLUME_NORM_COEFF_FOR_NDIMS * 2.0 * 64 * P.Mass[i] / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                double soft = All.ForceSoftening[P.Type[i]];
                if(soft != 0) {if((CellP.KernelRadiusDM[i] >1000.*soft)||(P.KernelRadius[i]<=0.01*soft)||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) {CellP.KernelRadiusDM[i] = soft;}}
            }
        }
    }
    if(ThisTask == 0) {printf("computing DM Vel_disp around gas particles.\n");}
    disp_density();
}
#endif
#endif


void test_id_uniqueness(void)
{
    double t0, t1;
#ifndef BOX_BND_PARTICLES
    int i;
    MyIDType *ids, *ids_first;
#endif

    if(ThisTask == 0)
    {
        printf("Testing ID uniqueness...\n");
    }

    if(NumPart == 0)
    {
        printf("need at least one particle per cpu\n");
        endrun(8);
    }

    t0 = my_second();

#ifndef BOX_BND_PARTICLES
    ids = (MyIDType *) mymalloc("ids", NumPart * sizeof(MyIDType));
    ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));

    for(i = 0; i < NumPart; i++)
        ids[i] = P.ID[i];

    parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);

    for(i = 1; i < NumPart; i++)
        if(ids[i] == ids[i - 1])
        {
            printf("non-unique ID=%llu found on task=%d   (i=%d NumPart=%d)\n", (unsigned long long) ids[i], ThisTask, i, NumPart);
            endrun(12);
        }

    MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

    if(ThisTask < NTask - 1)
        if(ids[NumPart - 1] == ids_first[ThisTask + 1])
        {
            printf("non-unique ID=%llu found on task=%d\n", (unsigned long long) ids[NumPart - 1], ThisTask);
            endrun(13);
        }

    myfree(ids_first);
    myfree(ids);
#endif

    t1 = my_second();

    if(ThisTask == 0)
    {
        printf("success.  took=%g sec\n", timediff(t0, t1));
    }
}

int compare_IDs(const void *a, const void *b)
{
    if(*((MyIDType *) a) < *((MyIDType *) b)) {return -1;}
    if(*((MyIDType *) a) > *((MyIDType *) b)) {return +1;}
    return 0;
}
