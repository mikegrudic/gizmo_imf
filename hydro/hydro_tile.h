#ifndef HYDRO_TILE_H
#define HYDRO_TILE_H

#define HYDRO_TILE_NGB_MAX 512

/* Tile arrays for neighbor data - gathered before the compute loop for cache efficiency.
   Instead of random SoA access like P.Mass[j] where j = ngblist[n], we gather all needed
   neighbor fields into contiguous arrays and access tile.mass[n] sequentially. */
struct HydroNeighborTile {
    int ngb_idx[HYDRO_TILE_NGB_MAX];          // original particle index j
    Vec3<MyDouble> pos[HYDRO_TILE_NGB_MAX];    // P.Pos[j]
    MyDouble mass[HYDRO_TILE_NGB_MAX];          // P.Mass[j]
    MyFloat h[HYDRO_TILE_NGB_MAX];              // P.KernelRadius[j]
    MyFloat drk_ngb_factor[HYDRO_TILE_NGB_MAX]; // P.DrkernNgbFactor[j]
    short int timebin[HYDRO_TILE_NGB_MAX];      // P.TimeBin[j]
    MyDouble density[HYDRO_TILE_NGB_MAX];       // CellP.Density[j]
    Vec3<MyDouble> velpred[HYDRO_TILE_NGB_MAX]; // CellP.VelPred[j]
    MyDouble pressure[HYDRO_TILE_NGB_MAX];      // CellP.Pressure[j]
    MyDouble ie_pred[HYDRO_TILE_NGB_MAX];       // CellP.InternalEnergyPred[j]
    MyDouble face_closure_err[HYDRO_TILE_NGB_MAX]; // CellP.FaceClosureError[j]
    SymmetricTensor2<MyDouble> nv_t[HYDRO_TILE_NGB_MAX]; // CellP.NV_T[j]
    MyDouble condition_number[HYDRO_TILE_NGB_MAX]; // CellP.ConditionNumber[j]
    Vec3<MyDouble> grad_density[HYDRO_TILE_NGB_MAX]; // CellP.Gradients.Density[j]
    Vec3<MyDouble> grad_pressure[HYDRO_TILE_NGB_MAX]; // CellP.Gradients.Pressure[j]
    Mat3<MyDouble> grad_velocity[HYDRO_TILE_NGB_MAX]; // CellP.Gradients.Velocity[j]
    MyFloat max_signal_vel[HYDRO_TILE_NGB_MAX]; // CellP.MaxSignalVel[j]
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
    MyDouble max_ke_ngb[HYDRO_TILE_NGB_MAX];   // CellP.MaxKineticEnergyNgb[j]
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
    Vec3<MyDouble> particle_vel[HYDRO_TILE_NGB_MAX]; // CellP.ParticleVel[j]
#endif
#ifdef MAGNETIC
    Vec3<MyDouble> bpred[HYDRO_TILE_NGB_MAX];  // CellP.BPred[j]
#ifdef DIVBCLEANING_DEDNER
    MyDouble phi_pred[HYDRO_TILE_NGB_MAX];      // CellP.PhiPred[j]
#endif
    MyFloat balpha[HYDRO_TILE_NGB_MAX];         // CellP.Balpha[j]
    Vec3<MyDouble> grad_b[HYDRO_TILE_NGB_MAX][3]; // CellP.Gradients.B[j]
#endif
#ifdef METALS
    MyFloat metallicity[HYDRO_TILE_NGB_MAX][NUM_METAL_SPECIES]; // P.Metallicity[j]
#endif
#ifdef COSMIC_RAY_FLUID
    MyFloat cosmic_ray_energy[HYDRO_TILE_NGB_MAX][N_CR_PARTICLE_BINS];
    MyFloat cosmic_ray_diffusion_coeff[HYDRO_TILE_NGB_MAX][N_CR_PARTICLE_BINS];
    Vec3<MyFloat> cosmic_ray_flux[HYDRO_TILE_NGB_MAX][N_CR_PARTICLE_BINS];
#endif
#ifdef GALSF_SUBGRID_WINDS
    MyFloat delay_time[HYDRO_TILE_NGB_MAX];
#endif
#ifdef HYDRO_SPH
    MyDouble drk_hydro_sum[HYDRO_TILE_NGB_MAX]; // CellP.DrkernHydroSumFactor[j]
#ifdef HYDRO_PRESSURE_SPH
    MyDouble egy_wt_density[HYDRO_TILE_NGB_MAX]; // CellP.EgyWtDensity[j]
#endif
#endif
    int count; // number of valid entries
    int use_tile; // 0 if numngb exceeded tile size
};

/* Gather neighbor data into tile */
static inline void hydro_tile_gather(HydroNeighborTile& tile, int* ngblist, int numngb) {
    tile.count = numngb;
    if(numngb > HYDRO_TILE_NGB_MAX) {tile.use_tile = 0; return;} /* fallback to direct SoA */
    tile.use_tile = 1;
    for(int n = 0; n < numngb; n++) {
        int j = ngblist[n];
        tile.ngb_idx[n] = j;
        tile.pos[n] = P.Pos[j];
        tile.mass[n] = P.Mass[j];
        tile.h[n] = P.KernelRadius[j];
        tile.drk_ngb_factor[n] = P.DrkernNgbFactor[j];
        tile.timebin[n] = P.TimeBin[j];
        tile.density[n] = CellP.Density[j];
        tile.velpred[n] = CellP.VelPred[j];
        tile.pressure[n] = CellP.Pressure[j];
        tile.ie_pred[n] = CellP.InternalEnergyPred[j];
        tile.face_closure_err[n] = CellP.FaceClosureError[j];
        tile.nv_t[n] = CellP.NV_T[j];
        tile.condition_number[n] = CellP.ConditionNumber[j];
        tile.grad_density[n] = CellP.Gradients.Density[j];
        tile.grad_pressure[n] = CellP.Gradients.Pressure[j];
        tile.grad_velocity[n] = CellP.Gradients.Velocity[j];
        tile.max_signal_vel[n] = CellP.MaxSignalVel[j];
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        tile.max_ke_ngb[n] = CellP.MaxKineticEnergyNgb[j];
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        tile.particle_vel[n] = CellP.ParticleVel[j];
#endif
#ifdef MAGNETIC
        for(int kk=0;kk<3;kk++) tile.bpred[n][kk] = Get_Gas_BField(j,kk);
#ifdef DIVBCLEANING_DEDNER
        tile.phi_pred[n] = Get_Gas_PhiField(j);
#endif
#endif
#ifdef METALS
        for(int kk=0;kk<NUM_METAL_SPECIES;kk++) tile.metallicity[n][kk] = P.Metallicity[j][kk];
#endif
#ifdef GALSF_SUBGRID_WINDS
        tile.delay_time[n] = CellP.DelayTime[j];
#endif
#ifdef HYDRO_SPH
        tile.drk_hydro_sum[n] = CellP.DrkernHydroSumFactor[j];
#ifdef HYDRO_PRESSURE_SPH
        tile.egy_wt_density[n] = CellP.EgyWtDensity[j];
#endif
#endif
    }
}

#endif // HYDRO_TILE_H
