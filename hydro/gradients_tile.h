#ifndef GRADIENTS_TILE_H
#define GRADIENTS_TILE_H

#define GRADIENTS_TILE_NGB_MAX 512

struct GradientsNeighborTile {
    int ngb_idx[GRADIENTS_TILE_NGB_MAX];
    Vec3<MyDouble> pos[GRADIENTS_TILE_NGB_MAX];
    MyDouble mass[GRADIENTS_TILE_NGB_MAX];
    MyFloat h[GRADIENTS_TILE_NGB_MAX];              // P.KernelRadius[j]
    MyDouble density[GRADIENTS_TILE_NGB_MAX];        // CellP.Density[j]
    MyDouble pressure[GRADIENTS_TILE_NGB_MAX];       // CellP.Pressure[j]
    Vec3<MyDouble> velpred[GRADIENTS_TILE_NGB_MAX];  // CellP.VelPred[j]
    MyDouble ie_pred[GRADIENTS_TILE_NGB_MAX];        // CellP.InternalEnergyPred[j]
    MyDouble condition_number[GRADIENTS_TILE_NGB_MAX]; // CellP.ConditionNumber[j]
#ifdef GALSF_SUBGRID_WINDS
    MyFloat delay_time[GRADIENTS_TILE_NGB_MAX];
#endif
#ifdef TURB_DIFF_DYNAMIC
    MyDouble norm_hat[GRADIENTS_TILE_NGB_MAX];       // CellP.Norm_hat[j]
    MyDouble velocity_bar[GRADIENTS_TILE_NGB_MAX][3]; // CellP.Velocity_bar[j]
#endif
    int count;
    int use_tile; /* 0 if numngb exceeded tile size, fall back to direct SoA */
};

static inline void gradients_tile_gather(GradientsNeighborTile& tile, int* ngblist, int numngb) {
    tile.count = numngb;
    if(numngb > GRADIENTS_TILE_NGB_MAX) {tile.use_tile = 0; return;} /* fallback to direct SoA */
    tile.use_tile = 1;
    for(int nn = 0; nn < numngb; nn++) {
        int j = ngblist[nn];
        tile.ngb_idx[nn] = j;
        tile.pos[nn] = P.Pos[j];
        tile.mass[nn] = P.Mass[j];
        tile.h[nn] = P.KernelRadius[j];
        tile.density[nn] = CellP.Density[j];
        tile.pressure[nn] = CellP.Pressure[j];
        tile.velpred[nn] = CellP.VelPred[j];
        tile.ie_pred[nn] = CellP.InternalEnergyPred[j];
        tile.condition_number[nn] = CellP.ConditionNumber[j];
#ifdef GALSF_SUBGRID_WINDS
        tile.delay_time[nn] = CellP.DelayTime[j];
#endif
#ifdef TURB_DIFF_DYNAMIC
        tile.norm_hat[nn] = CellP.Norm_hat[j];
        for(int kk = 0; kk < 3; kk++) {tile.velocity_bar[nn][kk] = CellP.Velocity_bar[j][kk];}
#endif
    }
}

#endif // GRADIENTS_TILE_H
