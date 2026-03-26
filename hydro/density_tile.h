#ifndef DENSITY_TILE_H
#define DENSITY_TILE_H

#define DENSITY_TILE_NGB_MAX 512

struct DensityNeighborTile {
    int ngb_idx[DENSITY_TILE_NGB_MAX];
    Vec3<MyDouble> pos[DENSITY_TILE_NGB_MAX];
    MyDouble mass[DENSITY_TILE_NGB_MAX];
    Vec3<MyDouble> velpred[DENSITY_TILE_NGB_MAX];
#ifdef HYDRO_PRESSURE_SPH
    MyDouble ie_pred[DENSITY_TILE_NGB_MAX];
#endif
#ifdef GALSF_SUBGRID_WINDS
    MyFloat delay_time[DENSITY_TILE_NGB_MAX];
#endif
    int count;
    int use_tile; /* 0 if numngb exceeded tile size, fall back to direct SoA */
};

static inline void density_tile_gather(DensityNeighborTile& tile, int* ngblist, int numngb) {
    tile.count = numngb;
    if(numngb > DENSITY_TILE_NGB_MAX) {tile.use_tile = 0; return;} /* fallback to direct SoA */
    tile.use_tile = 1;
    for(int nn = 0; nn < numngb; nn++) {
        int j = ngblist[nn];
        tile.ngb_idx[nn] = j;
        tile.pos[nn] = P.Pos[j];
        tile.mass[nn] = P.Mass[j];
        tile.velpred[nn] = CellP.VelPred[j];
#ifdef HYDRO_PRESSURE_SPH
        tile.ie_pred[nn] = CellP.InternalEnergyPred[j];
#endif
#ifdef GALSF_SUBGRID_WINDS
        tile.delay_time[nn] = CellP.DelayTime[j];
#endif
    }
}

/* Macros for tile access with fallback to direct SoA */
#define DTILE_POS(n,j) (dtile.use_tile ? dtile.pos[n] : P.Pos[j])
#define DTILE_MASS(n,j) (dtile.use_tile ? dtile.mass[n] : P.Mass[j])
#define DTILE_VELPRED(n,j) (dtile.use_tile ? dtile.velpred[n] : CellP.VelPred[j])

#endif
