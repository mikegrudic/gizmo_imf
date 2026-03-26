#ifndef GRAVTREE_TILE_H
#define GRAVTREE_TILE_H

/* Tile-based gather optimization for gravity tree walk.
 *
 * The gravity tree walk accesses particle SoA fields with random indices during
 * tree traversal. This header provides two optimizations:
 *
 * 1. GravExportTile: batched gather for filling GravDataIn export buffers.
 *    Instead of random SoA reads one-at-a-time in the export loop, we gather
 *    the core fields for a batch of particles into contiguous arrays first,
 *    then fill the export structs from the tile. This improves cache utilization
 *    when the export list is large (common in distributed runs).
 *
 * 2. GravLeafTile: accumulates leaf particle data encountered during the tree
 *    walk into a small contiguous buffer. When the buffer fills or a walk
 *    segment ends, the buffered interactions are flushed in a tight loop.
 *    This separates the random SoA gather from the force computation, letting
 *    the compiler vectorize the Newtonian force accumulation.
 */

#define GRAV_LEAF_TILE_SIZE 32 /* small enough to stay in L1 cache */
#define GRAV_EXPORT_TILE_SIZE 64

/* ---- Leaf particle interaction tile for the tree walk ---- */
struct GravLeafTile {
    double dx[GRAV_LEAF_TILE_SIZE];
    double dy[GRAV_LEAF_TILE_SIZE];
    double dz[GRAV_LEAF_TILE_SIZE];
    double mass[GRAV_LEAF_TILE_SIZE];
    double h_p[GRAV_LEAF_TILE_SIZE]; /* secondary softening */
    int ptype_sec[GRAV_LEAF_TILE_SIZE];
    int particle_index[GRAV_LEAF_TILE_SIZE]; /* original particle index, for optional physics that need it */
    double r2[GRAV_LEAF_TILE_SIZE]; /* pre-computed r^2 */
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    double zeta_sec[GRAV_LEAF_TILE_SIZE];
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
    double dvx[GRAV_LEAF_TILE_SIZE];
    double dvy[GRAV_LEAF_TILE_SIZE];
    double dvz[GRAV_LEAF_TILE_SIZE];
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    double r_source[GRAV_LEAF_TILE_SIZE];
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    double gasmass[GRAV_LEAF_TILE_SIZE];
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    SymmetricTensor2<MyFloat> tidal_tensorps_prevstep[GRAV_LEAF_TILE_SIZE];
#endif
    int count; /* number of valid entries */
};

static inline void grav_leaf_tile_reset(GravLeafTile& tile) {
    tile.count = 0;
}


/* ---- Export tile for batched GravDataIn filling ---- */
struct GravExportTile {
    int index[GRAV_EXPORT_TILE_SIZE]; /* particle indices (place) */
    Vec3<MyDouble> pos[GRAV_EXPORT_TILE_SIZE];
    MyDouble mass[GRAV_EXPORT_TILE_SIZE];
    int type[GRAV_EXPORT_TILE_SIZE];
    MyFloat oldacc[GRAV_EXPORT_TILE_SIZE];
    double soft[GRAV_EXPORT_TILE_SIZE];
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
    Vec3<MyDouble> vel[GRAV_EXPORT_TILE_SIZE];
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    double ags_zeta[GRAV_EXPORT_TILE_SIZE];
#endif
    int count;
};

/* Gather particle fields into export tile for a batch of exports.
 * Call this before filling GravDataIn to pre-load SoA data into contiguous memory. */
static inline void grav_export_tile_gather(GravExportTile& tile, struct data_index* DataIndexTable_ptr, int start, int count) {
    if(count > GRAV_EXPORT_TILE_SIZE) {count = GRAV_EXPORT_TILE_SIZE;}
    tile.count = count;
    for(int n = 0; n < count; n++) {
        int place = DataIndexTable_ptr[start + n].Index;
        tile.index[n] = place;
        tile.pos[n] = P.Pos[place];
        tile.mass[n] = P.Mass[place];
        tile.type[n] = P.Type[place];
        tile.oldacc[n] = P.OldAcc[place];
        tile.soft[n] = ForceSoftening_KernelRadius(place);
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
        tile.vel[n] = P.Vel[place];
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS)
        if((tile.type[n] == 0) && (P.KernelRadius[place] > All.ForceSoftening[tile.type[n]])) {tile.ags_zeta[n] = P.AGS_zeta[place];} else {tile.ags_zeta[n] = 0;}
#elif defined(ADAPTIVE_GRAVSOFT_FORALL)
        tile.ags_zeta[n] = P.AGS_zeta[place];
#endif
    }
}

/* Fill GravDataIn from a pre-gathered tile (avoids redundant SoA reads) */
static inline void grav_export_tile_fill(GravExportTile& tile, struct gravdata_in* GravDataIn_ptr, struct data_index* DataIndexTable_ptr, struct data_nodelist* DataNodeList_ptr, int start) {
    for(int n = 0; n < tile.count; n++) {
        int j = start + n;
        GravDataIn_ptr[j].Pos = tile.pos[n];
        GravDataIn_ptr[j].Type = tile.type[n];
        GravDataIn_ptr[j].Soft = tile.soft[n];
        GravDataIn_ptr[j].OldAcc = tile.oldacc[n];
        GravDataIn_ptr[j].Mass = tile.mass[n];
#if defined(SINGLE_STAR_TIMESTEPPING) || defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
        GravDataIn_ptr[j].Vel = tile.vel[n];
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        GravDataIn_ptr[j].AGS_zeta = tile.ags_zeta[n];
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        GravDataIn_ptr[j].Soft = P.AGS_KernelRadius[tile.index[n]]; /* override with adaptive value */
#endif
        /* Fields that require the original particle index and are not worth pre-caching
           (rarely accessed or behind ifdefs with complex logic) */
        int place = tile.index[n];
#if defined(SINK_DYNFRICTION_FROMTREE)
        if(tile.type[n]==5) {GravDataIn_ptr[j].Sink_Mass = P.Sink_Mass[place];}
#endif
#ifdef SINGLE_STAR_FIND_BINARIES
        if(tile.type[n] == 5)
        {
            GravDataIn_ptr[j].Min_Sink_OrbitalTime = P.Min_Sink_OrbitalTime[place];
            GravDataIn_ptr[j].comp_Mass = P.comp_Mass[place];
            GravDataIn_ptr[j].is_in_a_binary = P.is_in_a_binary[place];
            GravDataIn_ptr[j].comp_dx = P.comp_dx[place]; GravDataIn_ptr[j].comp_dv = P.comp_dv[place];
        }
        else {P.is_in_a_binary[place]=0;}
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
        GravDataIn_ptr[j].tidal_tensorps_prevstep=P.tidal_tensorps_prevstep[place];
#endif
        memcpy(GravDataIn_ptr[j].NodeList,DataNodeList_ptr[DataIndexTable_ptr[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
    }
}


/* ---- Leaf tile: gather a single leaf particle into the tile ---- */
static inline void grav_leaf_tile_add(GravLeafTile& tile, int no, double dx_in, double dy_in, double dz_in,
                                       double mass_in, double r2_in, double h_p_in, int ptype_sec_in
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                                       , double zeta_sec_in
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
                                       , double dvx_in, double dvy_in, double dvz_in
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
                                       , double r_source_in
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
                                       , double gasmass_in
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
                                       , SymmetricTensor2<MyFloat> tidal_tp
#endif
                                       ) {
    int c = tile.count;
    tile.dx[c] = dx_in;
    tile.dy[c] = dy_in;
    tile.dz[c] = dz_in;
    tile.mass[c] = mass_in;
    tile.r2[c] = r2_in;
    tile.h_p[c] = h_p_in;
    tile.ptype_sec[c] = ptype_sec_in;
    tile.particle_index[c] = no;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    tile.zeta_sec[c] = zeta_sec_in;
#endif
#if defined(COMPUTE_JERK_IN_GRAVTREE) || defined(SINK_DYNFRICTION_FROMTREE)
    tile.dvx[c] = dvx_in;
    tile.dvy[c] = dvy_in;
    tile.dvz[c] = dvz_in;
#endif
#ifdef GRAVITY_SPHERICAL_SYMMETRY
    tile.r_source[c] = r_source_in;
#endif
#ifdef GRAVTREE_CALCULATE_GAS_MASS_IN_NODE
    tile.gasmass[c] = gasmass_in;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FROM_TIDAL_CRITERION
    tile.tidal_tensorps_prevstep[c] = tidal_tp;
#endif
    tile.count = c + 1;
}


#endif /* GRAVTREE_TILE_H */
