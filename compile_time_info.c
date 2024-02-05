#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        BOX_SPATIAL_DIMENSION=3\n"
"        BOX_PERIODIC\n"
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        MAGNETIC\n"
"        GRAVITY_NOT_PERIODIC\n"
"        SINGLE_STAR_STARFORGE_DEFAULTS\n"
"        SINGLE_STAR_STARFORGE_PROTOSTELLAR_EVOLUTION=2\n"
"        SINGLE_STAR_FB_JETS\n"
"        SINGLE_STAR_FB_WINDS=0\n"
"        SINGLE_STAR_FB_SNE\n"
"        SINGLE_STAR_FB_RAD\n"
"        COOLING\n"
"        COSMIC_RAY_FLUID\n"
"        INPUT_COSMIC_RAY_ENERGY\n"
"        CRFLUID_M1=(300000.)\n"
"        CRFLUID_DIFFUSION_MODEL=0\n"
"        RT_SPEEDOFLIGHT_REDUCTION=(1.0e-3)\n"
"        CRFLUID_ION_ALFVEN_SPEED\n"
"        CRFLUID_ALT_RSOL_FORM\n"
"        OPENMP=2\n"
"        ADAPTIVE_TREEFORCE_UPDATE=0.0625\n"
"        SUBCYCLING_TEST\n"
"\n");
}
