#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"

void allocate_memory(void)
{
  size_t bytes;
  double bytes_tot = 0;
  int NTaskTimesThreads;
  NTaskTimesThreads = maxThreads * NTask;

  Exportflag = (int *) mymalloc("Exportflag", NTaskTimesThreads * sizeof(int));
  Exportindex = (int *) mymalloc("Exportindex", NTaskTimesThreads * sizeof(int));
  Exportnodecount = (int *) mymalloc("Exportnodecount", NTaskTimesThreads * sizeof(int));

  Send_count = (int *) mymalloc("Send_count", sizeof(int) * NTask);
  Send_offset = (int *) mymalloc("Send_offset", sizeof(int) * NTask);
  Recv_count = (int *) mymalloc("Recv_count", sizeof(int) * NTask);
  Recv_offset = (int *) mymalloc("Recv_offset", sizeof(int) * NTask);

  ProcessedFlag = (unsigned char *) mymalloc("ProcessedFlag", bytes = All.MaxPart * sizeof(unsigned char));
  bytes_tot += bytes;

  ActiveParticleList.reserve(All.MaxPart);
  NextInTimeBin.resize(All.MaxPart);
  PrevInTimeBin.resize(All.MaxPart);

  if(All.MaxPart > 0)
    {
      allocate_P(P, All.MaxPart, "P");
      if(ThisTask == 0) {printf("Allocated SoA particle data storage for %d particles.\n", All.MaxPart);}
    }

  if(All.MaxPartGas > 0)
    {
      allocate_CellP(CellP, All.MaxPartGas, "CellP");
      if(ThisTask == 0) {printf("Allocated SoA gas cell data storage for %d gas particles.\n", All.MaxPartGas);}

#ifdef CHIMES
      if (!(ChimesGasVars = (struct gasVariables *) mymalloc("gasVars", bytes = All.MaxPartGas * sizeof(struct gasVariables))))
	{
	  printf("failed to allocate memory for 'ChimesGasVars' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;
      if(ThisTask == 0) printf("Allocated %g MByte for storage of ChimesGasVars data.\n", bytes_tot / (1024.0 * 1024.0));
#endif
    }
}
