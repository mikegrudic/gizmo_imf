/* This is a generic code block designed for simple neighbor loops, so that they don't have to 
    be copy-pasted and can be generically optimized in a single place. specifically this is for
    the initial loop of particles on the local processor (and determination of which need passing) 

   Two blocks need to be defined or this will crash: 
   CONDITION_FOR_EVALUATION inserts the clause that actually determines
        whether or not to pass a particle to the main evaluation routine
   EVALUATION_CALL is the actual call, and needs to be written appropriately
 */
#if !defined(CONDITION_FOR_EVALUATION) || !defined(EVALUATION_CALL)
printf("Cannot compile the primary sub-loop without both CONDITION_FOR_EVALUATION and EVALUATION_CALL defined. Exiting. \n"); fflush(stdout); exit(995533);
#endif
/* variable assignment */
int i, j, *exportflag, *exportnodecount, *exportindex, *ngblist, thread_id = *(int *) p;
/* define the pointers needed for each thread to speak back regarding what needs processing */
ngblist = Ngblist.data() + thread_id * NumPart;
exportflag = Exportflag + thread_id * NTask;
exportnodecount = Exportnodecount + thread_id * NTask;
exportindex = Exportindex + thread_id * NTask;
/* Note: exportflag is local to each thread */
for(j = 0; j < NTask; j++) {exportflag[j] = -1;}
#ifdef _OPENMP
if(BufferCollisionFlag && thread_id) {return NULL;} /* force to serial for this subloop if threads simultaneously cross the Nexport bunchsize threshold */
#endif
/* now begin the actual loop */
#ifdef NGB_SEARCH_BOTH_WAYS
/* ---- Batched primary loop: grab NGB_BATCH_SIZE particles, one batched tree walk, process each ---- */
{
int _bi_indices[NGB_BATCH_SIZE];
int *_bi_ngblists[NGB_BATCH_SIZE];
int _bi_numngb[NGB_BATCH_SIZE];
MyDouble _bi_centers[NGB_BATCH_SIZE][3];
MyFloat _bi_rkerns[NGB_BATCH_SIZE];
int _bi_targets[NGB_BATCH_SIZE];
int _bi_valid[NGB_BATCH_SIZE];
int _bi_partition = NumPart / NGB_BATCH_SIZE;
for(j = 0; j < NGB_BATCH_SIZE; j++) {_bi_ngblists[j] = ngblist + j * _bi_partition;}

while(1)
{
    int _bi_count = 0, _bi_exit = 0;
#ifdef _OPENMP
#pragma omp critical(_nextlistprimblox_)
#endif
    {
        int _bi_attempts = 0, _bi_maxattempts = NGB_BATCH_SIZE * 4;
        while(_bi_count < NGB_BATCH_SIZE && _bi_attempts < _bi_maxattempts)
        {
            if(BufferFullFlag != 0 || NextParticle >= (int)ActiveParticleList.size()) {_bi_exit = 1; break;}
            int idx = ActiveParticleList[NextParticle]; NextParticle++;
            _bi_attempts++;
            if(ProcessedFlag[idx]) {continue;}
            _bi_indices[_bi_count++] = idx;
        }
    }
    if(_bi_count == 0) {break;}

    int _bi_nvalid = 0;
    for(int _bi = 0; _bi < _bi_count; _bi++)
    {
        i = _bi_indices[_bi];
        int _bi_passes = 0;
        CONDITION_FOR_EVALUATION { _bi_passes = 1; }
        _bi_valid[_bi] = _bi_passes;
        if(_bi_passes) {
            _bi_centers[_bi_nvalid][0] = P.Pos[i][0];
            _bi_centers[_bi_nvalid][1] = P.Pos[i][1];
            _bi_centers[_bi_nvalid][2] = P.Pos[i][2];
            _bi_rkerns[_bi_nvalid] = P.KernelRadius[i];
            _bi_targets[_bi_nvalid] = i;
            _bi_nvalid++;
        }
    }

    int _bi_did_batch = 0;
    if(_bi_nvalid >= 2)
    {
        int _bi_ret = ngb_treefind_pairs_threads_batched(_bi_nvalid, _bi_centers, _bi_rkerns,
                        _bi_targets, exportflag, exportnodecount, exportindex,
                        _bi_ngblists, _bi_numngb, NGB_SEARCH_BOTH_WAYS);
        if(_bi_ret == 0) {_bi_did_batch = 1;}
        else {break;}
    }

    int _bi_broken = 0, _bi_vi = 0;
    for(int _bi = 0; _bi < _bi_count && !_bi_broken; _bi++)
    {
        i = _bi_indices[_bi];
        if(ProcessedFlag[i]) {continue;}

        if(_bi_valid[_bi] && _bi_did_batch) {
            ngblist = _bi_ngblists[_bi_vi];
            ngb_set_prefound_list(_bi_ngblists[_bi_vi], _bi_numngb[_bi_vi]);
            _bi_vi++;
        } else if(_bi_valid[_bi]) {
            ngblist = _bi_ngblists[_bi_vi];
            _bi_vi++;
        }

        CONDITION_FOR_EVALUATION
        {
            if(EVALUATION_CALL < 0) {ngb_clear_prefound_list(); _bi_broken = 1; break;}
        }
        ngb_clear_prefound_list();
        ProcessedFlag[i] = 1;
    }
    if(_bi_broken || _bi_exit) {break;}
}
}
#else
/* now begin the actual loop, grabbing batches of particles to reduce critical section overhead */
#ifndef PRIMARY_LOOP_BATCH_SIZE
#define PRIMARY_LOOP_BATCH_SIZE 8
#endif
while(1)
{
    int batch[PRIMARY_LOOP_BATCH_SIZE], batch_count = 0;
#ifdef _OPENMP
#pragma omp critical(_nextlistprimblox_)
#endif
    {
        while(batch_count < PRIMARY_LOOP_BATCH_SIZE && BufferFullFlag == 0 && NextParticle < (int)ActiveParticleList.size())
        {
            int idx = ActiveParticleList[NextParticle];
            NextParticle++;
            if(!ProcessedFlag[idx]) {batch[batch_count++] = idx;}
        }
    }
    if(batch_count == 0) {break;}
    int buffer_full = 0;
    for(int b = 0; b < batch_count; b++)
    {
        i = batch[b];
        CONDITION_FOR_EVALUATION
        {
            if(EVALUATION_CALL < 0) {buffer_full = 1; break;} // export buffer has filled up //
        }
        ProcessedFlag[i] = 1; /* particle successfully finished */
    }
    if(buffer_full) {break;}
}
#endif /* NGB_SEARCH_BOTH_WAYS */
/* loop completed successfully */
return NULL;
