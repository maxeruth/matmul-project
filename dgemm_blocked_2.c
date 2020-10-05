const char* dgemm_desc = "Simple blocked dgemm.";

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 2)
#endif

/*
  A is M-by-K
  B is K-by-N
  C is M-by-N

  lda is the leading dimension of the matrix (the M of square_dgemm).
*/
void square_dgemm(const int M, const double* restrict A, const double* restrict B, double* restrict C)
{
    const int n_blocks = M / BLOCK_SIZE + (M%BLOCK_SIZE? 1 : 0);
    int kk, ii, j, k, i, kmin, imin;
    for (kk = 0; kk < n_blocks; ++kk) {
        for (ii = 0; ii < n_blocks; ++ii) {
            for (j = 0; j < M; ++j) {
		    for (k = kk*BLOCK_SIZE; k <(kk+1)*BLOCK_SIZE; ++k) {
			    const double r=B[j*M+k];
                	for (i=ii*BLOCK_SIZE; i< (ii+1)*BLOCK_SIZE;++i){
			       C[j*M+i]+=A[k*M+i]*r;	
            }
        }
    }
}
}
}

