const char* dgemm_desc = "Just the restrict keyword.";

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 8)
#endif

int min(int M, int N)
{
    return M < N ? M : N;
}

void square_dgemm(const int M, double r,
                  const double* restrict A, const double* restrict B, double* restrict C)
{
    int i, j, k, kk, jj;
    double sum;
    int en = BLOCK_SIZE * ((M + BLOCK_SIZE) /BLOCK_SIZE); /* Amount that fits evenly into blocks  */

    for (kk = 0; kk < en; kk += BLOCK_SIZE) {
        for (jj = 0; jj < en; jj += BLOCK_SIZE) {
            for (j = jj; j < min(jj+BLOCK_SIZE,M); j++) {
                for (k = kk; k < min(kk + BLOCK_SIZE, M); k++) {
		    r=B[j*M+k];
                    for (i = 0;i<M; i++) {
                        C[j*M+i] += A[k*M+i]*r;
                    }
                }
            }
        }
    }
}
