const char* dgemm_desc = "Adding a temporary variable for best.";


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdalign.h>
#include <math.h>
#include <omp.h>

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 32)
#endif

#ifndef ALIGN_SIZE
#define ALIGN_SIZE 32
#endif

int min(int M, int N)
{
    return M < N ? M : N;
}

void copy_block(const int M, int J, int K, const double* restrict B, double* restrict BBLOCK)
{
    int j;
    for (j = 0; j < J; ++j)
		memcpy(BBLOCK + j*BLOCK_SIZE, B + j*M, K*sizeof(double));
}

void multiply_block(int M, int jj, int kk, int J, int K, const double* restrict A, const double* restrict B_BLOCK, double* restrict C)
{
	int j,k,i;

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			const double r=B_BLOCK[j*BLOCK_SIZE+k]; 
			for (i = 0; i < M; i++) {
				C[j*M+i] += A[k*M+i]*r;
			}
		}
	}
}

void square_dgemm(const int M, 
                  const double* restrict A, const double* restrict B, double* restrict C)
{
    int kk, jj;
    double sum;
    //int en = BLOCK_SIZE * ((M + BLOCK_SIZE) /BLOCK_SIZE); /* Amount that fits evenly into blocks  */
	double* B_BLOCK = (double*) aligned_alloc(ALIGN_SIZE, BLOCK_SIZE * BLOCK_SIZE * sizeof(double));
	
    for (kk = 0; kk < M; kk += BLOCK_SIZE) {
        for (jj = 0; jj < M; jj += BLOCK_SIZE) {
			int J = min(BLOCK_SIZE, M - jj);
			int K = min(BLOCK_SIZE, M - kk);
			copy_block(M, J, K, B + jj*M + kk, B_BLOCK);
			multiply_block(M, jj, kk, J, K, A+kk*M, B_BLOCK, C + jj*M);
        }
    }
    
    free(B_BLOCK);
}



