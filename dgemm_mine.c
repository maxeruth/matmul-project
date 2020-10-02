#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdalign.h>
#include <math.h>
#include <omp.h>

#define MBLOCK 16 
typedef struct {
    alignas(32) double vec[MBLOCK * MBLOCK];
} fvec_t;

const char* dgemm_desc = "Just the restrict keyword.";

void copy_block_back(const int M, double* restrict A_OUT, double* restrict A_IN)
{
    int i, j;
    for (i = 0; i < MBLOCK; ++i)
        for (j = 0; j < MBLOCK; ++j)
            A_OUT[i * M + j] = A_IN[i * MBLOCK + j];
}

void copy_back(const int M, double* restrict A_OUT, fvec_t* restrict A_IN)
{
    int i, j;
    int N = M % MBLOCK == 0 ? M / MBLOCK : M / MBLOCK + 1;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            copy_block(M, A_OUT + MBLOCK * M * i + MBLOCK * j, A_IN[i * N + j].vec);
        }
    }
}

void copy_block(const int M, const double* restrict A_IN, double* restrict A_OUT)
{
    int i, j;
    for (i = 0; i < MBLOCK; ++i)
        for (j = 0; j < MBLOCK; ++j)
            A_OUT[i * MBLOCK + j] = A_IN[i * M + j];
}

void copy(const int M, const double* restrict A_IN, fvec_t* restrict A_OUT)
{
    int i, j;
    int N = M % MBLOCK == 0 ? M / MBLOCK : M / MBLOCK + 1;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            copy_block(M, A_IN + MBLOCK * M * i + MBLOCK * j, A_OUT[i * N + j].vec);
        }
    }
}

void multiply_block(double* restrict C, double* restrict A, double* restrict B)
{
    int j, k, i;
    for (j = 0; j < MBLOCK; ++j) {			
        for (k = 0; k < MBLOCK; ++k){
			for (i = 0; i < MBLOCK; ++i) {
                C[j*MBLOCK+i  ] += A[k*MBLOCK+i  ] * B[j*MBLOCK+k];
			}
        }
    }
}

void square_dgemm(const int M, 
                  const double* restrict A, const double* restrict B, double* restrict C)
{
    int i, j, k;
    int N = M % MBLOCK == 0 ? M / MBLOCK : M / MBLOCK + 1;
    fvec_t* A_BLOCK = (fvec_t*) aligned_alloc(32, N * N * sizeof(fvec_t));
    copy(M, A, A_BLOCK);

    fvec_t* B_BLOCK = (fvec_t*) aligned_alloc(32, N * N * sizeof(fvec_t));
    copy(M, B, B_BLOCK);

    fvec_t* C_BLOCK = (fvec_t*) aligned_alloc(32, N * N * sizeof(fvec_t));
    copy(M, C, C_BLOCK);

    for (j = 0; j < N; ++j) {
        for (k = 0; k < N; ++k) {
			for (i = 0; i < N; ++i) {
                multiply_block(C_BLOCK[j*N+i].vec, A_BLOCK[k*N+i].vec, B_BLOCK[j*N+k].vec);
			}
        }
    }

    copy_back(M, C, C_BLOCK);
}
