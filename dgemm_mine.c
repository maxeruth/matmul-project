const char* dgemm_desc = "Just the restrict keyword.";

void square_dgemm(const int M, 
                  const double* restrict A, const double* restrict B, double* restrict C)
{
    int i, j, k;
    double* restrict B_t = (double*) malloc(M * M * sizeof(double));
    
    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j) {
            B_t[i*M + j] = B[j*M + i];
        }
    }

    for (j = 0; j < M; ++j) {			
        for (k = 0; k < M; ++k){
			for (i = 0; i < M; ++i) {
                C[j*M+i  ] += A[k*M+i  ] * B_t[k*M+j];
			}
        }
    }
}
