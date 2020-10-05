const char* dgemm_desc = "Just the restrict keyword.";

void square_dgemm(const int M, 
                  const double* restrict A, const double* restrict B, double* restrict C)
{
    int i, j, k;
	for (k = 0; k < M; ++k){
		for (i = 0; i < M; ++i) {		
			for (j = 0; j < M; ++j) {	
                C[j*M+i  ] += A[k*M+i  ] * B[j*M+k];
			}
        }
    }
}
