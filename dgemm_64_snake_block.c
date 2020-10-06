const char* dgemm_desc = "Just the restrict keyword.";

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 64)
#endif

int min(int M, int N)
{
    return M < N ? M : N;
}

int iseven(int M)
{
	return (M+1)%2;
}

inline void do_block(const int M, const int jj, const int kk,
              const double* restrict A, const double* restrict B, double* restrict C)
{
	for (int j = jj; j < min(jj+BLOCK_SIZE,M); j++) {
		for (int k = kk; k < min(kk + BLOCK_SIZE, M); k++) {
		double r=B[j*M+k];
			for (int i = 0;i<M; i++) {
				C[j*M+i] += A[k*M+i]*r;
			}
		}
	}
}


void square_dgemm(const int M, 
                  const double* restrict A, const double* restrict B, double* restrict C)
{
	int J = M%BLOCK_SIZE==0 ? (M/BLOCK_SIZE-1)*BLOCK_SIZE : (M/BLOCK_SIZE)*BLOCK_SIZE; // Biggest value of jj or kk

    for (int kk = 0; kk <= J; kk += BLOCK_SIZE) {

        if (iseven(kk/BLOCK_SIZE)){
			for (int jj = 0; jj <= J; jj += BLOCK_SIZE) {
				for (int j = jj; j < min(jj+BLOCK_SIZE,M); j++) {
					for (int k = kk; k < min(kk + BLOCK_SIZE, M); k++) {
					double r=B[j*M+k];
						for (int i = 0;i<M; i++) {
							C[j*M+i] += A[k*M+i]*r;
						}
					}
				}
			}
		}
		else{
			for (int jj = J; jj >= 0; jj -= BLOCK_SIZE) {
				for (int j = jj; j < min(jj+BLOCK_SIZE,M); j++) {
					for (int k = kk; k < min(kk + BLOCK_SIZE, M); k++) {
					double r=B[j*M+k];
						for (int i = 0;i<M; i++) {
							C[j*M+i] += A[k*M+i]*r;
						}
					}
				}
			}
		}/**/
    }
}

/*
void square_dgemm(const int M,
                  const double* restrict A, const double* restrict B, double* restrict C)
{
    int i, j, k, kk, jj;
    double r;
    
	
	//printf("\nM = %d, J = %d\n",M,J);
	// Iterate over 
    for (kk = 0; kk < J; kk += BLOCK_SIZE) {
		for (jj = 0; jj <= J; jj += BLOCK_SIZE) {
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
    
}*/



