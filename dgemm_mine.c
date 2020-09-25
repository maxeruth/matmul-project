const char* dgemm_desc = "My awesome dgemm.";


void square_dgemm(const int M, const double *A, const double *B, double *C)
{
    int i, j, k;
    for (i = 0; i < M; i+=4) {
        for (j = 0; j < M; ++j) {
			if(M-i >= 4){
				int i1 = i + 1; 
				int i2 = i + 2; 
				int i3 = i + 3; 
				double cij0 = C[j*M+i];
				double cij1 = C[j*M+i1];
				double cij2 = C[j*M+i2];
				double cij3 = C[j*M+i3];
				for (k = 0; k < M; ++k){
					cij0 += A[k*M+i] * B[j*M+k];
					cij1 += A[k*M+i1] * B[j*M+k];
					cij2 += A[k*M+i2] * B[j*M+k];
					cij3 += A[k*M+i3] * B[j*M+k];
				}
				C[j*M+i] = cij0;
				C[j*M+i1] = cij1;
				C[j*M+i2] = cij2;
				C[j*M+i3] = cij3;
			}
			if(M-i == 3){
				int i1 = i + 1; 
				int i2 = i + 2; 
				double cij0 = C[j*M+i];
				double cij1 = C[j*M+i1];
				double cij2 = C[j*M+i2];
				for (k = 0; k < M; ++k){
					cij0 += A[k*M+i] * B[j*M+k];
					cij1 += A[k*M+i1] * B[j*M+k];
					cij2 += A[k*M+i2] * B[j*M+k];
				}
				C[j*M+i] = cij0;
				C[j*M+i1] = cij1;
				C[j*M+i2] = cij2;
			}if(M-i == 2){
				int i1 = i + 1; 
				double cij0 = C[j*M+i];
				double cij1 = C[j*M+i1];
				for (k = 0; k < M; ++k){
					cij0 += A[k*M+i] * B[j*M+k];
					cij1 += A[k*M+i1] * B[j*M+k];
				}
				C[j*M+i] = cij0;
				C[j*M+i1] = cij1;
			}
			if(M-i == 1){
				double cij0 = C[j*M+i];
				for (k = 0; k < M; ++k){
					cij0 += A[k*M+i] * B[j*M+k];
				}
				C[j*M+i] = cij0;
			}
			
        }
    }
}
