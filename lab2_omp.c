#include <malloc.h>
#include <omp.h>
#include <math.h>


// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */

int N,M;

void QR_decomposition(float* A,float* Q,float* R)
{
    float U[N][N];
    float U_square[N];
    float sqrt_U_square[N];

    for (int i = 0; i <N ; ++i) {
        U[i][0]=A[i][0];
    }

    for (int i = 1; i <N ; ++i) {

        for (int j = 0; j <N ; ++j) {
            U[j][i]=A[j][i];
        }

        for (int j = 0; j <i ; ++j) {
            float proj[N];
            float a=0;

            if (j==i-1){
                U_square[j]=0;
                for (int k = 0; k <N ; ++k) {
                    U_square[j]+=U[k][j]*U[k][j];
                }
            }

            for (int k = 0; k <N ; ++k) {
                a+=U[k][j]*A[k][i];
            }

            a=a/U_square[j];

            for (int k = 0; k <N ; ++k) {
                U[k][i]-=a*U[k][j];
            }

        }
    }

    U_square[N-1]=0;
    for (int k = 0; k <N ; ++k) {
        U_square[N-1]+=U[k][N-1]*U[k][N-1];
    }

    for (int i = 0; i <N ; ++i) {
        sqrt_U_square[i]=sqrt(U_square[i]);
    }

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            Q[j][i]=U[j][i]/sqrt_U_square[i];
        }
    }


    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            if (i>j) R[i][j]=0;
            else{
                float a=0;
                for (int k = 0; k <N ; ++k) {
                    a+=E[k][i]*A[k][j];
                }
                R[i][j]=a;
            }
        }
    }

}


void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
    ::M=M;
    ::N=N;

    float D_new[M][N];
    float U_new[M][M];
    float sigma[M][N];
    float V_new[N][N];
    float V_Tnew[N][N];

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            D_new[i][j]=D[N*i+j];
        }
    }

    float D_T[N][M];

    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <N ; ++j) {
            D_T[j][i]=D_new[i][j];
        }
    }

    float D_T_D[N][N];

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            float a=0;
            for (int k = 0; k <M ; ++k) {
                a+=D_T[i][k]*D_new[k][j];
            }
            D_T_D[i][j]=a;
        }
    }

    // Now we have to find the eigen values of D_T_D by using QR method.

    float D_not[N][N];
    float E[N][N];

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            D_not[i][j]=D_T_D[i][j];
            if (i==j) E[i][j]=1;
            else E[i][j]=0;
        }
    }

    // Until convergence-
    for (int i = 0; i <1000 ; ++i) {
        float Q[N][N];
        float R[N][N];

        QR_decomposition(&D_not,&Q,&R);

        for (int i = 0; i <N ; ++i) {
            for (int j = 0; j <N ; ++j) {
                float a=0;
                for (int k = 0; k <N; ++k) {
                    a+=R[i][k]*Q[k][j];
                }
                D_not[i][j]=a;
            }
        }

        for (int i = 0; i <N ; ++i) {
            for (int j = 0; j <N ; ++j) {
                float a=0;
                for (int k = 0; k <N; ++k) {
                    a+=E[i][k]*Q[k][j];
                }
                E[i][j]=a;
            }
        }

    }

    float eigen_values[N];

    for (int i = 0; i <N ; ++i) {
        eigen_values[i]=abs(D_not[i][i]);
    }

    float sqrt_eigen_values[N];

    for (int i = 0; i <N ; ++i) {
        sqrt_eigen_values[i]=sqrt(eigen_values[i][i]);
    }

    float sorted_singular_values[N];

    // Sorting the eigen values and computing the matrix V
    for (int i = 0; i <N ; ++i) {
        int max=0;
        int index=0;
        for (int j = 0; j <N ; ++j) {
            if (sqrt_eigen_values[j]>max){
                max=sqrt_eigen_values[j];
                index=j;
            }
        }
        sorted_singular_values[i]=max;
        for (int j = 0; j <N ; ++j) {
            V_new[j][i]=E[j][index];
        }
        sqrt_eigen_values[index]=0;
    }

    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <N ; ++j) {
            sigma[i][j]=0;
        }
    }

    // Computing matrix sigma
    for (int i = 0; i <N ; ++i) {
        sigma[i][i]=sorted_singular_values[i];
    }


    //Computing V transpose
    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            V_Tnew[j][i]=V_new[i][j];
        }
    }




}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
}
