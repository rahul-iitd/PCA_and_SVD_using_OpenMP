#include <malloc.h>
#include <omp.h>
#include <math.h>


int N,M;


void QR_decomposition(float** D_not,float** Q,float** R)
{
    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M; ++j) {
            Q[i][j]=0;
            R[i][j]=0;
        }
    }
    float U[M][M];
    float U_square[M];
    float sqrt_U_square[M];

    for (int i = 0; i <M ; ++i) {
        U[i][0]=D_not[i][0];
    }

    for (int i = 1; i <M ; ++i) {

        for (int j = 0; j <M ; ++j) {
            U[j][i]=D_not[j][i];
        }

        for (int j = 0; j <i ; ++j) {
            float proj[M];
            float a=0;

            if (j==i-1){
                U_square[j]=0;
                for (int k = 0; k <M ; ++k) {
                    U_square[j]+=U[k][j]*U[k][j];
                }
            }

            for (int k = 0; k <M ; ++k) {
                a+=U[k][j]*D_not[k][i];
            }

            a=a/U_square[j];

            for (int k = 0; k <M ; ++k) {
                U[k][i]-=a*U[k][j];
            }

        }
    }

    U_square[M-1]=0;
    for (int k = 0; k <M ; ++k) {
        U_square[M-1]+=U[k][M-1]*U[k][M-1];
    }

    for (int i = 0; i <M ; ++i) {
        sqrt_U_square[i]=sqrt(U_square[i]);
    }

    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
            Q[j][i]=U[j][i]/sqrt_U_square[i];
        }
    }


    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
            if (i>j) R[i][j]=0;
            else{
                float a=0;
                for (int k = 0; k <M ; ++k) {
                    a+=Q[k][i]*D_not[k][j];
                }
                R[i][j]=a;
            }
        }
    }

}


void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
    // We have to compute the SVD of D_Transpose.
    ::M=M;
    ::N=N;

    float D_new[M][N];
    float D_new_T[N][M];
    float U_new[N][N];
    float sigma[N][M];
    float V_new[M][M];
    float V_Tnew[M][M];

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            D_new[i][j]=D[N*i+j];
        }
    }

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <M ; ++j) {
            D_new_T[i][j]=D_new[j][i];
        }
    }

    //Now D is D transpose.


    float D_T_D[M][M];

    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
            float a=0;
            for (int k = 0; k <N ; ++k) {
                a+=D_new[i][k]*D_new_T[k][j];
            }
            D_T_D[i][j]=a;
        }
    }

    // Now we have to find the eigen values of D_T_D by using QR method.


    float E[M][M];
    float**  D_not;
    D_not = (float**) malloc(sizeof(float *) * M);
    for(int i = 0;i<M;i++){
        D_not[i]=(float*) malloc(sizeof(float) * M);
    }


    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
            D_not[i][j]=D_T_D[i][j];
            if (i==j) E[i][j]=1;
            else E[i][j]=0;
        }
    }

    // Until convergence-
    for(int i=0;i<1000;i++) {
//        printf("%d\n",i);
        i++;
        float**  Q;
        Q = (float**) malloc(sizeof(float *) * M);
        for(int i = 0;i<M;i++){
            Q[i]=(float*) malloc(sizeof(float) * M);
        }

        float**  R;
        R = (float**) malloc(sizeof(float *) * M);
        for(int i = 0;i<M;i++){
            R[i]=(float*) malloc(sizeof(float) * M);
        }

        QR_decomposition(D_not,Q,R);

        float D_not_new[M][M];

        for (int i = 0; i <M ; ++i) {
            for (int j = 0; j <M ; ++j) {
                float a=0;
                for (int k = 0; k <M; ++k) {
                    a+=R[i][k]*Q[k][j];
                }
                D_not[i][j]=a;
            }
        }

//        int count=0;
//        for (int i = 0; i <M ; ++i) {
//            for (int j = 0; j <M ; ++j) {
//                if (abs(D_not_new[i][j]-D_not[i][j])<=0.001) count++;
//                else break;
//            }
//        }
//        printf("%d\n",count);
//        if (count==M*M) break;

//        for (int i = 0; i <M ; ++i) {
//            for (int j = 0; j <M ; ++j) {
//                D_not[i][j]==D_not_new[i][j];
//            }
//        }

        float E_new[M][M];

        for (int i = 0; i <M ; ++i) {
            for (int j = 0; j <M; ++j) {
                float a=0;
                for (int k = 0; k <M; ++k) {
                    a+=E[i][k]*Q[k][j];
                }
                E_new[i][j]=a;
            }
        }

        for (int i = 0; i <M ; ++i) {
            for (int j = 0; j <M ; ++j) {
                E[i][j]=E_new[i][j];
            }
        }

    }

    float eigen_values[M];

    for (int i = 0; i <M ; ++i) {
        eigen_values[i]=abs(D_not[i][i]);
    }

    //printing eigen values to check
    for (int i = 0; i <M ; ++i) {
    }

    float sqrt_eigen_values[M];

    for (int i = 0; i <M ; ++i) {
        sqrt_eigen_values[i]=sqrt(eigen_values[i]);
    }

    float sorted_singular_values[M];    // taking only first N eigen values in decending order.

    // Sorting the eigen values and computing the matrix V
    for (int i = 0; i <M ; ++i) {
        float max=0;
        int index=0;
        for (int j = 0; j <M ; ++j) {
            if (sqrt_eigen_values[j]>max){
                max=sqrt_eigen_values[j];
                index=j;
            }
        }
        sorted_singular_values[i]=max;
        for (int j = 0; j <M ; ++j) {
            V_new[j][i]=E[j][index];
        }
        sqrt_eigen_values[index]=0;
    }

    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
        }
    }

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <M ; ++j) {
            sigma[i][j]=0;
        }
    }

    // Computing matrix sigma
    for (int i = 0; i <N ; ++i) {
        sigma[i][i]=sorted_singular_values[i];

    }


    //Computing V transpose
    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
            V_Tnew[j][i]=V_new[i][j];
        }
    }


    // Computing sigma inverse

    float sigma_inv[M][N];

    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <N ; ++j) {
            if (i!=j) sigma_inv[i][j]==0;
            else sigma_inv[i][j]=sorted_singular_values[i];
        }
    }

    // Computing U = D_new_T V Sigma_inv

    float intermediate[N][M];

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <M ; ++j) {
            float a=0;
            for (int k = 0; k <M ; ++k) {
                a+=D_new_T[i][k]*V_new[k][j];
            }
            intermediate[i][j]=a;
        }
    }

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            float a=0;
            for (int k = 0; k <M ; ++k) {
                a+=intermediate[i][k]*sigma_inv[k][j];
            }
            U_new[i][j]=a;
        }
    }

    for (int i = 0; i <N ; ++i) {
        SIGMA[0][i]=sigma[i][i];
        printf("%f\n",SIGMA[0][i]);
    }

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            U[0][N*i+j]=U_new[i][j];
        }
    }

    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
            V_T[0][M*i+j]=V_Tnew[i][j];
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
    float sum_eigen_values = 0;

    for (int i = 0; i <N ; ++i) {
        sum_eigen_values+=SIGMA[i];
    }

    int count=0;
    float sum=0;
    for (int i = 0; i <N ; ++i) {
        count+=1;
        sum+=SIGMA[i];
        if ((sum/sum_eigen_values)*100>=retention) break;
    }

    *K=count;

    float D_hat[M][count];

    float D_new[M][N];
    float U_new[N][N];

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            D_new[i][j]=D[N*i+j];
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            U_new[i][j]=U[N*i+j];
        }
    }

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j <count ; ++j) {
            float a=0;
            for (int k = 0; k <N ; ++k) {
                a+=D_new[i][k]*U_new[k][j];
            }
            D_hat[i][j]=a;
        }
    }

    D_HAT[0] = (float*) malloc(sizeof(float) * M*count);

    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <count ; ++j) {
            D_HAT[0][M*i+j]=D_hat[i][j];
        }
    }

}
