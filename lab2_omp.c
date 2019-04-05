#include <malloc.h>
#include <omp.h>
#include <math.h>


int N,M;


void QR_decomposition(double** D_not,double** Q,double** R)
{
    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M; ++j) {
            Q[i][j]=0;
            R[i][j]=0;
        }
    }
    double U[M][M];
    double U_square[M];
    double sqrt_U_square[M];

    for (int i = 0; i <M ; ++i) {
        U[i][0]=D_not[i][0];
    }

    for (int i = 1; i <M ; ++i) {

        for (int j = 0; j <M ; ++j) {
            U[j][i]=D_not[j][i];
        }

        for (int j = 0; j <i ; ++j) {
            double a=0;

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
                double a=0;
                for (int k = 0; k <M ; ++k) {
                    a+=Q[k][i]*D_not[k][j];
                }
                R[i][j]=a;
            }
        }
    }

}




void QR_decomposition_modified(double** D_not,double** Q,double** R)
{

    double V[N][N];

#pragma omp parallel for num_threads(8)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            V[i][j] = D_not[i][j];
        }
    }


    for (int i = 0; i <N ; ++i) {
        double a=0;
        for (int j = 0; j <N ; ++j) {
            a+=V[j][i]*V[j][i];
        }
        a=sqrt(a);

        R[i][i]=a;

        for (int j = 0; j <N ; ++j) {
            Q[j][i]=V[j][i]/a;
        }

        for (int j = i+1; j <N ; ++j) {
            double b=0;
            for (int k = 0; k <N ; ++k) {
                b+=Q[k][i]*V[k][j];
            }

            R[i][j]=b;

            for (int k = 0; k <N ; ++k) {
                V[k][j]=V[k][j]-b*Q[k][i];
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
    float U_new[M][M];
    float sigma[M][N];
    float V_new[N][N];
//    float V_Tnew[N][N];

#pragma omp parallel for num_threads(8)
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            D_new[i][j]=D[N*i+j];
        }
    }

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <M ; ++j) {
            D_new_T[i][j]=D_new[j][i];
        }
    }

    //Now D is D transpose.


    float D_T_D[N][N];

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            float a=0;
            for (int k = 0; k <M ; ++k) {
                a+=D_new_T[i][k]*D_new[k][j];
            }
            D_T_D[i][j]=a;
        }
    }

    // Now we have to find the eigen values of D_T_D by using QR method.


    double E[N][N];
    double**  D_not;
    double D_not_old[N][N];
    D_not = (double**) malloc(sizeof(double *) * N);
    for(int i = 0;i<N;i++){
        D_not[i]=(double*) malloc(sizeof(double) * N);
    }


    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            D_not[i][j]=D_T_D[i][j];
            if (i==j) E[i][j]=1;
            else E[i][j]=0;
        }
    }
    int c=0;
    // Until convergence-
    while(true){
        c++;
//        if (c>1) break;
//        printf("%d\n",i);
        double**  Q;
        Q = (double**) malloc(sizeof(double *) * N);
        for(int i = 0;i<N;i++){
            Q[i]=(double*) malloc(sizeof(double) * N);
        }

        double**  R;
        R = (double**) malloc(sizeof(double *) * N);
        for(int i = 0;i<N;i++){
            R[i]=(double*) malloc(sizeof(double) * N);
        }

        QR_decomposition_modified(D_not,Q,R);

//        if (c==1){
//            for (int i = 0; i <M ; ++i) {
//                for (int j = 0; j <M ; ++j) {
//                    printf("%f ",Q[i][j]);
//                }
//                printf("\n");
//            }
//            printf("\n");
//            for (int i = 0; i <M ; ++i) {
//                for (int j = 0; j <M ; ++j) {
//                    printf("%f ",R[i][j]);
//                }
//                printf("\n");
//            }
//            printf("\n");
//        }


#pragma omp parallel for num_threads(8)
        for (int i = 0; i <N ; ++i) {
            for (int j = 0; j <N ; ++j) {
                D_not_old[i][j]=D_not[i][j];
//                printf("%f\n",D_not_old[i][j]);
            }
        }

#pragma omp parallel for num_threads(8)
        for (int i = 0; i <N ; ++i) {
            for (int j = 0; j <N ; ++j) {
                double a=0;
                for (int k = 0; k <N; ++k) {
                    a+=R[i][k]*Q[k][j];
                }
                D_not[i][j]=a;
//                if (c%10==0) printf("%f\n",a);
            }
        }

        int count=0;
//        for (int i = 0; i <M ; ++i) {
//            for (int j = 0; j <M ; ++j) {
////                if (c%10==0){
////                    printf("%f\n",D_not_old[1][1]);
////                    printf("%f\n",D_not[1][1]);
////                }
//                if (abs(D_not_old[i][j]-D_not[i][j])<=0.001){
//                    count++;
////                    printf("%f\n",D_not_old[i][j]);
////                    printf("%f\n",D_not[i][j]);
//                }
////                else break;
//            }
//        }
#pragma omp parallel for num_threads(8)
        for (int i = 0; i <N ; ++i) {
            if (abs(D_not_old[i][i]-D_not[i][i])<=0.00001) count++;
        }

//        if (true){
//            for (int i = 0; i <M ; ++i) {
//                for (int j = 0; j <M ; ++j) {
//                    printf("%f ",D_not_old[i][j]);
//                }
//                printf("\n");
//            }
//            printf("\n");
//            for (int i = 0; i <M ; ++i) {
//                for (int j = 0; j <M ; ++j) {
//                    printf("%f ",D_not[i][j]);
//                }
//                printf("\n");
//            }
//            printf("\n");
//        }

//        printf("%d\n",count);
//        printf("%d\n",count2);
//        printf("%d\n",c);
        if (count==N) break;
        if (count>1000) break;

//        for (int i = 0; i <M ; ++i) {
//            for (int j = 0; j <M ; ++j) {
//                D_not[i][j]==D_not_new[i][j];
//            }
//        }

        double E_new[N][N];

#pragma omp parallel for num_threads(8)
        for (int i = 0; i <N ; ++i) {
            for (int j = 0; j <N; ++j) {
                double a=0;
                for (int k = 0; k <N; ++k) {
                    a+=E[i][k]*Q[k][j];
                }
                E_new[i][j]=a;
            }
        }

#pragma omp parallel for num_threads(8)
        for (int i = 0; i <N ; ++i) {
            for (int j = 0; j <N ; ++j) {
                E[i][j]=E_new[i][j];
            }
        }

    }

//    printf("%d\n",c);

    double eigen_values[N];

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <N ; ++i) {
        eigen_values[i]=abs(D_not[i][i]);
//        printf("%f\n",eigen_values[i]);
    }



    double sqrt_eigen_values[N];

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <N ; ++i) {
        sqrt_eigen_values[i]=sqrt(eigen_values[i]);
    }

    float sorted_singular_values[N];    // taking only first N eigen values in decending order.

    // Sorting the eigen values and computing the matrix V

    for (int i = 0; i <N ; ++i) {
        double max=0;
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



#pragma omp parallel for num_threads(8)
    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <N ; ++j) {
            sigma[i][j]=0;
        }
    }



    // Computing matrix sigma
#pragma omp parallel for num_threads(8)
    for (int i = 0; i <N ; ++i) {
        sigma[i][i]=sorted_singular_values[i];

    }


    //Computing V transpose
//#pragma omp parallel for num_threads(8)
//    for (int i = 0; i <N ; ++i) {
//        for (int j = 0; j <N ; ++j) {
//            V_Tnew[j][i]=V_new[i][j];
//        }
//    }


    // Computing sigma inverse

    float sigma_inv[N][M];

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <M ; ++j) {
            if (i!=j) sigma_inv[i][j]=0;
            else sigma_inv[i][j]=(1/sorted_singular_values[j]);
        }
    }

    // Computing U = D_new_T V Sigma_inv

    float intermediate[M][N];

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <N ; ++j) {
            float a=0;
            for (int k = 0; k <N ; ++k) {
                a+=D_new[i][k]*V_new[k][j];
            }
            intermediate[i][j]=a;
        }
    }

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
            float a=0;
            for (int k = 0; k <N ; ++k) {
                a+=intermediate[i][k]*sigma_inv[k][j];
            }
            U_new[i][j]=a;
        }
    }

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <N ; ++i) {
        SIGMA[0][i]=sigma[i][i];
    }

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            U[0][N*i+j]=V_new[i][j];
        }
    }

//    for (int i = 0; i <M ; ++i) {
//        for (int j = 0; j <M ; ++j) {
//            printf("%f ",intermediate[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//
//    for (int i = 0; i <M ; ++i) {
//        for (int j = 0; j <N ; ++j) {
//            printf("%f ",sigma_inv[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");


//    for (int i = 0; i <N ; ++i) {
//        for (int j = 0; j <N ; ++j) {
//            printf("%f ",U_new[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");

#pragma omp parallel for num_threads(8)
    for (int i = 0; i <M ; ++i) {
        for (int j = 0; j <M ; ++j) {
            V_T[0][M*i+j]=U_new[j][i];
        }
    }

//    float result[N][M];
//    float inter[N][M];
//    for (int i = 0; i <N ; ++i) {
//        for (int j = 0; j <M ; ++j) {
//            float a=0;
//            for (int k = 0; k <N ; ++k) {
//                a+=V_new[i][k]*sigma[j][k];
//            }
//            inter[i][j]=a;
//        }
//    }
//
//    for (int i = 0; i <N ; ++i) {
//        for (int j = 0; j <M ; ++j) {
//            float a=0;
//            for (int k = 0; k <M ; ++k) {
//                a+=inter[i][k]*U_new[j][k];
//            }
//            result[i][j]=a;
//        }
//    }
//
//    for (int i = 0; i <M ; ++i) {
//        for (int j = 0; j <N ; ++j) {
//            printf("%f ",result[j][i]);
//        }
//        printf("\n");
//    }
//    printf("\n");



}







void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    float sum_eigen_values = 0;

    for (int i = 0; i <N ; ++i) {
        sum_eigen_values+=SIGMA[i]*SIGMA[i];
    }

    int count=0;
    float sum=0;
    for (int i = 0; i <N ; ++i) {
        count+=1;
        sum+=SIGMA[i]*SIGMA[i];
        if ((sum/sum_eigen_values)*100>=retention) break;
    }

    *K=count;

    float D_hat[M][count];

    float D_new[M][N];
    float U_new[N][N];

#pragma omp parallel for num_threads(8)
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            D_new[i][j]=D[N*i+j];
        }
    }

#pragma omp parallel for num_threads(8)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            U_new[i][j]=U[N*i+j];
        }
    }

#pragma omp parallel for num_threads(8)
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
            D_HAT[0][count*i+j]=D_hat[i][j];
        }
    }

//    printf("K is %d\n",*K);
//    for (int i = 0; i <M ; ++i) {
//        for (int j = 0; j <count ; ++j) {
//            printf("%f ",D_hat[i][j]);
//        }
//        printf("\n");
//    }
}
