// вариант 00 10

#include <iostream>

#include <cmath>

using namespace std;

double ScalarProduce(int N, double *lambda, double *Y_m, double *Y_n)
{
    double produce = .0;
    for (int i = 0; i < N + 1; i++)
    {
        if (!isnan(lambda[i]) &&
            !isnan(Y_m[i]) &&
            !isnan(Y_n[i]) &&
            !isinf(lambda[i]) &&
            !isinf(Y_m[i]) &&
            !isinf(Y_n[i]))
        {
            produce += lambda[i] * Y_m[i] * Y_n[i];
        }
    }

    return produce;
}

double Lambda(int m, int N, double h)
{
    return 4.0 / (h * h) * pow(sin((M_PI * m) / (2 * N - 1)), 2);
}

void fullCoef(double *coef, int N)
{
    for (int i = 0; i < N + 1; i++)
    {
        if ((i == 0) || (i == N))
        {
            coef[i] = .0;
        }
        else
        {
            coef[i] = 1.0 / (2 * N - 1);
        }
    }
}

//

void fullY_m(double *Y_m, int N, int m)
{
    for (int k = 0; k < N + 1; ++k)
    {
        Y_m[k] = sin((M_PI * 2 * m * k) / (2 * N - 1));
    }
}

void printVec(double *vec, int N)
{
    for (int i = 0; i < N + 1; i++)
    {
        printf("%e,   ", vec[i]);
    }
    printf("\n\n");
}

double MaxScalarProduce(double *coef, int N)
{
    double *Y_m = new double[N + 1];
    double *Y_n = new double[N + 1];

    double maxProduce = .0;

    for (int m = 1; m < N; m++)
    {
        for (int n = 1; n < N; n++)
        {
            if (n != m)
            {
                fullY_m(Y_m, N, m);
                fullY_m(Y_n, N, n);
                printf("%d = m,     %d = n, ", m, n);
                printVec(Y_m, N);
                maxProduce = max(maxProduce, ScalarProduce(N, coef, Y_m, Y_n));
                // printf("Y_m: \n");
                // printVec(Y_m, N);
                // printf("Y_n: \n");
                // printVec(Y_n, N);
            }
        }
    }
    delete[] Y_m;
    delete[] Y_n;
    return maxProduce;
}

// len(tmp)=N+1
// len(Y_m)=N+1
// fulling q

double MaxNorm(int N, double *coefs)
{
    double eps = 1e-16;
    double h = 1.0 / N;
    double *tmp = new double[N + 1];
    double *Y_m = new double[N + 1];
    double maxNorm = .0;
    double norm_Ay_ly = .0;
    for (int m = 1; m < N; m++)
    {
        fullY_m(Y_m, N, m);

        // [0, ..., N-1] len=N
        for (int k = 0; k < N; k++)
        {
            if (k == 0)
            {
                tmp[k] = Y_m[k];
            }
            else
            {
                tmp[k] = (1. / (h * h) * (Y_m[k - 1] - 2 * Y_m[k] + Y_m[k + 1])) + (Lambda(m, N, h) * Y_m[k]);
            }
        }

        norm_Ay_ly = sqrt(ScalarProduce(N, coefs, tmp, tmp));
        if (Lambda(m, N, h) >= eps)
        {
            maxNorm = max(maxNorm, norm_Ay_ly / (Lambda(m, N, h)));
        }
        else
        {
            printf("Lambda is ZERO");
        }
        // printf("norm_Ay_ly :=%e  lamnda := %e \n", norm_Ay_ly, Lambda(m, N, h));
        maxNorm = max(maxNorm, norm_Ay_ly / (Lambda(m, N, h)));
    }

    delete[] Y_m;
    delete[] tmp;
    return maxNorm;
}

int main()
{
    int N;
    for (int k = 2; k < 3; ++k)
    {
        N = pow(2, k);

        double *coef = new double[N + 1];
        fullCoef(coef, N);
        // printf("coefs: \n");
        // printVec(coef, N);

        cout << "Max produce where k=" << k << endl;
        printf("%e \n", MaxScalarProduce(coef, N));

        cout << "Max norma where k=" << k << endl;
        printf("%e \n", MaxNorm(N, coef));
        cout << endl;
    }

    return 0;
}
