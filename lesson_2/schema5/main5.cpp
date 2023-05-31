#include <iostream>
#include <fstream>
#include <cmath>
#include <cfenv>
#include <float.h>

using namespace std;

double eps = 1e-15;
double MAX = 1e+50;

double exp_Ax(double A, double x)
{
    return exp(-A * x);
}

double MaxAbsDiff(double *y_k, double *y_xk, int k)
{
    double maxAbs = .0;
    double diff = .0;
    for (int i = 1; i < k; i++)
    {
        diff = abs(y_k[i] - y_xk[i]);
        if (diff > maxAbs)
        {
            maxAbs = diff;
        }
    }
    return maxAbs;
}

int main()
{
    feclearexcept(FE_ALL_EXCEPT);
    // if (fetestexcept(FE_INVALID))
    //     cout << "invalid result reported\n";
    // else
    //     cout << "invalid result not reported\n";

    ofstream outfile("/Users/olegkirillov/Desktop/Work/EVM_Vesna/github.com/lesson_2/schema5/output.txt");
    for (int n = 1; n <= 6; n++)
    {
        int A = 1;
        while (A <= 1000)
        {
            fenv_t env;
            feholdexcept(&env);
            int k = pow(10, n);

            double h = 1. / k;
            double *x_k = new double[k];
            x_k[0] = 0;
            x_k[1] = h;

            double *y_k = new double[k];
            double *y_xk = new double[k];
            y_k[0] = y_xk[0] = 1;
            y_k[1] = 1 - A * h;
            y_xk[1] = exp_Ax(A, x_k[1]);

            for (int i = 2; i < k; i++)
            {
                if (fetestexcept(/*FE_OVERFLOW | */ FE_UNDERFLOW))
                {
                    feclearexcept(/*FE_OVERFLOW | */ FE_UNDERFLOW);
                    throw runtime_error("FPE");
                }
                else
                {
                    x_k[i] = x_k[i - 1] + h;
                    if ((-2 * h * A * y_k[i - 1] + y_k[i - 2] > MAX) ||
                        (-2 * h * A * y_k[i - 1] + y_k[i - 2] < eps) ||
                        (exp_Ax(A, x_k[i]) > MAX) ||
                        (exp_Ax(A, x_k[i]) < eps))
                    {
                        continue;
                    }
                    else
                    {
                        y_k[i] = (2 * y_k[i - 1] - 0.5 * y_k[i - 2]) / (1.5 + A * h);
                        y_xk[i] = exp_Ax(A, x_k[i]);
                    }
                }
            }

            outfile << "Max abs where n = " << n << " and A = " << A << " is  " << MaxAbsDiff(y_k, y_xk, k) << endl;

            A *= 10;
        }
    }

    outfile.close();
}
