#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[]) {
    int m, n;
    double tol;

    int i, j, iter;

    m = atoi(argv[1]);
    n = atoi(argv[2]);
    tol = atof(argv[3]);

    double t[m + 2][n + 2], tnew[m + 2][n + 2], diff, difmax;

    printf("%d %d %lf\n", m, n, tol);

    // initialise temperature array
    for (i = 0; i <= m + 1; i++) {
        for (j = 0; j <= n + 1; j++) {
            t[i][j] = 30.0;
        }
    }

    // fix boundary conditions
    for (i = 1; i <= m; i++) {
        t[i][0] = 40.0;
        t[i][n + 1] = 50.0;
    }
    for (j = 1; j <= n; j++) {
        t[0][j] = 10.0;
        t[m + 1][j] = 30.0;
    }

    // main loop
    iter = 0;
    difmax = 100000.0; // Initialize difmax to a value larger than tol
    while (difmax > tol) {
        iter++;
        difmax = 0.0; // Reset difmax for each iteration

        // update temperature for next iteration
        #pragma omp parallel for private(i, j) shared(t, tnew) reduction(max : difmax)
        for (i = 1; i <= m; i++) {
            for (j = 1; j <= n; j++) {
                tnew[i][j] = (t[i - 1][j] + t[i + 1][j] + t[i][j - 1] + t[i][j + 1]) / 4.0;

                // work out maximum difference between old and new temperatures
                diff = fabs(tnew[i][j] - t[i][j]);
                if (diff > difmax) {
                    #pragma omp critical
                    if (diff > difmax) {
                        difmax = diff;
                    }
                }
                t[i][j] = tnew[i][j];
            }
        }
    }

    // print results
    printf("iter = %d  difmax = %9.11lf\n", iter, difmax);
    for (i = 0; i <= m + 1; i++) {
        for (j = 0; j <= n + 1; j++) {
            printf("%3.5lf ", t[i][j]);
        }
        printf("\n");
    }

    return 0;
}
