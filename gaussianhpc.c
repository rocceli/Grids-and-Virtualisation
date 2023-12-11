#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

void simulateHeatConduction(int m, int n, double tol) {
    double t[m + 2][n + 2], tnew[m + 2][n + 2], diff, difmax;

    for (int i = 0; i <= m + 1; i++) {
        for (int j = 0; j <= n + 1; j++) {
            if (i == 0)
                t[i][j] = 10.0; // Top boundary
            else if (i == m + 1)
                t[i][j] = 30.0; // Bottom boundary
            else if (j == 0)
                t[i][j] = 40.0; // Left boundary
            else if (j == n + 1)
                t[i][j] = 50.0; // Right boundary
            else
                t[i][j] = 30.0;
        }
    }

    double start_time = omp_get_wtime();

    int iter = 0;
    difmax = tol + 1.0;
    while (difmax > tol) {
        iter++;
        difmax = 0.0;

        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                tnew[i][j] = (t[i - 1][j] + t[i + 1][j] + t[i][j - 1] + t[i][j + 1]) / 4.0;

                diff = fabs(tnew[i][j] - t[i][j]);
                if (diff > difmax) {
                    difmax = diff;
                }
                t[i][j] = tnew[i][j];
            }
        }
    }

    double end_time = omp_get_wtime();
    printf("Problem Size: %dx%d, Execution Time: %lf seconds\n", m, n, end_time - start_time);
}

int main(int argc, char *argv[]) {
    int sizes[] = {100, 500, 1000};
    double tolerance = 0.001;

    for (int i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
        for (int j = 2; j <= 8; j++) { 
            omp_set_num_threads(j);
            printf("Simulation for size %dx%d, Threads: %d\n", sizes[i], sizes[i], j);
            simulateHeatConduction(sizes[i], sizes[i], tolerance);
        }
    }

    return 0;
}
