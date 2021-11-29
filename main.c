#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

void print_flatten_matrix(double *matrix_data, int n_x, int n_z) {
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_z; j++) {
            printf("%8.3f ", matrix_data[j * n_z + i]);
        }
        printf("\n");
    }
    printf("\n");
}

double boundary(double x, double z) {
    return 300;
}

double f(double x, double z) {
    return 1000 * exp(-pow(x - 5, 2) * pow(z - 5, 2));
}

int p_helper(int i, int j, int n_x) {
    return j * (n_x) + i;
}

double *solve(double x_bound, double z_bound, double k, double h_step) {
    int n_x = x_bound / h_step + 1;
    int n_z = z_bound / h_step + 1;
    int n = n_x * n_z;
    double *a_matrix_data = calloc(n * n, sizeof(double));
    double *b_vector = calloc(n, sizeof(double));
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_z; j++) {
            int p = p_helper(i, j, n_x);
            if (i == 0 || j == 0 || i == n_x - 1 || j == n_z - 1) {
                a_matrix_data[p * n + p] = 1;
                b_vector[p] = boundary(i * h_step, j * h_step);
            } else {
                a_matrix_data[p_helper(i - 1, j, n_x) * n + p] = 1;
                a_matrix_data[p_helper(i + 1, j, n_x) * n + p] = 1;
                a_matrix_data[p_helper(i, j - 1, n_x) * n + p] = 1;
                a_matrix_data[p_helper(i, j + 1, n_x) * n + p] = 1;
                a_matrix_data[p * n + p] = -4;

                b_vector[p] = -f(i * h_step, j * h_step) * h_step * h_step / k;
            }
        }
    }
    LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', n, n, 1, a_matrix_data, n, b_vector, n);
    free(a_matrix_data);
    return b_vector;
}

int main() {
    double a = 10;
    double b = 10;
    double h = 2;
    double k = 0.536;
    double *solution = solve(a, b, k, h);
    print_flatten_matrix(solution, a / h + 1, b / h + 1);
    return 0;
}
