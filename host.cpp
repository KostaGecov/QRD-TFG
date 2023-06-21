#include <fstream>
#include <iostream>

#include "kernel.cpp"

static index_t col_offset_geqrt = 0;    // offset to access the right column for GEQRT operation
static index_t col_offset_ttqrt = 0;    // offset to access the right column for TTQRT operation
static index_t n_iter_GEQRT = 4;  // 4 iteraciones iniciales, para controlar el
                                  // nº de operaciones GEQRT en cada paso
static index_t n_iter_TTQRT = 3;  // 3 iteraciones iniciales, para controlar el nº de operaciones TTQRT en
                                  // cada paso para i=2, 4; para i = 4, 3; para i = 6, 2; para i = 8, 1;

int main() {
    data_t A[TAM][TAM];
    data_t A_aux[TAM][TAM];

    data_t A_1[TAM_TILED][TAM];
    data_t A_2[TAM_TILED][TAM];
    data_t A_3[TAM_TILED][TAM];
    data_t A_4[TAM_TILED][TAM];
    std::ifstream data_in("data_in.dat");

    if (!data_in.is_open()) {
        // Error
    }

    for (index_t r = 0; r < TAM; r++) {
        for (int c = 0; c < TAM; c++) {
            data_in >> A[r][c];
        }
    }

divide_matrix_row_for:
    for (index_t r = 0; r < TAM; r++) {
    divide_matrix_col_for:
        for (index_t c = 0; c < TAM; c++) {
            if (r >= 0 && r < TAM_TILED) {
                A_1[r][c] = A[r][c];
            } else if (r >= TAM_TILED && r < TAM_TILED * 2) {
                A_2[r - TAM_TILED][c] = A[r][c];
            } else if (r >= TAM_TILED * 2 && r < TAM_TILED * 3) {
                A_3[r - TAM_TILED * 2][c] = A[r][c];
            } else if (r >= TAM_TILED * 3 && r < TAM_TILED * 4) {
                A_4[r - TAM_TILED * 3][c] = A[r][c];
            }
        }
    }

num_operations_for:
    for (index_t i = 2; i < NUM_OPERACIONES; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            switch (n_iter_GEQRT) {
                case 4:
                    krnl_givens_rotation(A_1, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_2, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_3, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 3:
                    krnl_givens_rotation(A_2, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_3, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 2:
                    krnl_givens_rotation(A_3, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 1:
                    krnl_givens_rotation(A_4, A_aux, GEQRT, col_offset_geqrt);
                    break;
                default:
                    break;
            }
            n_iter_GEQRT--;
            col_offset_geqrt += TAM_TILED;

        } else {
            switch (n_iter_TTQRT) {
                case 3:
                    krnl_givens_rotation(A_1, A_2, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_3, A_4, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_3, TTQRT, col_offset_ttqrt);
                    break;
                case 2:
                    krnl_givens_rotation(A_2, A_3, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_4, TTQRT, col_offset_ttqrt);
                    break;
                case 1:
                    krnl_givens_rotation(A_3, A_4, TTQRT, col_offset_ttqrt);
                    break;

                default:
                    break;
            }
            n_iter_TTQRT--;
            col_offset_ttqrt += TAM_TILED;
        }
    }

// Escritura de solución en matriz grande
write_sol_to_matrix_row_for:
    for (index_t r = 0; r < TAM; r++) {
    write_sol_to_matrix_col_for:
        for (index_t c = 0; c < TAM; c++) {
            if (r >= 0 && r < TAM_TILED) {
                A[r][c] = A_1[r][c];
            } else if (r >= TAM_TILED && r < TAM_TILED * 2) {
                A[r][c] = A_2[r - TAM_TILED][c];
            } else if (r >= TAM_TILED * 2 && r < TAM_TILED * 3) {
                A[r][c] = A_3[r - TAM_TILED * 2][c];
            } else if (r >= TAM_TILED * 3 && r < TAM_TILED * 4) {
                A[r][c] = A_4[r - TAM_TILED * 3][c];
            }
        }
    }

    // Print result matrix
    for (int i = 0; i < TAM; i++) {
        for (int j = 0; j < TAM; j++) {
            std::cout << A[i][j] << "  |  ";
        }
        std::cout << std::endl;
    }
    return 0;
}
