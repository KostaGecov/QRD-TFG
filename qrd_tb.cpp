#include <iostream>

#include "qrd.h"

int main() {
    data_t A[TAM_TILED][TAM] = {{3, 2, -1, 4, -7, 8, 3, 2, -1, 4, -7, 8,
                                 3, 2, -1, 4, -7, 8, 3, 2, -1, 4, -7, 8},
                                {2, 1, 9, 7, -1, 5, 2, 1, 9, 7, -1, 5,
                                 2, 1, 9, 7, -1, 5, 2, 1, 9, 7, -1, 5},
                                {-7, 5, 2, -6, 3, 4, -7, 5, 2, -6, 3, 4,
                                 -7, 5, 2, -6, 3, 4, -7, 5, 2, -6, 3, 4},
                                {-1, 2, 1, -4, 6, -9, -1, 2, 1, -4, 6, -9,
                                 -1, 2, 1, -4, 6, -9, -1, 2, 1, -4, 6, -9},
                                {-9, 16, 3, 1, 4, 11, -9, 16, 3, 1, 4, 11,
                                 -9, 16, 3, 1, 4, 11, -9, 16, 3, 1, 4, 11},
                                {4, -2, 3, 7, 2, 2, 4, -2, 3, 7, 2, 2,
                                 4, -2, 3, 7, 2, 2, 4, -2, 3, 7, 2, 2}};

    data_t A_rot[TAM_TILED][TAM];

    data_t A_1[TAM_TILED][TAM];
    data_t A_2[TAM_TILED][TAM];
    data_t A_3[TAM_TILED][TAM];
    data_t A_4[TAM_TILED][TAM];

    static index_t col_offset = 0;  // offset to access the right column for GEQRT operation
    static index_t n_iter = 4;      // 4 iteraciones iniciales, para controlar el nº de operaciones GEQRT en cada paso
                                    // para i=2, 4; para i = 4, 3; para i = 6, 2; para i = 8, 1;
    static index_t n_vez = 0;

divide_matrix_row_for:
    for (index_t r = 0; r < TAM; r++) {
    divide_matrix_col_for:
        for (index_t c = 0; c < TAM; c++) {
            if (r >= 0 && r <= 5) {
                A_1[r][c] = A[r][c];
            } else if (r >= 6 && r <= 11) {
                A_2[r][c] = A[r][c];
            } else if (r >= 12 && r <= 17) {
                A_3[r][c] = A[r][c];
            } else if (r >= 18 && r <= 23) {
                A_4[r][c] = A[r][c];
            }
        }
    }

    /*
     * ToDo: read input matrix and gold results from external file and implement
     * Tile algorithm
     */

num_operations_for:
    for (index_t i = 2; i < NUM_OPERACIONES; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            switch (n_iter) {
                case 4:
                    krnl_givens_rotation(A_1, A_rot, GEQRT, col_offset, n_iter);
                    krnl_givens_rotation(A_2, A_rot, GEQRT, col_offset, n_iter);
                    krnl_givens_rotation(A_3, A_rot, GEQRT, col_offset, n_iter);
                    krnl_givens_rotation(A_4, A_rot, GEQRT, col_offset, n_iter);
                    break;
                case 3:
                    krnl_givens_rotation(A_2, A_rot, GEQRT, col_offset, n_iter);
                    krnl_givens_rotation(A_3, A_rot, GEQRT, col_offset, n_iter);
                    krnl_givens_rotation(A_4, A_rot, GEQRT, col_offset, n_iter);
                    break;
                case 2:
                    krnl_givens_rotation(A_3, A_rot, GEQRT, col_offset, n_iter);
                    krnl_givens_rotation(A_4, A_rot, GEQRT, col_offset, n_iter);
                    break;
                case 1:
                    krnl_givens_rotation(A_4, A_rot, GEQRT, col_offset, n_iter);
                    break;
                default:
                    break;
            }

            n_iter--;
            col_offset += 6;

            /* for (index_t j = n_vez; j < TAM; j += 6) {  // bucle para controlar las submatrices que tengo que operar (4, 3, 2 o 1)
                krnl_givens_rotation(A, A_rot, GEQRT, col_offset, n_iter);
                col_offset += 6;
                n_iter--;
            }
            n_vez += 5; */

        } else {
            /*
            ToDo: estructurar como con operación GEQRT
            */
            krnl_givens_rotation(A, A_rot, TTQRT, 0, 0);
        }
    }

// Escritura de solución en matriz grande
write_sol_to_matrix_row_for:
    for (index_t r = 0; r < TAM; r++) {
    write_sol_to_matrix_col_for:
        for (index_t c = 0; c < TAM; c++) {
            if (r >= 0 && r <= 5) {
                A[r][c] = A_1[r][c];
            } else if (r >= 6 && r <= 11) {
                A[r][c] = A_2[r][c];
            } else if (r >= 12 && r <= 17) {
                A[r][c] = A_3[r][c];
            } else if (r >= 18 && r <= 23) {
                A[r][c] = A_4[r][c];
            }
        }
    }

    // Print result matrix
    for (int i = 0; i < TAM_TILED; i++) {
        for (int j = 0; j < TAM; j++) {
            std::cout << A_rot[i][j] << "  |  ";
        }
        std::cout << std::endl;
    }
    return 0;
}
