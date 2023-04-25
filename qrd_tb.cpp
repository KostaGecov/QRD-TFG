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
    static index_t n_iter = 3;      // 4 iteraciones iniciales
    static index_t n_vez = 0;

    for (index_t r = 0; r < TAM; r++) {
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
            for (index_t j = n_vez; j < 24; j += 6) {  // bucle para controlar las submatrices que tengo que operar (4, 3, 2 o 1)
                krnl_givens_rotation(A, A_rot, GEQRT, col_offset, n_iter);
                col_offset += 6;
                n_iter--;
            }
            n_vez += 5;
            /* ToDo: volver a escribir en matriz grande 24x24*/
        } else {
            krnl_givens_rotation(A, A_rot, TTQRT, 0, 0);
            /* ToDo: volver a escribir en matriz grande 24x24*/
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
