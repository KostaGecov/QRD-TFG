/**
 * @file host.cpp
 * @author Kosta Gecov (kostagecov@gmail.com)
 * @brief Vitis Hls host that implements the tiled QRD algorithm
 * @version 0.1
 * @date 2023-07-29
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <fstream>
#include <iostream>

#include "kernel.h"

void init_matrix(data_t matrix[TAM][TAM], std::fstream *file) {
    if (!file->is_open()) {
        std::cout << "Could not open file" << std::endl;
    } else {
        std::cout << "Opened file data file" << std::endl;
    }

initialize_matrix:
    for (index_t r = 0; r < TAM; r++) {
        for (index_t c = 0; c < TAM; c++) {
            *file >> matrix[r][c];
        }
    }
    file->close();
}

void init_matrix(float matrix[TAM][TAM], std::fstream *file) {
    if (!file->is_open()) {
        std::cout << "Could not open file" << std::endl;
    } else {
        std::cout << "Opened file data file" << std::endl;
    }

initialize_matrix:
    for (index_t r = 0; r < TAM; r++) {
        for (index_t c = 0; c < TAM; c++) {
            *file >> matrix[r][c];
        }
    }
    file->close();
}

float error(data_t A[TAM][TAM], float out_gold[TAM][TAM]) {
    float err = 0.0;
    data_t resA = 0.0;
    float resOut = 0.0;

    // to access just to the non zero elements (upper triangular matrix)
    for (index_t r = 0; r < TAM; r++) {
        for (index_t c = r; c < TAM; c++) {
            resA = A[r][c];
            resOut = out_gold[r][c];
            err += pow((abs((float)resA) - abs(resOut)), 2);
        }
    }

    // err / number of elements in the upper triangular matrix (including diagonal)
    return err / ((TAM * (TAM + 1)) / 2);
}

int main() {
    /**
     * offset to access the right column for GEQRT operation
     */
    static index_t col_offset_geqrt = 0;
    /**
     * offset to access the right column for TTQRT operation
     */
    static index_t col_offset_ttqrt = 0;
    /**
     * to control the GEQRT operations in each step
     */
    static index_t n_iter_GEQRT = 32;
    /**
     * to control the TTQRT operations in each step
     */
    static index_t n_iter_TTQRT = 31;

    /**
     * stores all 32 A matrices needed for tiled operations
     *
     */
    data_t A_tiled[NUM_TILED][TAM_TILED][TAM];
    data_t A[TAM][TAM];

    /**
     * stores all 32 Q matrices needed for tiled operations
     *
     */
    data_t Q_tiled[NUM_TILED][TAM_TILED][TAM];
    data_t Q[TAM][TAM];

    /**
     * to store the output data gold
     *
     */
    float out_gold[TAM][TAM];

    std::fstream data_in("data_in.dat", std::ios::in);
    init_matrix(A, &data_in);

    std::fstream data_out_gold("data_out_gold.dat", std::ios::in);
    init_matrix(out_gold, &data_out_gold);

initialize_Q:
    for (index_t r = 0; r < TAM; r++) {
        for (index_t c = 0; c < TAM; c++) {
            if (r == c) {
                Q[r][c] = 1;
            } else {
                Q[r][c] = 0;
            }
        }
    }

divide_matrices_row_for:
    for (index_t r = 0; r < TAM; r++) {
    	int tile_index = r / TAM_TILED;
    	int tile_offset = r % TAM_TILED;

    divide_matrices_col_for:
        for (index_t c = 0; c < TAM; c++) {
        	A_tiled[tile_index][tile_offset][c] = A[r][c];
			Q_tiled[tile_index][tile_offset][c] = Q[r][c];
        }
    }

num_operations_for:
    for (index_t i = 2; i < NUM_OPERACIONES; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            switch (n_iter_GEQRT) {
                case 32:
                    std::cout << "Entro en GEQRT 32" << std::endl;
                    for (index_t idx_mat_1 = 0; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_1, Q_1, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_2, Q_2, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_3, Q_3, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, Q_4, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5, Q_5, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, Q_6, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, Q_7, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, Q_8, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, Q_9, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, Q_10, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, Q_11, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, Q_12, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, Q_13, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, Q_14, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, Q_15, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, Q_16, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, Q_17, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, Q_18, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, Q_19, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, Q_20, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, Q_21, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, Q_22, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, Q_23, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, Q_24, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, Q_25, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, Q_26, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, Q_27, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, Q_28, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, Q_29, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, Q_30, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, Q_31, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, Q_32, GEQRT, col_offset_geqrt); */
                    break;
                case 31:
                    std::cout << "Entro en GEQRT 31" << std::endl;
                    for (index_t idx_mat_1 = 1; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_2,  Q_2,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_3,  Q_3,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4,  Q_4,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5,  Q_5,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6,  Q_6,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7,  Q_7,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8,  Q_8,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9,  Q_9,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 30:
                    std::cout << "Entro en GEQRT 30" << std::endl;
                    for (index_t idx_mat_1 = 2; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_3,  Q_3,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4,  Q_4,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5,  Q_5,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6,  Q_6,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7,  Q_7,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8,  Q_8,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9,  Q_9,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 29:
                    std::cout << "Entro en GEQRT 29" << std::endl;
                    for (index_t idx_mat_1 = 3; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_4,  Q_4,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5,  Q_5,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6,  Q_6,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7,  Q_7,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8,  Q_8,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9,  Q_9,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 28:
                    std::cout << "Entro en GEQRT 28" << std::endl;
                    for (index_t idx_mat_1 = 4; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_5,  Q_5,   GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6,  Q_6,   GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7,  Q_7,   GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8,  Q_8,   GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9,  Q_9,   GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 27:
                    std::cout << "Entro en GEQRT 27" << std::endl;
                    for (index_t idx_mat_1 = 5; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_6,  Q_6,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7,  Q_7,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8,  Q_8,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9,  Q_9,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 26:
                    std::cout << "Entro en GEQRT 26" << std::endl;
                    for (index_t idx_mat_1 = 6; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_7,  Q_7,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8,  Q_8,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9,  Q_9,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 25:
                    std::cout << "Entro en GEQRT 25" << std::endl;
                    for (index_t idx_mat_1 = 7; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_8,  Q_8,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9,  Q_9,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 24:
                    std::cout << "Entro en GEQRT 24" << std::endl;
                    for (index_t idx_mat_1 = 8; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_9,  Q_9,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 23:
                    std::cout << "Entro en GEQRT 23" << std::endl;
                    for (index_t idx_mat_1 = 9; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_10,  Q_10,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 22:
                    std::cout << "Entro en GEQRT 22" << std::endl;
                    for (index_t idx_mat_1 = 10; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_11,  Q_11,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 21:
                    std::cout << "Entro en GEQRT 21" << std::endl;
                    for (index_t idx_mat_1 = 11; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_12,  Q_12,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 20:
                    std::cout << "Entro en GEQRT 20" << std::endl;
                    for (index_t idx_mat_1 = 12; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_13,  Q_13,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 19:
                    std::cout << "Entro en GEQRT 19" << std::endl;
                    for (index_t idx_mat_1 = 13; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_14,  Q_14,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 18:
                    std::cout << "Entro en GEQRT 18" << std::endl;
                    for (index_t idx_mat_1 = 14; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_15,  Q_15,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 17:
                    std::cout << "Entro en GEQRT 17" << std::endl;
                    for (index_t idx_mat_1 = 15; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_16,  Q_16,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 16:
                    std::cout << "Entro en GEQRT 16" << std::endl;
                    for (index_t idx_mat_1 = 16; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_17,  Q_17,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 15:
                    std::cout << "Entro en GEQRT 15" << std::endl;
                    for (index_t idx_mat_1 = 17; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_18,  Q_18,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 14:
                    std::cout << "Entro en GEQRT 14" << std::endl;
                    for (index_t idx_mat_1 = 18; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_19,  Q_19,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 13:
                    std::cout << "Entro en GEQRT 13" << std::endl;
                    for (index_t idx_mat_1 = 19; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_20,  Q_20,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 12:
                    std::cout << "Entro en GEQRT 12" << std::endl;
                    for (index_t idx_mat_1 = 20; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_21,  Q_21,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 11:
                    std::cout << "Entro en GEQRT 11" << std::endl;
                    for (index_t idx_mat_1 = 21; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_22,  Q_22,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 10:
                    std::cout << "Entro en GEQRT 10" << std::endl;
                    for (index_t idx_mat_1 = 22; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_23,  Q_23,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 9:
                    std::cout << "Entro en GEQRT 9" << std::endl;
                    for (index_t idx_mat_1 = 23; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_24,  Q_24,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 8:
                    std::cout << "Entro en GEQRT 8" << std::endl;
                    for (index_t idx_mat_1 = 24; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_25,  Q_25,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 7:
                    std::cout << "Entro en GEQRT 7" << std::endl;
                    for (index_t idx_mat_1 = 25; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_26,  Q_26,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 6:
                    std::cout << "Entro en GEQRT 6" << std::endl;
                    for (index_t idx_mat_1 = 26; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_27,  Q_27,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 5:
                    std::cout << "Entro en GEQRT 5" << std::endl;
                    for (index_t idx_mat_1 = 27; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_28,  Q_28,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 4:
                    std::cout << "Entro en GEQRT 4" << std::endl;
                    for (index_t idx_mat_1 = 28; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_29,  Q_29,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 3:
                    std::cout << "Entro en GEQRT 3" << std::endl;
                    for (index_t idx_mat_1 = 29; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_30,  Q_30,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 2:
                    std::cout << "Entro en GEQRT 2" << std::endl;
                    for (index_t idx_mat_1 = 30; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }
                    /* krnl_givens_rotation(A_31,  Q_31,  GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                case 1:
                    std::cout << "Entro en GEQRT 1" << std::endl;
                    krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, 31, 0);

                    /* krnl_givens_rotation(A_32,  Q_32,  GEQRT, col_offset_geqrt); */
                    break;
                default:
                    break;
            }
            n_iter_GEQRT--;
            col_offset_geqrt += TAM_TILED;

        } else {
            switch (n_iter_TTQRT) {
                case 31:
                    std::cout << "Entro en TTQRT 31" << std::endl;
                    for (index_t idx_mat_1 = 0, idx_mat_2 = 1; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_1, A_2, Q_1, Q_2, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_3, A_4, Q_3, Q_4, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_5, A_6, Q_5, Q_6TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_8, Q_7, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_10, Q_9, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, Q_11, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, Q_13, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 0, idx_mat_2 = 2; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_1, A_3, Q_1, Q_3, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_5, A_7, Q_5, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_11, Q_9, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, Q_13, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 0, idx_mat_2 = 4; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_1, A_5, Q_1, Q_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_13, Q_9, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_21, Q_17, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, Q_25, Q_29, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 0, idx_mat_2 = 8; idx_mat_2 < 25; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_1, A_9, Q_1, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_25, Q_17, Q_25, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 0, 16);

                    /* krnl_givens_rotation(A_1, A_17, Q_1, Q_17, TTQRT, col_offset_ttqrt); */
                    break;
                case 30:
                    std::cout << "Entro en TTQRT 30" << std::endl;
                    for (index_t idx_mat_1 = 1, idx_mat_2 = 2; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_2, A_3, Q_2, Q_3, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_4, A_5, Q_4, Q_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_6, A_7, Q_6, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_9, Q_8, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_11, Q_10, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, Q_12, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 1, idx_mat_2 = 3; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_2, A_4, Q_2, Q_4, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_6, A_8, Q_6, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_12, Q_10, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, Q_14, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 1, idx_mat_2 = 5; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_2, A_6, Q_2, Q_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_14, Q_10, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_22, Q_18, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, Q_26, Q_30, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 1, idx_mat_2 = 9; idx_mat_2 < 26; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_2, A_10, Q_2, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_26, Q_18, Q_26, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 1, 17);

                    /* krnl_givens_rotation(A_2, A_18, Q_2, Q_18, TTQRT, col_offset_ttqrt); */
                    break;
                case 29:
                    std::cout << "Entro en TTQRT 29" << std::endl;
                    for (index_t idx_mat_1 = 2, idx_mat_2 = 3; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_3, A_4, Q_3, Q_4, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_5, A_6, Q_5, Q_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_8, Q_7, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_10, Q_9, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, Q_11, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, Q_13, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 2, idx_mat_2 = 4; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_3, A_5, Q_3, Q_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_9, Q_7, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_13, Q_11, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, Q_15, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 2, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_3, A_7, Q_3, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_15, Q_11, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_23, Q_19, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_31, Q_27, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 2, idx_mat_2 = 10; idx_mat_2 < 27; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_3, A_11, Q_3, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_27, Q_19, Q_27, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 2, 18);

                    /* krnl_givens_rotation(A_3, A_19, Q_3, Q_19, TTQRT, col_offset_ttqrt); */
                    break;
                case 28:
                    std::cout << "Entro en TTQRT 28" << std::endl;
                    for (index_t idx_mat_1 = 3, idx_mat_2 = 4; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_4, A_5, Q_4, Q_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_6, A_7, Q_6, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_9, Q_8, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_11, Q_10, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, Q_12, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 3, idx_mat_2 = 5; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_4, A_6, Q_4, Q_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_10, Q_8, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_14, Q_12, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, Q_16, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 3, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_4, A_8, Q_4, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_16, Q_12, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_24, Q_20, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, Q_28, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 3, idx_mat_2 = 11; idx_mat_2 < 28; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_4, A_12, Q_4, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_28, Q_20, Q_28, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 3, 19);

                    /* krnl_givens_rotation(A_4, A_20, Q_4, Q_20, TTQRT, col_offset_ttqrt); */
                    break;
                case 27:
                    std::cout << "Entro en TTQRT 27" << std::endl;
                    for (index_t idx_mat_1 = 4, idx_mat_2 = 5; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_5, A_6, Q_5, Q_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_8, Q_7, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_10, Q_9, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, Q_11, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, Q_13, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 4, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_5, A_7, Q_5, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_11, Q_9, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, Q_13, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 4, idx_mat_2 = 8; idx_mat_2 < 25; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_5, A_9, Q_5, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_17, Q_13, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_25, Q_21, Q_25, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 4, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_5, A_13, Q_5, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_29, Q_21, Q_29, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 4, 20);

                    /* krnl_givens_rotation(A_5, A_21, Q_5, Q_21, TTQRT, col_offset_ttqrt); */
                    break;
                case 26:
                    std::cout << "Entro en TTQRT 26" << std::endl;
                    for (index_t idx_mat_1 = 5, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_6, A_7, Q_6, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_9, Q_8, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_11, Q_10, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, Q_12, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 5, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_6, A_8, Q_6, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_12, Q_10, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, Q_14, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 5, idx_mat_2 = 9; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_6, A_10, Q_6, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_18, Q_14, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_26, Q_22, Q_26, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 5, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_6, A_14, Q_6, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_30, Q_22, Q_30, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 5, 21);

                    /* krnl_givens_rotation(A_6, A_22, Q_6, Q_22, TTQRT, col_offset_ttqrt); */
                    break;
                case 25:
                    std::cout << "Entro en TTQRT 25" << std::endl;
                    for (index_t idx_mat_1 = 6, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_7, A_8, Q_7, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_10, Q_9, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, Q_11, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, Q_13, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 6, idx_mat_2 = 8; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_7, A_9, Q_7, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_13, Q_11, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, Q_15, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 6, idx_mat_2 = 10; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_7, A_11, Q_7, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_19, Q_15, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_27, Q_23, Q_27, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 6, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_7, A_15, Q_7, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_31, Q_23, Q_31, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 6, 22);
                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 6, 26);
                    /* krnl_givens_rotation(A_7, A_23, Q_7, Q_23, TTQRT, col_offset_ttqrt); */
                    break;
                case 24:
                    std::cout << "Entro en TTQRT 24" << std::endl;
                    for (index_t idx_mat_1 = 7, idx_mat_2 = 8; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_8, A_9, Q_8, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_11, Q_10, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, Q_12, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 7, idx_mat_2 = 9; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_8, A_10, Q_8, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_14, Q_12, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, Q_16, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 7, idx_mat_2 = 11; idx_mat_2 < 28; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_8, A_12, Q_8, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_20, Q_16, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_28, Q_24, Q_28, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 7, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_8, A_16, Q_8, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_32, Q_24, Q_32, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 7, 23);

                    /* krnl_givens_rotation(A_8, A_24, Q_8, Q_24, TTQRT, col_offset_ttqrt); */
                    break;
                case 23:
                    std::cout << "Entro en TTQRT 23" << std::endl;
                    for (index_t idx_mat_1 = 8, idx_mat_2 = 9; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_9, A_10, Q_9, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, Q_11, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, Q_13, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 8, idx_mat_2 = 10; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_9, A_11, Q_9, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, Q_13, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 8, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_9, A_13, Q_9, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_21, Q_17, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, Q_25, Q_29, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 8, idx_mat_2 = 16; idx_mat_2 < 25; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_9, A_17, Q_9, Q_17, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_25, Q_9, Q_25, TTQRT, col_offset_ttqrt); */
                    break;
                case 22:
                    std::cout << "Entro en TTQRT 22" << std::endl;
                    for (index_t idx_mat_1 = 9, idx_mat_2 = 10; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_10, A_11, Q_10, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, Q_12, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 9, idx_mat_2 = 11; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_10, A_12, Q_10, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, Q_14, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 9, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_10, A_14, Q_10, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_22, Q_18, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, Q_26, Q_30, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 9, idx_mat_2 = 17; idx_mat_2 < 26; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_10, A_18, Q_10, Q_18, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_26, Q_10, Q_26, TTQRT, col_offset_ttqrt); */
                    break;
                case 21:
                    for (index_t idx_mat_1 = 10, idx_mat_2 = 11; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_11, A_12, Q_11, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, Q_13, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 10, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_11, A_13, Q_11, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, Q_15, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 10, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_11, A_15, Q_11, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_23, Q_19, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_31, Q_27, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 10, idx_mat_2 = 18; idx_mat_2 < 27; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_11, A_19, Q_11, Q_19, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_27, Q_11, Q_27, TTQRT, col_offset_ttqrt); */
                    break;
                case 20:
                    for (index_t idx_mat_1 = 11, idx_mat_2 = 12; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_12, A_13, Q_12, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 11, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_12, A_14, Q_12, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, Q_16, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 11, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_12, A_16, Q_12, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_24, Q_20, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, Q_28, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 11, idx_mat_2 = 19; idx_mat_2 < 28; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_12, A_20, Q_12, Q_20, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_28, Q_12, Q_28, TTQRT, col_offset_ttqrt); */
                    break;
                case 19:
                    for (index_t idx_mat_1 = 12, idx_mat_2 = 13; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_13, A_14, Q_13, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 12, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_13, A_15, Q_13, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 12, idx_mat_2 = 16; idx_mat_2 < 25; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_13, A_17, Q_13, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_25, Q_21, Q_25, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 12, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_13, A_21, Q_13, Q_21, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_29, Q_13, Q_29, TTQRT, col_offset_ttqrt); */
                    break;
                case 18:
                    for (index_t idx_mat_1 = 13, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 13, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_14, A_16, Q_14, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 13, idx_mat_2 = 17; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_14, A_18, Q_14, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_26, Q_22, Q_26, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 13, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_14, A_22, Q_14, Q_22, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_30, Q_14, Q_30, TTQRT, col_offset_ttqrt); */
                    break;
                case 17:
                    for (index_t idx_mat_1 = 14, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 14, idx_mat_2 = 16; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_15, A_17, Q_15, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 14, idx_mat_2 = 18; idx_mat_2 < 27; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_15, A_19, Q_15, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_27, Q_23, Q_27, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 14, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_15, A_23, Q_15, Q_23, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_31, Q_15, Q_31, TTQRT, col_offset_ttqrt); */
                    break;
                case 16:
                    for (index_t idx_mat_1 = 15, idx_mat_2 = 16; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 15, idx_mat_2 = 17; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_16, A_18, Q_16, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 15, idx_mat_2 = 19; idx_mat_2 < 28; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_16, A_20, Q_16, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_28, Q_24, Q_28, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 15, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_16, A_24, Q_16, Q_24, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_32, Q_16, Q_32, TTQRT, col_offset_ttqrt); */
                    break;
                case 15:
                    for (index_t idx_mat_1 = 16, idx_mat_2 = 17; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 16, idx_mat_2 = 18; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 16, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_17, A_21, Q_17, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, Q_25, Q_29, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 16, 24);

                    /* krnl_givens_rotation(A_17, A_25, Q_17, Q_25, TTQRT, col_offset_ttqrt); */
                    break;
                case 14:
                    for (index_t idx_mat_1 = 17, idx_mat_2 = 18; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }
                    /* krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 17, idx_mat_2 = 19; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 17, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_18, A_22, Q_18, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, Q_26, Q_30, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 17, 25);

                    /* krnl_givens_rotation(A_18, A_26, Q_18, Q_26, TTQRT, col_offset_ttqrt); */
                    break;
                case 13:
                    for (index_t idx_mat_1 = 18, idx_mat_2 = 19; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 18, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 18, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_19, A_23, Q_19, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_31, Q_27, Q_31, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 18, 26);

                    /* krnl_givens_rotation(A_19, A_27, Q_19, Q_27, TTQRT, col_offset_ttqrt); */
                    break;
                case 12:
                    for (index_t idx_mat_1 = 19, idx_mat_2 = 20; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 19, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 19, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_20, A_24, Q_20, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, Q_28, Q_32, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 19, 27);

                    /* krnl_givens_rotation(A_20, A_28, Q_20, Q_28, TTQRT, col_offset_ttqrt); */
                    break;
                case 11:
                    for (index_t idx_mat_1 = 20, idx_mat_2 = 21; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 20, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 20, idx_mat_2 = 24; idx_mat_2 < 29; idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_21, A_25, Q_21, Q_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_21, A_29, Q_21, Q_29, TTQRT, col_offset_ttqrt); */
                    break;
                case 10:
                    for (index_t idx_mat_1 = 21, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 21, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 21, idx_mat_2 = 25; idx_mat_2 < 30; idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_22, A_26, Q_22, Q_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_22, A_30, Q_22, Q_30, TTQRT, col_offset_ttqrt); */
                    break;
                case 9:
                    for (index_t idx_mat_1 = 22, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 22, idx_mat_2 = 24; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 22, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_23, A_27, Q_23, Q_27, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_23, A_31, Q_23, Q_31, TTQRT, col_offset_ttqrt); */
                    break;
                case 8:
                    for (index_t idx_mat_1 = 23, idx_mat_2 = 24; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 23, idx_mat_2 = 25; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 23, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_24, A_28, Q_24, Q_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_24, A_32, Q_24, Q_32, TTQRT, col_offset_ttqrt); */
                    break;
                case 7:
                    for (index_t idx_mat_1 = 24, idx_mat_2 = 25; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 24, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 24, 28);

                    /* krnl_givens_rotation(A_25, A_29, Q_25, Q_29, TTQRT, col_offset_ttqrt); */
                    break;
                case 6:
                    for (index_t idx_mat_1 = 25, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 25, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 25, 29);

                    /* krnl_givens_rotation(A_26, A_30, Q_26, Q_30, TTQRT, col_offset_ttqrt); */
                    break;
                case 5:
                    for (index_t idx_mat_1 = 26, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 26, idx_mat_2 = 28; idx_mat_2 < 31; idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_27, A_31, Q_27, Q_31, TTQRT, col_offset_ttqrt); */
                    break;
                case 4:
                    for (index_t idx_mat_1 = 27, idx_mat_2 = 28; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt); */

                    for (index_t idx_mat_1 = 27, idx_mat_2 = 29; idx_mat_2 < 32; idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_28, A_32, Q_28, Q_32, TTQRT, col_offset_ttqrt); */
                    break;
                case 3:
                    for (index_t idx_mat_1 = 28, idx_mat_2 = 29; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 28, 30);

                    /* krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt); */
                    break;
                case 2:
                    for (index_t idx_mat_1 = 29, idx_mat_2 = 30; idx_mat_2 < 32; idx_mat_2++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    /* krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt); */
                    break;
                case 1:
                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 30, 31);

                    /* krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt); */
                    break;
                default:
                    break;
            }
            n_iter_TTQRT--;
            col_offset_ttqrt += TAM_TILED;
        }
    }

write_sol_to_matrix_row_for:
    for (index_t r = 0; r < TAM; r++) {
    	int tile_index = r / TAM_TILED;
		int tile_offset = r % TAM_TILED;
    write_sol_to_matrix_col_for:
        for (index_t c = 0; c < TAM; c++) {
        	A[r][c] = A_tiled[tile_index][tile_offset][c];
        	Q[r][c] = Q_tiled[tile_index][tile_offset][c];
        }
    }

    // Print R matrix
    std::cout << "R Matrix: " << std::endl;
    for (int i = 0; i < TAM; i++) {
        for (int j = 0; j < TAM; j++) {
            std::cout << A[i][j] << "  |  ";
        }
        std::cout << std::endl;
    }

    std::cout << "ECM = " << error(A, out_gold);

    return 0;
}
