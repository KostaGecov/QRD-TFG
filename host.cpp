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
    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = 0; c < TAM; c++) {
            *file >> matrix[r][c];
        }
    }
    file->close();
}

void init_matrix(float matrix[TAM][TAM], std::fstream *file) {
    if (!file->is_open()) {
        std::cerr << "Error opening file" << std::endl;
    } else {
        std::cout << "Opened file" << std::endl;
    }

initialize_matrix:
    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = 0; c < TAM; c++) {
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
    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = r; c < TAM; c++) {
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
    static uint16_t col_offset_geqrt = 0;
    /**
     * offset to access the right column for TTQRT operation
     */
    static uint16_t col_offset_ttqrt = 0;
    /**
     * to control the GEQRT operations in each step
     */
    static uint16_t n_iter_GEQRT = 32;
    /**
     * to control the TTQRT operations in each step
     */
    static uint16_t n_iter_TTQRT = 31;

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
    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = 0; c < TAM; c++) {
            if (r == c) {
                Q[r][c] = 1;
            } else {
                Q[r][c] = 0;
            }
        }
    }

divide_matrices_row_for:
    for (uint16_t r = 0; r < TAM; r++) {
        int tile_index = r / TAM_TILED;
        int tile_offset = r % TAM_TILED;

    divide_matrices_col_for:
        for (uint16_t c = 0; c < TAM; c++) {
            A_tiled[tile_index][tile_offset][c] = A[r][c];
            Q_tiled[tile_index][tile_offset][c] = Q[r][c];
        }
    }

num_operations_for:
    for (uint16_t i = 2; i < NUM_OPERACIONES; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            switch (n_iter_GEQRT) {
                case 32:
                    for (uint16_t idx_mat_1 = 0; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 31:
                    for (uint16_t idx_mat_1 = 1; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 30:
                    for (uint16_t idx_mat_1 = 2; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 29:
                    for (uint16_t idx_mat_1 = 3; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 28:
                    for (uint16_t idx_mat_1 = 4; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 27:
                    for (uint16_t idx_mat_1 = 5; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 26:
                    for (uint16_t idx_mat_1 = 6; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 25:
                    for (uint16_t idx_mat_1 = 7; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 24:
                    for (uint16_t idx_mat_1 = 8; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 23:
                    for (uint16_t idx_mat_1 = 9; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 22:
                    for (uint16_t idx_mat_1 = 10; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 21:
                    for (uint16_t idx_mat_1 = 11; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 20:
                    for (uint16_t idx_mat_1 = 12; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 19:
                    for (uint16_t idx_mat_1 = 13; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 18:
                    for (uint16_t idx_mat_1 = 14; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 17:
                    for (uint16_t idx_mat_1 = 15; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 16:
                    for (uint16_t idx_mat_1 = 16; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 15:
                    for (uint16_t idx_mat_1 = 17; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 14:
                    for (uint16_t idx_mat_1 = 18; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 13:
                    for (uint16_t idx_mat_1 = 19; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 12:
                    for (uint16_t idx_mat_1 = 20; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 11:
                    for (uint16_t idx_mat_1 = 21; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 10:
                    for (uint16_t idx_mat_1 = 22; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 9:
                    for (uint16_t idx_mat_1 = 23; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 8:
                    for (uint16_t idx_mat_1 = 24; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 7:
                    for (uint16_t idx_mat_1 = 25; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 6:
                    for (uint16_t idx_mat_1 = 26; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 5:
                    for (uint16_t idx_mat_1 = 27; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 4:
                    for (uint16_t idx_mat_1 = 28; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 3:
                    for (uint16_t idx_mat_1 = 29; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 2:
                    for (uint16_t idx_mat_1 = 30; idx_mat_1 < 32; idx_mat_1++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0);
                    }

                    break;
                case 1:
                    krnl_givens_rotation(A_tiled, Q_tiled, GEQRT, col_offset_geqrt, 31, 0);

                    break;
                default:
                    break;
            }
            n_iter_GEQRT--;
            col_offset_geqrt += TAM_TILED;

        } else {
            switch (n_iter_TTQRT) {
                case 31:
                    for (uint16_t idx_mat_1 = 0, idx_mat_2 = 1; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 0, idx_mat_2 = 2; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 0, idx_mat_2 = 4; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 0, idx_mat_2 = 8; idx_mat_2 < 25; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 0, 16);

                    break;
                case 30:
                    for (uint16_t idx_mat_1 = 1, idx_mat_2 = 2; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 1, idx_mat_2 = 3; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 1, idx_mat_2 = 5; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 1, idx_mat_2 = 9; idx_mat_2 < 26; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 1, 17);

                    break;
                case 29:
                    for (uint16_t idx_mat_1 = 2, idx_mat_2 = 3; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 2, idx_mat_2 = 4; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 2, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 2, idx_mat_2 = 10; idx_mat_2 < 27; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 2, 18);

                    break;
                case 28:
                    for (uint16_t idx_mat_1 = 3, idx_mat_2 = 4; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 3, idx_mat_2 = 5; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 3, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 3, idx_mat_2 = 11; idx_mat_2 < 28; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 3, 19);

                    break;
                case 27:
                    for (uint16_t idx_mat_1 = 4, idx_mat_2 = 5; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 4, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 4, idx_mat_2 = 8; idx_mat_2 < 25; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 4, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 4, 20);

                    break;
                case 26:
                    for (uint16_t idx_mat_1 = 5, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 5, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 5, idx_mat_2 = 9; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 5, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 5, 21);

                    break;
                case 25:
                    for (uint16_t idx_mat_1 = 6, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 6, idx_mat_2 = 8; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 6, idx_mat_2 = 10; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 6, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 6, 22);
                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 6, 26);

                    break;
                case 24:
                    for (uint16_t idx_mat_1 = 7, idx_mat_2 = 8; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 7, idx_mat_2 = 9; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 7, idx_mat_2 = 11; idx_mat_2 < 28; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 7, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 16, idx_mat_2 += 16) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 7, 23);

                    break;
                case 23:
                    for (uint16_t idx_mat_1 = 8, idx_mat_2 = 9; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 8, idx_mat_2 = 10; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 8, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 8, idx_mat_2 = 16; idx_mat_2 < 25; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 22:
                    for (uint16_t idx_mat_1 = 9, idx_mat_2 = 10; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 9, idx_mat_2 = 11; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 9, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 9, idx_mat_2 = 17; idx_mat_2 < 26; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 21:
                    for (uint16_t idx_mat_1 = 10, idx_mat_2 = 11; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 10, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 10, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 10, idx_mat_2 = 18; idx_mat_2 < 27; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 20:
                    for (uint16_t idx_mat_1 = 11, idx_mat_2 = 12; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 11, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 11, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 11, idx_mat_2 = 19; idx_mat_2 < 28; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 19:
                    for (uint16_t idx_mat_1 = 12, idx_mat_2 = 13; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 12, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 12, idx_mat_2 = 16; idx_mat_2 < 25; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 12, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 18:
                    for (uint16_t idx_mat_1 = 13, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 13, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 13, idx_mat_2 = 17; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 13, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 17:
                    for (uint16_t idx_mat_1 = 14, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 14, idx_mat_2 = 16; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 14, idx_mat_2 = 18; idx_mat_2 < 27; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 14, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 16:
                    for (uint16_t idx_mat_1 = 15, idx_mat_2 = 16; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 15, idx_mat_2 = 17; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 15, idx_mat_2 = 19; idx_mat_2 < 28; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 15, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 15:
                    for (uint16_t idx_mat_1 = 16, idx_mat_2 = 17; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 16, idx_mat_2 = 18; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 16, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 16, 24);

                    break;
                case 14:
                    for (uint16_t idx_mat_1 = 17, idx_mat_2 = 18; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 17, idx_mat_2 = 19; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 17, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 17, 25);

                    break;
                case 13:
                    for (uint16_t idx_mat_1 = 18, idx_mat_2 = 19; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 18, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 18, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 18, 26);

                    break;
                case 12:
                    for (uint16_t idx_mat_1 = 19, idx_mat_2 = 20; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 19, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 19, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 19, 27);

                    break;
                case 11:
                    for (uint16_t idx_mat_1 = 20, idx_mat_2 = 21; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 20, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 20, idx_mat_2 = 24; idx_mat_2 < 29; idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 10:
                    for (uint16_t idx_mat_1 = 21, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 21, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 21, idx_mat_2 = 25; idx_mat_2 < 30; idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 9:
                    for (uint16_t idx_mat_1 = 22, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 22, idx_mat_2 = 24; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 22, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 8:
                    for (uint16_t idx_mat_1 = 23, idx_mat_2 = 24; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 23, idx_mat_2 = 25; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 23, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 7:
                    for (uint16_t idx_mat_1 = 24, idx_mat_2 = 25; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 24, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 24, 28);

                    break;
                case 6:
                    for (uint16_t idx_mat_1 = 25, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 25, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 25, 29);

                    break;
                case 5:
                    for (uint16_t idx_mat_1 = 26, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 26, idx_mat_2 = 28; idx_mat_2 < 31; idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 4:
                    for (uint16_t idx_mat_1 = 27, idx_mat_2 = 28; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    for (uint16_t idx_mat_1 = 27, idx_mat_2 = 29; idx_mat_2 < 32; idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 3:
                    for (uint16_t idx_mat_1 = 28, idx_mat_2 = 29; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 28, 30);

                    break;
                case 2:
                    for (uint16_t idx_mat_1 = 29, idx_mat_2 = 30; idx_mat_2 < 32; idx_mat_2++) {
                        krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2);
                    }

                    break;
                case 1:
                    krnl_givens_rotation(A_tiled, Q_tiled, TTQRT, col_offset_ttqrt, 30, 31);

                    break;
                default:
                    break;
            }
            n_iter_TTQRT--;
            col_offset_ttqrt += TAM_TILED;
        }
    }

write_sol_to_matrix_row_for:
    for (uint16_t r = 0; r < TAM; r++) {
        int tile_index = r / TAM_TILED;
        int tile_offset = r % TAM_TILED;
    write_sol_to_matrix_col_for:
        for (uint16_t c = 0; c < TAM; c++) {
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
