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

/* void traspose(Q[TAM_TILED][TAM]) {
}

void multiply_matrices(A[TAM_TILED][TAM], Q[TAM_TILED][TAM]) {
}

void update(A[TAM_TILED][TAM], Q[TAM_TILED][TAM]) {
    traspose(Q);
    multiply(A, Q);
} */

/* void matrix_mul(data_t A[TAM][TAM], data_t Q[TAM][TAM], data_t A_res[TAM][TAM]) {
    for (index_t i = 0; i < n; i++) {
        for (index_t j = 0; j < n; j++) {
            for (index_t k = 0; k < n; k++) {
                A_res[i][j] += A[i][k] * Q[k][j];
            }
        }
    }
} */

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

    bool sign[N_ITER];

    /**
     * stores all 32 A matrices needed for tiled operations
     *
     */
    data_t A_tiled[NUM_TILED][TAM_TILED][TAM];

    data_t A[TAM][TAM];
    data_t A_res[TAM][TAM];
    data_t A_aux[NUM_TILED][TAM_TILED][TAM];  // Not used
/*  data_t A_1[TAM_TILED][TAM];
    data_t A_2[TAM_TILED][TAM];
    data_t A_3[TAM_TILED][TAM];
    data_t A_4[TAM_TILED][TAM];
    data_t A_5[TAM_TILED][TAM];
    data_t A_6[TAM_TILED][TAM];
    data_t A_7[TAM_TILED][TAM];
    data_t A_8[TAM_TILED][TAM];
    data_t A_9[TAM_TILED][TAM];
    data_t A_10[TAM_TILED][TAM];
    data_t A_11[TAM_TILED][TAM];
    data_t A_12[TAM_TILED][TAM];
    data_t A_13[TAM_TILED][TAM];
    data_t A_14[TAM_TILED][TAM];
    data_t A_15[TAM_TILED][TAM];
    data_t A_16[TAM_TILED][TAM];
    data_t A_17[TAM_TILED][TAM];
    data_t A_18[TAM_TILED][TAM];
    data_t A_19[TAM_TILED][TAM];
    data_t A_20[TAM_TILED][TAM];
    data_t A_21[TAM_TILED][TAM];
    data_t A_22[TAM_TILED][TAM];
    data_t A_23[TAM_TILED][TAM];
    data_t A_24[TAM_TILED][TAM];
    data_t A_25[TAM_TILED][TAM];
    data_t A_26[TAM_TILED][TAM];
    data_t A_27[TAM_TILED][TAM];
    data_t A_28[TAM_TILED][TAM];
    data_t A_29[TAM_TILED][TAM];
    data_t A_30[TAM_TILED][TAM];
    data_t A_31[TAM_TILED][TAM];
    data_t A_32[TAM_TILED][TAM]; */

    /**
     * stores all 32 Q matrices needed for tiled operations
     *
     */
    data_t Q_tiled[NUM_TILED][TAM_TILED][TAM];

    data_t Q[TAM][TAM];
    data_t Q_aux[NUM_TILED][TAM_TILED][TAM];  // Not used
/*  data_t Q_1[TAM_TILED][TAM];
    data_t Q_2[TAM_TILED][TAM];
    data_t Q_3[TAM_TILED][TAM];
    data_t Q_4[TAM_TILED][TAM];
    data_t Q_5[TAM_TILED][TAM];
    data_t Q_6[TAM_TILED][TAM];
    data_t Q_7[TAM_TILED][TAM];
    data_t Q_8[TAM_TILED][TAM];
    data_t Q_9[TAM_TILED][TAM];
    data_t Q_10[TAM_TILED][TAM];
    data_t Q_11[TAM_TILED][TAM];
    data_t Q_12[TAM_TILED][TAM];
    data_t Q_13[TAM_TILED][TAM];
    data_t Q_14[TAM_TILED][TAM];
    data_t Q_15[TAM_TILED][TAM];
    data_t Q_16[TAM_TILED][TAM];
    data_t Q_17[TAM_TILED][TAM];
    data_t Q_18[TAM_TILED][TAM];
    data_t Q_19[TAM_TILED][TAM];
    data_t Q_20[TAM_TILED][TAM];
    data_t Q_21[TAM_TILED][TAM];
    data_t Q_22[TAM_TILED][TAM];
    data_t Q_23[TAM_TILED][TAM];
    data_t Q_24[TAM_TILED][TAM];
    data_t Q_25[TAM_TILED][TAM];
    data_t Q_26[TAM_TILED][TAM];
    data_t Q_27[TAM_TILED][TAM];
    data_t Q_28[TAM_TILED][TAM];
    data_t Q_29[TAM_TILED][TAM];
    data_t Q_30[TAM_TILED][TAM];
    data_t Q_31[TAM_TILED][TAM];
    data_t Q_32[TAM_TILED][TAM]; */

    std::ifstream data_in("data_in.dat");

    if (!data_in.is_open()) {
        std::cerr << "Could not open data_in.dat" << std::endl;
        return -1;
    }

initialize_matrices:
    for (index_t r = 0; r < TAM; r++) {
        for (index_t c = 0; c < TAM; c++) {
            data_in >> A[r][c];
            if (r == c) {
                Q[r][c] = 1;
            } else {
                Q[r][c] = 0;
            }
        }
    }
    data_in.close();

divide_matrices_row_for:
    for (index_t r = 0; r < TAM; r++) {
    divide_matrices_col_for:
        for (index_t c = 0; c < TAM; c++) {
            if (r >= 0 && r < TAM_TILED) {
                A_tiled[1 - 1][r][c] = A[r][c];
                Q_tiled[1 - 1][r][c] = Q[r][c];
            } else if (r >= TAM_TILED && r < TAM_TILED * 2) {
                A_tiled[2 - 1][r - TAM_TILED][c] = A[r][c];
                Q_tiled[2 - 1][r - TAM_TILED][c] = Q[r][c];
            } else if (r >= TAM_TILED * 2 && r < TAM_TILED * 3) {
                A_tiled[3 - 1][r - TAM_TILED * 2][c] = A[r][c];
                Q_tiled[3 - 1][r - TAM_TILED * 2][c] = Q[r][c];
            } else if (r >= TAM_TILED * 3 && r < TAM_TILED * 4) {
                A_tiled[4 - 1][r - TAM_TILED * 3][c] = A[r][c];
                Q_tiled[4 - 1][r - TAM_TILED * 3][c] = Q[r][c];
            } else if (r >= TAM_TILED * 4 && r < TAM_TILED * 5) {
                A_tiled[5 - 1][r - TAM_TILED * 4][c] = A[r][c];
                Q_tiled[5 - 1][r - TAM_TILED * 4][c] = Q[r][c];
            } else if (r >= TAM_TILED * 5 && r < TAM_TILED * 6) {
                A_tiled[6 - 1][r - TAM_TILED * 5][c] = A[r][c];
                Q_tiled[6 - 1][r - TAM_TILED * 5][c] = Q[r][c];
            } else if (r >= TAM_TILED * 6 && r < TAM_TILED * 7) {
                A_tiled[7 - 1][r - TAM_TILED * 6][c] = A[r][c];
                Q_tiled[7 - 1][r - TAM_TILED * 6][c] = Q[r][c];
            } else if (r >= TAM_TILED * 7 && r < TAM_TILED * 8) {
                A_tiled[8 - 1][r - TAM_TILED * 7][c] = A[r][c];
                Q_tiled[8 - 1][r - TAM_TILED * 7][c] = Q[r][c];
            } else if (r >= TAM_TILED * 8 && r < TAM_TILED * 9) {
                A_tiled[9 - 1][r - TAM_TILED * 8][c] = A[r][c];
                Q_tiled[9 - 1][r - TAM_TILED * 8][c] = Q[r][c];
            } else if (r >= TAM_TILED * 9 && r < TAM_TILED * 10) {
                A_tiled[10 - 1][r - TAM_TILED * 9][c] = A[r][c];
                Q_tiled[10 - 1][r - TAM_TILED * 9][c] = Q[r][c];
            } else if (r >= TAM_TILED * 10 && r < TAM_TILED * 11) {
                A_tiled[11 - 1][r - TAM_TILED * 10][c] = A[r][c];
                Q_tiled[11 - 1][r - TAM_TILED * 10][c] = Q[r][c];
            } else if (r >= TAM_TILED * 11 && r < TAM_TILED * 12) {
                A_tiled[12 - 1][r - TAM_TILED * 11][c] = A[r][c];
                Q_tiled[12 - 1][r - TAM_TILED * 11][c] = Q[r][c];
            } else if (r >= TAM_TILED * 12 && r < TAM_TILED * 13) {
                A_tiled[13 - 1][r - TAM_TILED * 12][c] = A[r][c];
                Q_tiled[13 - 1][r - TAM_TILED * 12][c] = Q[r][c];
            } else if (r >= TAM_TILED * 13 && r < TAM_TILED * 14) {
                A_tiled[14 - 1][r - TAM_TILED * 13][c] = A[r][c];
                Q_tiled[14 - 1][r - TAM_TILED * 13][c] = Q[r][c];
            } else if (r >= TAM_TILED * 14 && r < TAM_TILED * 15) {
                A_tiled[15 - 1][r - TAM_TILED * 14][c] = A[r][c];
                Q_tiled[15 - 1][r - TAM_TILED * 14][c] = Q[r][c];
            } else if (r >= TAM_TILED * 15 && r < TAM_TILED * 16) {
                A_tiled[16 - 1][r - TAM_TILED * 15][c] = A[r][c];
                Q_tiled[16 - 1][r - TAM_TILED * 15][c] = Q[r][c];
            } else if (r >= TAM_TILED * 16 && r < TAM_TILED * 17) {
                A_tiled[17 - 1][r - TAM_TILED * 16][c] = A[r][c];
                Q_tiled[17 - 1][r - TAM_TILED * 16][c] = Q[r][c];
            } else if (r >= TAM_TILED * 17 && r < TAM_TILED * 18) {
                A_tiled[18 - 1][r - TAM_TILED * 17][c] = A[r][c];
                Q_tiled[18 - 1][r - TAM_TILED * 17][c] = Q[r][c];
            } else if (r >= TAM_TILED * 18 && r < TAM_TILED * 19) {
                A_tiled[19 - 1][r - TAM_TILED * 18][c] = A[r][c];
                Q_tiled[19 - 1][r - TAM_TILED * 18][c] = Q[r][c];
            } else if (r >= TAM_TILED * 19 && r < TAM_TILED * 20) {
                A_tiled[20 - 1][r - TAM_TILED * 19][c] = A[r][c];
                Q_tiled[20 - 1][r - TAM_TILED * 19][c] = Q[r][c];
            } else if (r >= TAM_TILED * 20 && r < TAM_TILED * 21) {
                A_tiled[21 - 1][r - TAM_TILED * 20][c] = A[r][c];
                Q_tiled[21 - 1][r - TAM_TILED * 20][c] = Q[r][c];
            } else if (r >= TAM_TILED * 21 && r < TAM_TILED * 22) {
                A_tiled[22 - 1][r - TAM_TILED * 21][c] = A[r][c];
                Q_tiled[22 - 1][r - TAM_TILED * 21][c] = Q[r][c];
            } else if (r >= TAM_TILED * 22 && r < TAM_TILED * 23) {
                A_tiled[23 - 1][r - TAM_TILED * 22][c] = A[r][c];
                Q_tiled[23 - 1][r - TAM_TILED * 22][c] = Q[r][c];
            } else if (r >= TAM_TILED * 23 && r < TAM_TILED * 24) {
                A_tiled[24 - 1][r - TAM_TILED * 23][c] = A[r][c];
                Q_tiled[24 - 1][r - TAM_TILED * 23][c] = Q[r][c];
            } else if (r >= TAM_TILED * 24 && r < TAM_TILED * 25) {
                A_tiled[25 - 1][r - TAM_TILED * 24][c] = A[r][c];
                Q_tiled[25 - 1][r - TAM_TILED * 24][c] = Q[r][c];
            } else if (r >= TAM_TILED * 25 && r < TAM_TILED * 26) {
                A_tiled[26 - 1][r - TAM_TILED * 25][c] = A[r][c];
                Q_tiled[26 - 1][r - TAM_TILED * 25][c] = Q[r][c];
            } else if (r >= TAM_TILED * 26 && r < TAM_TILED * 27) {
                A_tiled[27 - 1][r - TAM_TILED * 26][c] = A[r][c];
                Q_tiled[27 - 1][r - TAM_TILED * 26][c] = Q[r][c];
            } else if (r >= TAM_TILED * 27 && r < TAM_TILED * 28) {
                A_tiled[28 - 1][r - TAM_TILED * 27][c] = A[r][c];
                Q_tiled[28 - 1][r - TAM_TILED * 27][c] = Q[r][c];
            } else if (r >= TAM_TILED * 28 && r < TAM_TILED * 29) {
                A_tiled[29 - 1][r - TAM_TILED * 28][c] = A[r][c];
                Q_tiled[29 - 1][r - TAM_TILED * 28][c] = Q[r][c];
            } else if (r >= TAM_TILED * 29 && r < TAM_TILED * 30) {
                A_tiled[30 - 1][r - TAM_TILED * 29][c] = A[r][c];
                Q_tiled[30 - 1][r - TAM_TILED * 29][c] = Q[r][c];
            } else if (r >= TAM_TILED * 30 && r < TAM_TILED * 31) {
                A_tiled[31 - 1][r - TAM_TILED * 30][c] = A[r][c];
                Q_tiled[31 - 1][r - TAM_TILED * 30][c] = Q[r][c];
            } else /*  if (r >= TAM_TILED * 31 && r < TAM_TILED * 32) */ {
                A_tiled[32 - 1][r - TAM_TILED * 31][c] = A[r][c];
                Q_tiled[32 - 1][r - TAM_TILED * 31][c] = Q[r][c];
            }
        }
    }

num_operations_for:
    for (index_t i = 2; i < NUM_OPERACIONES; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            switch (n_iter_GEQRT) {
                case 32:
                    for (index_t n_mat = 1; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /*                     krnl_givens_rotation(A_1, A_aux, Q_1, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_2, A_aux, Q_2, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_3, A_aux, Q_3, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_4, A_aux, Q_4, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_5, A_aux, Q_5, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_6, A_aux, Q_6, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_7, A_aux, Q_7, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_8, A_aux, Q_8, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_9, A_aux, Q_9, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                                        krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 31:
                    for (index_t n_mat = 2; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_2, A_aux, Q_2, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_3, A_aux, Q_3, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, A_aux, Q_4, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5, A_aux, Q_5, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, Q_6, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, Q_7, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, Q_8, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, Q_9, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 30:
                    for (index_t n_mat = 3; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_3, A_aux, Q_3, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, A_aux, Q_4, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5, A_aux, Q_5, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, Q_6, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, Q_7, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, Q_8, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, Q_9, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 29:
                    for (index_t n_mat = 4; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_4, A_aux, Q_4, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5, A_aux, Q_5, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, Q_6, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, Q_7, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, Q_8, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, Q_9, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 28:
                    for (index_t n_mat = 5; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_5, A_aux, Q_5, A_aux, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, Q_6, A_aux, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, Q_7, A_aux, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, Q_8, A_aux, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, Q_9, A_aux, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 27:
                    for (index_t n_mat = 6; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_6, A_aux, Q_6, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, Q_7, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, Q_8, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, Q_9, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 26:
                    for (index_t n_mat = 7; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_7, A_aux, Q_7, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, Q_8, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, Q_9, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 25:
                    for (index_t n_mat = 8; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_8, A_aux, Q_8, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, Q_9, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 24:
                    for (index_t n_mat = 9; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_9, A_aux, Q_9, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 23:
                    for (index_t n_mat = 10; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_10, A_aux, Q_10, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 22:
                    for (index_t n_mat = 11; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_11, A_aux, Q_11, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 21:
                    for (index_t n_mat = 12; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_12, A_aux, Q_12, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 20:
                    for (index_t n_mat = 13; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_13, A_aux, Q_13, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 19:
                    for (index_t n_mat = 14; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_14, A_aux, Q_14, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 18:
                    for (index_t n_mat = 15; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_15, A_aux, Q_15, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 17:
                    for (index_t n_mat = 16; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_16, A_aux, Q_16, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 16:
                    for (index_t n_mat = 17; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_17, A_aux, Q_17, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 15:
                    for (index_t n_mat = 18; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_18, A_aux, Q_18, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 14:
                    for (index_t n_mat = 19; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_19, A_aux, Q_19, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 13:
                    for (index_t n_mat = 20; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_20, A_aux, Q_20, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 12:
                    for (index_t n_mat = 21; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_21, A_aux, Q_21, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 11:
                    for (index_t n_mat = 22; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_22, A_aux, Q_22, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 10:
                    for (index_t n_mat = 23; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_23, A_aux, Q_23, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 9:
                    for (index_t n_mat = 24; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_24, A_aux, Q_24, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 8:
                    for (index_t n_mat = 25; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_25, A_aux, Q_25, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 7:
                    for (index_t n_mat = 26; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_26, A_aux, Q_26, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 6:
                    for (index_t n_mat = 27; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_27, A_aux, Q_27, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 5:
                    for (index_t n_mat = 28; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_28, A_aux, Q_28, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 4:
                    for (index_t n_mat = 29; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_29, A_aux, Q_29, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 3:
                    for (index_t n_mat = 30; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_30, A_aux, Q_30, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 2:
                    for (index_t n_mat = 31; i < 33; i++) {
                        krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);
                    }
                    /* krnl_givens_rotation(A_31, A_aux, Q_31, Q_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                case 1:
                    n_mat = 32;
                    krnl_givens_rotation(A_tiled, A_aux, Q_tiled, Q_aux, GEQRT, col_offset_geqrt, n_mat);

                    /* krnl_givens_rotation(A_32, A_aux, Q_32, Q_aux, GEQRT, col_offset_geqrt); */
                    break;
                default:
                    break;
            }
            n_iter_GEQRT--;
            col_offset_geqrt += TAM_TILED;

        } else {
            switch (n_iter_TTQRT) {
                case 31:
                    krnl_givens_rotation(A_1, A_2, Q_1, Q_2, TTQRT, col_offset_ttqrt);
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
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_3, Q_1, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_5, A_7, Q_5, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_11, Q_9, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, Q_13, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_5, Q_1, Q_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_13, Q_9, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_21, Q_17, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, Q_25, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_9, Q_1, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_25, Q_17, Q_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_17, Q_1, Q_17, TTQRT, col_offset_ttqrt);
                    break;
                case 30:
                    krnl_givens_rotation(A_2, A_3, Q_2, Q_3, TTQRT, col_offset_ttqrt);
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
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_4, Q_2, Q_4, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_6, A_8, Q_6, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_12, Q_10, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, Q_14, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_6, Q_2, Q_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_14, Q_10, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_22, Q_18, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, Q_26, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_10, Q_2, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_26, Q_18, Q_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_18, Q_2, Q_18, TTQRT, col_offset_ttqrt);
                    break;
                case 29:
                    krnl_givens_rotation(A_3, A_4, Q_3, Q_4, TTQRT, col_offset_ttqrt);
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
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_3, A_5, Q_3, Q_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_9, Q_7, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_13, Q_11, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, Q_15, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_3, A_7, Q_3, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_15, Q_11, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_23, Q_19, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_31, Q_27, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_3, A_11, Q_3, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_27, Q_19, Q_27, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_3, A_19, Q_3, Q_19, TTQRT, col_offset_ttqrt);
                    break;
                case 28:
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
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_4, A_6, Q_4, Q_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_10, Q_8, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_14, Q_12, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, Q_16, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_4, A_8, Q_4, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_16, Q_12, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_24, Q_20, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, Q_28, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_4, A_12, Q_4, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_28, Q_20, Q_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_4, A_20, Q_4, Q_20, TTQRT, col_offset_ttqrt);
                    break;
                case 27:
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
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_5, A_7, Q_5, Q_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_11, Q_9, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, Q_13, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_5, A_9, Q_5, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_17, Q_13, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_25, Q_21, Q_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_5, A_13, Q_5, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_29, Q_21, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_5, A_21, Q_5, Q_21, TTQRT, col_offset_ttqrt);
                    break;
                case 26:
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
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_6, A_8, Q_6, Q_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_12, Q_10, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, Q_14, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_6, A_10, Q_6, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_18, Q_14, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_26, Q_22, Q_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_6, A_14, Q_6, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_30, Q_22, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_6, A_22, Q_6, Q_22, TTQRT, col_offset_ttqrt);
                    break;
                case 25:
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
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_7, A_9, Q_7, Q_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_13, Q_11, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, Q_15, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_7, A_11, Q_7, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_19, Q_15, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_27, Q_23, Q_27, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_7, A_15, Q_7, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_31, Q_23, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_7, A_23, Q_7, Q_23, TTQRT, col_offset_ttqrt);
                    break;
                case 24:
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
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_8, A_10, Q_8, Q_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_14, Q_12, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, Q_16, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_8, A_12, Q_8, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_20, Q_16, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_28, Q_24, Q_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_8, A_16, Q_8, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_32, Q_24, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_8, A_24, Q_8, Q_24, TTQRT, col_offset_ttqrt);
                    break;
                case 23:
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
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_11, Q_9, Q_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, Q_13, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_13, Q_9, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_21, Q_17, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, Q_25, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_17, Q_9, Q_17, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_25, Q_9, Q_25, TTQRT, col_offset_ttqrt);
                    break;
                case 22:
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
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_12, Q_10, Q_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, Q_14, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_14, Q_10, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_22, Q_18, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, Q_26, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_18, Q_10, Q_18, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_26, Q_10, Q_26, TTQRT, col_offset_ttqrt);
                    break;
                case 21:
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
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_13, Q_11, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, Q_15, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_15, Q_11, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_23, Q_19, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_31, Q_27, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_19, Q_11, Q_19, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_27, Q_11, Q_27, TTQRT, col_offset_ttqrt);
                    break;
                case 20:
                    krnl_givens_rotation(A_12, A_13, Q_12, Q_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_14, Q_12, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, Q_16, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_16, Q_12, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_24, Q_20, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, Q_28, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_20, Q_12, Q_20, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_28, Q_12, Q_28, TTQRT, col_offset_ttqrt);
                    break;
                case 19:
                    krnl_givens_rotation(A_13, A_14, Q_13, Q_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_15, Q_13, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_17, Q_13, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_25, Q_21, Q_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_21, Q_13, Q_21, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_29, Q_13, Q_29, TTQRT, col_offset_ttqrt);
                    break;
                case 18:
                    krnl_givens_rotation(A_14, A_15, Q_14, Q_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_16, Q_14, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_18, Q_14, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_26, Q_22, Q_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_22, Q_14, Q_22, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_30, Q_14, Q_30, TTQRT, col_offset_ttqrt);
                    break;
                case 17:
                    krnl_givens_rotation(A_15, A_16, Q_15, Q_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_17, Q_15, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_19, Q_15, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_27, Q_23, Q_27, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_23, Q_15, Q_23, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_31, Q_15, Q_31, TTQRT, col_offset_ttqrt);
                    break;
                case 16:
                    krnl_givens_rotation(A_16, A_17, Q_16, Q_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_18, Q_16, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_20, Q_16, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_28, Q_24, Q_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_24, Q_16, Q_24, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_32, Q_16, Q_32, TTQRT, col_offset_ttqrt);
                    break;
                case 15:
                    krnl_givens_rotation(A_17, A_18, Q_17, Q_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_17, A_19, Q_17, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_17, A_21, Q_17, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, Q_25, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_17, A_25, Q_17, Q_25, TTQRT, col_offset_ttqrt);
                    break;
                case 14:
                    krnl_givens_rotation(A_18, A_19, Q_18, Q_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_18, A_20, Q_18, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_18, A_22, Q_18, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, Q_26, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_18, A_26, Q_18, Q_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_18, A_32, Q_18, Q_32, TTQRT, col_offset_ttqrt);
                    break;
                case 13:
                    krnl_givens_rotation(A_19, A_20, Q_19, Q_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_19, A_21, Q_19, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_19, A_23, Q_19, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_31, Q_27, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_19, A_27, Q_19, Q_27, TTQRT, col_offset_ttqrt);
                    break;
                case 12:
                    krnl_givens_rotation(A_20, A_21, Q_20, Q_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_20, A_22, Q_20, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_20, A_24, Q_20, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, Q_28, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_20, A_28, Q_20, Q_28, TTQRT, col_offset_ttqrt);
                    break;
                case 11:
                    krnl_givens_rotation(A_21, A_22, Q_21, Q_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_21, A_23, Q_21, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_21, A_25, Q_21, Q_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_21, A_29, Q_21, Q_29, TTQRT, col_offset_ttqrt);
                    break;
                case 10:
                    krnl_givens_rotation(A_22, A_23, Q_22, Q_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_22, A_24, Q_22, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_22, A_26, Q_22, Q_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_22, A_30, Q_22, Q_30, TTQRT, col_offset_ttqrt);
                    break;
                case 9:
                    krnl_givens_rotation(A_23, A_24, Q_23, Q_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_23, A_25, Q_23, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_23, A_27, Q_23, Q_27, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_23, A_31, Q_23, Q_31, TTQRT, col_offset_ttqrt);
                    break;
                case 8:
                    krnl_givens_rotation(A_24, A_25, Q_24, Q_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_24, A_26, Q_24, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_24, A_28, Q_24, Q_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_24, A_32, Q_24, Q_32, TTQRT, col_offset_ttqrt);
                    break;
                case 7:
                    krnl_givens_rotation(A_25, A_26, Q_25, Q_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_25, A_27, Q_25, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_25, A_29, Q_25, Q_29, TTQRT, col_offset_ttqrt);
                    break;
                case 6:
                    krnl_givens_rotation(A_26, A_27, Q_26, Q_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_26, A_28, Q_26, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_26, A_30, Q_26, Q_30, TTQRT, col_offset_ttqrt);
                    break;
                case 5:
                    krnl_givens_rotation(A_27, A_28, Q_27, Q_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_27, A_29, Q_27, Q_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_27, A_31, Q_27, Q_31, TTQRT, col_offset_ttqrt);
                    break;
                case 4:
                    krnl_givens_rotation(A_28, A_29, Q_28, Q_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_28, A_30, Q_28, Q_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_28, A_32, Q_28, Q_32, TTQRT, col_offset_ttqrt);
                    break;
                case 3:
                    krnl_givens_rotation(A_29, A_30, Q_29, Q_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_29, A_31, Q_29, Q_31, TTQRT, col_offset_ttqrt);
                    break;
                case 2:
                    krnl_givens_rotation(A_30, A_31, Q_30, Q_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_30, A_32, Q_30, Q_32, TTQRT, col_offset_ttqrt);
                    break;
                case 1:
                    krnl_givens_rotation(A_31, A_32, Q_31, Q_32, TTQRT, col_offset_ttqrt);
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
                A[r][c] = A_tiled[1 - 1][r][c];
                Q[r][c] = Q_tiled[1 - 1][r][c];
            } else if (r >= TAM_TILED && r < TAM_TILED * 2) {
                A[r][c] = A_tiled[2 - 1][r - TAM_TILED][c];
                Q[r][c] = Q_tiled[2 - 1][r - TAM_TILED][c];
            } else if (r >= TAM_TILED * 2 && r < TAM_TILED * 3) {
                A[r][c] = A_tiled[3 - 1][r - TAM_TILED * 2][c];
                Q[r][c] = Q_tiled[3 - 1][r - TAM_TILED * 2][c];
            } else if (r >= TAM_TILED * 3 && r < TAM_TILED * 4) {
                A[r][c] = A_tiled[4 - 1][r - TAM_TILED * 3][c];
                Q[r][c] = Q_tiled[4 - 1][r - TAM_TILED * 3][c];
            } else if (r >= TAM_TILED * 4 && r < TAM_TILED * 5) {
                A[r][c] = A_tiled[5 - 1][r - TAM_TILED * 4][c];
                Q[r][c] = Q_tiled[5 - 1][r - TAM_TILED * 4][c];
            } else if (r >= TAM_TILED * 5 && r < TAM_TILED * 6) {
                A[r][c] = A_tiled[6 - 1][r - TAM_TILED * 5][c];
                Q[r][c] = Q_tiled[6 - 1][r - TAM_TILED * 5][c];
            } else if (r >= TAM_TILED * 6 && r < TAM_TILED * 7) {
                A[r][c] = A_tiled[7 - 1][r - TAM_TILED * 6][c];
                Q[r][c] = Q_tiled[7 - 1][r - TAM_TILED * 6][c];
            } else if (r >= TAM_TILED * 7 && r < TAM_TILED * 8) {
                A[r][c] = A_tiled[8 - 1][r - TAM_TILED * 7][c];
                Q[r][c] = Q_tiled[8 - 1][r - TAM_TILED * 7][c];
            } else if (r >= TAM_TILED * 8 && r < TAM_TILED * 9) {
                A[r][c] = A_tiled[9 - 1][r - TAM_TILED * 8][c];
                Q[r][c] = Q_tiled[9 - 1][r - TAM_TILED * 8][c];
            } else if (r >= TAM_TILED * 9 && r < TAM_TILED * 10) {
                A[r][c] = A_tiled[10 - 1][r - TAM_TILED * 9][c];
                Q[r][c] = Q_tiled[10 - 1][r - TAM_TILED * 9][c];
            } else if (r >= TAM_TILED * 10 && r < TAM_TILED * 11) {
                A[r][c] = A_tiled[11 - 1][r - TAM_TILED * 10][c];
                Q[r][c] = Q_tiled[11 - 1][r - TAM_TILED * 10][c];
            } else if (r >= TAM_TILED * 11 && r < TAM_TILED * 12) {
                A[r][c] = A_tiled[12 - 1][r - TAM_TILED * 11][c];
                Q[r][c] = Q_tiled[12 - 1][r - TAM_TILED * 11][c];
            } else if (r >= TAM_TILED * 12 && r < TAM_TILED * 13) {
                A[r][c] = A_tiled[13 - 1][r - TAM_TILED * 12][c];
                Q[r][c] = Q_tiled[13 - 1][r - TAM_TILED * 12][c];
            } else if (r >= TAM_TILED * 13 && r < TAM_TILED * 14) {
                A[r][c] = A_tiled[14 - 1][r - TAM_TILED * 13][c];
                Q[r][c] = Q_tiled[14 - 1][r - TAM_TILED * 13][c];
            } else if (r >= TAM_TILED * 14 && r < TAM_TILED * 15) {
                A[r][c] = A_tiled[15 - 1][r - TAM_TILED * 14][c];
                Q[r][c] = Q_tiled[15 - 1][r - TAM_TILED * 14][c];
            } else if (r >= TAM_TILED * 15 && r < TAM_TILED * 16) {
                A[r][c] = A_tiled[16 - 1][r - TAM_TILED * 15][c];
                Q[r][c] = Q_tiled[16 - 1][r - TAM_TILED * 15][c];
            } else if (r >= TAM_TILED * 16 && r < TAM_TILED * 17) {
                A[r][c] = A_tiled[17 - 1][r - TAM_TILED * 16][c];
                Q[r][c] = Q_tiled[17 - 1][r - TAM_TILED * 16][c];
            } else if (r >= TAM_TILED * 17 && r < TAM_TILED * 18) {
                A[r][c] = A_tiled[18 - 1][r - TAM_TILED * 17][c];
                Q[r][c] = Q_tiled[18 - 1][r - TAM_TILED * 17][c];
            } else if (r >= TAM_TILED * 18 && r < TAM_TILED * 19) {
                A[r][c] = A_tiled[19 - 1][r - TAM_TILED * 18][c];
                Q[r][c] = Q_tiled[19 - 1][r - TAM_TILED * 18][c];
            } else if (r >= TAM_TILED * 19 && r < TAM_TILED * 20) {
                A[r][c] = A_tiled[20 - 1][r - TAM_TILED * 19][c];
                Q[r][c] = Q_tiled[20 - 1][r - TAM_TILED * 19][c];
            } else if (r >= TAM_TILED * 20 && r < TAM_TILED * 21) {
                A[r][c] = A_tiled[21 - 1][r - TAM_TILED * 20][c];
                Q[r][c] = Q_tiled[21 - 1][r - TAM_TILED * 20][c];
            } else if (r >= TAM_TILED * 21 && r < TAM_TILED * 22) {
                A[r][c] = A_tiled[22 - 1][r - TAM_TILED * 21][c];
                Q[r][c] = Q_tiled[22 - 1][r - TAM_TILED * 21][c];
            } else if (r >= TAM_TILED * 22 && r < TAM_TILED * 23) {
                A[r][c] = A_tiled[23 - 1][r - TAM_TILED * 22][c];
                Q[r][c] = Q_tiled[23 - 1][r - TAM_TILED * 22][c];
            } else if (r >= TAM_TILED * 23 && r < TAM_TILED * 24) {
                A[r][c] = A_tiled[24 - 1][r - TAM_TILED * 23][c];
                Q[r][c] = Q_tiled[24 - 1][r - TAM_TILED * 23][c];
            } else if (r >= TAM_TILED * 24 && r < TAM_TILED * 25) {
                A[r][c] = A_tiled[25 - 1][r - TAM_TILED * 24][c];
                Q[r][c] = Q_tiled[25 - 1][r - TAM_TILED * 24][c];
            } else if (r >= TAM_TILED * 25 && r < TAM_TILED * 26) {
                A[r][c] = A_tiled[26 - 1][r - TAM_TILED * 25][c];
                Q[r][c] = Q_tiled[26 - 1][r - TAM_TILED * 25][c];
            } else if (r >= TAM_TILED * 26 && r < TAM_TILED * 27) {
                A[r][c] = A_tiled[27 - 1][r - TAM_TILED * 26][c];
                Q[r][c] = Q_tiled[27 - 1][r - TAM_TILED * 26][c];
            } else if (r >= TAM_TILED * 27 && r < TAM_TILED * 28) {
                A[r][c] = A_tiled[28 - 1][r - TAM_TILED * 27][c];
                Q[r][c] = Q_tiled[28 - 1][r - TAM_TILED * 27][c];
            } else if (r >= TAM_TILED * 28 && r < TAM_TILED * 29) {
                A[r][c] = A_tiled[29 - 1][r - TAM_TILED * 28][c];
                Q[r][c] = Q_tiled[29 - 1][r - TAM_TILED * 28][c];
            } else if (r >= TAM_TILED * 29 && r < TAM_TILED * 30) {
                A[r][c] = A_tiled[30 - 1][r - TAM_TILED * 29][c];
                Q[r][c] = Q_tiled[30 - 1][r - TAM_TILED * 29][c];
            } else if (r >= TAM_TILED * 30 && r < TAM_TILED * 31) {
                A[r][c] = A_tiled[31 - 1][r - TAM_TILED * 30][c];
                Q[r][c] = Q_tiled[31 - 1][r - TAM_TILED * 30][c];
            } else /*  if (r >= TAM_TILED * 31 && r < TAM_TILED * 32) */ {
                A[r][c] = A_tiled[32 - 1][r - TAM_TILED * 31][c];
                Q[r][c] = Q_tiled[32 - 1][r - TAM_TILED * 31][c];
            }
        }
    }

    // Print R matrix
    /*     for (int i = 0; i < TAM; i++) {
            for (int j = 0; j < TAM; j++) {
                std::cout << A[i][j] << "  |  ";
            }
            std::cout << std::endl;
        }
        // Print Q matrix
        for (int i = 0; i < TAM; i++) {
            for (int j = 0; j < TAM; j++) {
                std::cout << Q[i][j] << "  |  ";
            }
            std::cout << std::endl;
        } */

    // matrix_mul(A, Q, A_res);

    // Print A_res matrix
    for (int i = 0; i < TAM; i++) {
        for (int j = 0; j < TAM; j++) {
            std::cout << A_res[i][j] << "  |  ";
        }
        std::cout << std::endl;
    }

    std::ofstream qrd_out("qrd_out.dat");

    if (!qrd_out.is_open()) {
        std::cerr << "Could not open qrd_out.dat" << std::endl;
        return -1;
    }

    for (index_t r = 0; r < TAM; r++) {
        for (index_t c = 0; c < TAM; c++) {
            qrd_out << A_res[r][c];
        }
    }
    qrd_out.close();

    return 0;
}
