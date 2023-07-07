#include <fstream>
#include <iostream>

#include "kernel.cpp"

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
     * 4 iteraciones iniciales, para controlar el nº de operaciones GEQRT en cada paso
     */
    static index_t n_iter_GEQRT = 32;
    /**
     * 3 iteraciones iniciales, para controlar el nº de operaciones TTQRT en cada paso
     *  para i=2, 4; para i = 4, 3; para i = 6, 2; para i = 8, 1;
     */
    static index_t n_iter_TTQRT = 31;

    data_t A[TAM][TAM];
    data_t A_aux[TAM][TAM];

    data_t A_1[TAM_TILED][TAM];
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
    data_t A_32[TAM_TILED][TAM];

    std::ifstream data_in("data_in.dat");

    if (!data_in.is_open()) {
        std::cerr << "Could not open data_in.dat" << std : endl;
        return -1;
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
            } else if (r >= TAM_TILED * 4 && r < TAM_TILED * 5) {
                A_5[r - TAM_TILED * 4][c] = A[r][c];
            } else if (r >= TAM_TILED * 5 && r < TAM_TILED * 6) {
                A_6[r - TAM_TILED * 5][c] = A[r][c];
            } else if (r >= TAM_TILED * 6 && r < TAM_TILED * 7) {
                A_7[r - TAM_TILED * 6][c] = A[r][c];
            } else if (r >= TAM_TILED * 7 && r < TAM_TILED * 8) {
                A_8[r - TAM_TILED * 7][c] = A[r][c];
            } else if (r >= TAM_TILED * 8 && r < TAM_TILED * 9) {
                A_9[r - TAM_TILED * 8][c] = A[r][c];
            } else if (r >= TAM_TILED * 9 && r < TAM_TILED * 10) {
                A_10[r - TAM_TILED * 9][c] = A[r][c];
            } else if (r >= TAM_TILED * 10 && r < TAM_TILED * 11) {
                A_11[r - TAM_TILED * 10][c] = A[r][c];
            } else if (r >= TAM_TILED * 11 && r < TAM_TILED * 12) {
                A_12[r - TAM_TILED * 11][c] = A[r][c];
            } else if (r >= TAM_TILED * 12 && r < TAM_TILED * 13) {
                A_13[r - TAM_TILED * 12][c] = A[r][c];
            } else if (r >= TAM_TILED * 13 && r < TAM_TILED * 14) {
                A_14[r - TAM_TILED * 13][c] = A[r][c];
            } else if (r >= TAM_TILED * 14 && r < TAM_TILED * 15) {
                A_15[r - TAM_TILED * 14][c] = A[r][c];
            } else if (r >= TAM_TILED * 15 && r < TAM_TILED * 16) {
                A_16[r - TAM_TILED * 15][c] = A[r][c];
            } else if (r >= TAM_TILED * 16 && r < TAM_TILED * 17) {
                A_17[r - TAM_TILED * 16][c] = A[r][c];
            } else if (r >= TAM_TILED * 17 && r < TAM_TILED * 18) {
                A_18[r - TAM_TILED * 17][c] = A[r][c];
            } else if (r >= TAM_TILED * 18 && r < TAM_TILED * 19) {
                A_19[r - TAM_TILED * 18][c] = A[r][c];
            } else if (r >= TAM_TILED * 19 && r < TAM_TILED * 20) {
                A_20[r - TAM_TILED * 19][c] = A[r][c];
            } else if (r >= TAM_TILED * 20 && r < TAM_TILED * 21) {
                A_21[r - TAM_TILED * 20][c] = A[r][c];
            } else if (r >= TAM_TILED * 21 && r < TAM_TILED * 22) {
                A_22[r - TAM_TILED * 21][c] = A[r][c];
            } else if (r >= TAM_TILED * 22 && r < TAM_TILED * 23) {
                A_23[r - TAM_TILED * 22][c] = A[r][c];
            } else if (r >= TAM_TILED * 23 && r < TAM_TILED * 24) {
                A_24[r - TAM_TILED * 23][c] = A[r][c];
            } else if (r >= TAM_TILED * 24 && r < TAM_TILED * 25) {
                A_25[r - TAM_TILED * 24][c] = A[r][c];
            } else if (r >= TAM_TILED * 25 && r < TAM_TILED * 26) {
                A_26[r - TAM_TILED * 25][c] = A[r][c];
            } else if (r >= TAM_TILED * 26 && r < TAM_TILED * 27) {
                A_27[r - TAM_TILED * 26][c] = A[r][c];
            } else if (r >= TAM_TILED * 27 && r < TAM_TILED * 28) {
                A_28[r - TAM_TILED * 27][c] = A[r][c];
            } else if (r >= TAM_TILED * 28 && r < TAM_TILED * 29) {
                A_29[r - TAM_TILED * 28][c] = A[r][c];
            } else if (r >= TAM_TILED * 29 && r < TAM_TILED * 30) {
                A_30[r - TAM_TILED * 29][c] = A[r][c];
            } else if (r >= TAM_TILED * 30 && r < TAM_TILED * 31) {
                A_31[r - TAM_TILED * 30][c] = A[r][c];
            } else /*  if (r >= TAM_TILED * 31 && r < TAM_TILED * 32) */ {
                A_32[r - TAM_TILED * 31][c] = A[r][c];
            }
        }
    }

num_operations_for:
    for (index_t i = 2; i < NUM_OPERACIONES; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            switch (n_iter_GEQRT) {
                case 32:
                    krnl_givens_rotation(A_1, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_2, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_3, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 31:
                    krnl_givens_rotation(A_2, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_3, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 30:
                    krnl_givens_rotation(A_3, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_4, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 29:
                    krnl_givens_rotation(A_4, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_5, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 28:
                    krnl_givens_rotation(A_5, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_6, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 27:
                    krnl_givens_rotation(A_6, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_7, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 26:
                    krnl_givens_rotation(A_7, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_8, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 25:
                    krnl_givens_rotation(A_8, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 24:
                    krnl_givens_rotation(A_9, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 23:
                    krnl_givens_rotation(A_10, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 22:
                    krnl_givens_rotation(A_11, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 21:
                    krnl_givens_rotation(A_12, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 20:
                    krnl_givens_rotation(A_13, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 19:
                    krnl_givens_rotation(A_14, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 18:
                    krnl_givens_rotation(A_15, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 17:
                    krnl_givens_rotation(A_16, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 16:
                    krnl_givens_rotation(A_17, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 15:
                    krnl_givens_rotation(A_18, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 14:
                    krnl_givens_rotation(A_19, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 13:
                    krnl_givens_rotation(A_20, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 12:
                    krnl_givens_rotation(A_21, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 11:
                    krnl_givens_rotation(A_22, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 10:
                    krnl_givens_rotation(A_23, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 9:
                    krnl_givens_rotation(A_24, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 8:
                    krnl_givens_rotation(A_25, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 7:
                    krnl_givens_rotation(A_26, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 6:
                    krnl_givens_rotation(A_27, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 5:
                    krnl_givens_rotation(A_28, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 4:
                    krnl_givens_rotation(A_29, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 3:
                    krnl_givens_rotation(A_30, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 2:
                    krnl_givens_rotation(A_31, A_aux, GEQRT, col_offset_geqrt);
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                case 1:
                    krnl_givens_rotation(A_32, A_aux, GEQRT, col_offset_geqrt);
                    break;
                default:
                    break;
            }
            n_iter_GEQRT--;
            col_offset_geqrt += TAM_TILED;

        } else {
            switch (n_iter_TTQRT) {
                case 31:
                    krnl_givens_rotation(A_1, A_2, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_3, A_4, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_5, A_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_3, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_5, A_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_1, A_17, TTQRT, col_offset_ttqrt);
                    break;
                case 30:
                    krnl_givens_rotation(A_2, A_3, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_4, A_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_6, A_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_4, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_6, A_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_2, A_18, TTQRT, col_offset_ttqrt);
                    break;
                case 29:
                    krnl_givens_rotation(A_3, A_4, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_5, A_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_3, A_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_3, A_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_3, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_27, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_3, A_19, TTQRT, col_offset_ttqrt);
                    break;
                case 28:
                    krnl_givens_rotation(A_4, A_5, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_6, A_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_4, A_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_4, A_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_4, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_4, A_20, TTQRT, col_offset_ttqrt);
                    break;
                case 27:
                    krnl_givens_rotation(A_5, A_6, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_7, A_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_5, A_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_5, A_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_5, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_5, A_21, TTQRT, col_offset_ttqrt);
                    break;
                case 26:
                    krnl_givens_rotation(A_6, A_7, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_8, A_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_6, A_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_6, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_6, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_6, A_22, TTQRT, col_offset_ttqrt);
                    break;
                case 25:
                    krnl_givens_rotation(A_7, A_8, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_9, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_7, A_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_7, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_37, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_7, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_7, A_23, TTQRT, col_offset_ttqrt);
                    break;
                case 24:
                    krnl_givens_rotation(A_8, A_9, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_10, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_8, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_8, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_8, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_8, A_24, TTQRT, col_offset_ttqrt);
                    break;
                case 23:
                    krnl_givens_rotation(A_9, A_10, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_11, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_17, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_9, A_25, TTQRT, col_offset_ttqrt);
                    break;
                case 22:
                    krnl_givens_rotation(A_10, A_11, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_12, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_18, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_10, A_26, TTQRT, col_offset_ttqrt);
                    break;
                case 21:
                    krnl_givens_rotation(A_11, A_12, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_13, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_23, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_19, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_11, A_27, TTQRT, col_offset_ttqrt);
                    break;
                case 20:
                    krnl_givens_rotation(A_12, A_13, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_14, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_20, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_12, A_28, TTQRT, col_offset_ttqrt);
                    break;
                case 19:
                    krnl_givens_rotation(A_13, A_14, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_15, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_21, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_13, A_29, TTQRT, col_offset_ttqrt);
                    break;
                case 18:
                    krnl_givens_rotation(A_14, A_15, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_16, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_22, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_14, A_30, TTQRT, col_offset_ttqrt);
                    break;
                case 17:
                    krnl_givens_rotation(A_15, A_16, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_27, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_23, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_15, A_31, TTQRT, col_offset_ttqrt);
                    break;
                case 16:
                    krnl_givens_rotation(A_16, A_17, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_24, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_16, A_32, TTQRT, col_offset_ttqrt);
                    break;
                case 15:
                    krnl_givens_rotation(A_17, A_18, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_17, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_17, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_17, A_25, TTQRT, col_offset_ttqrt);
                    break;
                case 14:
                    krnl_givens_rotation(A_18, A_19, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_18, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_18, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_18, A_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_18, A_32, TTQRT, col_offset_ttqrt);
                    break;
                case 13:
                    krnl_givens_rotation(A_19, A_20, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_19, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_19, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_19, A_27, TTQRT, col_offset_ttqrt);
                    break;
                case 12:
                    krnl_givens_rotation(A_20, A_21, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_20, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_20, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_20, A_28, TTQRT, col_offset_ttqrt);
                    break;
                case 11:
                    krnl_givens_rotation(A_21, A_22, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_21, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_21, A_25, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_21, A_29, TTQRT, col_offset_ttqrt);
                    break;
                case 10:
                    krnl_givens_rotation(A_22, A_23, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_22, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_22, A_26, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_22, A_30, TTQRT, col_offset_ttqrt);
                    break;
                case 9:
                    krnl_givens_rotation(A_23, A_24, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_23, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_23, A_27, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_23, A_31, TTQRT, col_offset_ttqrt);
                    break;
                case 8:
                    krnl_givens_rotation(A_24, A_25, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_24, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_24, A_28, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_24, A_32, TTQRT, col_offset_ttqrt);
                    break;
                case 7:
                    krnl_givens_rotation(A_25, A_26, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_25, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_25, A_29, TTQRT, col_offset_ttqrt);
                    break;
                case 6:
                    krnl_givens_rotation(A_26, A_27, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_26, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_26, A_30, TTQRT, col_offset_ttqrt);
                    break;
                case 5:
                    krnl_givens_rotation(A_27, A_28, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_27, A_29, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_27, A_31, TTQRT, col_offset_ttqrt);
                    break;
                case 4:
                    krnl_givens_rotation(A_28, A_29, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_28, A_30, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_28, A_32, TTQRT, col_offset_ttqrt);
                    break;
                case 3:
                    krnl_givens_rotation(A_29, A_30, TTQRT, col_offset_ttqrt);
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_29, A_31, TTQRT, col_offset_ttqrt);
                    break;
                case 2:
                    krnl_givens_rotation(A_30, A_31, TTQRT, col_offset_ttqrt);

                    krnl_givens_rotation(A_30, A_32, TTQRT, col_offset_ttqrt);
                    break;
                case 1:
                    krnl_givens_rotation(A_31, A_32, TTQRT, col_offset_ttqrt);
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
            } else if (r >= TAM_TILED * 4 && r < TAM_TILED * 5) {
                A[r][c] = A_5[r - TAM_TILED * 4][c];
            } else if (r >= TAM_TILED * 5 && r < TAM_TILED * 6) {
                A[r][c] = A_6[r - TAM_TILED * 5][c];
            } else if (r >= TAM_TILED * 6 && r < TAM_TILED * 7) {
                A[r][c] = A_7[r - TAM_TILED * 6][c];
            } else if (r >= TAM_TILED * 7 && r < TAM_TILED * 8) {
                A[r][c] = A_8[r - TAM_TILED * 7][c];
            } else if (r >= TAM_TILED * 8 && r < TAM_TILED * 9) {
                A[r][c] = A_9[r - TAM_TILED * 8][c];
            } else if (r >= TAM_TILED * 9 && r < TAM_TILED * 10) {
                A[r][c] = A_10[r - TAM_TILED * 9][c];
            } else if (r >= TAM_TILED * 10 && r < TAM_TILED * 11) {
                A[r][c] = A_11[r - TAM_TILED * 10][c];
            } else if (r >= TAM_TILED * 11 && r < TAM_TILED * 12) {
                A[r][c] = A_12[r - TAM_TILED * 11][c];
            } else if (r >= TAM_TILED * 12 && r < TAM_TILED * 13) {
                A[r][c] = A_13[r - TAM_TILED * 12][c];
            } else if (r >= TAM_TILED * 13 && r < TAM_TILED * 14) {
                A[r][c] = A_14[r - TAM_TILED * 13][c];
            } else if (r >= TAM_TILED * 14 && r < TAM_TILED * 15) {
                A[r][c] = A_15[r - TAM_TILED * 14][c];
            } else if (r >= TAM_TILED * 15 && r < TAM_TILED * 16) {
                A[r][c] = A_16[r - TAM_TILED * 15][c];
            } else if (r >= TAM_TILED * 16 && r < TAM_TILED * 17) {
                A[r][c] = A_17[r - TAM_TILED * 16][c];
            } else if (r >= TAM_TILED * 17 && r < TAM_TILED * 18) {
                A[r][c] = A_18[r - TAM_TILED * 17][c];
            } else if (r >= TAM_TILED * 18 && r < TAM_TILED * 19) {
                A[r][c] = A_19[r - TAM_TILED * 18][c];
            } else if (r >= TAM_TILED * 19 && r < TAM_TILED * 20) {
                A[r][c] = A_20[r - TAM_TILED * 19][c];
            } else if (r >= TAM_TILED * 20 && r < TAM_TILED * 21) {
                A[r][c] = A_21[r - TAM_TILED * 20][c];
            } else if (r >= TAM_TILED * 21 && r < TAM_TILED * 22) {
                A[r][c] = A_22[r - TAM_TILED * 21][c];
            } else if (r >= TAM_TILED * 22 && r < TAM_TILED * 23) {
                A[r][c] = A_23[r - TAM_TILED * 22][c];
            } else if (r >= TAM_TILED * 23 && r < TAM_TILED * 24) {
                A[r][c] = A_24[r - TAM_TILED * 23][c];
            } else if (r >= TAM_TILED * 24 && r < TAM_TILED * 25) {
                A[r][c] = A_25[r - TAM_TILED * 24][c];
            } else if (r >= TAM_TILED * 25 && r < TAM_TILED * 26) {
                A[r][c] = A_26[r - TAM_TILED * 25][c];
            } else if (r >= TAM_TILED * 26 && r < TAM_TILED * 27) {
                A[r][c] = A_27[r - TAM_TILED * 26][c];
            } else if (r >= TAM_TILED * 27 && r < TAM_TILED * 28) {
                A[r][c] = A_28[r - TAM_TILED * 27][c];
            } else if (r >= TAM_TILED * 28 && r < TAM_TILED * 29) {
                A[r][c] = A_29[r - TAM_TILED * 28][c];
            } else if (r >= TAM_TILED * 29 && r < TAM_TILED * 30) {
                A[r][c] = A_30[r - TAM_TILED * 29][c];
            } else if (r >= TAM_TILED * 30 && r < TAM_TILED * 31) {
                A[r][c] = A_31[r - TAM_TILED * 30][c];
            } else /*  if (r >= TAM_TILED * 31 && r < TAM_TILED * 32) */ {
                A[r][c] = A_32[r - TAM_TILED * 31][c];
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
