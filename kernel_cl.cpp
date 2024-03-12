/**
 * @file kernel.cpp
 * @author Kosta Gecov (kostagecov@gmail.com)
 * @brief Vitis Hls kernel that implements the tiled QRD
 * @version 0.1
 * @date 2023-07-29
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_stream.h>

#include <iostream>

#define FIXED_POINT 24
#define FX_POINT_INT 6

#define TAM_TILED 8
#define TAM 256

#define N_ITER (FIXED_POINT - 1)  // word_lenght - 1

#define GEQRT 0
#define TTQRT 1

/**
 * @brief 24 bits fixed point data, 6 for integer value and 18 for decimals.
 * The more bits it has, more shifts can be performed later, so the approximation to 0 will be more precise.
 * For bigger matrices, we need bigger data formats to be able to calculate the right result.
 *
 */
typedef ap_fixed<FIXED_POINT, FX_POINT_INT, AP_RND> data_t;

// Scale factor to compensate rotations
const data_t SCALE_FACTOR = 0.607252935008881;

// TRIPCOUNT identifiers
const unsigned int N_ELEM_ROW = TAM;
const unsigned int ITER = N_ITER;
const unsigned int TILED_SIZE = TAM_TILED;

class Rotator {
   public:
    // after data is read from an hls::stream<>, it cannot be read again
    // Stream is FIFO type
    hls::stream<data_t, TAM> row_x_in;
    hls::stream<data_t, TAM> row_y_in;
    hls::stream<data_t, TAM> row_x_out;
    hls::stream<data_t, TAM> row_y_out;

    // Position of rows to be rotated
    unsigned int row_x, row_y, col;

   public:
    /**
     * @brief Constructor to initialize rotator with positions
     *
     * @param x
     * @param y
     * @param c
     */
    Rotator(unsigned int x, unsigned int y, unsigned int c);

    /**
     * @brief Perform a givens rotation to two input streams and stores the results in other two streams
     *
     * @param row_x_in
     * @param row_y_in
     * @param row_x_out
     * @param row_y_out
     * @param col_rotator From where the rotation should start
     */
    void givens_rotation(hls::stream<data_t, TAM>& row_x_in,
                         hls::stream<data_t, TAM>& row_y_in,
                         hls::stream<data_t, TAM>& row_x_out,
                         hls::stream<data_t, TAM>& row_y_out,
                         unsigned int col_rotator);
};

/**
 * @brief Read input rows using blocking write to streams
 *
 * @param input pointer to the flattened input matrix
 * @param row_in_1
 * @param row_in_2
 * @param row_in_3
 * @param row_in_4
 * @param row_in_5
 * @param row_in_6
 * @param row_in_7
 * @param row_in_8
 */

void read_input_rows(data_t* input,
                     hls::stream<data_t, TAM>& row_in_1,
                     hls::stream<data_t, TAM>& row_in_2,
                     hls::stream<data_t, TAM>& row_in_3,
                     hls::stream<data_t, TAM>& row_in_4,
                     hls::stream<data_t, TAM>& row_in_5,
                     hls::stream<data_t, TAM>& row_in_6,
                     hls::stream<data_t, TAM>& row_in_7,
                     hls::stream<data_t, TAM>& row_in_8);

void cordic(data_t x[TAM], data_t y[TAM], data_t x_aux[TAM], bool sign, uint8_t n_iter, unsigned int col);

void write_output_rows(data_t* output,
                       hls::stream<data_t, TAM>& row_out_1,
                       hls::stream<data_t, TAM>& row_out_2,
                       hls::stream<data_t, TAM>& row_out_3,
                       hls::stream<data_t, TAM>& row_out_4,
                       hls::stream<data_t, TAM>& row_out_5,
                       hls::stream<data_t, TAM>& row_out_6,
                       hls::stream<data_t, TAM>& row_out_7,
                       hls::stream<data_t, TAM>& row_out_8);

void kernel_givens_rotation_GE(data_t* input_tile_1, data_t* output_tile_1, uint8_t col_offset);

void kernel_givens_rotation_TT(data_t* input_tile_1, data_t* input_tile_2, data_t* output_tile_1, data_t* output_tile_2, uint8_t col_offset);

extern "C" {
/**
 * @brief Performs the givens rotation of the tiles. It is the top function
 *
 * @param A_tile A_tiled matrix
 * @param type_op It can be GEQRT or TTQRT
 * @param col_offset Offset used to avoid reading the positions that had already become 0
 * @param idx_mat_1 To access the right A and Q matrices
 * @param idx_mat_2 To access the right A and Q matrices. In case the operation is GEQRT, this value is 0
 */
void kernel_givens_rotation(data_t* input_tile_1, data_t* input_tile_2,
                            data_t* output_tile_1, data_t* output_tile_2,
                            uint8_t type_op, uint8_t col_offset);
}

Rotator::Rotator(unsigned int x, unsigned int y, unsigned int c) {
#pragma HLS INLINE off
    // actually, row_x and row_y are not used, they have an indicative role
    // while declaring rotators objects
    Rotator::row_x = x;
    Rotator::row_y = y;
    Rotator::col = c;
}

void read_input_rows(data_t* input,
                     hls::stream<data_t, TAM>& row_in_1,
                     hls::stream<data_t, TAM>& row_in_2,
                     hls::stream<data_t, TAM>& row_in_3,
                     hls::stream<data_t, TAM>& row_in_4,
                     hls::stream<data_t, TAM>& row_in_5,
                     hls::stream<data_t, TAM>& row_in_6,
                     hls::stream<data_t, TAM>& row_in_7,
                     hls::stream<data_t, TAM>& row_in_8) {
#pragma HLS INLINE off

read_input_rows_for:
    for (uint16_t j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT avg = N_ELEM_ROW max = N_ELEM_ROW min = N_ELEM_ROW
        row_in_1.write(input[j]);
        row_in_2.write(input[j + 256]);
        row_in_3.write(input[j + 256 * 2]);
        row_in_4.write(input[j + 256 * 3]);
        row_in_5.write(input[j + 256 * 4]);
        row_in_6.write(input[j + 256 * 5]);
        row_in_7.write(input[j + 256 * 6]);
        row_in_8.write(input[j + 256 * 7]);
    }
}

void cordic(data_t x[TAM], data_t y[TAM], data_t x_aux[TAM], bool sign, uint8_t n_iter, unsigned int col) {
#pragma HLS INLINE off
aux_var_for:
	for (uint16_t i = col; i < TAM; i++) {
#pragma HLS LOOP_TRIPCOUNT max = N_ELEM_ROW min = N_ELEM_ROW
		x_aux[i] = x[i];
	}

    // If Y is negative, we need to add to it so that it gets closer to zero
    // and to the contrary with X coordinate
    if (sign) {
    column_rotation_pos_for:
        for (uint16_t j = col; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = N_ELEM_ROW min = TILED_SIZE
            x[j] = x[j] - (y[j] >> n_iter);
            y[j] = y[j] + (x_aux[j] >> n_iter);
        }
    } else {
    column_rotation_neg_for:
        for (uint16_t j = col; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = N_ELEM_ROW min = TILED_SIZE
            x[j] = x[j] + (y[j] >> n_iter);
            y[j] = y[j] - (x_aux[j] >> n_iter);
        }
    }
}

void Rotator::givens_rotation(hls::stream<data_t, TAM>& row_x_in,
                              hls::stream<data_t, TAM>& row_y_in,
                              hls::stream<data_t, TAM>& row_x_out,
                              hls::stream<data_t, TAM>& row_y_out,
                              unsigned int col_rotator) {
#pragma HLS INLINE off
    bool sign = false;
    uint16_t i = 0, j = 0, k = 0, s = 0;
    data_t x[TAM] = {0}, y[TAM] = {0}, x_aux[TAM] = {0};

read_input_data:
    for (j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT avg = N_ELEM_ROW max = N_ELEM_ROW min = N_ELEM_ROW
        row_x_in.read(x[j]);
        row_y_in.read(y[j]);
    }

    // Choose the right sign for the rotation,
    // taking into account the coordinates' quadrants
    if (x[col_rotator] < 0) {
    sign_for:
        for (s = col_rotator; s < TAM; s++) {
#pragma HLS LOOP_TRIPCOUNT max = N_ELEM_ROW min = TILED_SIZE
            x[s] = -x[s];
            y[s] = -y[s];
        }
    }

iterations_for:
    for (k = 0; k < N_ITER; k++) {
#pragma HLS LOOP_TRIPCOUNT max = ITER min = ITER
        sign = (y[col_rotator] < 0);

        cordic(x, y, x_aux, sign, k, col_rotator);
    }

    if ((y[col_rotator] < 0.001) && (y[col_rotator] > -0.001)) {
        y[col_rotator] = 0;
    }

scale_factor_for:
    for (j = col_rotator; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = N_ELEM_ROW min = TILED_SIZE
        x[j] = x[j] * SCALE_FACTOR;
        y[j] = y[j] * SCALE_FACTOR;
    }

write_output_data:
    for (j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = N_ELEM_ROW min = N_ELEM_ROW
        row_x_out.write(x[j]);
        row_y_out.write(y[j]);
    }
}

void write_output_rows(data_t* output,
                       hls::stream<data_t, TAM>& row_out_1,
                       hls::stream<data_t, TAM>& row_out_2,
                       hls::stream<data_t, TAM>& row_out_3,
                       hls::stream<data_t, TAM>& row_out_4,
                       hls::stream<data_t, TAM>& row_out_5,
                       hls::stream<data_t, TAM>& row_out_6,
                       hls::stream<data_t, TAM>& row_out_7,
                       hls::stream<data_t, TAM>& row_out_8) {
#pragma HLS INLINE off

write_output_rows_for:
    for (uint16_t j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT avg = N_ELEM_ROW max = N_ELEM_ROW min = N_ELEM_ROW
        row_out_1.read(output[j]);
        row_out_2.read(output[j + 256]);
        row_out_3.read(output[j + 256 * 2]);
        row_out_4.read(output[j + 256 * 3]);
        row_out_5.read(output[j + 256 * 4]);
        row_out_6.read(output[j + 256 * 5]);
        row_out_7.read(output[j + 256 * 6]);
        row_out_8.read(output[j + 256 * 7]);
    }
}

void kernel_givens_rotation_GE(data_t* input_tile_1, data_t* output_tile_1, uint8_t col_offset) {
    // Rotators for GEQRT operation
    Rotator Rot1_GE(0, 1, 0 /* + col_offset */);
    Rotator Rot2_GE(2, 3, 0);
    Rotator Rot3_GE(4, 5, 0);
    Rotator Rot4_GE(6, 7, 0);

    Rotator Rot5_GE(0, 2, 0);
    Rotator Rot6_GE(4, 6, 0);
    Rotator Rot7_GE(1, 3, 1);
    Rotator Rot8_GE(5, 7, 1);

    Rotator Rot9_GE(0, 4, 0);
    Rotator Rot10_GE(1, 2, 1);
    Rotator Rot11_GE(3, 7, 2);
    Rotator Rot12_GE(5, 6, 1);

    Rotator Rot13_GE(2, 3, 2);
    Rotator Rot14_GE(4, 5, 1);

    Rotator Rot15_GE(1, 4, 1);
    Rotator Rot16_GE(2, 5, 2);
    Rotator Rot17_GE(3, 7, 3);

    Rotator Rot18_GE(3, 5, 3);
    Rotator Rot19_GE(4, 6, 2);

    Rotator Rot20_GE(2, 4, 2);
    Rotator Rot21_GE(3, 6, 3);
    Rotator Rot22_GE(5, 7, 4);

    Rotator Rot23_GE(3, 4, 3);
    Rotator Rot24_GE(5, 6, 4);

    Rotator Rot25_GE(4, 5, 4);
    Rotator Rot26_GE(6, 7, 5);

    Rotator Rot27_GE(5, 6, 5);
    Rotator Rot28_GE(6, 7, 6);

    uint8_t column1 = Rot1_GE.col + col_offset;
    uint8_t column2 = Rot2_GE.col + col_offset;
    uint8_t column3 = Rot3_GE.col + col_offset;
    uint8_t column4 = Rot4_GE.col + col_offset;
    uint8_t column5 = Rot5_GE.col + col_offset;
    uint8_t column6 = Rot6_GE.col + col_offset;
    uint8_t column7 = Rot7_GE.col + col_offset;
    uint8_t column8 = Rot8_GE.col + col_offset;
    uint8_t column9 = Rot9_GE.col + col_offset;
    uint8_t column10 = Rot10_GE.col + col_offset;
    uint8_t column11 = Rot11_GE.col + col_offset;
    uint8_t column12 = Rot12_GE.col + col_offset;
    uint8_t column13 = Rot13_GE.col + col_offset;
    uint8_t column14 = Rot14_GE.col + col_offset;
    uint8_t column15 = Rot15_GE.col + col_offset;
    uint8_t column16 = Rot16_GE.col + col_offset;
    uint8_t column17 = Rot17_GE.col + col_offset;
    uint8_t column18 = Rot18_GE.col + col_offset;
    uint8_t column19 = Rot19_GE.col + col_offset;
    uint8_t column20 = Rot20_GE.col + col_offset;
    uint8_t column21 = Rot21_GE.col + col_offset;
    uint8_t column22 = Rot22_GE.col + col_offset;
    uint8_t column23 = Rot23_GE.col + col_offset;
    uint8_t column24 = Rot24_GE.col + col_offset;
    uint8_t column25 = Rot25_GE.col + col_offset;
    uint8_t column26 = Rot26_GE.col + col_offset;
    uint8_t column27 = Rot27_GE.col + col_offset;
    uint8_t column28 = Rot28_GE.col + col_offset;

#pragma HLS DATAFLOW

    read_input_rows(input_tile_1, Rot1_GE.row_x_in, Rot1_GE.row_y_in, Rot2_GE.row_x_in,
                    Rot2_GE.row_y_in, Rot3_GE.row_x_in, Rot3_GE.row_y_in, Rot4_GE.row_x_in, Rot4_GE.row_y_in);

    // Warning dataflow: rot.col  + col_offset

    Rot1_GE.givens_rotation(Rot1_GE.row_x_in, Rot1_GE.row_y_in, Rot1_GE.row_x_out,
                            Rot1_GE.row_y_out, column1);

    Rot2_GE.givens_rotation(Rot2_GE.row_x_in, Rot2_GE.row_y_in, Rot2_GE.row_x_out,
                            Rot2_GE.row_y_out, column2);

    Rot3_GE.givens_rotation(Rot3_GE.row_x_in, Rot3_GE.row_y_in, Rot3_GE.row_x_out,
                            Rot3_GE.row_y_out, column3);

    Rot4_GE.givens_rotation(Rot4_GE.row_x_in, Rot4_GE.row_y_in, Rot4_GE.row_x_out,
                            Rot4_GE.row_y_out, column4);

    Rot5_GE.givens_rotation(Rot1_GE.row_x_out, Rot2_GE.row_x_out, Rot5_GE.row_x_out,
                            Rot5_GE.row_y_out, column5);

    Rot6_GE.givens_rotation(Rot3_GE.row_x_out, Rot4_GE.row_x_out, Rot6_GE.row_x_out,
                            Rot6_GE.row_y_out, column6);

    Rot7_GE.givens_rotation(Rot1_GE.row_y_out, Rot2_GE.row_y_out, Rot7_GE.row_x_out,
                            Rot7_GE.row_y_out, column7);

    Rot8_GE.givens_rotation(Rot3_GE.row_y_out, Rot4_GE.row_y_out, Rot8_GE.row_x_out,
                            Rot8_GE.row_y_out, column8);

    Rot9_GE.givens_rotation(Rot5_GE.row_x_out, Rot6_GE.row_x_out, Rot9_GE.row_x_out,
                            Rot9_GE.row_y_out, column9);

    Rot10_GE.givens_rotation(Rot7_GE.row_x_out, Rot5_GE.row_y_out, Rot10_GE.row_x_out,
                             Rot10_GE.row_y_out, column10);

    Rot11_GE.givens_rotation(Rot7_GE.row_y_out, Rot8_GE.row_y_out, Rot11_GE.row_x_out,
                             Rot11_GE.row_y_out, column11);

    Rot12_GE.givens_rotation(Rot8_GE.row_x_out, Rot6_GE.row_y_out, Rot12_GE.row_x_out,
                             Rot12_GE.row_y_out, column12);

    Rot13_GE.givens_rotation(Rot10_GE.row_y_out, Rot11_GE.row_x_out, Rot13_GE.row_x_out,
                             Rot13_GE.row_y_out, column13);

    Rot14_GE.givens_rotation(Rot9_GE.row_y_out, Rot12_GE.row_x_out, Rot14_GE.row_x_out,
                             Rot14_GE.row_y_out, column14);

    Rot15_GE.givens_rotation(Rot10_GE.row_x_out, Rot14_GE.row_x_out, Rot15_GE.row_x_out,
                             Rot15_GE.row_y_out, column15);

    Rot16_GE.givens_rotation(Rot13_GE.row_x_out, Rot14_GE.row_y_out, Rot16_GE.row_x_out,
                             Rot16_GE.row_y_out, column16);

    Rot17_GE.givens_rotation(Rot13_GE.row_y_out, Rot11_GE.row_y_out, Rot17_GE.row_x_out,
                             Rot17_GE.row_y_out, column17);

    Rot18_GE.givens_rotation(Rot17_GE.row_x_out, Rot16_GE.row_y_out, Rot18_GE.row_x_out,
                             Rot18_GE.row_y_out, column18);

    Rot19_GE.givens_rotation(Rot15_GE.row_y_out, Rot12_GE.row_y_out, Rot19_GE.row_x_out,
                             Rot19_GE.row_y_out, column19);

    Rot20_GE.givens_rotation(Rot16_GE.row_x_out, Rot19_GE.row_x_out, Rot20_GE.row_x_out,
                             Rot20_GE.row_y_out, column20);

    Rot21_GE.givens_rotation(Rot18_GE.row_x_out, Rot19_GE.row_y_out, Rot21_GE.row_x_out,
                             Rot21_GE.row_y_out, column21);

    Rot22_GE.givens_rotation(Rot18_GE.row_y_out, Rot17_GE.row_y_out, Rot22_GE.row_x_out,
                             Rot22_GE.row_y_out, column22);

    Rot23_GE.givens_rotation(Rot21_GE.row_x_out, Rot20_GE.row_y_out, Rot23_GE.row_x_out,
                             Rot23_GE.row_y_out, column23);

    Rot24_GE.givens_rotation(Rot22_GE.row_x_out, Rot21_GE.row_y_out, Rot24_GE.row_x_out,
                             Rot24_GE.row_y_out, column24);

    Rot25_GE.givens_rotation(Rot23_GE.row_y_out, Rot24_GE.row_x_out, Rot25_GE.row_x_out,
                             Rot25_GE.row_y_out, column25);

    Rot26_GE.givens_rotation(Rot24_GE.row_y_out, Rot22_GE.row_y_out, Rot26_GE.row_x_out,
                             Rot26_GE.row_y_out, column26);

    Rot27_GE.givens_rotation(Rot25_GE.row_y_out, Rot26_GE.row_x_out, Rot27_GE.row_x_out,
                             Rot27_GE.row_y_out, column27);

    Rot28_GE.givens_rotation(Rot27_GE.row_y_out, Rot26_GE.row_y_out, Rot28_GE.row_x_out,
                             Rot28_GE.row_y_out, column28);

    write_output_rows(output_tile_1, Rot9_GE.row_x_out, Rot15_GE.row_x_out, Rot20_GE.row_x_out,
                      Rot23_GE.row_x_out, Rot25_GE.row_x_out, Rot27_GE.row_x_out, Rot28_GE.row_x_out, Rot28_GE.row_y_out);
}

void kernel_givens_rotation_TT(data_t* input_tile_1, data_t* input_tile_2,
                               data_t* output_tile_1, data_t* output_tile_2,
                               uint8_t col_offset) {
    // Rotators for TTQRT operation
    Rotator Rot1_TT(0, 0, 0);
    Rotator Rot2_TT(1, 1, 1);
    Rotator Rot3_TT(2, 2, 2);
    Rotator Rot4_TT(3, 3, 3);
    Rotator Rot5_TT(4, 4, 4);
    Rotator Rot6_TT(5, 5, 5);
    Rotator Rot7_TT(6, 6, 6);
    Rotator Rot8_TT(7, 7, 7);

    Rotator Rot9_TT(1, 0, 1);
    Rotator Rot10_TT(2, 1, 2);
    Rotator Rot11_TT(3, 2, 3);
    Rotator Rot12_TT(4, 3, 4);
    Rotator Rot13_TT(5, 4, 5);
    Rotator Rot14_TT(6, 5, 6);
    Rotator Rot15_TT(7, 6, 7);

    Rotator Rot16_TT(2, 0, 2);
    Rotator Rot17_TT(3, 1, 3);
    Rotator Rot18_TT(4, 2, 4);
    Rotator Rot19_TT(5, 3, 5);
    Rotator Rot20_TT(6, 4, 6);
    Rotator Rot21_TT(7, 5, 7);

    Rotator Rot22_TT(3, 0, 3);
    Rotator Rot23_TT(4, 1, 4);
    Rotator Rot24_TT(5, 2, 5);
    Rotator Rot25_TT(6, 3, 6);
    Rotator Rot26_TT(7, 4, 7);

    Rotator Rot27_TT(4, 0, 4);
    Rotator Rot28_TT(5, 1, 5);
    Rotator Rot29_TT(6, 2, 6);
    Rotator Rot30_TT(7, 3, 7);

    Rotator Rot31_TT(5, 0, 5);
    Rotator Rot32_TT(6, 1, 6);
    Rotator Rot33_TT(7, 2, 7);

    Rotator Rot34_TT(6, 0, 6);
    Rotator Rot35_TT(7, 1, 7);

    Rotator Rot36_TT(7, 0, 7);

    uint8_t column1 = Rot1_TT.col + col_offset;
    uint8_t column2 = Rot2_TT.col + col_offset;
    uint8_t column3 = Rot3_TT.col + col_offset;
    uint8_t column4 = Rot4_TT.col + col_offset;
    uint8_t column5 = Rot5_TT.col + col_offset;
    uint8_t column6 = Rot6_TT.col + col_offset;
    uint8_t column7 = Rot7_TT.col + col_offset;
    uint8_t column8 = Rot8_TT.col + col_offset;
    uint8_t column9 = Rot9_TT.col + col_offset;
    uint8_t column10 = Rot10_TT.col + col_offset;
    uint8_t column11 = Rot11_TT.col + col_offset;
    uint8_t column12 = Rot12_TT.col + col_offset;
    uint8_t column13 = Rot13_TT.col + col_offset;
    uint8_t column14 = Rot14_TT.col + col_offset;
    uint8_t column15 = Rot15_TT.col + col_offset;
    uint8_t column16 = Rot16_TT.col + col_offset;
    uint8_t column17 = Rot17_TT.col + col_offset;
    uint8_t column18 = Rot18_TT.col + col_offset;
    uint8_t column19 = Rot19_TT.col + col_offset;
    uint8_t column20 = Rot20_TT.col + col_offset;
    uint8_t column21 = Rot21_TT.col + col_offset;
    uint8_t column22 = Rot22_TT.col + col_offset;
    uint8_t column23 = Rot23_TT.col + col_offset;
    uint8_t column24 = Rot24_TT.col + col_offset;
    uint8_t column25 = Rot25_TT.col + col_offset;
    uint8_t column26 = Rot26_TT.col + col_offset;
    uint8_t column27 = Rot27_TT.col + col_offset;
    uint8_t column28 = Rot28_TT.col + col_offset;
    uint8_t column29 = Rot29_TT.col + col_offset;
    uint8_t column30 = Rot30_TT.col + col_offset;
    uint8_t column31 = Rot31_TT.col + col_offset;
    uint8_t column32 = Rot32_TT.col + col_offset;
    uint8_t column33 = Rot33_TT.col + col_offset;
    uint8_t column34 = Rot34_TT.col + col_offset;
    uint8_t column35 = Rot35_TT.col + col_offset;
    uint8_t column36 = Rot36_TT.col + col_offset;

#pragma HLS DATAFLOW
    // Read X coordinates rows from first matrix
    read_input_rows(input_tile_1, Rot1_TT.row_x_in, Rot2_TT.row_x_in,
                    Rot3_TT.row_x_in, Rot4_TT.row_x_in, Rot5_TT.row_x_in,
                    Rot6_TT.row_x_in, Rot7_TT.row_x_in, Rot8_TT.row_x_in);

    // Read Y coordinates rows from second matrix
    read_input_rows(input_tile_2, Rot1_TT.row_y_in, Rot2_TT.row_y_in,
                    Rot3_TT.row_y_in, Rot4_TT.row_y_in, Rot5_TT.row_y_in,
                    Rot6_TT.row_y_in, Rot7_TT.row_y_in, Rot8_TT.row_y_in);

    Rot1_TT.givens_rotation(Rot1_TT.row_x_in, Rot1_TT.row_y_in,
                            Rot1_TT.row_x_out, Rot1_TT.row_y_out,
                            column1);

    Rot2_TT.givens_rotation(Rot2_TT.row_x_in, Rot2_TT.row_y_in,
                            Rot2_TT.row_x_out, Rot2_TT.row_y_out,
                            column2);

    Rot3_TT.givens_rotation(Rot3_TT.row_x_in, Rot3_TT.row_y_in,
                            Rot3_TT.row_x_out, Rot3_TT.row_y_out,
                            column3);

    Rot4_TT.givens_rotation(Rot4_TT.row_x_in, Rot4_TT.row_y_in,
                            Rot4_TT.row_x_out, Rot4_TT.row_y_out,
                            column4);

    Rot5_TT.givens_rotation(Rot5_TT.row_x_in, Rot5_TT.row_y_in,
                            Rot5_TT.row_x_out, Rot5_TT.row_y_out,
                            column5);

    Rot6_TT.givens_rotation(Rot6_TT.row_x_in, Rot6_TT.row_y_in,
                            Rot6_TT.row_x_out, Rot6_TT.row_y_out,
                            column6);

    Rot7_TT.givens_rotation(Rot7_TT.row_x_in, Rot7_TT.row_y_in,
                            Rot7_TT.row_x_out, Rot7_TT.row_y_out,
                            column7);

    Rot8_TT.givens_rotation(Rot8_TT.row_x_in, Rot8_TT.row_y_in,
                            Rot8_TT.row_x_out, Rot8_TT.row_y_out,
                            column8);

    Rot9_TT.givens_rotation(Rot2_TT.row_x_out, Rot1_TT.row_y_out,
                            Rot9_TT.row_x_out, Rot9_TT.row_y_out,
                            column9);

    Rot10_TT.givens_rotation(Rot3_TT.row_x_out, Rot2_TT.row_y_out,
                             Rot10_TT.row_x_out, Rot10_TT.row_y_out,
                             column10);

    Rot11_TT.givens_rotation(Rot4_TT.row_x_out, Rot3_TT.row_y_out,
                             Rot11_TT.row_x_out, Rot11_TT.row_y_out,
                             column11);

    Rot12_TT.givens_rotation(Rot5_TT.row_x_out, Rot4_TT.row_y_out,
                             Rot12_TT.row_x_out, Rot12_TT.row_y_out,
                             column12);

    Rot13_TT.givens_rotation(Rot6_TT.row_x_out, Rot5_TT.row_y_out,
                             Rot13_TT.row_x_out, Rot13_TT.row_y_out,
                             column13);

    Rot14_TT.givens_rotation(Rot7_TT.row_x_out, Rot6_TT.row_y_out,
                             Rot14_TT.row_x_out, Rot14_TT.row_y_out,
                             column14);

    Rot15_TT.givens_rotation(Rot8_TT.row_x_out, Rot7_TT.row_y_out,
                             Rot15_TT.row_x_out, Rot15_TT.row_y_out,
                             column15);

    Rot16_TT.givens_rotation(Rot10_TT.row_x_out, Rot9_TT.row_y_out,
                             Rot16_TT.row_x_out, Rot16_TT.row_y_out,
                             column16);

    Rot17_TT.givens_rotation(Rot11_TT.row_x_out, Rot10_TT.row_y_out,
                             Rot17_TT.row_x_out, Rot17_TT.row_y_out,
                             column17);

    Rot18_TT.givens_rotation(Rot12_TT.row_x_out, Rot11_TT.row_y_out,
                             Rot18_TT.row_x_out, Rot18_TT.row_y_out,
                             column18);

    Rot19_TT.givens_rotation(Rot13_TT.row_x_out, Rot12_TT.row_y_out,
                             Rot19_TT.row_x_out, Rot19_TT.row_y_out,
                             column19);

    Rot20_TT.givens_rotation(Rot14_TT.row_x_out, Rot13_TT.row_y_out,
                             Rot20_TT.row_x_out, Rot20_TT.row_y_out,
                             column20);

    Rot21_TT.givens_rotation(Rot15_TT.row_x_out, Rot14_TT.row_y_out,
                             Rot21_TT.row_x_out, Rot21_TT.row_y_out,
                             column21);

    Rot22_TT.givens_rotation(Rot17_TT.row_x_out, Rot16_TT.row_y_out,
                             Rot22_TT.row_x_out, Rot22_TT.row_y_out,
                             column22);

    Rot23_TT.givens_rotation(Rot18_TT.row_x_out, Rot17_TT.row_y_out,
                             Rot23_TT.row_x_out, Rot23_TT.row_y_out,
                             column23);

    Rot24_TT.givens_rotation(Rot19_TT.row_x_out, Rot18_TT.row_y_out,
                             Rot24_TT.row_x_out, Rot24_TT.row_y_out,
                             column24);

    Rot25_TT.givens_rotation(Rot20_TT.row_x_out, Rot19_TT.row_y_out,
                             Rot25_TT.row_x_out, Rot25_TT.row_y_out,
                             column25);

    Rot26_TT.givens_rotation(Rot21_TT.row_x_out, Rot20_TT.row_y_out,
                             Rot26_TT.row_x_out, Rot26_TT.row_y_out,
                             column26);

    Rot27_TT.givens_rotation(Rot23_TT.row_x_out, Rot22_TT.row_y_out,
                             Rot27_TT.row_x_out, Rot27_TT.row_y_out,
                             column27);

    Rot28_TT.givens_rotation(Rot24_TT.row_x_out, Rot23_TT.row_y_out,
                             Rot28_TT.row_x_out, Rot28_TT.row_y_out,
                             column28);

    Rot29_TT.givens_rotation(Rot25_TT.row_x_out, Rot24_TT.row_y_out,
                             Rot29_TT.row_x_out, Rot29_TT.row_y_out,
                             column29);

    Rot30_TT.givens_rotation(Rot26_TT.row_x_out, Rot25_TT.row_y_out,
                             Rot30_TT.row_x_out, Rot30_TT.row_y_out,
                             column30);

    Rot31_TT.givens_rotation(Rot28_TT.row_x_out, Rot27_TT.row_y_out,
                             Rot31_TT.row_x_out, Rot31_TT.row_y_out,
                             column31);

    Rot32_TT.givens_rotation(Rot29_TT.row_x_out, Rot28_TT.row_y_out,
                             Rot32_TT.row_x_out, Rot32_TT.row_y_out,
                             column32);

    Rot33_TT.givens_rotation(Rot30_TT.row_x_out, Rot29_TT.row_y_out,
                             Rot33_TT.row_x_out, Rot33_TT.row_y_out,
                             column33);

    Rot34_TT.givens_rotation(Rot32_TT.row_x_out, Rot31_TT.row_y_out,
                             Rot34_TT.row_x_out, Rot34_TT.row_y_out,
                             column34);

    Rot35_TT.givens_rotation(Rot33_TT.row_x_out, Rot32_TT.row_y_out,
                             Rot35_TT.row_x_out, Rot35_TT.row_y_out,
                             column35);

    Rot36_TT.givens_rotation(Rot35_TT.row_x_out, Rot34_TT.row_y_out,
                             Rot36_TT.row_x_out, Rot36_TT.row_y_out,
                             column36);

    write_output_rows(output_tile_1, Rot1_TT.row_x_out, Rot9_TT.row_x_out, Rot16_TT.row_x_out,
                      Rot22_TT.row_x_out, Rot27_TT.row_x_out, Rot31_TT.row_x_out, Rot34_TT.row_x_out, Rot36_TT.row_x_out);

    write_output_rows(output_tile_2, Rot36_TT.row_y_out, Rot35_TT.row_y_out, Rot33_TT.row_y_out,
                      Rot30_TT.row_y_out, Rot26_TT.row_y_out, Rot21_TT.row_y_out, Rot15_TT.row_y_out, Rot8_TT.row_y_out);
}

extern "C" {
void kernel_givens_rotation(data_t* input_tile_1, data_t* input_tile_2,
                            data_t* output_tile_1, data_t* output_tile_2,
                            uint8_t type_op, uint8_t col_offset) {
#pragma HLS INTERFACE mode = m_axi bundle = amem0 port = input_tile_1 offset = slave
#pragma HLS INTERFACE mode = m_axi bundle = amem1 port = input_tile_2 offset = slave
#pragma HLS INTERFACE mode = m_axi bundle = amem2 port = output_tile_1 offset = slave
#pragma HLS INTERFACE mode = m_axi bundle = amem3 port = output_tile_2 offset = slave

    if (type_op == GEQRT) {
        kernel_givens_rotation_GE(input_tile_1, output_tile_1, col_offset);
    } else if (type_op == TTQRT) {
        kernel_givens_rotation_TT(input_tile_1, input_tile_2, output_tile_1, output_tile_2, col_offset);
    }
}
}

