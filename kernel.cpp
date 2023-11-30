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
#include "kernel.h"

Rotator::Rotator(int x, int y, int c) {
    // actually, row_x and row_y are not used, they have an indicative rol
    // when declaring de rotator objects
    Rotator::row_x = x;
    Rotator::row_y = y;
    Rotator::col = c;
}

void read_input_rows(data_t Matrix[NUM_TILED][TAM_TILED][TAM],
                     index_t idx_mat,
                     hls::stream<data_t, TAM>& row_in_1,
                     hls::stream<data_t, TAM>& row_in_2,
                     hls::stream<data_t, TAM>& row_in_3,
                     hls::stream<data_t, TAM>& row_in_4,
                     hls::stream<data_t, TAM>& row_in_5,
                     hls::stream<data_t, TAM>& row_in_6,
                     hls::stream<data_t, TAM>& row_in_7,
                     hls::stream<data_t, TAM>& row_in_8) {
#pragma HLS INLINE off

// Read the rows from the input array and write them to the streams
read_input_rows_for:
    for (index_t j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT avg = 256 max = 256 min = 256
        row_in_1.write(Matrix[idx_mat][0][j]);
        row_in_2.write(Matrix[idx_mat][1][j]);
        row_in_3.write(Matrix[idx_mat][2][j]);
        row_in_4.write(Matrix[idx_mat][3][j]);
        row_in_5.write(Matrix[idx_mat][4][j]);
        row_in_6.write(Matrix[idx_mat][5][j]);
        row_in_7.write(Matrix[idx_mat][6][j]);
        row_in_8.write(Matrix[idx_mat][7][j]);
    }
}

void Rotator::givens_rotation(hls::stream<data_t, TAM>& row_x_in,
                              hls::stream<data_t, TAM>& row_y_in,
                              hls::stream<data_t, TAM>& row_x_out,
                              hls::stream<data_t, TAM>& row_y_out,
                              hls::stream<data_t, TAM>& q_u_in,
                              hls::stream<data_t, TAM>& q_v_in,
                              hls::stream<data_t, TAM>& q_u_out,
                              hls::stream<data_t, TAM>& q_v_out,
                              int col_rotator) {
#pragma HLS INLINE off
    data_t x[TAM] = {0}, y[TAM] = {0};
    data_t u[TAM] = {0}, v[TAM] = {0};
#pragma HLS ARRAY_PARTITION dim = 1 factor = 32 type = block variable = x
#pragma HLS ARRAY_PARTITION dim = 1 factor = 32 type = block variable = y
    data_t aux = 0;

read_input_data:
    for (index_t j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT avg = 256 max = 256 min = 256
        // if (!row_x_in.empty() && !row_y_in.empty() && !q_u_in.empty() && !q_v_in.empty()) {
        //        x[j] = row_x_in.read();
        row_x_in.read(x[j]);
        //        y[j] = row_y_in.read();
        row_y_in.read(y[j]);
        //        u[j] = q_u_in.read();
        q_u_in.read(u[j]);
        //        v[j] = q_v_in.read();
        q_v_in.read(v[j]);
        // } else {
        // std::cout << "Empty stream" << std::endl;
        // }
    }

    // Choose the right sign for the rotation, taking into account the quadrants
    // of the coordinates
    if (x[col_rotator] < 0) {
    sign_for:
        for (index_t s = col_rotator; s < TAM; s++) {
#pragma HLS LOOP_TRIPCOUNT max = 256 min = 8
            /*if (y[s] >= 0) {
                aux = x[s];
                x[s] = y[s];
                y[s] = -aux;

                aux = u[s];
                u[s] = v[s];
                v[s] = -aux;
            } else {
                aux = x[s];
                x[s] = -y[s];
                y[s] = aux;

                aux = u[s];
                u[s] = -v[s];
                v[s] = aux;
            }*/
            x[s] = -x[s];
            y[s] = -y[s];
            u[s] = -u[s];
            v[s] = -v[s];
        }
    }

iterations_for:
    for (index_t k = 0; k < N_ITER; k++) {
#pragma HLS LOOP_TRIPCOUNT max = 20 min = 20
        data_t x_prev = x[col_rotator];
        data_t u_prev = u[col_rotator];
    column_rotation_for:
        for (index_t j = col_rotator; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = 256 min = 8
#pragma HLS UNROLL factor = 4
            // If Y is negative, we need to add to it so that it gets closer to zero
            // and to the contrary with X coordinate

            if (y[col_rotator] < 0) {
                x[j] = x[j] - (y[j] >> k);
                y[j] = y[j] + (x_prev >> k);

                u[j] = u[j] - (v[j] >> k);
                v[j] = v[j] + (u_prev >> k);

            } else {
                x[j] = x[j] + (y[j] >> k);
                y[j] = y[j] - (x_prev >> k);

                u[j] = u[j] + (v[j] >> k);
                v[j] = v[j] - (u_prev >> k);
            }
        }
    }

    y[col_rotator] = 0;

scale_factor_for:
    for (index_t j = col_rotator; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = 256 min = 8
        x[j] = x[j] * SCALE_FACTOR;
        y[j] = y[j] * SCALE_FACTOR;

        u[j] = u[j] * SCALE_FACTOR;
        v[j] = v[j] * SCALE_FACTOR;
    }

write_output_data:
    for (index_t j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = 256 min = 256
        row_x_out.write(x[j]);
        row_y_out.write(y[j]);

        q_u_out.write(u[j]);
        q_v_out.write(v[j]);
    }
}

extern "C" {
void krnl_givens_rotation(data_t A_tile[NUM_TILED][TAM_TILED][TAM],
                          data_t Q_tile[NUM_TILED][TAM_TILED][TAM],
                          index_t type_op, index_t col_offset,
                          index_t idx_mat_1, index_t idx_mat_2) {
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = A_tile
#pragma HLS ARRAY_PARTITION dim = 2 type = complete variable = A_tile
#pragma HLS ARRAY_PARTITION dim = 3 factor = 32 type = block variable = A_tile
#pragma HLS ARRAY_PARTITION dim = 1 type = complete variable = Q_tile
#pragma HLS ARRAY_PARTITION dim = 2 type = complete variable = Q_tile
#pragma HLS ARRAY_PARTITION dim = 3 factor = 32 type = block variable = Q_tile
#pragma HLS DATAFLOW
#pragma HLS TOP name = krnl_givens_rotation

    if (type_op == GEQRT) {
    GEQRT_OPERATION:
        // Rotators for GEQRT operation
        Rotator Rot1_GE(0, 1, 0);
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

        read_input_rows(A_tile, idx_mat_1, Rot1_GE.row_x_in, Rot1_GE.row_y_in, Rot2_GE.row_x_in,
                        Rot2_GE.row_y_in, Rot3_GE.row_x_in, Rot3_GE.row_y_in, Rot4_GE.row_x_in, Rot4_GE.row_y_in);

        read_input_rows(Q_tile, idx_mat_1, Rot1_GE.q_u_in, Rot1_GE.q_v_in, Rot2_GE.q_u_in,
                        Rot2_GE.q_v_in, Rot3_GE.q_u_in, Rot3_GE.q_v_in, Rot4_GE.q_u_in, Rot4_GE.q_v_in);

        Rot1_GE.givens_rotation(Rot1_GE.row_x_in, Rot1_GE.row_y_in, Rot1_GE.row_x_out,
                                Rot1_GE.row_y_out, Rot1_GE.q_u_in, Rot1_GE.q_v_in, Rot1_GE.q_u_out, Rot1_GE.q_v_out, Rot1_GE.col + col_offset);
        Rot2_GE.givens_rotation(Rot2_GE.row_x_in, Rot2_GE.row_y_in, Rot2_GE.row_x_out,
                                Rot2_GE.row_y_out, Rot2_GE.q_u_in, Rot2_GE.q_v_in, Rot2_GE.q_u_out, Rot2_GE.q_v_out, Rot2_GE.col + col_offset);
        Rot3_GE.givens_rotation(Rot3_GE.row_x_in, Rot3_GE.row_y_in, Rot3_GE.row_x_out,
                                Rot3_GE.row_y_out, Rot3_GE.q_u_in, Rot3_GE.q_v_in, Rot3_GE.q_u_out, Rot3_GE.q_v_out, Rot3_GE.col + col_offset);
        Rot4_GE.givens_rotation(Rot4_GE.row_x_in, Rot4_GE.row_y_in, Rot4_GE.row_x_out,
                                Rot4_GE.row_y_out, Rot4_GE.q_u_in, Rot4_GE.q_v_in, Rot4_GE.q_u_out, Rot4_GE.q_v_out, Rot4_GE.col + col_offset);
        Rot5_GE.givens_rotation(Rot1_GE.row_x_out, Rot2_GE.row_x_out, Rot5_GE.row_x_out,
                                Rot5_GE.row_y_out, Rot1_GE.q_u_out, Rot2_GE.q_u_out, Rot5_GE.q_u_out, Rot5_GE.q_v_out, Rot5_GE.col + col_offset);
        Rot6_GE.givens_rotation(Rot3_GE.row_x_out, Rot4_GE.row_x_out, Rot6_GE.row_x_out,
                                Rot6_GE.row_y_out, Rot3_GE.q_u_out, Rot4_GE.q_u_out, Rot6_GE.q_u_out, Rot6_GE.q_v_out, Rot6_GE.col + col_offset);
        Rot7_GE.givens_rotation(Rot1_GE.row_y_out, Rot2_GE.row_y_out, Rot7_GE.row_x_out,
                                Rot7_GE.row_y_out, Rot1_GE.q_v_out, Rot2_GE.q_v_out, Rot7_GE.q_u_out, Rot7_GE.q_v_out, Rot7_GE.col + col_offset);
        Rot8_GE.givens_rotation(Rot3_GE.row_y_out, Rot4_GE.row_y_out, Rot8_GE.row_x_out,
                                Rot8_GE.row_y_out, Rot3_GE.q_v_out, Rot4_GE.q_v_out, Rot8_GE.q_u_out, Rot8_GE.q_v_out, Rot8_GE.col + col_offset);
        Rot9_GE.givens_rotation(Rot5_GE.row_x_out, Rot6_GE.row_x_out, Rot9_GE.row_x_out,
                                Rot9_GE.row_y_out, Rot5_GE.q_u_out, Rot6_GE.q_u_out, Rot9_GE.q_u_out, Rot9_GE.q_v_out, Rot9_GE.col + col_offset);
        Rot10_GE.givens_rotation(Rot7_GE.row_x_out, Rot5_GE.row_y_out, Rot10_GE.row_x_out,
                                 Rot10_GE.row_y_out, Rot7_GE.q_u_out, Rot5_GE.q_v_out, Rot10_GE.q_u_out, Rot10_GE.q_v_out, Rot10_GE.col + col_offset);
        Rot11_GE.givens_rotation(Rot7_GE.row_y_out, Rot8_GE.row_y_out, Rot11_GE.row_x_out,
                                 Rot11_GE.row_y_out, Rot7_GE.q_v_out, Rot8_GE.q_v_out, Rot11_GE.q_u_out, Rot11_GE.q_v_out, Rot11_GE.col + col_offset);
        Rot12_GE.givens_rotation(Rot8_GE.row_x_out, Rot6_GE.row_y_out, Rot12_GE.row_x_out,
                                 Rot12_GE.row_y_out, Rot8_GE.q_u_out, Rot6_GE.q_v_out, Rot12_GE.q_u_out, Rot12_GE.q_v_out, Rot12_GE.col + col_offset);
        Rot13_GE.givens_rotation(Rot10_GE.row_y_out, Rot11_GE.row_x_out, Rot13_GE.row_x_out,
                                 Rot13_GE.row_y_out, Rot10_GE.q_v_out, Rot11_GE.q_u_out, Rot13_GE.q_u_out, Rot13_GE.q_v_out, Rot13_GE.col + col_offset);
        Rot14_GE.givens_rotation(Rot9_GE.row_y_out, Rot12_GE.row_x_out, Rot14_GE.row_x_out,
                                 Rot14_GE.row_y_out, Rot9_GE.q_v_out, Rot12_GE.q_u_out, Rot14_GE.q_u_out, Rot14_GE.q_v_out, Rot14_GE.col + col_offset);
        Rot15_GE.givens_rotation(Rot10_GE.row_x_out, Rot14_GE.row_x_out, Rot15_GE.row_x_out,
                                 Rot15_GE.row_y_out, Rot10_GE.q_u_out, Rot14_GE.q_u_out, Rot15_GE.q_u_out, Rot15_GE.q_v_out, Rot15_GE.col + col_offset);
        Rot16_GE.givens_rotation(Rot13_GE.row_x_out, Rot14_GE.row_y_out, Rot16_GE.row_x_out,
                                 Rot16_GE.row_y_out, Rot13_GE.q_u_out, Rot14_GE.q_v_out, Rot16_GE.q_u_out, Rot16_GE.q_v_out, Rot16_GE.col + col_offset);
        Rot17_GE.givens_rotation(Rot13_GE.row_y_out, Rot11_GE.row_y_out, Rot17_GE.row_x_out,
                                 Rot17_GE.row_y_out, Rot13_GE.q_v_out, Rot11_GE.q_v_out, Rot17_GE.q_u_out, Rot17_GE.q_v_out, Rot17_GE.col + col_offset);
        Rot18_GE.givens_rotation(Rot17_GE.row_x_out, Rot16_GE.row_y_out, Rot18_GE.row_x_out,
                                 Rot18_GE.row_y_out, Rot17_GE.q_u_out, Rot16_GE.q_v_out, Rot18_GE.q_u_out, Rot18_GE.q_v_out, Rot18_GE.col + col_offset);
        Rot19_GE.givens_rotation(Rot15_GE.row_y_out, Rot12_GE.row_y_out, Rot19_GE.row_x_out,
                                 Rot19_GE.row_y_out, Rot15_GE.q_v_out, Rot12_GE.q_v_out, Rot19_GE.q_u_out, Rot19_GE.q_v_out, Rot19_GE.col + col_offset);
        Rot20_GE.givens_rotation(Rot16_GE.row_x_out, Rot19_GE.row_x_out, Rot20_GE.row_x_out,
                                 Rot20_GE.row_y_out, Rot16_GE.q_u_out, Rot19_GE.q_u_out, Rot20_GE.q_u_out, Rot20_GE.q_v_out, Rot20_GE.col + col_offset);
        Rot21_GE.givens_rotation(Rot18_GE.row_x_out, Rot19_GE.row_y_out, Rot21_GE.row_x_out,
                                 Rot21_GE.row_y_out, Rot18_GE.q_u_out, Rot19_GE.q_v_out, Rot21_GE.q_u_out, Rot21_GE.q_v_out, Rot21_GE.col + col_offset);
        Rot22_GE.givens_rotation(Rot18_GE.row_y_out, Rot17_GE.row_y_out, Rot22_GE.row_x_out,
                                 Rot22_GE.row_y_out, Rot18_GE.q_v_out, Rot17_GE.q_v_out, Rot22_GE.q_u_out, Rot22_GE.q_v_out, Rot22_GE.col + col_offset);
        Rot23_GE.givens_rotation(Rot21_GE.row_x_out, Rot20_GE.row_y_out, Rot23_GE.row_x_out,
                                 Rot23_GE.row_y_out, Rot21_GE.q_u_out, Rot20_GE.q_v_out, Rot23_GE.q_u_out, Rot23_GE.q_v_out, Rot23_GE.col + col_offset);
        Rot24_GE.givens_rotation(Rot22_GE.row_x_out, Rot21_GE.row_y_out, Rot24_GE.row_x_out,
                                 Rot24_GE.row_y_out, Rot22_GE.q_u_out, Rot21_GE.q_v_out, Rot24_GE.q_u_out, Rot24_GE.q_v_out, Rot24_GE.col + col_offset);
        Rot25_GE.givens_rotation(Rot23_GE.row_y_out, Rot24_GE.row_x_out, Rot25_GE.row_x_out,
                                 Rot25_GE.row_y_out, Rot23_GE.q_v_out, Rot24_GE.q_u_out, Rot25_GE.q_u_out, Rot25_GE.q_v_out, Rot25_GE.col + col_offset);
        Rot26_GE.givens_rotation(Rot24_GE.row_y_out, Rot22_GE.row_y_out, Rot26_GE.row_x_out,
                                 Rot26_GE.row_y_out, Rot24_GE.q_v_out, Rot22_GE.q_v_out, Rot26_GE.q_u_out, Rot26_GE.q_v_out, Rot26_GE.col + col_offset);
        Rot27_GE.givens_rotation(Rot25_GE.row_y_out, Rot26_GE.row_x_out, Rot27_GE.row_x_out,
                                 Rot27_GE.row_y_out, Rot25_GE.q_v_out, Rot26_GE.q_u_out, Rot27_GE.q_u_out, Rot27_GE.q_v_out, Rot27_GE.col + col_offset);
        Rot28_GE.givens_rotation(Rot27_GE.row_y_out, Rot26_GE.row_y_out, Rot28_GE.row_x_out,
                                 Rot28_GE.row_y_out, Rot27_GE.q_v_out, Rot26_GE.q_v_out, Rot28_GE.q_u_out, Rot28_GE.q_v_out, Rot28_GE.col + col_offset);

    // Write output streams to matrix A_tile
    write_output_streams_col_for:
        for (index_t c = 0; c < TAM; c++) {
#pragma HLS LOOP_TRIPCOUNT max = 256 min = 256
        write_output_streams_row_for:
            for (index_t r = 0; r < TAM_TILED; r++) {
#pragma HLS LOOP_TRIPCOUNT max = 8 min = 8
                if (r == 0) {
                    Rot9_GE.row_x_out.read(A_tile[idx_mat_1][r][c]);
                    Rot9_GE.q_u_out.read(Q_tile[idx_mat_1][r][c]);
                } else if (r == 1) {
                    Rot15_GE.row_x_out.read(A_tile[idx_mat_1][r][c]);
                    Rot15_GE.q_u_out.read(Q_tile[idx_mat_1][r][c]);
                } else if (r == 2) {
                    Rot20_GE.row_x_out.read(A_tile[idx_mat_1][r][c]);
                    Rot20_GE.q_u_out.read(Q_tile[idx_mat_1][r][c]);
                } else if (r == 3) {
                    Rot23_GE.row_x_out.read(A_tile[idx_mat_1][r][c]);
                    Rot23_GE.q_u_out.read(Q_tile[idx_mat_1][r][c]);
                } else if (r == 4) {
                    Rot25_GE.row_x_out.read(A_tile[idx_mat_1][r][c]);
                    Rot25_GE.q_u_out.read(Q_tile[idx_mat_1][r][c]);
                } else if (r == 5) {
                    Rot27_GE.row_x_out.read(A_tile[idx_mat_1][r][c]);
                    Rot27_GE.q_u_out.read(Q_tile[idx_mat_1][r][c]);
                } else if (r == 6) {
                    Rot28_GE.row_x_out.read(A_tile[idx_mat_1][r][c]);
                    Rot28_GE.q_u_out.read(Q_tile[idx_mat_1][r][c]);
                } else {
                    Rot28_GE.row_y_out.read(A_tile[idx_mat_1][r][c]);
                    Rot28_GE.q_v_out.read(Q_tile[idx_mat_1][r][c]);
                }
            }
        }

    } else if (type_op == TTQRT) {
    TTQRT_OPERATION:
        // Rotators for GEQRT operation
        Rotator Rot1_TT(0, 0, 0);
        Rotator Rot2_TT(1, 1, 1);
        Rotator Rot3_TT(2, 2, 2);
        Rotator Rot4_TT(3, 3, 3);
        Rotator Rot5_TT(4, 4, 4);
        Rotator Rot6_TT(5, 5, 5);
        Rotator Rot7_TT(6, 6, 6);
        Rotator Rot8_TT(7, 7, 7);

        Rotator Rot9_TT(0, 0, 1);
        Rotator Rot10_TT(1, 1, 2);
        Rotator Rot11_TT(2, 2, 3);
        Rotator Rot12_TT(3, 3, 4);
        Rotator Rot13_TT(4, 4, 5);
        Rotator Rot14_TT(5, 5, 6);
        Rotator Rot15_TT(6, 6, 7);

        Rotator Rot16_TT(0, 0, 2);
        Rotator Rot17_TT(1, 1, 3);
        Rotator Rot18_TT(2, 2, 4);
        Rotator Rot19_TT(3, 3, 5);
        Rotator Rot20_TT(4, 4, 6);
        Rotator Rot21_TT(5, 5, 7);

        Rotator Rot22_TT(0, 0, 3);
        Rotator Rot23_TT(1, 1, 4);
        Rotator Rot24_TT(2, 2, 5);
        Rotator Rot25_TT(3, 3, 6);
        Rotator Rot26_TT(4, 4, 7);

        Rotator Rot27_TT(0, 0, 4);
        Rotator Rot28_TT(1, 1, 5);
        Rotator Rot29_TT(2, 2, 6);
        Rotator Rot30_TT(3, 3, 7);

        Rotator Rot31_TT(0, 0, 5);
        Rotator Rot32_TT(1, 1, 6);
        Rotator Rot33_TT(2, 2, 7);

        Rotator Rot34_TT(0, 0, 6);
        Rotator Rot35_TT(1, 1, 7);

        Rotator Rot36_TT(0, 0, 7);

        // Read X coordinates rows from first matrix
        read_input_rows(A_tile, idx_mat_1, Rot1_TT.row_x_in, Rot2_TT.row_x_in,
                        Rot3_TT.row_x_in, Rot4_TT.row_x_in, Rot5_TT.row_x_in,
                        Rot6_TT.row_x_in, Rot7_TT.row_x_in, Rot8_TT.row_x_in);

        // Read Y coordinates rows from second matrix
        read_input_rows(A_tile, idx_mat_2, Rot1_TT.row_y_in, Rot2_TT.row_y_in,
                        Rot3_TT.row_y_in, Rot4_TT.row_y_in, Rot5_TT.row_y_in,
                        Rot6_TT.row_y_in, Rot7_TT.row_y_in, Rot8_TT.row_y_in);

        read_input_rows(Q_tile, idx_mat_1, Rot1_TT.q_u_in, Rot2_TT.q_u_in,
                        Rot3_TT.q_u_in, Rot4_TT.q_u_in, Rot5_TT.q_u_in,
                        Rot6_TT.q_u_in, Rot7_TT.q_u_in, Rot8_TT.q_u_in);

        read_input_rows(Q_tile, idx_mat_2, Rot1_TT.q_v_in, Rot2_TT.q_v_in,
                        Rot3_TT.q_v_in, Rot4_TT.q_v_in, Rot5_TT.q_v_in,
                        Rot6_TT.q_v_in, Rot7_TT.q_v_in, Rot8_TT.q_v_in);

        Rot1_TT.givens_rotation(Rot1_TT.row_x_in, Rot1_TT.row_y_in,
                                Rot1_TT.row_x_out, Rot1_TT.row_y_out,
                                Rot1_TT.q_u_in, Rot1_TT.q_v_in,
                                Rot1_TT.q_u_out, Rot1_TT.q_v_out,
                                Rot1_TT.col + col_offset);

        Rot2_TT.givens_rotation(Rot2_TT.row_x_in, Rot2_TT.row_y_in,
                                Rot2_TT.row_x_out, Rot2_TT.row_y_out,
                                Rot2_TT.q_u_in, Rot2_TT.q_v_in,
                                Rot2_TT.q_u_out, Rot2_TT.q_v_out,
                                Rot2_TT.col + col_offset);

        Rot3_TT.givens_rotation(Rot3_TT.row_x_in, Rot3_TT.row_y_in,
                                Rot3_TT.row_x_out, Rot3_TT.row_y_out,
                                Rot3_TT.q_u_in, Rot3_TT.q_v_in,
                                Rot3_TT.q_u_out, Rot3_TT.q_v_out,
                                Rot3_TT.col + col_offset);
        Rot4_TT.givens_rotation(Rot4_TT.row_x_in, Rot4_TT.row_y_in,
                                Rot4_TT.row_x_out, Rot4_TT.row_y_out,
                                Rot4_TT.q_u_in, Rot4_TT.q_v_in,
                                Rot4_TT.q_u_out, Rot4_TT.q_v_out,
                                Rot4_TT.col + col_offset);
        Rot5_TT.givens_rotation(Rot5_TT.row_x_in, Rot5_TT.row_y_in,
                                Rot5_TT.row_x_out, Rot5_TT.row_y_out,
                                Rot5_TT.q_u_in, Rot5_TT.q_v_in,
                                Rot5_TT.q_u_out, Rot5_TT.q_v_out,
                                Rot5_TT.col + col_offset);
        Rot6_TT.givens_rotation(Rot6_TT.row_x_in, Rot6_TT.row_y_in,
                                Rot6_TT.row_x_out, Rot6_TT.row_y_out,
                                Rot6_TT.q_u_in, Rot6_TT.q_v_in,
                                Rot6_TT.q_u_out, Rot6_TT.q_v_out,
                                Rot6_TT.col + col_offset);
        Rot7_TT.givens_rotation(Rot7_TT.row_x_in, Rot7_TT.row_y_in,
                                Rot7_TT.row_x_out, Rot7_TT.row_y_out,
                                Rot7_TT.q_u_in, Rot7_TT.q_v_in,
                                Rot7_TT.q_u_out, Rot7_TT.q_v_out,
                                Rot7_TT.col + col_offset);
        Rot8_TT.givens_rotation(Rot8_TT.row_x_in, Rot8_TT.row_y_in,
                                Rot8_TT.row_x_out, Rot8_TT.row_y_out,
                                Rot8_TT.q_u_in, Rot8_TT.q_v_in,
                                Rot8_TT.q_u_out, Rot8_TT.q_v_out,
                                Rot8_TT.col + col_offset);
        Rot9_TT.givens_rotation(Rot1_TT.row_x_out, Rot1_TT.row_y_out,
                                Rot9_TT.row_x_out, Rot9_TT.row_y_out,
                                Rot1_TT.q_u_out, Rot1_TT.q_v_out,
                                Rot9_TT.q_u_out, Rot9_TT.q_v_out,
                                Rot9_TT.col + col_offset);
        Rot10_TT.givens_rotation(Rot2_TT.row_x_out, Rot2_TT.row_y_out,
                                 Rot10_TT.row_x_out, Rot10_TT.row_y_out,
                                 Rot2_TT.q_u_out, Rot2_TT.q_v_out,
                                 Rot10_TT.q_u_out, Rot10_TT.q_v_out,
                                 Rot10_TT.col + col_offset);
        Rot11_TT.givens_rotation(Rot3_TT.row_x_out, Rot3_TT.row_y_out,
                                 Rot11_TT.row_x_out, Rot11_TT.row_y_out,
                                 Rot3_TT.q_u_out, Rot3_TT.q_v_out,
                                 Rot11_TT.q_u_out, Rot11_TT.q_v_out,
                                 Rot11_TT.col + col_offset);
        Rot12_TT.givens_rotation(Rot4_TT.row_x_out, Rot4_TT.row_y_out,
                                 Rot12_TT.row_x_out, Rot12_TT.row_y_out,
                                 Rot4_TT.q_u_out, Rot4_TT.q_v_out,
                                 Rot12_TT.q_u_out, Rot12_TT.q_v_out,
                                 Rot12_TT.col + col_offset);
        Rot13_TT.givens_rotation(Rot5_TT.row_x_out, Rot5_TT.row_y_out,
                                 Rot13_TT.row_x_out, Rot13_TT.row_y_out,
                                 Rot5_TT.q_u_out, Rot5_TT.q_v_out,
                                 Rot13_TT.q_u_out, Rot13_TT.q_v_out,
                                 Rot13_TT.col + col_offset);
        Rot14_TT.givens_rotation(Rot6_TT.row_x_out, Rot6_TT.row_y_out,
                                 Rot14_TT.row_x_out, Rot14_TT.row_y_out,
                                 Rot6_TT.q_u_out, Rot6_TT.q_v_out,
                                 Rot14_TT.q_u_out, Rot14_TT.q_v_out,
                                 Rot14_TT.col + col_offset);
        Rot15_TT.givens_rotation(Rot7_TT.row_x_out, Rot7_TT.row_y_out,
                                 Rot15_TT.row_x_out, Rot15_TT.row_y_out,
                                 Rot7_TT.q_u_out, Rot7_TT.q_v_out,
                                 Rot15_TT.q_u_out, Rot15_TT.q_v_out,
                                 Rot15_TT.col + col_offset);
        Rot16_TT.givens_rotation(Rot9_TT.row_x_out, Rot9_TT.row_y_out,
                                 Rot16_TT.row_x_out, Rot16_TT.row_y_out,
                                 Rot9_TT.q_u_out, Rot9_TT.q_v_out,
                                 Rot16_TT.q_u_out, Rot16_TT.q_v_out,
                                 Rot16_TT.col + col_offset);
        Rot17_TT.givens_rotation(Rot10_TT.row_x_out, Rot10_TT.row_y_out,
                                 Rot17_TT.row_x_out, Rot17_TT.row_y_out,
                                 Rot10_TT.q_u_out, Rot10_TT.q_v_out,
                                 Rot17_TT.q_u_out, Rot17_TT.q_v_out,
                                 Rot17_TT.col + col_offset);
        Rot18_TT.givens_rotation(Rot11_TT.row_x_out, Rot11_TT.row_y_out,
                                 Rot18_TT.row_x_out, Rot18_TT.row_y_out,
                                 Rot11_TT.q_u_out, Rot11_TT.q_v_out,
                                 Rot18_TT.q_u_out, Rot18_TT.q_v_out,
                                 Rot18_TT.col + col_offset);
        Rot19_TT.givens_rotation(Rot12_TT.row_x_out, Rot12_TT.row_y_out,
                                 Rot19_TT.row_x_out, Rot19_TT.row_y_out,
                                 Rot12_TT.q_u_out, Rot12_TT.q_v_out,
                                 Rot19_TT.q_u_out, Rot19_TT.q_v_out,
                                 Rot19_TT.col + col_offset);
        Rot20_TT.givens_rotation(Rot13_TT.row_x_out, Rot13_TT.row_y_out,
                                 Rot20_TT.row_x_out, Rot20_TT.row_y_out,
                                 Rot13_TT.q_u_out, Rot13_TT.q_v_out,
                                 Rot20_TT.q_u_out, Rot20_TT.q_v_out,
                                 Rot20_TT.col + col_offset);
        Rot21_TT.givens_rotation(Rot14_TT.row_x_out, Rot14_TT.row_y_out,
                                 Rot21_TT.row_x_out, Rot21_TT.row_y_out,
                                 Rot14_TT.q_u_out, Rot14_TT.q_v_out,
                                 Rot21_TT.q_u_out, Rot21_TT.q_v_out,
                                 Rot21_TT.col + col_offset);
        Rot22_TT.givens_rotation(Rot16_TT.row_x_out, Rot16_TT.row_y_out,
                                 Rot22_TT.row_x_out, Rot22_TT.row_y_out,
                                 Rot16_TT.q_u_out, Rot16_TT.q_v_out,
                                 Rot22_TT.q_u_out, Rot22_TT.q_v_out,
                                 Rot22_TT.col + col_offset);
        Rot23_TT.givens_rotation(Rot17_TT.row_x_out, Rot17_TT.row_y_out,
                                 Rot23_TT.row_x_out, Rot23_TT.row_y_out,
                                 Rot17_TT.q_u_out, Rot17_TT.q_v_out,
                                 Rot23_TT.q_u_out, Rot23_TT.q_v_out,
                                 Rot23_TT.col + col_offset);
        Rot24_TT.givens_rotation(Rot18_TT.row_x_out, Rot18_TT.row_y_out,
                                 Rot24_TT.row_x_out, Rot24_TT.row_y_out,
                                 Rot18_TT.q_u_out, Rot18_TT.q_v_out,
                                 Rot24_TT.q_u_out, Rot24_TT.q_v_out,
                                 Rot24_TT.col + col_offset);
        Rot25_TT.givens_rotation(Rot19_TT.row_x_out, Rot19_TT.row_y_out,
                                 Rot25_TT.row_x_out, Rot25_TT.row_y_out,
                                 Rot19_TT.q_u_out, Rot19_TT.q_v_out,
                                 Rot25_TT.q_u_out, Rot25_TT.q_v_out,
                                 Rot25_TT.col + col_offset);
        Rot26_TT.givens_rotation(Rot20_TT.row_x_out, Rot20_TT.row_y_out,
                                 Rot26_TT.row_x_out, Rot26_TT.row_y_out,
                                 Rot20_TT.q_u_out, Rot20_TT.q_v_out,
                                 Rot26_TT.q_u_out, Rot26_TT.q_v_out,
                                 Rot26_TT.col + col_offset);
        Rot27_TT.givens_rotation(Rot22_TT.row_x_out, Rot22_TT.row_y_out,
                                 Rot27_TT.row_x_out, Rot27_TT.row_y_out,
                                 Rot22_TT.q_u_out, Rot22_TT.q_v_out,
                                 Rot27_TT.q_u_out, Rot27_TT.q_v_out,
                                 Rot27_TT.col + col_offset);
        Rot28_TT.givens_rotation(Rot23_TT.row_x_out, Rot23_TT.row_y_out,
                                 Rot28_TT.row_x_out, Rot28_TT.row_y_out,
                                 Rot23_TT.q_u_out, Rot23_TT.q_v_out,
                                 Rot28_TT.q_u_out, Rot28_TT.q_v_out,
                                 Rot28_TT.col + col_offset);
        Rot29_TT.givens_rotation(Rot24_TT.row_x_out, Rot24_TT.row_y_out,
                                 Rot29_TT.row_x_out, Rot29_TT.row_y_out,
                                 Rot24_TT.q_u_out, Rot24_TT.q_v_out,
                                 Rot29_TT.q_u_out, Rot29_TT.q_v_out,
                                 Rot29_TT.col + col_offset);
        Rot30_TT.givens_rotation(Rot25_TT.row_x_out, Rot25_TT.row_y_out,
                                 Rot30_TT.row_x_out, Rot30_TT.row_y_out,
                                 Rot25_TT.q_u_out, Rot25_TT.q_v_out,
                                 Rot30_TT.q_u_out, Rot30_TT.q_v_out,
                                 Rot30_TT.col + col_offset);
        Rot31_TT.givens_rotation(Rot27_TT.row_x_out, Rot27_TT.row_y_out,
                                 Rot31_TT.row_x_out, Rot31_TT.row_y_out,
                                 Rot27_TT.q_u_out, Rot27_TT.q_v_out,
                                 Rot31_TT.q_u_out, Rot31_TT.q_v_out,
                                 Rot31_TT.col + col_offset);
        Rot32_TT.givens_rotation(Rot28_TT.row_x_out, Rot28_TT.row_y_out,
                                 Rot32_TT.row_x_out, Rot32_TT.row_y_out,
                                 Rot28_TT.q_u_out, Rot28_TT.q_v_out,
                                 Rot32_TT.q_u_out, Rot32_TT.q_v_out,
                                 Rot32_TT.col + col_offset);
        Rot33_TT.givens_rotation(Rot29_TT.row_x_out, Rot29_TT.row_y_out,
                                 Rot33_TT.row_x_out, Rot33_TT.row_y_out,
                                 Rot29_TT.q_u_out, Rot29_TT.q_v_out,
                                 Rot33_TT.q_u_out, Rot33_TT.q_v_out,
                                 Rot33_TT.col + col_offset);
        Rot34_TT.givens_rotation(Rot31_TT.row_x_out, Rot31_TT.row_y_out,
                                 Rot34_TT.row_x_out, Rot34_TT.row_y_out,
                                 Rot31_TT.q_u_out, Rot31_TT.q_v_out,
                                 Rot34_TT.q_u_out, Rot34_TT.q_v_out,
                                 Rot34_TT.col + col_offset);
        Rot35_TT.givens_rotation(Rot32_TT.row_x_out, Rot32_TT.row_y_out,
                                 Rot35_TT.row_x_out, Rot35_TT.row_y_out,
                                 Rot32_TT.q_u_out, Rot32_TT.q_v_out,
                                 Rot35_TT.q_u_out, Rot35_TT.q_v_out,
                                 Rot35_TT.col + col_offset);
        Rot36_TT.givens_rotation(Rot34_TT.row_x_out, Rot34_TT.row_y_out,
                                 Rot36_TT.row_x_out, Rot36_TT.row_y_out,
                                 Rot34_TT.q_u_out, Rot34_TT.q_v_out,
                                 Rot36_TT.q_u_out, Rot36_TT.q_v_out,
                                 Rot36_TT.col + col_offset);

    // Write output streams to matrix A_tile and 0s to A_tiled_2
    write_output_streams_col_TTQRT_for:
        for (index_t c = 0; c < TAM; c++) {
#pragma HLS LOOP_TRIPCOUNT max = 256 min = 256
        write_output_streams_row_TTQRT_for:
            for (index_t r = 0; r < TAM_TILED; r++) {
#pragma HLS LOOP_TRIPCOUNT max = 8 min = 8
                if (r == 0) {
                    A_tile[idx_mat_1][r][c] = Rot36_TT.row_x_out.read();
                    A_tile[idx_mat_2][r][c] = Rot36_TT.row_y_out.read();
                    Q_tile[idx_mat_1][r][c] = Rot36_TT.q_u_out.read();
                    Q_tile[idx_mat_2][r][c] = Rot36_TT.q_v_out.read();
                } else if (r == 1) {
                    A_tile[idx_mat_1][r][c] = Rot35_TT.row_x_out.read();
                    A_tile[idx_mat_2][r][c] = Rot35_TT.row_y_out.read();
                    Q_tile[idx_mat_1][r][c] = Rot35_TT.q_u_out.read();
                    Q_tile[idx_mat_2][r][c] = Rot35_TT.q_v_out.read();
                } else if (r == 2) {
                    A_tile[idx_mat_1][r][c] = Rot33_TT.row_x_out.read();
                    A_tile[idx_mat_2][r][c] = Rot33_TT.row_y_out.read();
                    Q_tile[idx_mat_1][r][c] = Rot33_TT.q_u_out.read();
                    Q_tile[idx_mat_2][r][c] = Rot33_TT.q_v_out.read();
                } else if (r == 3) {
                    A_tile[idx_mat_1][r][c] = Rot30_TT.row_x_out.read();
                    A_tile[idx_mat_2][r][c] = Rot30_TT.row_y_out.read();
                    Q_tile[idx_mat_1][r][c] = Rot30_TT.q_u_out.read();
                    Q_tile[idx_mat_2][r][c] = Rot30_TT.q_v_out.read();
                } else if (r == 4) {
                    A_tile[idx_mat_1][r][c] = Rot26_TT.row_x_out.read();
                    A_tile[idx_mat_2][r][c] = Rot26_TT.row_y_out.read();
                    Q_tile[idx_mat_1][r][c] = Rot26_TT.q_u_out.read();
                    Q_tile[idx_mat_2][r][c] = Rot26_TT.q_v_out.read();
                } else if (r == 5) {
                    A_tile[idx_mat_1][r][c] = Rot21_TT.row_x_out.read();
                    A_tile[idx_mat_2][r][c] = Rot21_TT.row_y_out.read();
                    Q_tile[idx_mat_1][r][c] = Rot21_TT.q_u_out.read();
                    Q_tile[idx_mat_2][r][c] = Rot21_TT.q_v_out.read();
                } else if (r == 6) {
                    A_tile[idx_mat_1][r][c] = Rot15_TT.row_x_out.read();
                    A_tile[idx_mat_2][r][c] = Rot15_TT.row_y_out.read();
                    Q_tile[idx_mat_1][r][c] = Rot15_TT.q_u_out.read();
                    Q_tile[idx_mat_2][r][c] = Rot15_TT.q_v_out.read();
                } else {
                    A_tile[idx_mat_1][r][c] = Rot8_TT.row_x_out.read();
                    A_tile[idx_mat_2][r][c] = Rot8_TT.row_y_out.read();
                    Q_tile[idx_mat_1][r][c] = Rot8_TT.q_u_out.read();
                    Q_tile[idx_mat_2][r][c] = Rot8_TT.q_v_out.read();
                }
            }
        }
    }
}
}
