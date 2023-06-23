#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_stream.h>
#include <iostream>

#define TAM_INDEX 6
#define FIXED_POINT 28
#define FX_POINT_INT 10

#define TAM_TILED 32
#define TAM 256
#define N_ITER 40
#define NUM_OPERACIONES 17  //(15 + 2(offset))

#define GEQRT 0
#define TTQRT 1

// The more bits it has, more shifts can be performed later, so the
// approximation to 0 will be more precise For bigger matrices, we need bigger
// data formats to be able to calculate the right result
typedef ap_fixed<FIXED_POINT, FX_POINT_INT, AP_RND> data_t;  // 24 bits fixed point data, 10 for integer value and 14 for
                                                             // decimals

// data type used for indexes variables in for loops
typedef ap_uint<TAM_INDEX> index_t;  // Max value inside code is 15 (n_iter)

class Rotator {
   public:
    // after data is read from an hls::stream<>, it cannot be read again
    // Stream is FIFO type
    hls::stream<data_t, TAM> row_x_in, row_y_in;
    hls::stream<data_t, TAM> row_x_out, row_y_out;

    int row_x, row_y, col;  // posiciones de las filas a rotar. En teor�a las
                            // columnas son las mismas en ambas filas

   public:
    // Rotator();
    Rotator(int x, int y, int c);
    void givens_rotation(hls::stream<data_t, TAM>& row_x_in,
                         hls::stream<data_t, TAM>& row_y_in,
                         hls::stream<data_t, TAM>& row_x_out,
                         hls::stream<data_t, TAM>& row_y_out, int col_rotator);
};

void read_input_rows(data_t A[TAM][TAM],
                     hls::stream<data_t, TAM>& row_in_1,
                     hls::stream<data_t, TAM>& row_in_2,
                     hls::stream<data_t, TAM>& row_in_3,
                     hls::stream<data_t, TAM>& row_in_4,
                     hls::stream<data_t, TAM>& row_in_5,
                     hls::stream<data_t, TAM>& row_in_6);
extern "C" {
void krnl_givens_rotation(data_t A_tiled_1[TAM_TILED][TAM], data_t A_tiled_2[TAM_TILED][TAM],
                          index_t type_op, index_t col_offset);
}

const data_t SCALE_FACTOR = 0.6072529;

Rotator::Rotator(int x, int y, int c) {
    Rotator::row_x = x;  // Realmente no hace falta
    Rotator::row_y = y;  // Realmente no hace falta
    Rotator::col = c;
}

// Read input rows using blocking write to streams
void read_input_rows(data_t Matrix[TAM_TILED][TAM],
                     hls::stream<data_t, TAM>& row_in_1,
                     hls::stream<data_t, TAM>& row_in_2,
                     hls::stream<data_t, TAM>& row_in_3,
                     hls::stream<data_t, TAM>& row_in_4,
                     hls::stream<data_t, TAM>& row_in_5,
                     hls::stream<data_t, TAM>& row_in_6) {
#pragma HLS INLINE off

// Read the rows from the input array and write them to the streams
read_input_rows_for:
    for (index_t j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT avg = 24 max = 24 min = 0
        row_in_1.write(Matrix[0][j]);
        row_in_2.write(Matrix[1][j]);
        row_in_3.write(Matrix[2][j]);
        row_in_4.write(Matrix[3][j]);
        row_in_5.write(Matrix[4][j]);
        row_in_6.write(Matrix[5][j]);
    }
}

// Performs the givens rotation of two streams and writes the outputs in other
// two streams
void Rotator::givens_rotation(hls::stream<data_t, TAM>& row_x_in,
                              hls::stream<data_t, TAM>& row_y_in,
                              hls::stream<data_t, TAM>& row_x_out,
                              hls::stream<data_t, TAM>& row_y_out,
                              int col_rotator) {
#pragma HLS INLINE off
    bool sign;
    data_t x[TAM], y[TAM];
    data_t aux;

read_input_data:
    for (index_t j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT avg = 24 max = 24 min = 0
        x[j] = row_x_in.read();
        y[j] = row_y_in.read();
    }

    // Choose the right sign for the rotation, taking into account the quadrants
    // of the coordinates
    if (x[col_rotator] < 0) {
    sign_for:
        for (index_t s = col_rotator; s < TAM; s++) {
#pragma HLS LOOP_TRIPCOUNT max = 24 min = 0
            // Debo cambiar el signo a los 24 elementos de la fila
            if (y[s] >= 0) {
                aux = x[s];
                x[s] = y[s];
                y[s] = -aux;
            } else {
                aux = x[s];
                x[s] = -y[s];
                y[s] = aux;
            }
        }
    }

iterations_for:
    for (index_t k = 0; k < N_ITER; k++) {
#pragma HLS LOOP_TRIPCOUNT max = 30 min = 0
    column_rotation_for:
        for (index_t j = col_rotator; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = 24 min = 0
            if (y[col_rotator] < 0) {
                sign = true;
            } else {
                sign = false;
            }
            // If Y is negative, we need to add to it so that it gets closer to zero
            // and to the contrary with X coordinate
            if (sign) {
                x[j] = x[j] - (y[j] >> k);
                y[j] = y[j] + (x[j] >> k);

            } else {
                x[j] = x[j] + (y[j] >> k);
                y[j] = y[j] - (x[j] >> k);
            }
            // if j >= TAM_TILED
            // do: operaciones opuestas? No
        }
    }

scale_factor_for:
    for (index_t j = col_rotator; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = 24 min = 0
        x[j] = x[j] * SCALE_FACTOR;
        y[j] = y[j] * SCALE_FACTOR;
    }

write_output_data:
    for (index_t j = 0; j < TAM; j++) {
#pragma HLS LOOP_TRIPCOUNT max = 24 min = 0
        row_x_out.write(x[j]);
        row_y_out.write(y[j]);
    }
}

extern "C" {
// Dataflow function
void krnl_givens_rotation(data_t A_tiled_1[TAM_TILED][TAM],
                          data_t A_tiled_2[TAM_TILED][TAM],
                          index_t type_op, index_t col_offset) {
#pragma HLS ARRAY_PARTITION dim = 2 factor = 6 type = block variable = A_tiled_1
#pragma HLS ARRAY_PARTITION dim = 2 factor = 6 type = block variable = A_tiled_2
#pragma HLS DATAFLOW
#pragma HLS TOP name = krnl_givens_rotation

    if (type_op == GEQRT) {
    GEQRT_OPERATION:
        // Rotators for GEQRT operation
        Rotator Rot1_GE(0, 1, 0);
        Rotator Rot2_GE(2, 3, 0);
        Rotator Rot3_GE(4, 5, 0);

        Rotator Rot4_GE(2, 4, 0);
        Rotator Rot5_GE(3, 5, 1);

        Rotator Rot6_GE(0, 2, 0);
        Rotator Rot7_GE(1, 3, 1);

        Rotator Rot8_GE(2, 4, 1);
        Rotator Rot9_GE(3, 5, 2);

        Rotator Rot10_GE(1, 2, 1);
        Rotator Rot11_GE(3, 4, 2);

        Rotator Rot12_GE(2, 3, 2);
        Rotator Rot13_GE(4, 5, 3);

        Rotator Rot14_GE(3, 4, 3);

        Rotator Rot15_GE(4, 5, 4);

        read_input_rows(A_tiled_1, Rot1_GE.row_x_in, Rot1_GE.row_y_in, Rot2_GE.row_x_in,
                        Rot2_GE.row_y_in, Rot3_GE.row_x_in, Rot3_GE.row_y_in);

        Rot1_GE.givens_rotation(Rot1_GE.row_x_in, Rot1_GE.row_y_in, Rot1_GE.row_x_out,
                             Rot1_GE.row_y_out, Rot1_GE.col + col_offset);
        Rot2_GE.givens_rotation(Rot2_GE.row_x_in, Rot2_GE.row_y_in, Rot2_GE.row_x_out,
                             Rot2_GE.row_y_out, Rot2_GE.col + col_offset);
        Rot3_GE.givens_rotation(Rot3_GE.row_x_in, Rot3_GE.row_y_in, Rot3_GE.row_x_out,
                             Rot3_GE.row_y_out, Rot3_GE.col + col_offset);
        Rot4_GE.givens_rotation(Rot2_GE.row_x_out, Rot3_GE.row_x_out, Rot4_GE.row_x_out,
                             Rot4_GE.row_y_out, Rot4_GE.col + col_offset);
        Rot5_GE.givens_rotation(Rot2_GE.row_y_out, Rot3_GE.row_y_out, Rot5_GE.row_x_out,
                             Rot5_GE.row_y_out, Rot5_GE.col + col_offset);
        Rot6_GE.givens_rotation(Rot1_GE.row_x_out, Rot4_GE.row_x_out, Rot6_GE.row_x_out,
                             Rot6_GE.row_y_out, Rot6_GE.col + col_offset);
        Rot7_GE.givens_rotation(Rot1_GE.row_y_out, Rot5_GE.row_x_out, Rot7_GE.row_x_out,
                             Rot7_GE.row_y_out, Rot7_GE.col + col_offset);
        Rot8_GE.givens_rotation(Rot6_GE.row_y_out, Rot4_GE.row_y_out, Rot8_GE.row_x_out,
                             Rot8_GE.row_y_out, Rot8_GE.col + col_offset);
        Rot9_GE.givens_rotation(Rot7_GE.row_y_out, Rot5_GE.row_y_out, Rot9_GE.row_x_out,
                             Rot9_GE.row_y_out, Rot9_GE.col + col_offset);
        Rot10_GE.givens_rotation(Rot7_GE.row_x_out, Rot8_GE.row_x_out, Rot10_GE.row_x_out,
                              Rot10_GE.row_y_out, Rot10_GE.col + col_offset);
        Rot11_GE.givens_rotation(Rot9_GE.row_x_out, Rot8_GE.row_y_out, Rot11_GE.row_x_out,
                              Rot11_GE.row_y_out, Rot11_GE.col + col_offset);
        Rot12_GE.givens_rotation(Rot10_GE.row_y_out, Rot11_GE.row_x_out, Rot12_GE.row_x_out,
                              Rot12_GE.row_y_out, Rot12_GE.col + col_offset);
        Rot13_GE.givens_rotation(Rot11_GE.row_y_out, Rot9_GE.row_y_out, Rot13_GE.row_x_out,
                              Rot13_GE.row_y_out, Rot13_GE.col + col_offset);
        Rot14_GE.givens_rotation(Rot12_GE.row_y_out, Rot13_GE.row_x_out, Rot14_GE.row_x_out,
                              Rot14_GE.row_y_out, Rot14_GE.col + col_offset);
        Rot15_GE.givens_rotation(Rot14_GE.row_y_out, Rot13_GE.row_y_out, Rot15_GE.row_x_out,
                              Rot15_GE.row_y_out, Rot15_GE.col + col_offset);

    // Write output streams to matrix A_tiled_1
    write_output_streams_col_for:
        for (index_t c = 0; c < TAM; c++) {
#pragma HLS LOOP_TRIPCOUNT max = 24 min = 0

        write_output_streams_row_for:
            for (index_t r = 0; r < TAM_TILED; r++) {
#pragma HLS LOOP_TRIPCOUNT max = 6 min = 0

                if (r == 0)
                    A_tiled_1[r][c] = Rot6_GE.row_x_out.read();
                else if (r == 1)
                    A_tiled_1[r][c] = Rot10_GE.row_x_out.read();
                else if (r == 2)
                    A_tiled_1[r][c] = Rot12_GE.row_x_out.read();
                else if (r == 3)
                    A_tiled_1[r][c] = Rot14_GE.row_x_out.read();
                else if (r == 4)
                    A_tiled_1[r][c] = Rot15_GE.row_x_out.read();
                else
                    A_tiled_1[r][c] = Rot15_GE.row_y_out.read();
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

        Rotator Rot7_TT(0, 0, 1);
        Rotator Rot8_TT(1, 1, 2);
        Rotator Rot9_TT(2, 2, 3);
        Rotator Rot10_TT(3, 3, 4);
        Rotator Rot11_TT(4, 4, 5);

        Rotator Rot12_TT(0, 0, 2);
        Rotator Rot13_TT(1, 1, 3);
        Rotator Rot14_TT(2, 2, 4);
        Rotator Rot15_TT(3, 3, 5);

        Rotator Rot16_TT(0, 0, 3);
        Rotator Rot17_TT(1, 1, 4);
        Rotator Rot18_TT(2, 2, 5);

        Rotator Rot19_TT(0, 0, 4);
        Rotator Rot20_TT(1, 1, 5);

        Rotator Rot21_TT(0, 0, 5);

        // Read X coordinates rows from first matrix
        read_input_rows(A_tiled_1, Rot1_TT.row_x_in, Rot2_TT.row_x_in,
                        Rot3_TT.row_x_in, Rot4_TT.row_x_in, Rot5_TT.row_x_in,
                        Rot6_TT.row_x_in);
        // Read Y coordinates rows from second matrix
        read_input_rows(A_tiled_2, Rot1_TT.row_y_in, Rot2_TT.row_y_in,
                        Rot3_TT.row_y_in, Rot4_TT.row_y_in, Rot5_TT.row_y_in,
                        Rot6_TT.row_y_in);

        Rot1_TT.givens_rotation(Rot1_TT.row_x_in, Rot1_TT.row_y_in,
                                Rot1_TT.row_x_out, Rot1_TT.row_y_out,
                                Rot1_TT.col + col_offset);
        Rot2_TT.givens_rotation(Rot2_TT.row_x_in, Rot2_TT.row_y_in,
                                Rot2_TT.row_x_out, Rot2_TT.row_y_out,
                                Rot2_TT.col + col_offset);
        Rot3_TT.givens_rotation(Rot3_TT.row_x_in, Rot3_TT.row_y_in,
                                Rot3_TT.row_x_out, Rot3_TT.row_y_out,
                                Rot3_TT.col + col_offset);
        Rot4_TT.givens_rotation(Rot4_TT.row_x_in, Rot4_TT.row_y_in,
                                Rot4_TT.row_x_out, Rot4_TT.row_y_out,
                                Rot4_TT.col + col_offset);
        Rot5_TT.givens_rotation(Rot5_TT.row_x_in, Rot5_TT.row_y_in,
                                Rot5_TT.row_x_out, Rot5_TT.row_y_out,
                                Rot5_TT.col + col_offset);
        Rot6_TT.givens_rotation(Rot6_TT.row_x_in, Rot6_TT.row_y_in,
                                Rot6_TT.row_x_out, Rot6_TT.row_y_out,
                                Rot6_TT.col + col_offset);
        Rot7_TT.givens_rotation(Rot1_TT.row_x_out, Rot1_TT.row_y_out,
                                Rot7_TT.row_x_out, Rot7_TT.row_y_out,
                                Rot7_TT.col + col_offset);
        Rot8_TT.givens_rotation(Rot2_TT.row_x_out, Rot2_TT.row_y_out,
                                Rot8_TT.row_x_out, Rot8_TT.row_y_out,
                                Rot8_TT.col + col_offset);
        Rot9_TT.givens_rotation(Rot3_TT.row_x_out, Rot3_TT.row_y_out,
                                Rot9_TT.row_x_out, Rot9_TT.row_y_out,
                                Rot9_TT.col + col_offset);
        Rot10_TT.givens_rotation(Rot4_TT.row_x_out, Rot4_TT.row_y_out,
                                 Rot10_TT.row_x_out, Rot10_TT.row_y_out,
                                 Rot10_TT.col + col_offset);
        Rot11_TT.givens_rotation(Rot5_TT.row_x_out, Rot5_TT.row_y_out,
                                 Rot11_TT.row_x_out, Rot11_TT.row_y_out,
                                 Rot11_TT.col + col_offset);
        Rot12_TT.givens_rotation(Rot7_TT.row_x_out, Rot7_TT.row_y_out,
                                 Rot12_TT.row_x_out, Rot12_TT.row_y_out,
                                 Rot12_TT.col + col_offset);
        Rot13_TT.givens_rotation(Rot8_TT.row_x_out, Rot8_TT.row_y_out,
                                 Rot13_TT.row_x_out, Rot13_TT.row_y_out,
                                 Rot13_TT.col + col_offset);
        Rot14_TT.givens_rotation(Rot9_TT.row_x_out, Rot9_TT.row_y_out,
                                 Rot14_TT.row_x_out, Rot14_TT.row_y_out,
                                 Rot14_TT.col + col_offset);
        Rot15_TT.givens_rotation(Rot10_TT.row_x_out, Rot10_TT.row_y_out,
                                 Rot15_TT.row_x_out, Rot15_TT.row_y_out,
                                 Rot15_TT.col + col_offset);
        Rot16_TT.givens_rotation(Rot12_TT.row_x_out, Rot12_TT.row_y_out,
                                 Rot16_TT.row_x_out, Rot16_TT.row_y_out,
                                 Rot16_TT.col + col_offset);
        Rot17_TT.givens_rotation(Rot13_TT.row_x_out, Rot13_TT.row_y_out,
                                 Rot17_TT.row_x_out, Rot17_TT.row_y_out,
                                 Rot17_TT.col + col_offset);
        Rot18_TT.givens_rotation(Rot14_TT.row_x_out, Rot14_TT.row_y_out,
                                 Rot18_TT.row_x_out, Rot18_TT.row_y_out,
                                 Rot18_TT.col + col_offset);
        Rot19_TT.givens_rotation(Rot16_TT.row_x_out, Rot16_TT.row_y_out,
                                 Rot19_TT.row_x_out, Rot19_TT.row_y_out,
                                 Rot19_TT.col + col_offset);
        Rot20_TT.givens_rotation(Rot17_TT.row_x_out, Rot17_TT.row_y_out,
                                 Rot20_TT.row_x_out, Rot20_TT.row_y_out,
                                 Rot20_TT.col + col_offset);
        Rot21_TT.givens_rotation(Rot19_TT.row_x_out, Rot19_TT.row_y_out,
                                 Rot21_TT.row_x_out, Rot21_TT.row_y_out,
                                 Rot21_TT.col + col_offset);

    /*
        ToDo: Write 6 output rows to first matrix and Write 6 output rows to
        second matrix
    */
    // Write output streams to matrix A_tiled_1
    write_output_streams_col_TTQRT_for:
        for (index_t c = 0; c < TAM; c++) {
#pragma HLS LOOP_TRIPCOUNT max = 24 min = 0

        write_output_streams_row_TTQRT_for:
            for (index_t r = 0; r < TAM_TILED; r++) {
#pragma HLS LOOP_TRIPCOUNT max = 6 min = 0

                if (r == 0)
                    A_tiled_1[r][c] = Rot21_TT.row_x_out.read();
                else if (r == 1)
                    A_tiled_1[r][c] = Rot20_TT.row_x_out.read();
                else if (r == 2)
                    A_tiled_1[r][c] = Rot18_TT.row_x_out.read();
                else if (r == 3)
                    A_tiled_1[r][c] = Rot15_TT.row_x_out.read();
                else if (r == 4)
                    A_tiled_1[r][c] = Rot11_TT.row_x_out.read();
                else
                    A_tiled_1[r][c] = Rot6_TT.row_x_out.read();
            }
        }
    write_output_streams_col_TTQRT_y_for:
        for (index_t c = 0; c < TAM; c++) {
#pragma HLS LOOP_TRIPCOUNT max = 24 min = 0

        write_output_streams_row_TTQRT_y_for:
            for (index_t r = 0; r < TAM_TILED; r++) {
#pragma HLS LOOP_TRIPCOUNT max = 6 min = 0

                if (r == 0)
                    A_tiled_2[r][c] = Rot21_TT.row_y_out.read();
                else if (r == 1)
                    A_tiled_2[r][c] = Rot20_TT.row_y_out.read();
                else if (r == 2)
                    A_tiled_2[r][c] = Rot18_TT.row_y_out.read();
                else if (r == 3)
                    A_tiled_2[r][c] = Rot15_TT.row_y_out.read();
                else if (r == 4)
                    A_tiled_2[r][c] = Rot11_TT.row_y_out.read();
                else
                    A_tiled_2[r][c] = Rot6_TT.row_y_out.read();
            }
        }
    }
}
}
