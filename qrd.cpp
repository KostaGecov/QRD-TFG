#include <iostream>

#include "qrd.h"

/*
 * Bloque funciones de la clase Rotator
 *
 */

Rotator::Rotator(int x, int y, int c) {
    Rotator::row_x = x;
    Rotator::row_y = y;
    Rotator::col = c;
}

// Read input rows using blocking write to streams
void Rotator::read_input_rows(data_t A[TAM][TAM], hls::stream<data_t, 4> & row_x_in, hls::stream<data_t, 4> & row_y_in) {
    // Read the rows from the input array and write them to the streams
    read_input_rows_for: for (index_t j = 0; j < TAM; j++) {
        row_x_in.write(A[Rotator::row_x][j]);
        row_y_in.write(A[Rotator::row_y][j]);
    }
}

// Performs the givens rotation of two streams and writes the outputs in other two streams
void Rotator::givens_rotation(hls::stream<data_t, 4> & row_x_in, hls::stream<data_t, 4> & row_y_in, hls::stream<data_t, 4> & row_x_out, hls::stream<data_t, 4> & row_y_out) {
    bool sign;
    data_t x[TAM], y[TAM];
    data_t aux;

    read_input_data:
        for (index_t j = 0; j < TAM; j++) {
            #pragma HLS pipeline
            x[j] = row_x_in.read();
            y[j] = row_y_in.read();
        }

    // Choose the right sign for the rotation, taking into account the quadrants of the coordinates
    if (x[Rotator::col] < 0) {
        sign_for: for (index_t s = Rotator::col; s < TAM; s++) {
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
            column_rotation_for: for (index_t j = Rotator::col; j < TAM; j++) {
                if (y[Rotator::col] < 0) {
                    sign = true;
                } else {
                    sign = false;
                }
                // Si la  Y es negativa, hay que sumarle para que se vaya acercando a cero
                // y la X es al contrario
                if (sign) {
                    x[j] = x[j] - (y[j] >> k);
                    y[j] = y[j] + (x[j] >> k);

                } else {
                    x[j] = x[j] + (y[j] >> k);
                    y[j] = y[j] - (x[j] >> k);
                }
            }
        }

    write_output_data:
        for (index_t j = 0; j < TAM; j++) {
            #pragma HLS pipeline
            row_x_out.write(x[j]);
            row_y_out.write(y[j]);
        }
}

// Writes the streams that contain the rotated rows to a matrix
void write_streams_to_matrix(data_t A_rot[TAM][TAM], hls::stream<data_t, 4> & row_x_out_1, hls::stream<data_t, 4> & row_y_out_2, hls::stream<data_t, 4> & row_x_out_3, hls::stream<data_t, 4> & row_y_out_4) {
    write_streams_to_matrix_for: for (index_t j = 0; j < TAM; j++) {
        A_rot[0][j] = row_x_out_1.read();
        A_rot[1][j] = row_y_out_2.read();
        A_rot[2][j] = row_x_out_3.read();
        A_rot[3][j] = row_y_out_4.read();
    }
}

// Dataflow function
void krnl_givens_rotation(data_t A[TAM][TAM], data_t A_rot[TAM][TAM]) {
    // Dataflow region
    Rotator rot1(0, 1, 0);
    Rotator rot2(2, 3, 0);
    Rotator rot3(0, 2, 0);
    Rotator rot4(1, 3, 1);
    Rotator rot5(1, 2, 1);
    Rotator rot6(2, 3, 2);

    // Read input rows and perform Givens rotations
    rot1.read_input_rows(A, rot1.row_x_in, rot1.row_y_in);
    rot1.givens_rotation(rot1.row_x_in, rot1.row_y_in, rot1.row_x_out, rot1.row_y_out);

    rot2.read_input_rows(A, rot2.row_x_in, rot2.row_y_in);
    rot2.givens_rotation(rot2.row_x_in, rot2.row_y_in, rot2.row_x_out, rot2.row_y_out);

    rot3.givens_rotation(rot1.row_x_out, rot2.row_x_out, rot3.row_x_out, rot3.row_y_out);
    rot4.givens_rotation(rot1.row_y_out, rot2.row_y_out, rot4.row_x_out, rot4.row_y_out);
    rot5.givens_rotation(rot4.row_x_out, rot3.row_y_out, rot5.row_x_out, rot5.row_y_out);
    rot6.givens_rotation(rot5.row_y_out, rot4.row_y_out, rot6.row_x_out, rot6.row_y_out);

    // Write output streams to matrix A_rot
    for (index_t j = 0; j < TAM; j++) {
        for (index_t i = 0; i < 4; i++) {
            if (i == 0) A_rot[i][j] = rot3.row_x_out.read();
            else if (i == 1) A_rot[i][j] = rot5.row_x_out.read();
            else if (i == 2) A_rot[i][j] = rot6.row_x_out.read();
            else A_rot[i][j] = rot6.row_y_out.read();
        }
    }
    //    write_streams_to_matrix(A_rot, rot3.row_x_out, rot5.row_x_out, rot6.row_x_out, rot6.row_y_out); // Produces II violation
}
