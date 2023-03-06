#ifndef _QRD_
#define _QRD_

#define TAM 6
#define N_ITER 30

#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_stream.h>

// The more bits it has, more shifts can be performed later, so the approximation to 0 will be more precise
// For bigger matrices, we need bigger data formats to be able to calculate the right result
typedef ap_fixed<24, 10, AP_RND> data_t; // 24 bits fixed point data, 10 for integer value and 14 for decimals
// data type used for indexes variables in for loops
typedef ap_uint<5> index_t; // Max value inside code is 15 (n_iter)

class Rotator{
	public:
		// after data is read from an hls::stream<>, it cannot be read again
		// Stream is FIFO type
		hls::stream<data_t, 6> row_x_in, row_y_in;
		hls::stream<data_t, 6> row_x_out, row_y_out;
		int row_x, row_y, col; // posiciones de las filas a rotar. En teoría las columnas son las mismas en ambas filas

	public:
		Rotator();
		Rotator(int x, int y, int c);
//		void read_input_rows(data_t A[TAM][TAM], hls::stream<data_t, 6> &row_x_in, hls::stream<data_t, 6> &row_y_in);
		void givens_rotation(hls::stream<data_t, 6> &row_x_in, hls::stream<data_t, 6> &row_y_in, hls::stream<data_t, 6> &row_x_out, hls::stream<data_t, 6> &row_y_out);
		void write_output_rows(data_t A_rot[TAM][TAM], hls::stream<data_t, 6> &row_x_out, hls::stream<data_t, 6> &row_y_out);
};

//void write_streams_to_matrix(data_t A_rot[TAM][TAM], hls::stream < data_t, 4 > & row_x_out_1, hls::stream < data_t, 4 > & row_y_out_2, hls::stream < data_t, 4 > & row_x_out_3, hls::stream < data_t, 4 > & row_y_out_4);
void read_input_rows(data_t A[TAM][TAM], hls::stream<data_t, 6> & row_x_in_1, hls::stream<data_t, 6> & row_y_in_1, hls::stream<data_t, 6> & row_x_in_2, hls::stream<data_t, 6> & row_y_in_2, hls::stream<data_t, 6> & row_x_in_3, hls::stream<data_t, 6> & row_y_in_3);
void krnl_givens_rotation(data_t A[TAM][TAM], data_t A_rot[TAM][TAM]);

#endif
