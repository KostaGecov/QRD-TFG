#ifndef _QRD_
#define _QRD_

#define TAM_INDEX 5
#define FIXED_POINT 24
#define FX_POINT_INT 10

#define TAM_TILED 6
#define TAM 24
#define N_ITER 30
#define NUM_OPERACIONES 7

#define GEQRT 0
#define TTQRT 1

#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_stream.h>

// The more bits it has, more shifts can be performed later, so the approximation to 0 will be more precise
// For bigger matrices, we need bigger data formats to be able to calculate the right result
typedef ap_fixed<FIXED_POINT, FX_POINT_INT, AP_RND> data_t; // 24 bits fixed point data, 10 for integer value and 14 for decimals

// data type used for indexes variables in for loops
typedef ap_uint<TAM_INDEX> index_t; // Max value inside code is 15 (n_iter)

class Rotator{
	public:
		// after data is read from an hls::stream<>, it cannot be read again
		// Stream is FIFO type
		hls::stream<data_t, TAM> row_x_in, row_y_in;
		hls::stream<data_t, TAM> row_x_out, row_y_out;
//		hls::stream<data_t, TAM> full_row_x_in, full_row_y_in;
//		hls::stream<data_t, TAM> full_row_x_out, full_row_y_out;
		int row_x, row_y, col; // posiciones de las filas a rotar. En teor�a las columnas son las mismas en ambas filas

	public:
		Rotator();
		Rotator(int x, int y, int c);
		void givens_rotation(hls::stream<data_t, TAM> &row_x_in, hls::stream<data_t, TAM> &row_y_in, hls::stream<data_t, TAM> &row_x_out, hls::stream<data_t, TAM> &row_y_out, int col_rotator);
		void write_output_rows(data_t A_rot[TAM][TAM], hls::stream<data_t, TAM> &row_x_out, hls::stream<data_t, TAM> &row_y_out);
};

void read_input_rows(data_t A[TAM][TAM], hls::stream<data_t, TAM> & row_x_in_1, hls::stream<data_t, TAM> & row_y_in_1, hls::stream<data_t, TAM> & row_x_in_2, hls::stream<data_t, TAM> & row_y_in_2, hls::stream<data_t, TAM> & row_x_in_3, hls::stream<data_t, TAM> & row_y_in_3);
void krnl_givens_rotation(data_t A[TAM][TAM], data_t A_rot[TAM][TAM], index_t type_op, index_t col_offset);

#endif
