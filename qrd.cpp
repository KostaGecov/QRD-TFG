#include <iostream>
#include "qrd.h"
/*
class Rotator{
	public:
		hls::stream<data_t> row_x_in, row_y_in;
		hls::stream<data_t> row_x_out, row_y_out;

		Rotator(){
			// Empty constructor
		}

		/*
		 * Member functions
		 *
		 */
		/*void rotation_giv(hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in, hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out){
			Rotator rotator[7];
			data_t x[TAM], y[TAM];
			data_t x_sh[TAM], y_sh[TAM];
			int n_iter = 15;


			rotation:
			for(int n = 0; n < n_iter; n++){
				for(int m = 0; m < TAM; m++){
					x_sh[m] = x[m] >> n;
					y_sh[m] = y[m] >> n;

					if (y[0] < 0) {
						x[m] = x[m] + y_sh[m];
						y[m] = y[m] - x_sh[m];
					} else {
						x[m] = x[m] - y_sh[m];
						y[m] = y[m] + x_sh[m];
					}
				}
			}
		}
};
*/
void read_input(data_t M[TAM][TAM], hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in, int i){
	for(int j = 0; j < TAM; j++){
		row_x_in << M[i-1][j];
		row_y_in << M[i][j];
	}
}


// after data is read from an hls::stream<>, it cannot be read again

// Stream es una FIFO, tenerlo en cuenta para los índices de lectura/escritura

void rot_givens(hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in, hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out) {
	int n_iter = 15;
	bool sign;
	data_t x[TAM], y[TAM], x_sh[TAM], y_sh[TAM];

	for(int i = 0; i < TAM; i++){
		row_x_in >> x[i];
		row_y_in >> y[i];
	}

	if(y[0] < 0){
		sign = true;
	}else{
		sign = false;
	}

	for(int k = 0; k < n_iter; k++){
		for(int i = 0; i < TAM; i++){
			if(sign){
				x[i] = x[i] - (y[i] >> k);
				y[i] = y[i] + (x[i] >> k);
			}else{
				x[i] = x[i] + (y[i] >> k);
				y[i] = y[i] - (x[i] >> k);
			}
		}
	}

	for(int i = 0; i < TAM; i++){
		row_x_out << x[i];
		row_y_out << y[i];
	}
}

void write_output(data_t M[TAM][TAM], hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out, int i){
	for(int j = 0; j < TAM; j++){
		row_x_out >> M[i-1][j];
		row_y_out >> M[i][j];
	}
}

void krnl_givens_rotation(data_t M[TAM][TAM], int i){
	hls::stream<data_t> row_x_in;
	hls::stream<data_t> row_y_in;
	hls::stream<data_t> row_x_out;
	hls::stream<data_t> row_y_out;

	read_input(M, row_x_in, row_y_in, i);
	rot_givens(row_x_in, row_y_in, row_x_out, row_y_out);
	write_output(M, row_x_out, row_y_out, i);

}
