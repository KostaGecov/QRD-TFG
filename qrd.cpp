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

		void rotation_giv(hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in, hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out){
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

void read_input(data_t A[TAM][TAM], hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in, int i){
	rows_to_stream:
	for(int j = 0; j < TAM; j++){
		row_x_in << A[i-1][j];
		row_y_in << A[i][j];
//		row_x_in.write(A[i-1][j]);
//		row_y_in.write(A[i][j]);
	}
	// los datos se guardan en otro sentido [3, 2, 1, 0]
}


// after data is read from an hls::stream<>, it cannot be read again

// Stream es una FIFO, tenerlo en cuenta para los índices de lectura/escritura

void rot_givens(hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in, hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out, int i) {
	int n_iter = 15;
	bool sign;
	data_t x[TAM], y[TAM], x_sh[TAM], y_sh[TAM];
	data_t aux;

	read_input_data:
	for(int j = 0; j < TAM; j++){
		row_x_in >> x[j];
		row_y_in >> y[j];
//		std::cout << "COORD X leida: " << x[j] << std::endl;
//		std::cout << "COORD Y leida: " << y[j] << std::endl;
	}


	check_sign:
	for(int s = 0; s < TAM; s++){
		if(x[0] < 0){
			if(y[s] >= 0){
				aux = x[s];
				x[s] = y[s];
				y[s] = -aux;
			}else {
				aux = x[s];
				x[s] = -y[s];
				y[s] = aux;
			}
		}
//		if(y[s] < 0){
//			sign[s] = true;
//		}else{
//			sign[s] = false;
//		}
	}

// solamente roto a 0 primer elemento de cada fila, tengo que seleccionar de alguna forma el par de coordenadas para hacer la Y = 0



		// con este for pretendo cambiar las coordenadas X/Y que quiero hacer 0, de las dos filas que me entran como argumentos
	change_pivot_for:
	for(int t = 0; t < i; t++){
//		if(y[t] < 0){
//			sign = true;
//		}else{
//			sign = false;
//		}
		rotation_for:
		for(int k = 0; k < n_iter; k++){
			column_rotation_for:
			for(int j = 0; j < TAM; j++){
				if((t+j) < TAM){
					if(y[t] < 0){
						sign = true;
					}else{
						sign = false;
					}
					if(sign){
						x[t+j] = x[t+j] - (y[t+j] >> k);
						y[t+j] = y[t+j] + (x[t+j] >> k);

					}else{
						x[t+j] = x[t+j] + (y[t+j] >> k);
						y[t+j] = y[t+j] - (x[t+j] >> k);
					}
				}
			}
		}
	}

	write_output_data:
	for(int j = 0; j < TAM; j++){
		row_x_out << x[j];
		row_y_out << y[j];
	}
}

void write_output(data_t A_rot[TAM][TAM], hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out, int i){
	stream_to_rows:
	for(int j = 0; j < TAM; j++){
//			row_x_out >> A_rot[i-1][j];
//			row_y_out >> A_rot[i][j];
		A_rot[i-1][j] = row_x_out.read();
		A_rot[i][j] = row_y_out.read();
	}
}

void krnl_givens_rotation(data_t A[TAM][TAM], data_t A_rot[TAM][TAM]){
	hls::stream<data_t> row_x_in;
	hls::stream<data_t> row_y_in;
	hls::stream<data_t> row_x_out;
	hls::stream<data_t> row_y_out;

	for(int i = TAM - 1; i > 0; --i){
		read_input(A, row_x_in, row_y_in, i);
		rot_givens(row_x_in, row_y_in, row_x_out, row_y_out, i);
		write_output(A_rot, row_x_out, row_y_out, i);
	}
}
