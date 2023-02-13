#include <iostream>
#include "qrd.h"

/*
 * Bloque funciones de la clase Rotator
 *
 */

Rotator::Rotator(int x, int y, int c){
	Rotator::row_x = x;
	Rotator::row_y = y;
	Rotator::col = c;
}

void Rotator::read_input_rows(data_t A[TAM][TAM], hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in){
	row_to_stream:
	for(int j = 0; j < TAM; j++){
		row_x_in << A[Rotator::row_x][j];
		row_y_in << A[Rotator::row_y][j];
	}
}

void Rotator::givens_rotation(hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in, hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out){
	int n_iter = 15;
	bool sign;
	data_t x[TAM], y[TAM], x_sh[TAM], y_sh[TAM];
	data_t aux;

	read_input_data:
	for(int j = 0; j < TAM; j++){
		row_x_in >> x[j];
		row_y_in >> y[j];
	}

	// Bucle para seleccionar el signo adecuado, teniendo en cuenta los cuadrantes de las coordenadas
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
	}

	rotation_for:
	for(int k = 0; k < n_iter; k++){
		column_rotation_for:
		for(int j = Rotator::col; j < TAM; j++){
			// if(j < TAM){
			if(y[Rotator::col] < 0){
				sign = true;
			}else{
				sign = false;
			}
			if(sign){
				x[j] = x[j] - (y[j] >> k);
				y[j] = y[j] + (x[j] >> k);

			}else{
				x[j] = x[j] + (y[j] >> k);
				y[j] = y[j] - (x[j] >> k);
			}
			// }
		}
	}

	write_output_data:
	for(int j = 0; j < TAM; j++){
		row_x_out << x[j];
		row_y_out << y[j];
	}
}

void Rotator::write_output_rows(data_t A_rot[TAM][TAM], hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out){
	stream_to_rows:
	for(int j = 0; j < TAM; j++){
		row_x_out >> A_rot[Rotator::row_x][j];
		row_y_out >> A_rot[Rotator::row_y][j];
	}
}

void Rotator::krnl_gr(data_t A[TAM][TAM], data_t A_rot[TAM][TAM]){
//	hls::stream<data_t> row_x_in;
//	hls::stream<data_t> row_y_in;
//	hls::stream<data_t> row_x_out;
//	hls::stream<data_t> row_y_out;

	Rotator::read_input_rows(A, row_x_in, row_y_in);
	//givens_rotation();
	Rotator::write_output_rows(A_rot, row_x_out, row_y_out);
}

/*
 * Bloque funciones Dataflow Antiguo
 *
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

void rot_givens(hls::stream<data_t> &row_x_in, hls::stream<data_t> &row_y_in, hls::stream<data_t> &row_x_out, hls::stream<data_t> &row_y_out, int i) {
	int n_iter = 30;
	bool sign;
	data_t x[TAM], y[TAM], x_sh[TAM], y_sh[TAM];
	data_t aux;

	read_input_data:
	for(int j = 0; j < TAM; j++){
		row_x_in >> x[j];
		row_y_in >> y[j];
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

	// con este for pretendo cambiar las coordenadas X/Y que quiero hacer 0, de las dos filas que me entran como argumentos
	change_pivot_for:
	for(int t = 0; t < i; t++){
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
		row_x_out >> A_rot[i-1][j];
		row_y_out >> A_rot[i][j];
	}
}

void krnl_givens_rotation(data_t A[TAM][TAM], data_t A_rot[TAM][TAM]){
	krnl_for:
	for(int i = TAM - 1; i > 0; --i){
		hls::stream<data_t> row_x_in;
		hls::stream<data_t> row_y_in;
		hls::stream<data_t> row_x_out;
		hls::stream<data_t> row_y_out;


		/*
		 * ToDo: ver cómo hacer que se hagan de forma paralela
		 *
		 */

		Rotator rot1(0, 1, 0);
		Rotator rot2(2, 3, 0);

		Rotator rot3(0, 2, 0);
		Rotator rot4(1, 3, 1);

		Rotator rot5(1, 2, 1);

		Rotator rot6(2, 3, 2);

		read_input(A, row_x_in, row_y_in, i);
		rot_givens(row_x_in, row_y_in, row_x_out, row_y_out, i);
		write_output(A_rot, row_x_out, row_y_out, i);
	}
}
