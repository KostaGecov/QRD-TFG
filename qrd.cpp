#include <iostream>

#include "qrd.h"

Rotator::Rotator(int x, int y, int c) {
    Rotator::row_x = x;	// Realmente no hace falta
    Rotator::row_y = y; // Realmente no hace falta
    Rotator::col = c;
}

// Read input rows using blocking write to streams
void read_input_rows(data_t A[TAM_TILED][TAM], hls::stream<data_t, TAM> & row_x_in_1, hls::stream<data_t, TAM> & row_y_in_1, hls::stream<data_t, TAM> & row_x_in_2, hls::stream<data_t, TAM> & row_y_in_2, hls::stream<data_t, TAM> & row_x_in_3, hls::stream<data_t, TAM> & row_y_in_3) {
    // Read the rows from the input array and write them to the streams
    read_input_rows_for:
	for (index_t j = 0; j < TAM; j++) {
    	row_x_in_1.write(A[0][j]);
    	row_y_in_1.write(A[1][j]);
        row_x_in_2.write(A[2][j]);
        row_y_in_2.write(A[3][j]);
        row_x_in_3.write(A[4][j]);
        row_y_in_3.write(A[5][j]);
    }
}

// Performs the givens rotation of two streams and writes the outputs in other two streams
void Rotator::givens_rotation(hls::stream<data_t, TAM> & row_x_in, hls::stream<data_t, TAM> & row_y_in, hls::stream<data_t, TAM> & row_x_out, hls::stream<data_t, TAM> & row_y_out, int col_rotator) {
    bool sign;
    data_t x[TAM], y[TAM];
    data_t aux;

    read_input_data:
	for (index_t j = 0; j < TAM; j++) {
		x[j] = row_x_in.read();
		y[j] = row_y_in.read();
	}

    // Choose the right sign for the rotation, taking into account the quadrants of the coordinates
    if (x[col_rotator] < 0) {
        sign_for:
		for (index_t s = col_rotator; s < TAM; s++) { // Debo cambiar el signo a los 24 elementos de la fila
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
        column_rotation_for:
		for (index_t j = col_rotator; j < TAM; j++) {
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
			// do: operaciones opuestas?
		}
	}

    write_output_data:
	for (index_t j = 0; j < TAM; j++) {
		row_x_out.write(x[j]);
		row_y_out.write(y[j]);
	}
}

// Dataflow function
void krnl_givens_rotation(data_t A[TAM_TILED][TAM], data_t A_rot[TAM_TILED][TAM], index_t type_op, index_t col_offset) {
	if(type_op == GEQRT){
		// Rotators for GEQRT operation
		Rotator rot1(0, 1, 0);
		Rotator rot2(2, 3, 0);
		Rotator rot3(4, 5, 0);

		Rotator rot4(2, 4, 0);
		Rotator rot5(3, 5, 1);

		Rotator rot6(0, 2, 0);
		Rotator rot7(1, 3, 1);

		Rotator rot8(2, 4, 1);
		Rotator rot9(3, 5, 2);

		Rotator rot10(1, 2, 1);
		Rotator rot11(3, 4, 2);

		Rotator rot12(2, 3, 2);
		Rotator rot13(4, 5, 3);

		Rotator rot14(3, 4, 3);

		Rotator rot15(4, 5, 4);

		for(index_t j = 0; j < TAM / TAM_TILED; j++){
			read_input_rows(A, rot1.row_x_in, rot1.row_y_in, rot2.row_x_in, rot2.row_y_in, rot3.row_x_in, rot3.row_y_in);

			rot1.givens_rotation(rot1.row_x_in, rot1.row_y_in, rot1.row_x_out, rot1.row_y_out, rot1.col + col_offset);
			rot2.givens_rotation(rot2.row_x_in, rot2.row_y_in, rot2.row_x_out, rot2.row_y_out, rot2.col + col_offset);
			rot3.givens_rotation(rot3.row_x_in, rot3.row_y_in, rot3.row_x_out, rot3.row_y_out, rot3.col + col_offset);
			rot4.givens_rotation(rot2.row_x_out, rot3.row_x_out, rot4.row_x_out, rot4.row_y_out, rot4.col + col_offset);
			rot5.givens_rotation(rot2.row_y_out, rot3.row_y_out, rot5.row_x_out, rot5.row_y_out, rot5.col + col_offset);
			rot6.givens_rotation(rot1.row_x_out, rot4.row_x_out, rot6.row_x_out, rot6.row_y_out, rot6.col + col_offset);
			rot7.givens_rotation(rot1.row_y_out, rot5.row_x_out, rot7.row_x_out, rot7.row_y_out, rot7.col + col_offset);
			rot8.givens_rotation(rot6.row_y_out, rot4.row_y_out, rot8.row_x_out, rot8.row_y_out, rot8.col + col_offset);
			rot9.givens_rotation(rot7.row_y_out, rot5.row_y_out, rot9.row_x_out, rot9.row_y_out, rot9.col + col_offset);
			rot10.givens_rotation(rot7.row_x_out, rot8.row_x_out, rot10.row_x_out, rot10.row_y_out, rot10.col + col_offset);
			rot11.givens_rotation(rot9.row_x_out, rot8.row_y_out, rot11.row_x_out, rot11.row_y_out, rot11.col + col_offset);
			rot12.givens_rotation(rot10.row_y_out, rot11.row_x_out, rot12.row_x_out, rot12.row_y_out, rot12.col + col_offset);
			rot13.givens_rotation(rot11.row_y_out, rot9.row_y_out, rot13.row_x_out, rot13.row_y_out, rot13.col + col_offset);
			rot14.givens_rotation(rot12.row_y_out, rot13.row_x_out, rot14.row_x_out, rot14.row_y_out, rot14.col + col_offset);
			rot15.givens_rotation(rot14.row_y_out, rot13.row_y_out, rot15.row_x_out, rot15.row_y_out, rot15.col + col_offset);
		}
		// Write output streams to matrix A_rot
		for (index_t c = 0; c < TAM; c++) {
			write_output_streams_row_for:
			for (index_t r = 0; r < TAM_TILED; r++) {
				if (r == 0) A_rot[r][c] = rot6.row_x_out.read();
				else if (r == 1) A_rot[r][c] = rot10.row_x_out.read();
				else if (r == 2) A_rot[r][c] = rot12.row_x_out.read();
				else if (r == 3) A_rot[r][c] = rot14.row_x_out.read();
				else if (r == 4) A_rot[r][c] = rot15.row_x_out.read();
				else A_rot[r][c] = rot15.row_y_out.read();
			}
		}
	} else if(type_op == TTQRT){
		// Rotators for GEQRT operation
		Rotator rot1(0, 0, 0);
		Rotator rot2(1, 1, 1);
		Rotator rot3(2, 2, 2);
		Rotator rot4(3, 3, 3);
		Rotator rot5(4, 4, 4);
		Rotator rot6(5, 5, 5);

		Rotator rot7(0, 0, 1);
		Rotator rot8(1, 1, 2);
		Rotator rot9(2, 2, 3);
		Rotator rot10(3, 3, 4);
		Rotator rot11(4, 4, 5);
		
		Rotator rot12(0, 0, 2);
		Rotator rot13(1, 1, 3);
		Rotator rot14(2, 2, 4);
		Rotator rot15(3, 3, 5);
		
		Rotator rot16(0, 0, 3);
		Rotator rot17(1, 1, 4);
		Rotator rot18(2, 2, 5);
		
		Rotator rot19(0, 0, 4);
		Rotator rot20(1, 1, 5);
		
		Rotator rot21(0, 0, 5);

		// Read 6 input rows from first matrix
		read_input_rows(A, rot1.row_x_in, rot1.row_y_in, rot2.row_x_in, rot2.row_y_in, rot3.row_x_in, rot3.row_y_in);
		// Read 6 input rows from second matrix
		read_input_rows(A, rot1.row_x_in, rot1.row_y_in, rot2.row_x_in, rot2.row_y_in, rot3.row_x_in, rot3.row_y_in);

		/*
			ToDo: Complete with correct rows
		*/
		rot1.givens_rotation(rot1.row_x_in, rot1.row_y_in, rot1.row_x_out, rot1.row_y_out, rot1.col + col_offset);
		rot2.givens_rotation(rot2.row_x_in, rot2.row_y_in, rot2.row_x_out, rot2.row_y_out, rot2.col + col_offset);
		rot3.givens_rotation(rot3.row_x_in, rot3.row_y_in, rot3.row_x_out, rot3.row_y_out, rot3.col + col_offset);
		rot4.givens_rotation(rot2.row_x_out, rot3.row_x_out, rot4.row_x_out, rot4.row_y_out, rot4.col + col_offset);
		rot5.givens_rotation(rot2.row_y_out, rot3.row_y_out, rot5.row_x_out, rot5.row_y_out, rot5.col + col_offset);
		rot6.givens_rotation(rot1.row_x_out, rot4.row_x_out, rot6.row_x_out, rot6.row_y_out, rot6.col + col_offset);
		rot7.givens_rotation(rot1.row_y_out, rot5.row_x_out, rot7.row_x_out, rot7.row_y_out, rot7.col + col_offset);
		rot8.givens_rotation(rot6.row_y_out, rot4.row_y_out, rot8.row_x_out, rot8.row_y_out, rot8.col + col_offset);
		rot9.givens_rotation(rot7.row_y_out, rot5.row_y_out, rot9.row_x_out, rot9.row_y_out, rot9.col + col_offset);
		rot10.givens_rotation(rot7.row_x_out, rot8.row_x_out, rot10.row_x_out, rot10.row_y_out, rot10.col + col_offset);
		rot11.givens_rotation(rot9.row_x_out, rot8.row_y_out, rot11.row_x_out, rot11.row_y_out, rot11.col + col_offset);
		rot12.givens_rotation(rot10.row_y_out, rot11.row_x_out, rot12.row_x_out, rot12.row_y_out, rot12.col + col_offset);
		rot13.givens_rotation(rot11.row_y_out, rot9.row_y_out, rot13.row_x_out, rot13.row_y_out, rot13.col + col_offset);
		rot14.givens_rotation(rot12.row_y_out, rot13.row_x_out, rot14.row_x_out, rot14.row_y_out, rot14.col + col_offset);
		rot15.givens_rotation(rot14.row_y_out, rot13.row_y_out, rot15.row_x_out, rot15.row_y_out, rot15.col + col_offset);
		rot16.givens_rotation(rot14.row_y_out, rot13.row_y_out, rot15.row_x_out, rot15.row_y_out, rot15.col + col_offset);
		rot17.givens_rotation(rot14.row_y_out, rot13.row_y_out, rot15.row_x_out, rot15.row_y_out, rot15.col + col_offset);
		rot18.givens_rotation(rot14.row_y_out, rot13.row_y_out, rot15.row_x_out, rot15.row_y_out, rot15.col + col_offset);
		rot19.givens_rotation(rot14.row_y_out, rot13.row_y_out, rot15.row_x_out, rot15.row_y_out, rot15.col + col_offset);
		rot20.givens_rotation(rot14.row_y_out, rot13.row_y_out, rot15.row_x_out, rot15.row_y_out, rot15.col + col_offset);
		rot21.givens_rotation(rot14.row_y_out, rot13.row_y_out, rot15.row_x_out, rot15.row_y_out, rot15.col + col_offset);

		// Write 6 output rows to first matrix
		// Write 6 output rows to second matrix
	
	}
}
