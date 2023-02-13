#include <iostream>
#include "qrd.h"

int main(){
/*	data_t A[TAM][TAM] = {
							{3, 2, -1, 4, 1},
							{2, 1, 5, 7, 3},
							{8, 5, 2, -6, 5},
							{-1, 2, 1, -4, 7},
							{-9, 16, 3, 1, 4}
						 };

*/
	data_t A[TAM][TAM] = {
							{3, 2, -1, 4},
							{2, 1, 9, 7},
							{-7, 5, 2, -6},
							{-1, 2, 1, -4},
						 };

	data_t A_rot[TAM][TAM];

	Rotator rot1(0, 1, 0);
	Rotator rot2(2, 3, 0);
	rot1.read_input_rows(A, rot1.row_x_in, rot1.row_y_in);
	rot2.read_input_rows(A, rot2.row_x_in, rot2.row_y_in);
	rot1.givens_rotation(rot1.row_x_in, rot1.row_y_in, rot1.row_x_out, rot1.row_y_out);
	rot2.givens_rotation(rot2.row_x_in, rot2.row_y_in, rot2.row_x_out, rot2.row_y_out);
	rot1.write_output_rows(A, rot1.row_x_out, rot1.row_y_out);
	rot2.write_output_rows(A, rot2.row_x_out, rot2.row_y_out);


	Rotator rot3(0, 2, 0);
	Rotator rot4(1, 3, 1);
	rot3.read_input_rows(A, rot3.row_x_in, rot3.row_y_in);
	rot4.read_input_rows(A, rot4.row_x_in, rot4.row_y_in);
	rot3.givens_rotation(rot3.row_x_in, rot3.row_y_in, rot3.row_x_out, rot3.row_y_out);
	rot4.givens_rotation(rot4.row_x_in, rot4.row_y_in, rot4.row_x_out, rot4.row_y_out);
	rot3.write_output_rows(A, rot3.row_x_out, rot3.row_y_out);
	rot4.write_output_rows(A, rot4.row_x_out, rot4.row_y_out);

	Rotator rot5(1, 2, 1);
	rot5.read_input_rows(A, rot5.row_x_in, rot5.row_y_in);
	rot5.givens_rotation(rot5.row_x_in, rot5.row_y_in, rot5.row_x_out, rot5.row_y_out);
	rot5.write_output_rows(A, rot5.row_x_out, rot5.row_y_out);

	Rotator rot6(2, 3, 2);
	rot6.read_input_rows(A, rot6.row_x_in, rot6.row_y_in);
	rot6.givens_rotation(rot6.row_x_in, rot6.row_y_in, rot6.row_x_out, rot6.row_y_out);
	rot6.write_output_rows(A, rot6.row_x_out, rot6.row_y_out);

//	krnl_givens_rotation(A, A_rot);

    // Print result matrix
    for(int i = 0; i < TAM; i++){
		for(int j = 0; j < TAM; j++){
			std::cout << A[i][j] << "  |  ";
		}
		std:: cout << std::endl;
	}
    return 0;
}
