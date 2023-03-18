#include <iostream>

#include "qrd.h"

int main(){
	data_t A[TAM][TAM] = {
							{3, 2, -1, 4, -7, 8},
							{2, 1, 9, 7, -1, 5},
							{-7, 5, 2, -6, 3, 4},
							{-1, 2, 1, -4, 6, -9},
							{-9, 16, 3, 1, 4, 11},
							{4, -2, 3, 7, 2, 2}
						 };

	data_t A_rot[TAM][TAM];

	/*
	 * ToDo: read input matrix and gold results from external file and implement Tile algorithm
	 */

	krnl_givens_rotation(A, A_rot);

    // Print result matrix
    for(int i = 0; i < TAM; i++){
		for(int j = 0; j < TAM; j++){
			std::cout << A_rot[i][j] << "  |  ";
		}
		std:: cout << std::endl;
	}
    return 0;
}
