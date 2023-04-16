#include <iostream>

#include "qrd.h"

int main(){
	data_t A[TAM_TILED][TAM] = {
							{3, 2, -1, 4, -7, 8, 3, 2, -1, 4, -7, 8, 3, 2, -1, 4, -7, 8, 3, 2, -1, 4, -7, 8},
							{2, 1, 9, 7, -1, 5, 2, 1, 9, 7, -1, 5, 2, 1, 9, 7, -1, 5, 2, 1, 9, 7, -1, 5},
							{-7, 5, 2, -6, 3, 4, -7, 5, 2, -6, 3, 4, -7, 5, 2, -6, 3, 4, -7, 5, 2, -6, 3, 4},
							{-1, 2, 1, -4, 6, -9, -1, 2, 1, -4, 6, -9, -1, 2, 1, -4, 6, -9, -1, 2, 1, -4, 6, -9},
							{-9, 16, 3, 1, 4, 11, -9, 16, 3, 1, 4, 11, -9, 16, 3, 1, 4, 11, -9, 16, 3, 1, 4, 11},
							{4, -2, 3, 7, 2, 2, 4, -2, 3, 7, 2, 2, 4, -2, 3, 7, 2, 2, 4, -2, 3, 7, 2, 2}
						 };

	data_t A_rot[TAM_TILED][TAM];
	static index_t col_offset = 0; // offset to access the right column for GEQRT operation

	/*
	 * ToDo: read input matrix and gold results from external file and implement Tile algorithm
	 */


	for(index_t i = 2; i < NUM_OPERACIONES + 2; i++){
		// GEQRT operation
		if(i % 2 == 1){
			krnl_givens_rotation(A, A_rot, GEQRT, col_offset);
			col_offset += 6;
		} else {
			krnl_givens_rotation(A, A_rot, TTQRT, 0);
		}
	}

    // Print result matrix
    for(int i = 0; i < TAM_TILED; i++){
		for(int j = 0; j < TAM; j++){
			std::cout << A_rot[i][j] << "  |  ";
		}
		std:: cout << std::endl;
	}
    return 0;
}
