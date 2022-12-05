#include <iostream>
#include "qrd.h"

int main(){

	static int i = TAM -1;

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
							{2, 1, 5, 7},
							{0, 5, 2, -6},
							{-1, 2, 1, 0}
						 };

	data_t A_rot[TAM][TAM];

	/*
	 * ToDo: Pasar i como puntero para que dentro del dataflow pueda ser modificada y así ir iterando las filas de la matriz
	 *	krnl_givens_rotation(A, A_rot, i);
	 */

	for(i; i > 0; i--){
		krnl_givens_rotation(A, A_rot, i);
	}

    // Print result matrix
    for(int f = 0; f < TAM; f++){
		for(int j = 0; j < TAM; j++){

			std::cout << A_rot[f][j] << "  |  ";
		}
		std:: cout << std::endl;
	}
    return 0;
}
