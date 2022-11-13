#include <iostream>
#include "qrd.h"

int main(){

	data_t A[TAM][TAM] = {
							{3, 2, -1, 4, 1},
							{2, 1, 5, 7, 3},
							{8, 5, 2, -6, 5},
							{-1, 2, 1, -4, 7},
							{-9, 16, 3, 1, 4}
						 };

/*
	data_t A[TAM][TAM] = {
							{3, 2, -1, 4},
							{2, 1, 5, 7},
							{0, 5, 2, -6},
							{-1, 2, 1, 0}
						 };
*/
    rot_givens(A);

    // Print result matrix
    for(int i = 0; i < TAM; i++){
		for(int j = 0; j < TAM; j++){

			std::cout << A[i][j] << "  |  ";
		}
		std:: cout << std::endl;
	}

    return 0;
}
