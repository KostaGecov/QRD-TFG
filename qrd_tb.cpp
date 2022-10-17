#include <iostream>
#include "qrd.h"

int main(){

	data_t A[TAM][TAM] = {
							{3, 2, -1, 4, 1},
							{2, 1, 5, 7, 3},
							{0, 5, 2, -6, 5},
							{-1, 2, 1, 0, 7},
							{-9, 16, 8, -1, 4}
						 };

/*
	data_t A[TAM][TAM] = {
							{3, 2, -1, 4},
							{2, 1, 5, 7},
							{0, 5, 2, -6},
							{-1, 2, 1, 0}
						 };

    bool sign;
    static int n_iter = 15;
*/

    rot_givens(A);

/*
    // Primera columna

    rot_givens_succ(A, A[2][0], A[3][0], sign, n_iter, 2, 3, 0);

    rot_givens_succ(A, A[0][0], A[1][0], sign, n_iter, 0, 1, 0);

    rot_givens_succ(A, A[0][0], A[2][0], sign, n_iter, 0, 2, 0);

    // Segunda columna

    rot_givens_succ(A, A[2][1], A[3][1], sign, n_iter, 2, 3, 1);

    rot_givens_succ(A, A[1][1], A[2][1], sign, n_iter, 1, 2, 1);

    // Tercera columna

    rot_givens_succ(A, A[2][2], A[3][2], sign, n_iter, 2, 3, 2);

*/

    // Print result matrix
    for(int i = 0; i < TAM; i++){
		for(int j = 0; j < TAM; j++){

			std::cout << A[i][j] << "  |  ";
		}
		std:: cout << std::endl;
	}

    return 0;
}