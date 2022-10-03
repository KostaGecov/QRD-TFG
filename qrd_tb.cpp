#include <iostream>
#include "qrd.h"


int main(){

	//typedef ap_fixed<16, 7, AP_RND> data_t; // 16 bits fixed point data, 7 for integer value and 3 for decimals
	Matrix A = {{3, 2, -1, 4},
			    {2, 1, 5, 7},
			    {0, 5, 2, -6},
			    {-1, 2, 1, 0}};

	Matrix R = {{1, 0, 0, 0},
			    {0, 1, 0, 0},
			    {0, 0, 1, 0},
			    {0, 0, 0, 1}};

	Matrix res(TAM, std::vector<data_t>(TAM, 0)); // Upper-Triangular Matrix R

	for(int i = 0; i < TAM; i++){
		for(int j = 0; j < TAM; j++){
			res[i][j] = A[i][j];
		}
	}


    bool sign = (A[TAM-1][0] >= 0);

    rot_givens(/*A, R, */res, sign);


    // Print result matrix
    for(int i = 0; i < TAM; i++){
		for(int j = 0; j < TAM; j++){
			std::cout << res[i][j] << " ";
		}
		std:: cout << std::endl;
	}
    return 0;
}
