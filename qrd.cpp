#include <iostream>
#include "qrd.h"

void rot_givens(/*const Matrix &A, const Matrix &R, */Matrix &res, bool sign){
	data_t coord_X, coord_Y;
	data_t coord_X_shif, coord_Y_shif;

	int n_iter = 1;
	
    // Accesing coordinates X and Y

	// The access to elements X and Y is incorrect, got to fix it
     for(int j = 0; j < TAM-1; j++){ // Columns
	    for(int i = TAM-1; i >= 0; i--){ // Rows
            if(i > j){

                coord_X = res[i-1][j];
                coord_Y = res[i][j];

                coord_X_shif = coord_X >> n_iter;
                coord_Y_shif = coord_Y >> n_iter;

                if(sign){
                	res[i-1][j] = coord_X + coord_Y_shif;
                	res[i][j] = coord_Y - coord_X_shif;
                }else{
                	res[i-1][j] = coord_X - coord_Y_shif;
                    res[i][j] = coord_Y + coord_X_shif;
                }

                sign = !sign;
                n_iter++;
            }
	    }
    }
}
