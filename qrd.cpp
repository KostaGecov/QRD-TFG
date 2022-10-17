#include <iostream>
#include "qrd.h"


void rot_givens(data_t A[TAM][TAM]) {
    data_t coord_X, coord_Y;
    data_t coord_X_shif, coord_Y_shif;

    int n_iter = 15;
    bool sign;

    // Accesing coordinates X and Y
    rot_columns_for:
    for (int j = 0; j < TAM - 1; j++) { // Columns
    	rot_rows_for:
        for (int i = TAM - 1; i >= 0; i--) { // Rows
            if ((i > j) && (A[i][j] != 0)) {

                coord_X = A[i - 1][j];
                coord_Y = A[i][j];

                // Vectorization (Coordinate Y = 0)
                for (int k = 0; k < n_iter; k++) {
                    if (coord_Y < 0)
                        sign = true;
                    else
                        sign = false;

                    coord_X_shif = coord_X >> k;
                    coord_Y_shif = coord_Y >> k;

                    if (sign) {
                        coord_X = coord_X - coord_Y_shif;
                        coord_Y = coord_Y + coord_X_shif;
                    } else {
                        coord_X = coord_X + coord_Y_shif;
                        coord_Y = coord_Y - coord_X_shif;
                    }
                }
                A[i - 1][j] = coord_X;
                A[i][j] = coord_Y;

                coord_X = 0;
                coord_Y = 0;
                coord_X_shif = 0;
                coord_Y_shif = 0;
            }

        }
    }
}

void rot_givens_succ(Matrix & A, data_t X, data_t Y, bool & sign, int n_iter, int row_X, int row_Y, int col) {
    data_t coord_X, coord_Y;
    data_t coord_X_shif, coord_Y_shif;

    coord_X = X;
    coord_Y = Y;

    // Vectorization (Coordinate Y = 0)
    for (int i = 0; i < n_iter; i++) {
        if (coord_Y < 0)
            sign = false;
        else
            sign = true;

        coord_X_shif = coord_X >> i;
        coord_Y_shif = coord_Y >> i;

        if (sign) {
            coord_X = coord_X - coord_Y_shif;
            coord_Y = coord_Y + coord_X_shif;
        } else {
            coord_X = coord_X + coord_Y_shif;
            coord_Y = coord_Y - coord_X_shif;
        }
    }

    A[row_X][col] = coord_X;
    A[row_Y][col] = coord_Y;

}
