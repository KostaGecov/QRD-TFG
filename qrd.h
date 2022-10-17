#ifndef _QRD_
#define _QRD_

#define TAM 5

#include <ap_fixed.h>
#include <vector>


// The more bits it has, more shifts can be performed later, so the approximation to 0 will be more precise
// For bigger matrices, we need bigger data formats to be able to calculate the right result
typedef ap_fixed<24, 10, AP_RND> data_t; // 24 bits fixed point data, 10 for integer value and 14 for decimals

typedef std::vector<std::vector<data_t>> Matrix;

void rot_givens(data_t A[TAM][TAM]);

void rot_givens_succ(Matrix &A, data_t X, data_t Y, bool &sign, int n_iter, int row_X, int row_Y, int col);

#endif
