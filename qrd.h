#ifndef _QRD_
#define _QRD_

#define TAM 4

#include <ap_fixed.h>
#include <vector>

typedef ap_fixed<16, 7, AP_RND> data_t; // 16 bits fixed point data, 7 for integer value and 3 for decimals

typedef std::vector<std::vector<data_t>> Matrix;

void rot_givens(/*const Matrix &A, const Matrix &R, */Matrix &res, bool sign);

#endif
