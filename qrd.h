#ifndef _QRD_
#define _QRD_

<<<<<<< HEAD
#define TAM 5
=======
#define TAM 4
>>>>>>> f28c6b479ae2575b504ce2e2ebd85abf00dc5715

#include <ap_fixed.h>
#include <vector>

<<<<<<< HEAD

// The more bits it has, more shifts can be performed later, so the approximation to 0 will be more precise
// For bigger matrices, we need bigger data formats to be able to calculate the right result
typedef ap_fixed<24, 10, AP_RND> data_t; // 24 bits fixed point data, 10 for integer value and 14 for decimals

typedef std::vector<std::vector<data_t>> Matrix;

void rot_givens(data_t A[TAM][TAM]);

void rot_givens_succ(Matrix &A, data_t X, data_t Y, bool &sign, int n_iter, int row_X, int row_Y, int col);
=======
typedef ap_fixed<16, 7, AP_RND> data_t; // 16 bits fixed point data, 7 for integer value and 3 for decimals

typedef std::vector<std::vector<data_t>> Matrix;

void rot_givens(/*const Matrix &A, const Matrix &R, */Matrix &res, bool sign);
>>>>>>> f28c6b479ae2575b504ce2e2ebd85abf00dc5715

#endif
