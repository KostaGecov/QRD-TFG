#define TAM 1024

#include <cmath>
#include <iostream>
#include <math.h>
#include <chrono>
#include <fstream>

const float EPSILON = 1e-10;

float A[TAM][TAM] = {0.0};

void init_matrix(float matrix[TAM][TAM], std::fstream* file) {
    if (!file->is_open()) {
        std::cerr << "Error opening file" << std::endl;
    } else {
        std::cout << "Opened file" << std::endl;
    }

    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = 0; c < TAM; c++) {
            *file >> matrix[r][c];
        }
    }
    file->close();
}

float mse(float A[TAM][TAM]) {
    float err = 0.0;
    float resA = 0.0;
    float resOut = 0.0;

//    std::fstream data_out_gold("data_out_gold.dat", std::ios::in);    
    std::fstream data_out_gold("random_matrix_1024_gold.dat", std::ios::in);
    if (!data_out_gold.is_open()) {
        std::cerr << "Error opening data_out_gold" << std::endl;
    } else {
        std::cout << "Opened file" << std::endl;
    }

    // to access just to the non zero elements (upper triangular matrix)
    for (uint16_t r = 0; r < TAM; r++) {
        //for (uint16_t c = r; c < TAM; c++) {
        for (uint16_t c = 0; c < TAM; c++) {
            resA = A[r][c];
            data_out_gold >> resOut;
            // resOut = out_gold[r][c];
            err += pow((abs(resA) - abs(resOut)), 2);
            resOut = 0;
        }
    }

    // err / number of elements in the upper triangular matrix (including diagonal)
    return err / ((TAM * (TAM + 1)) / 2);
}

void givens_rotation(float a, float b, float& c, float& s) {
    if (std::abs(b) < EPSILON) {
        c = 1;
        s = 0;
    } else {
        float tau;
        if (std::abs(b) > std::abs(a)) {
            tau = -a / b;
            s = 1 / std::sqrt(1 + tau * tau);
            c = s * tau;
        } else {
            tau = -b / a;
            c = 1 / std::sqrt(1 + tau * tau);
            s = c * tau;
        }
    }
}

void qr_decomposition_givens(float matrix[TAM][TAM], int n) {
    //float Q[TAM][TAM];
    //for (int i = 0; i < n; ++i) {
    //    for (int j = 0; j < n; ++j) {
    //        Q[i][j] = (i == j) ? 1 : 0;  // Inicializar matriz Q como la matriz identidad
    //   }
    //}

    for (int j = 0; j < n; ++j) {
        for (int i = n - 1; i > j; --i) {
            float c, s;
            givens_rotation(matrix[j][j], matrix[i][j], c, s);

            // Actualizar R
            float temp1, temp2;
            for (int k = 0; k < n; ++k) {
                temp1 = matrix[j][k];
                temp2 = matrix[i][k];
                matrix[j][k] = c * temp1 - s * temp2;
                matrix[i][k] = s * temp1 + c * temp2;
            }

            // Actualizar Q
            //for (int k = 0; k < n; ++k) {
            //    temp1 = Q[j][k];
            //    temp2 = Q[i][k];
            //    Q[j][k] = c * temp1 - s * temp2;
            //    Q[i][k] = s * temp1 + c * temp2;
            //}
        }
    }

//    std::cout << "Q:\n";
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            std::cout << Q[i][j] << " ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "R:\n";
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            std::cout << matrix[i][j] << " ";
//        }
//        std::cout << "\n";
//    }
}

int main() {
    std::cout << "Start of program" << std::endl;
    
    float err = 0.0f;
    
    auto start = std::chrono::high_resolution_clock::now();

    std::fstream data_in("random_matrix_1024.dat", std::ios::in);
//    std::fstream data_in("data_in.dat", std::ios::in);
    init_matrix(A, &data_in);

//    std::fstream data_out_gold("random_matrix_1024_1_gold.dat", std::ios::in);
//    std::fstream data_out_gold("data_out_gold.dat", std::ios::in);
//    init_matrix(out_gold, &data_out_gold);

    qr_decomposition_givens(A, TAM);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = end - start;
    err = mse(A);
    std::cout << "Mean Squared Error = " << err << std::endl;
    std::cout << "SW Execution time: " << std::chrono::duration<double, std::milli>(diff).count()<< std::endl;

    return 0;
}
