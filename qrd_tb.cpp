#include <fstream>
#include <iostream>

// OpenCL utility layer include
#include "qrd.h"
#include "xcl2.hpp"

static index_t col_offset = 0;    // offset to access the right column for GEQRT operation
static index_t n_iter_GEQRT = 4;  // 4 iteraciones iniciales, para controlar el
                                  // nº de operaciones GEQRT en cada paso
static index_t n_iter_TTQRT = 3;  // 3 iteraciones iniciales, para controlar el nº de operaciones TTQRT en
                                  // cada paso para i=2, 4; para i = 4, 3; para i = 6, 2; para i = 8, 1;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }

    auto binaryFile = argv[1];

    cl_int err;
    cl::Context context;
    cl::CommandQueue q;
    cl::Kernel krnl_givens_rotation;

    auto devices = xcl::get_xil_devices();

    // Create Program and Kernel
    auto fileBuf = xcl::read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    bool valid_device = false;
    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));

        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            OCL_CHECK(err, krnl_givens_rotation = cl::Kernel(program, "givens_rotation_kernel", &err));
            valid_device = true;
            break;  // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    data_t A_file[TAM][TAM];

    std::ifstream data_in("data_in.txt");

    if (!data_in.is_open()) {
        // Error
    }

    for (index_t r = 0; r < TAM; r++) {
        for (int c = 0; c < TAM; c++) {
            data_in >> A_file[r][c];
        }
    }

    data_t A[TAM][TAM] = {
        {-1, 5, -6, 0, -2, 8, -1, -2, 2, -6, -4, 9, -2, -3, -7, -8, -8, 8, -6, -8, 5, -8, -9, 5},
        {-1, 2, -1, -9, -7, 6, -9, 5, -4, -7, 6, -7, 8, -4, -1, 1, -1, -4, 7, 5, 1, -8, 5, 7},
        {5, 0, 2, -8, 8, 7, 9, -9, -1, -6, 4, -6, 3, 9, 2, -4, 2, 7, -2, -7, 7, -6, -1, -1},
        {-8, -7, 6, 9, -6, 0, -3, -4, -9, -3, -4, -1, -3, 6, 2, -3, 2, -4, -3, 8, 3, -2, -5, -7},
        {3, 0, 6, -7, 7, -9, -5, 7, -4, 9, 5, 0, 0, -1, 7, -2, 2, 0, 3, 4, -6, -6, -6, 5},
        {-6, 1, 4, 0, 4, 2, -10, 9, 1, 8, 3, 7, 1, -7, -3, 1, -6, 9, -10, 4, 3, -6, 3, -6},
        {-2, -1, 5, -5, -4, -8, -3, 4, -5, -8, -2, -4, -8, 0, -1, 8, 8, 6, -2, 4, -10, -1, -4, 0},
        {1, 2, 5, 0, -3, -2, 4, 6, -1, -10, -3, -5, -3, -9, 7, -3, 7, -4, -1, 6, -7, 1, 1, -5},
        {-7, 1, 7, 1, -5, 1, 7, -2, 2, -8, 6, 7, 6, 1, -8, -3, 7, 8, -5, -6, 1, -10, 9, 1},
        {-9, 4, -5, -5, -5, -10, -9, 3, -9, 3, -9, 0, -9, 5, -9, -1, 9, 8, -6, 4, -3, 7, 0, 1},
        {-5, 0, 5, -6, 4, 8, 6, 5, 3, 3, -9, -3, -7, -2, -6, 0, 6, -3, -4, -4, 8, -4, -4, 6},
        {2, 3, -3, 4, -9, -3, 3, 7, -6, -3, -9, -1, 9, -3, 4, -9, -10, 2, 4, 4, 3, 9, 3, 8},
        {-2, 1, -9, 3, 0, -5, 0, -3, -4, -7, -10, -4, -5, -9, 6, -9, 7, -4, 4, 5, -7, 8, 0, 1},
        {7, 0, -9, -4, 5, 9, 4, 5, 8, 4, 2, -10, -9, 0, 6, -1, -2, -10, 5, 9, -1, -5, -7, -1},
        {-2, -4, 5, -6, 6, -10, 2, -7, -10, -10, -7, -1, -9, -6, -7, 6, 6, -3, 1, 7, 4, 7, -3, -3},
        {3, -1, -1, -9, -1, 3, -10, -4, 4, -9, -7, -10, -5, -8, -6, -4, 0, -4, -6, 2, 0, 2, 6, -9},
        {-8, -5, 6, -8, 3, 8, -7, 4, 3, 3, -4, -6, 6, -2, -3, 9, -3, 0, 2, -10, -4, -7, 2, -5},
        {-4, 8, -9, 2, -7, -10, 9, -2, 9, 9, 3, -3, 2, -6, 8, 5, 5, -8, -7, -8, -6, -10, -4, -9},
        {4, 0, 5, 8, 0, -7, -7, -7, 3, 4, -3, -6, -3, 6, 3, -4, -8, 0, 5, -6, -5, 5, -1, -6},
        {4, 0, 9, -7, -6, 8, -8, -5, 2, 5, -8, 5, -5, 8, 7, -4, -7, -1, 7, -5, -5, -3, -10, 5},
        {-2, -5, 7, 0, -7, 9, 2, -8, 5, -10, 8, 7, -9, -10, -9, -5, -7, -2, -1, -6, 5, 0, 9, 4},
        {-9, -8, -4, -10, -9, 3, -7, 4, -1, -10, -7, -1, -3, -3, -6, 2, 0, 9, -8, -8, -5, -9, -9, 0},
        {-6, -2, 5, 4, 8, 5, -9, -9, -10, 2, -9, -7, -1, -10, -5, 3, -3, -8, 9, -3, -8, 7, -3, -1},
        {-5, 9, 1, -10, -3, -5, -7, -1, -10, -5, -9, -8, 4, -2, -1, -3, 9, 1, -6, -2, -9, -3, 7, -4}};

    data_t A_aux[TAM][TAM];

    data_t A_1[TAM_TILED][TAM];
    data_t A_2[TAM_TILED][TAM];
    data_t A_3[TAM_TILED][TAM];
    data_t A_4[TAM_TILED][TAM];

divide_matrix_row_for:
    for (index_t r = 0; r < TAM; r++) {
    divide_matrix_col_for:
        for (index_t c = 0; c < TAM; c++) {
            if (r >= 0 && r <= 5) {
                A_1[r][c] = A[r][c];
            } else if (r >= 6 && r <= 11) {
                A_2[r - 6][c] = A[r][c];
            } else if (r >= 12 && r <= 17) {
                A_3[r - 12][c] = A[r][c];
            } else if (r >= 18 && r <= 23) {
                A_4[r - 18][c] = A[r][c];
            }
        }
    }

num_operations_for:
    for (index_t i = 2; i < NUM_OPERACIONES; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            switch (n_iter_GEQRT) {
                case 4:
                    krnl_givens_rotation.setArg(0, A_1);
                    krnl_givens_rotation.setArg(1, A_aux);
                    krnl_givens_rotation.setArg(2, GEQRT);
                    krnl_givens_rotation.setArg(3, col_offset);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_2);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_3);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_4);
                    q.enqueueTask(krnl_givens_rotation);
                    break;
                case 3:
                    krnl_givens_rotation.setArg(0, A_2);
                    krnl_givens_rotation.setArg(1, A_aux);
                    krnl_givens_rotation.setArg(2, GEQRT);
                    krnl_givens_rotation.setArg(3, col_offset);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_3);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_4);
                    q.enqueueTask(krnl_givens_rotation);
                    break;
                case 2:
                    krnl_givens_rotation.setArg(0, A_3);
                    krnl_givens_rotation.setArg(1, A_aux);
                    krnl_givens_rotation.setArg(2, GEQRT);
                    krnl_givens_rotation.setArg(3, col_offset);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_4);
                    q.enqueueTask(krnl_givens_rotation);
                    break;
                case 1:
                    krnl_givens_rotation.setArg(0, A_4);
                    krnl_givens_rotation.setArg(1, A_aux);
                    krnl_givens_rotation.setArg(2, GEQRT);
                    krnl_givens_rotation.setArg(3, col_offset);
                    q.enqueueTask(krnl_givens_rotation);
                    break;
                default:
                    break;
            }
            n_iter_GEQRT--;
            col_offset += 6;
        } else {
            switch (n_iter_TTQRT) {
                case 3:
                    krnl_givens_rotation.setArg(0, A_1);
                    krnl_givens_rotation.setArg(1, A_2);
                    krnl_givens_rotation.setArg(2, TTQRT);
                    krnl_givens_rotation.setArg(3, 0);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_3);
                    krnl_givens_rotation.setArg(1, A_4);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_1);
                    krnl_givens_rotation.setArg(1, A_3);
                    q.enqueueTask(krnl_givens_rotation);
                    break;
                case 2:
                    krnl_givens_rotation.setArg(0, A_2);
                    krnl_givens_rotation.setArg(1, A_3);
                    q.enqueueTask(krnl_givens_rotation);
                    krnl_givens_rotation.setArg(0, A_2);
                    krnl_givens_rotation.setArg(1, A_4);
                    q.enqueueTask(krnl_givens_rotation);
                    break;
                case 1:
                    krnl_givens_rotation.setArg(0, A_3);
                    krnl_givens_rotation.setArg(1, A_4);
                    q.enqueueTask(krnl_givens_rotation);
                    break;
                default:
                    break;
            }
            n_iter_TTQRT--;
        }
    }

    // Wait for all the tasks in the queue to complete
    q.finish();

// Escritura de solución en matriz grande
write_sol_to_matrix_row_for:
    for (index_t r = 0; r < TAM; r++) {
    write_sol_to_matrix_col_for:
        for (index_t c = 0; c < TAM; c++) {
            if (r >= 0 && r <= 5) {
                A[r][c] = A_1[r][c];
            } else if (r >= 6 && r <= 11) {
                A[r][c] = A_2[r - 6][c];
            } else if (r >= 12 && r <= 17) {
                A[r][c] = A_3[r - 12][c];
            } else if (r >= 18 && r <= 23) {
                A[r][c] = A_4[r - 18][c];
            }
        }
    }

    // Print result matrix
    for (int i = 0; i < TAM; i++) {
        for (int j = 0; j < TAM; j++) {
            std::cout << A[i][j] << "  |  ";
        }
        std::cout << std::endl;
    }
    return 0;
}
