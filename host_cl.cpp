#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_stream.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "xcl2.hpp"

#define FIXED_POINT 32
#define FX_POINT_INT 6

#define TAM_TILED 8
#define TAM 256
#define NUM_TILED 32
#define FLATTEN_SIZE TAM *TAM_TILED

#define NUM_OPERACIONES 65  // (63 + 2(offset))

#define GEQRT 0
#define TTQRT 1

/**
 * @brief 32 bits fixed point data, 6 for integer value and 26 for decimals.
 * The more bits it has, more shifts can be performed later, so the approximation to 0 will be more precise.
 * For bigger matrices, we need bigger data formats to be able to calculate the right result.
 *
 */
typedef ap_fixed<FIXED_POINT, FX_POINT_INT, AP_RND> data_t;

/**
 * @brief initialize data_t type matrix with values from input file
 *
 * @param matrix data_t type values
 * @param file input file
 */
bool init_matrix(data_t matrix[TAM][TAM], std::fstream *file);

/**
 * @brief initialize float type matrix with values from input file
 *
 * @param matrix float type values
 * @param file input file
 */
bool init_matrix(float matrix[TAM][TAM], std::fstream *file);

/**
 * @brief calculate mean squared error of the result
 *
 * @param A
 * @param out_gold
 * @return float error
 */
float mse(data_t A[TAM][TAM], float out_gold[TAM][TAM]);

/**
 * @brief prior to kernel execution
 *
 * @param qrd_kernel kernel to execute
 * @param q command queue for the device
 * @param context context of the device
 */
bool tiled_qr_decomposition(cl::Kernel qrd_kernel, cl::CommandQueue q, cl::Context context);

/**
 * @brief executes kernel
 *
 * @param A_matrix_ptr
 * @param qrd_kernel
 * @param q
 */
void kernel_execute(data_t A_tiled[NUM_TILED][TAM_TILED][TAM], uint8_t type_op, uint8_t col_offset,
                    uint8_t idx_mat_1, uint8_t idx_mat_2, cl::Kernel qrd_kernel, cl::CommandQueue q, cl::Context context);

void flatten_matrix(data_t matrix[NUM_TILED][TAM_TILED][TAM], data_t fl_matrix[FLATTEN_SIZE], uint8_t idx_mat);

void unflatten_matrix(data_t matrix[NUM_TILED][TAM_TILED][TAM], data_t fl_matrix[FLATTEN_SIZE], uint8_t idx_mat);

/**
 * offset to access the right column for GEQRT operation
 */
static uint16_t col_offset_geqrt = 0;
/**
 * offset to access the right column for TTQRT operation
 */
static uint16_t col_offset_ttqrt = 0;
/**
 * to control the GEQRT operations in each step
 */
static uint16_t n_iter_GEQRT = 32;
/**
 * to control the TTQRT operations in each step
 */
static uint16_t n_iter_TTQRT = 31;

/**
 * stores all 32 A matrices needed for tiled operations
 *
 */
data_t A_tiled[NUM_TILED][TAM_TILED][TAM];
data_t A[TAM][TAM];

/**
 * to store the output data gold
 *
 */
float out_gold[TAM][TAM];

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }

    // XCLBIN file to program the FPGA
    std::string binaryFile = argv[1];
    cl_int err;
    cl::Context context;
    cl::CommandQueue q;
    cl::Program program;
    cl::Kernel qrd_kernel;

    // OPENCL HOST CODE AREA START
    auto start1 = std::chrono::high_resolution_clock::now();

    // get_xil_devices() is a utility API which will find the Xilinx
    // platforms and will return list of devices connected to Xilinx platform
    auto devices = xcl::get_xil_devices();

    // read_binary_file() is a utility API which will load the binaryFile
    // and will return the pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binaryFile);

    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};

    bool valid_device = false;

    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));

        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        OCL_CHECK(err, program = cl::Program(context, {device}, bins, nullptr, &err));

        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            OCL_CHECK(err, qrd_kernel = cl::Kernel(program, "kernel_givens_rotation", &err));
            valid_device = true;
            break;  // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    // Finished creating context, command queue and kernel
    auto end1 = std::chrono::high_resolution_clock::now();
    auto diff1 = end1 - start1;
    std::cout << "FPGA programming time: " << std::chrono::duration<double, std::milli>(diff1).count() << std::endl;

    // Now, we start to execute the kernel
    auto start = std::chrono::high_resolution_clock::now();
    if (!tiled_qr_decomposition(qrd_kernel, q, context))
        return -1;
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = end - start;
    std::cout << "Execution time: " << std::chrono::duration<double, std::milli>(diff).count();

    return 0;
}

bool tiled_qr_decomposition(cl::Kernel qrd_kernel, cl::CommandQueue q, cl::Context context) {
    // TODO: pasar path ficheros .dat como argumento entrada
    std::fstream data_in("/home/kgecov/workspace_vitis/qrd_system/data_in.dat", std::ios::in);
    if (!init_matrix(A, &data_in))
        return false;

    std::fstream data_out_gold("/home/kgecov/workspace_vitis/qrd_system/data_out_gold.dat", std::ios::in);
    if (!init_matrix(out_gold, &data_out_gold))
        return false;

divide_matrices_row_for:
    for (uint16_t r = 0; r < TAM; r++) {
        int tile_index = r / TAM_TILED;
        int tile_offset = r % TAM_TILED;

    divide_matrices_col_for:
        for (uint16_t c = 0; c < TAM; c++) {
            A_tiled[tile_index][tile_offset][c] = A[r][c];
        }
    }

num_operations_for:
    for (uint16_t i = 2; i < NUM_OPERACIONES; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            switch (n_iter_GEQRT) {
                case 32:
                    for (uint16_t idx_mat_1 = 0; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 31:
                    for (uint16_t idx_mat_1 = 1; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 30:
                    for (uint16_t idx_mat_1 = 2; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 29:
                    for (uint16_t idx_mat_1 = 3; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 28:
                    for (uint16_t idx_mat_1 = 4; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 27:
                    for (uint16_t idx_mat_1 = 5; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 26:
                    for (uint16_t idx_mat_1 = 6; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 25:
                    for (uint16_t idx_mat_1 = 7; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 24:
                    for (uint16_t idx_mat_1 = 8; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 23:
                    for (uint16_t idx_mat_1 = 9; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 22:
                    for (uint16_t idx_mat_1 = 10; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 21:
                    for (uint16_t idx_mat_1 = 11; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 20:
                    for (uint16_t idx_mat_1 = 12; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 19:
                    for (uint16_t idx_mat_1 = 13; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 18:
                    for (uint16_t idx_mat_1 = 14; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 17:
                    for (uint16_t idx_mat_1 = 15; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 16:
                    for (uint16_t idx_mat_1 = 16; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 15:
                    for (uint16_t idx_mat_1 = 17; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 14:
                    for (uint16_t idx_mat_1 = 18; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 13:
                    for (uint16_t idx_mat_1 = 19; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 12:
                    for (uint16_t idx_mat_1 = 20; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 11:
                    for (uint16_t idx_mat_1 = 21; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 10:
                    for (uint16_t idx_mat_1 = 22; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 9:
                    for (uint16_t idx_mat_1 = 23; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 8:
                    for (uint16_t idx_mat_1 = 24; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 7:
                    for (uint16_t idx_mat_1 = 25; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 6:
                    for (uint16_t idx_mat_1 = 26; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 5:
                    for (uint16_t idx_mat_1 = 27; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 4:
                    for (uint16_t idx_mat_1 = 28; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 3:
                    for (uint16_t idx_mat_1 = 29; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 2:
                    for (uint16_t idx_mat_1 = 30; idx_mat_1 < 32; idx_mat_1++) {
                        kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context);
                    }

                    break;
                case 1:
                    kernel_execute(A_tiled, GEQRT, col_offset_geqrt, 31, 0, qrd_kernel, q, context);

                    break;
                default:
                    break;
            }
            n_iter_GEQRT--;
            col_offset_geqrt += TAM_TILED;

        } else {
            switch (n_iter_TTQRT) {
                case 31:
                    for (uint16_t idx_mat_1 = 0, idx_mat_2 = 1; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 0, idx_mat_2 = 2; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 0, idx_mat_2 = 4; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 0, idx_mat_2 = 8; idx_mat_2 < 25; idx_mat_1 += 16, idx_mat_2 += 16) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 0, 16, qrd_kernel, q, context);

                    break;
                case 30:
                    for (uint16_t idx_mat_1 = 1, idx_mat_2 = 2; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 1, idx_mat_2 = 3; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 1, idx_mat_2 = 5; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 1, idx_mat_2 = 9; idx_mat_2 < 26; idx_mat_1 += 16, idx_mat_2 += 16) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 1, 17, qrd_kernel, q, context);

                    break;
                case 29:
                    for (uint16_t idx_mat_1 = 2, idx_mat_2 = 3; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 2, idx_mat_2 = 4; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 2, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 2, idx_mat_2 = 10; idx_mat_2 < 27; idx_mat_1 += 16, idx_mat_2 += 16) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 2, 18, qrd_kernel, q, context);

                    break;
                case 28:
                    for (uint16_t idx_mat_1 = 3, idx_mat_2 = 4; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 3, idx_mat_2 = 5; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 3, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 3, idx_mat_2 = 11; idx_mat_2 < 28; idx_mat_1 += 16, idx_mat_2 += 16) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 3, 19, qrd_kernel, q, context);

                    break;
                case 27:
                    for (uint16_t idx_mat_1 = 4, idx_mat_2 = 5; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 4, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 4, idx_mat_2 = 8; idx_mat_2 < 25; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 4, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 16, idx_mat_2 += 16) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 4, 20, qrd_kernel, q, context);

                    break;
                case 26:
                    for (uint16_t idx_mat_1 = 5, idx_mat_2 = 6; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 5, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 5, idx_mat_2 = 9; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 5, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 16, idx_mat_2 += 16) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 5, 21, qrd_kernel, q, context);

                    break;
                case 25:
                    for (uint16_t idx_mat_1 = 6, idx_mat_2 = 7; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 6, idx_mat_2 = 8; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 6, idx_mat_2 = 10; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 6, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 16, idx_mat_2 += 16) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 6, 22, qrd_kernel, q, context);

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 6, 26, qrd_kernel, q, context);

                    break;
                case 24:
                    for (uint16_t idx_mat_1 = 7, idx_mat_2 = 8; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 7, idx_mat_2 = 9; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 7, idx_mat_2 = 11; idx_mat_2 < 28; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 7, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 16, idx_mat_2 += 16) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 7, 23, qrd_kernel, q, context);

                    break;
                case 23:
                    for (uint16_t idx_mat_1 = 8, idx_mat_2 = 9; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 8, idx_mat_2 = 10; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 8, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 8, idx_mat_2 = 16; idx_mat_2 < 25; idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 22:
                    for (uint16_t idx_mat_1 = 9, idx_mat_2 = 10; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 9, idx_mat_2 = 11; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 9, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 9, idx_mat_2 = 17; idx_mat_2 < 26; idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 21:
                    for (uint16_t idx_mat_1 = 10, idx_mat_2 = 11; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 10, idx_mat_2 = 12; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 10, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 10, idx_mat_2 = 18; idx_mat_2 < 27; idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 20:
                    for (uint16_t idx_mat_1 = 11, idx_mat_2 = 12; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 11, idx_mat_2 = 13; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 11, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 11, idx_mat_2 = 19; idx_mat_2 < 28; idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 19:
                    for (uint16_t idx_mat_1 = 12, idx_mat_2 = 13; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 12, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 12, idx_mat_2 = 16; idx_mat_2 < 25; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 12, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 18:
                    for (uint16_t idx_mat_1 = 13, idx_mat_2 = 14; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 13, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 13, idx_mat_2 = 17; idx_mat_2 < 26; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 13, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 17:
                    for (uint16_t idx_mat_1 = 14, idx_mat_2 = 15; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 14, idx_mat_2 = 16; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 14, idx_mat_2 = 18; idx_mat_2 < 27; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 14, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 16:
                    for (uint16_t idx_mat_1 = 15, idx_mat_2 = 16; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 15, idx_mat_2 = 17; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 15, idx_mat_2 = 19; idx_mat_2 < 28; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 15, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 15:
                    for (uint16_t idx_mat_1 = 16, idx_mat_2 = 17; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 16, idx_mat_2 = 18; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 16, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 16, 24, qrd_kernel, q, context);

                    break;
                case 14:
                    for (uint16_t idx_mat_1 = 17, idx_mat_2 = 18; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 17, idx_mat_2 = 19; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 17, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 17, 25, qrd_kernel, q, context);

                    break;
                case 13:
                    for (uint16_t idx_mat_1 = 18, idx_mat_2 = 19; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 18, idx_mat_2 = 20; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 18, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 18, 26, qrd_kernel, q, context);

                    break;
                case 12:
                    for (uint16_t idx_mat_1 = 19, idx_mat_2 = 20; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 19, idx_mat_2 = 21; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 19, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 8, idx_mat_2 += 8) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 19, 27, qrd_kernel, q, context);

                    break;
                case 11:
                    for (uint16_t idx_mat_1 = 20, idx_mat_2 = 21; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 20, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 20, idx_mat_2 = 24; idx_mat_2 < 29; idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 10:
                    for (uint16_t idx_mat_1 = 21, idx_mat_2 = 22; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 21, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 21, idx_mat_2 = 25; idx_mat_2 < 30; idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 9:
                    for (uint16_t idx_mat_1 = 22, idx_mat_2 = 23; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 22, idx_mat_2 = 24; idx_mat_2 < 29; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 22, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 8:
                    for (uint16_t idx_mat_1 = 23, idx_mat_2 = 24; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 23, idx_mat_2 = 25; idx_mat_2 < 30; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 23, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 7:
                    for (uint16_t idx_mat_1 = 24, idx_mat_2 = 25; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 24, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 24, 28, qrd_kernel, q, context);

                    break;
                case 6:
                    for (uint16_t idx_mat_1 = 25, idx_mat_2 = 26; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 25, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_1 += 4, idx_mat_2 += 4) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 25, 29, qrd_kernel, q, context);

                    break;
                case 5:
                    for (uint16_t idx_mat_1 = 26, idx_mat_2 = 27; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 26, idx_mat_2 = 28; idx_mat_2 < 31; idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 4:
                    for (uint16_t idx_mat_1 = 27, idx_mat_2 = 28; idx_mat_2 < 31; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    for (uint16_t idx_mat_1 = 27, idx_mat_2 = 29; idx_mat_2 < 32; idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 3:
                    for (uint16_t idx_mat_1 = 28, idx_mat_2 = 29; idx_mat_2 < 32; idx_mat_1 += 2, idx_mat_2 += 2) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 28, 30, qrd_kernel, q, context);

                    break;
                case 2:
                    for (uint16_t idx_mat_1 = 29, idx_mat_2 = 30; idx_mat_2 < 32; idx_mat_2++) {
                        kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, idx_mat_1, idx_mat_2, qrd_kernel, q, context);
                    }

                    break;
                case 1:
                    kernel_execute(A_tiled, TTQRT, col_offset_ttqrt, 30, 31, qrd_kernel, q, context);

                    break;
                default:
                    break;
            }
            n_iter_TTQRT--;
            col_offset_ttqrt += TAM_TILED;
        }
    }

write_sol_to_matrix_row_for:
    for (uint16_t r = 0; r < TAM; r++) {
        int tile_index = r / TAM_TILED;
        int tile_offset = r % TAM_TILED;
    write_sol_to_matrix_col_for:
        for (uint16_t c = 0; c < TAM; c++) {
            A[r][c] = A_tiled[tile_index][tile_offset][c];
        }
    }

    // Print R matrix
    std::cout << "R Matrix: " << std::endl;
    for (int i = 0; i < TAM; i++) {
        for (int j = 0; j < TAM; j++) {
            std::cout << A[i][j] << "  |  ";
        }
        std::cout << std::endl;
    }

    std::cout << "ECM = " << mse(A, out_gold) << std::endl;
    return true;
}

void kernel_execute(data_t A_tiled[NUM_TILED][TAM_TILED][TAM], uint8_t type_op, uint8_t col_offset, uint8_t idx_mat_1, uint8_t idx_mat_2, cl::Kernel qrd_kernel, cl::CommandQueue q, cl::Context context) {
    cl_int error;
    uint8_t aux;
    data_t flattened_matrix_1[FLATTEN_SIZE];
    data_t flattened_matrix_2[FLATTEN_SIZE];

    if (type_op == GEQRT) {
        // Flatten input matrix
        flatten_matrix(A_tiled, flattened_matrix_1, idx_mat_1);

        // Read and write pointers from kernel to host
        OCL_CHECK(error, cl::Buffer geqrt_input_matrix_ptr(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));
        // need to create unused buffer to pass it to kernel as argument because ttqrt uses 2 input buffers
        OCL_CHECK(error, cl::Buffer geqrt_input_unused_ptr(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, sizeof(uint8_t), NULL, &error));
        OCL_CHECK(error, cl::Buffer geqrt_output_matrix_ptr(context, CL_MEM_WRITE_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));
        // need to create unused buffer to pass it to kernel as argument because ttqrt uses 2 output buffers
        OCL_CHECK(error, cl::Buffer geqrt_output_unused_ptr(context, CL_MEM_READ_ONLY, sizeof(uint8_t), NULL, &error));

        // Set kernel arguments
        OCL_CHECK(error, error = qrd_kernel.setArg(0, sizeof(cl_mem), &geqrt_input_matrix_ptr));
        OCL_CHECK(error, error = qrd_kernel.setArg(1, sizeof(cl_mem), &geqrt_input_unused_ptr));
        OCL_CHECK(error, error = qrd_kernel.setArg(2, sizeof(cl_mem), &geqrt_output_matrix_ptr));
        OCL_CHECK(error, error = qrd_kernel.setArg(3, sizeof(cl_mem), &geqrt_output_unused_ptr));
        OCL_CHECK(error, error = qrd_kernel.setArg(4, type_op));
        OCL_CHECK(error, error = qrd_kernel.setArg(5, col_offset));

        // Send input buffer to kernel
        OCL_CHECK(error, error = q.enqueueWriteBuffer({geqrt_input_matrix_ptr},
                                                      CL_FALSE,
                                                      0,
                                                      FLATTEN_SIZE * sizeof(data_t),
                                                      flattened_matrix_1,
                                                      nullptr, nullptr));

        // Execute kernel
        OCL_CHECK(error, error = q.enqueueTask(qrd_kernel));

        // Wait for kernel to finish
        OCL_CHECK(error, error = q.finish());

        // Write output buffer data back to A_tiled
        OCL_CHECK(error, error = q.enqueueReadBuffer({geqrt_output_matrix_ptr},
                                                     CL_TRUE,
                                                     0,
                                                     FLATTEN_SIZE * sizeof(data_t),
                                                     flattened_matrix_1,
                                                     nullptr, nullptr));

        unflatten_matrix(A_tiled, flattened_matrix_1, idx_mat_1);

    } else if (type_op == TTQRT) {
        // Flatten input matrices
        flatten_matrix(A_tiled, flattened_matrix_1, idx_mat_1);
        flatten_matrix(A_tiled, flattened_matrix_2, idx_mat_2);

        // Create two kernel input buffers and two output buffer
        OCL_CHECK(error, cl::Buffer ttqrt_input_matrix_1_ptr(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));
        OCL_CHECK(error, cl::Buffer ttqrt_input_matrix_2_ptr(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));
        OCL_CHECK(error, cl::Buffer ttqrt_output_matrix_1_ptr(context, CL_MEM_WRITE_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));
        OCL_CHECK(error, cl::Buffer ttqrt_output_matrix_2_ptr(context, CL_MEM_WRITE_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));

        // Set kernel arguments
        OCL_CHECK(error, error = qrd_kernel.setArg(0, sizeof(cl_mem), &ttqrt_input_matrix_1_ptr));
        OCL_CHECK(error, error = qrd_kernel.setArg(1, sizeof(cl_mem), &ttqrt_input_matrix_2_ptr));
        OCL_CHECK(error, error = qrd_kernel.setArg(2, sizeof(cl_mem), &ttqrt_output_matrix_1_ptr));
        OCL_CHECK(error, error = qrd_kernel.setArg(3, sizeof(cl_mem), &ttqrt_output_matrix_2_ptr));
        OCL_CHECK(error, error = qrd_kernel.setArg(4, type_op));
        OCL_CHECK(error, error = qrd_kernel.setArg(5, col_offset));

        // Send input buffer to kernel
        OCL_CHECK(error, error = q.enqueueWriteBuffer({ttqrt_input_matrix_1_ptr},
                                                      CL_FALSE,
                                                      0,
                                                      FLATTEN_SIZE * sizeof(data_t),
                                                      flattened_matrix_1,
                                                      nullptr, nullptr));

        OCL_CHECK(error, error = q.enqueueWriteBuffer({ttqrt_input_matrix_2_ptr},
                                                      CL_FALSE,
                                                      0,
                                                      FLATTEN_SIZE * sizeof(data_t),
                                                      flattened_matrix_2,
                                                      nullptr, nullptr));

        // Execute kernel
        OCL_CHECK(error, error = q.enqueueTask(qrd_kernel));

        // Wait for kernel to finish
        OCL_CHECK(error, error = q.finish());

        // Write output buffer data back to A_tiled
        OCL_CHECK(error, error = q.enqueueReadBuffer({ttqrt_output_matrix_1_ptr},
                                                     CL_TRUE,
                                                     0,
                                                     FLATTEN_SIZE * sizeof(data_t),
                                                     flattened_matrix_1,
                                                     nullptr, nullptr));
        OCL_CHECK(error, error = q.enqueueReadBuffer({ttqrt_output_matrix_2_ptr},
                                                     CL_TRUE,
                                                     0,
                                                     FLATTEN_SIZE * sizeof(data_t),
                                                     flattened_matrix_2,
                                                     nullptr, nullptr));

        unflatten_matrix(A_tiled, flattened_matrix_1, idx_mat_1);
        unflatten_matrix(A_tiled, flattened_matrix_2, idx_mat_2);
    }
}

bool init_matrix(data_t matrix[TAM][TAM], std::fstream *file) {
    if (!file->is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return false;
    } else {
        std::cout << "Opened file" << std::endl;
    }

initialize_matrix:
    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = 0; c < TAM; c++) {
            *file >> matrix[r][c];
        }
    }
    file->close();
    return true;
}

bool init_matrix(float matrix[TAM][TAM], std::fstream *file) {
    if (!file->is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return false;
    } else {
        std::cout << "Opened file" << std::endl;
    }

initialize_matrix:
    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = 0; c < TAM; c++) {
            *file >> matrix[r][c];
        }
    }
    file->close();
    return true;
}

float mse(data_t A[TAM][TAM], float out_gold[TAM][TAM]) {
    float err = 0.0;
    data_t resA = 0.0;
    float resOut = 0.0;

    // to access just to the non zero elements (upper triangular matrix)
    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = r; c < TAM; c++) {
            resA = A[r][c];
            resOut = out_gold[r][c];
            err += pow((abs((float)resA) - abs(resOut)), 2);
        }
    }

    // err / number of elements in the upper triangular matrix (including diagonal)
    return err / ((TAM * (TAM + 1)) / 2);
}

void flatten_matrix(data_t matrix[NUM_TILED][TAM_TILED][TAM], data_t fl_matrix[FLATTEN_SIZE], uint8_t idx_mat) {
    // Flatten the 3D array into 1D array
    for (int j = 0; j < TAM_TILED; j++) {
        for (int k = 0; k < TAM; k++) {
            fl_matrix[j * TAM + k] = matrix[idx_mat][j][k];
        }
    }
}

void unflatten_matrix(data_t matrix[NUM_TILED][TAM_TILED][TAM], data_t fl_matrix[FLATTEN_SIZE], uint8_t idx_mat) {
    // Unflatten the 1D array into 3D matrix
    for (int j = 0; j < TAM_TILED; j++) {
        for (int k = 0; k < TAM; k++) {
            matrix[idx_mat][j][k] = fl_matrix[j * TAM + k];
        }
    }
}
