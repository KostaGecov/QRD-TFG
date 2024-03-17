#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_stream.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "xcl2.hpp"

#define FIXED_POINT 24
#define FX_POINT_INT 6

#define TAM_TILED 8
#define TAM 256
#define NUM_TILED (TAM / TAM_TILED)
#define FLATTEN_SIZE (TAM * TAM_TILED)

#define NUM_OPERATIONS 63
#define NUM_OP_GEQRT 32
#define NUM_OP_TTQRT 31

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
void tiled_qr_decomposition(data_t A_tiled[NUM_TILED][TAM_TILED][TAM], cl::Kernel qrd_kernel, cl::CommandQueue q, cl::Context context);

/**
 * @brief executes kernel
 *
 * @param A_matrix_ptr
 * @param qrd_kernel
 * @param q
 */
void kernel_execute(data_t A_tiled[NUM_TILED][TAM_TILED][TAM], uint8_t type_op, uint8_t col_offset, uint8_t idx_mat_1, uint8_t idx_mat_2, cl::Kernel qrd_kernel, cl::CommandQueue q, cl::Context context,
                    cl::Buffer input_matrix_1_ptr, cl::Buffer input_matrix_2_ptr, cl::Buffer output_matrix_1_ptr, cl::Buffer output_matrix_2_ptr);

/**
 * @brief Writes elements of multidimensiona matrix in order into an array
 *
 * @param matrix
 * @param fl_matrix
 * @param idx_mat
 */
void flatten_matrix(data_t matrix[NUM_TILED][TAM_TILED][TAM], data_t fl_matrix[FLATTEN_SIZE], uint8_t idx_mat);

/**
 * @brief Writes elements of result array from kernel in order back into the multidimensional matrix
 *
 * @param matrix
 * @param fl_matrix
 * @param idx_mat
 */
void unflatten_matrix(data_t matrix[NUM_TILED][TAM_TILED][TAM], data_t fl_matrix[FLATTEN_SIZE], uint8_t idx_mat);

// Offset to access the right column for GEQRT operation
static uint16_t col_offset_geqrt = 0;

// Offset to access the right column for TTQRT operation
static uint16_t col_offset_ttqrt = 0;

// To control the GEQRT operations in each step
static uint16_t n_iter_GEQRT = 32;

// To control the TTQRT operations in each step
static uint16_t n_iter_TTQRT = 31;

// Stores all 32 A matrices needed for tiled operations
data_t A_tiled[NUM_TILED][TAM_TILED][TAM];

// Stores input matrix
data_t A[TAM][TAM];

// Stores output data gold matrix
float out_gold[TAM][TAM];

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }

    // XCLBIN file to program the FPGA
    std::string binaryFile = argv[3];
    cl_int err;
    cl::Context context;
    cl::CommandQueue q;
    cl::Program program;
    cl::Kernel qrd_kernel;

    unsigned int tile_index = 0;
    unsigned int tile_offset = 0;

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

    std::fstream data_in(argv[1], std::ios::in);

    if (!init_matrix(A, &data_in))
        return -1;

    std::fstream data_out_gold(argv[2], std::ios::in);
    if (!init_matrix(out_gold, &data_out_gold))
        return -1;

    for (uint16_t r = 0; r < TAM; r++) {
        tile_index = r / TAM_TILED;
        tile_offset = r % TAM_TILED;

        for (uint16_t c = 0; c < TAM; c++) {
            A_tiled[tile_index][tile_offset][c] = A[r][c];
        }
    }

    // Now, we start to execute the kernel
    auto start = std::chrono::high_resolution_clock::now();
    tiled_qr_decomposition(A_tiled, qrd_kernel, q, context);
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = end - start;

    for (uint16_t r = 0; r < TAM; r++) {
        tile_index = r / TAM_TILED;
        tile_offset = r % TAM_TILED;

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

    std::cout << "Mean Squared Error = " << mse(A, out_gold) << std::endl;
    std::cout << "FPGA Execution time: " << std::chrono::duration<double, std::milli>(diff).count();

    return 0;
}

void tiled_qr_decomposition(data_t A_tiled[NUM_TILED][TAM_TILED][TAM], cl::Kernel qrd_kernel, cl::CommandQueue q, cl::Context context) {
    cl_int error = 0;

    // Create two kernel input buffers and two output buffer
    OCL_CHECK(error, cl::Buffer input_matrix_1_ptr(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_HOST_WRITE_ONLY | CL_MEM_READ_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));
    OCL_CHECK(error, cl::Buffer input_matrix_2_ptr(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_HOST_WRITE_ONLY | CL_MEM_READ_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));
    OCL_CHECK(error, cl::Buffer output_matrix_1_ptr(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_HOST_READ_ONLY | CL_MEM_WRITE_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));
    OCL_CHECK(error, cl::Buffer output_matrix_2_ptr(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_HOST_READ_ONLY | CL_MEM_WRITE_ONLY, FLATTEN_SIZE * sizeof(data_t), NULL, &error));

    // Set kernel arguments
    OCL_CHECK(error, error = qrd_kernel.setArg(0, sizeof(cl_mem), &input_matrix_1_ptr));
    OCL_CHECK(error, error = qrd_kernel.setArg(1, sizeof(cl_mem), &input_matrix_2_ptr));
    OCL_CHECK(error, error = qrd_kernel.setArg(2, sizeof(cl_mem), &output_matrix_1_ptr));
    OCL_CHECK(error, error = qrd_kernel.setArg(3, sizeof(cl_mem), &output_matrix_2_ptr));

    for (uint16_t i = 0; i < NUM_OPERATIONS; i++) {
        // GEQRT operation
        if (i % 2 == 0) {
            for (uint16_t idx_mat_1 = NUM_OP_GEQRT - n_iter_GEQRT; idx_mat_1 < NUM_TILED; idx_mat_1++) {
                kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context,
                               input_matrix_1_ptr, input_matrix_2_ptr, output_matrix_1_ptr, output_matrix_2_ptr);
            }

            // if (n_iter_GEQRT > 1) {
            //     kernel_execute(A_tiled, GEQRT, col_offset_geqrt + TAM_TILED, 0, 0, qrd_kernel, q, context,
            //                    input_matrix_1_ptr, input_matrix_2_ptr, output_matrix_1_ptr, output_matrix_2_ptr);
            // }
            n_iter_GEQRT--;
            col_offset_geqrt += TAM_TILED;

        } else {
            for (uint16_t idx_mat_1 = (NUM_OP_TTQRT - n_iter_TTQRT), idx_mat_2 = ((NUM_OP_TTQRT - n_iter_TTQRT) + 1); idx_mat_2 < NUM_TILED; idx_mat_2++) {
                kernel_execute(A_tiled, GEQRT, col_offset_geqrt, idx_mat_1, 0, qrd_kernel, q, context,
                               input_matrix_1_ptr, input_matrix_2_ptr, output_matrix_1_ptr, output_matrix_2_ptr);
            }

            n_iter_TTQRT--;
            col_offset_ttqrt += TAM_TILED;
        }
    }
}

void kernel_execute(data_t A_tiled[NUM_TILED][TAM_TILED][TAM], uint8_t type_op, uint8_t col_offset, uint8_t idx_mat_1, uint8_t idx_mat_2, cl::Kernel qrd_kernel, cl::CommandQueue q, cl::Context context,
                    cl::Buffer input_matrix_1_ptr, cl::Buffer input_matrix_2_ptr, cl::Buffer output_matrix_1_ptr, cl::Buffer output_matrix_2_ptr) {
    cl_int error;
    data_t *flattened_matrix_write_1 = (data_t *)q.enqueueMapBuffer(input_matrix_1_ptr, CL_TRUE, CL_MAP_WRITE, 0, FLATTEN_SIZE * sizeof(data_t), nullptr, nullptr, &error);
    data_t *flattened_matrix_write_2 = (data_t *)q.enqueueMapBuffer(input_matrix_2_ptr, CL_TRUE, CL_MAP_WRITE, 0, FLATTEN_SIZE * sizeof(data_t), nullptr, nullptr, &error);
    data_t *flattened_matrix_read_1 = (data_t *)q.enqueueMapBuffer(output_matrix_1_ptr, CL_TRUE, CL_MAP_READ, 0, FLATTEN_SIZE * sizeof(data_t), nullptr, nullptr, &error);
    data_t *flattened_matrix_read_2 = (data_t *)q.enqueueMapBuffer(output_matrix_2_ptr, CL_TRUE, CL_MAP_READ, 0, FLATTEN_SIZE * sizeof(data_t), nullptr, nullptr, &error);

    OCL_CHECK(error, error = qrd_kernel.setArg(4, type_op));
    OCL_CHECK(error, error = qrd_kernel.setArg(5, col_offset));

    if (type_op == GEQRT) {
        // Flatten input matrix
        flatten_matrix(A_tiled, flattened_matrix_write_1, idx_mat_1);

        // Send input buffer to kernel
        OCL_CHECK(error, error = q.enqueueMigrateMemObjects({input_matrix_1_ptr}, 0));

        // Execute kernel
        OCL_CHECK(error, error = q.enqueueTask(qrd_kernel));

        // Wait for kernel to finish
        OCL_CHECK(error, error = q.finish());

        // Write output buffer data back to A_tiled
        OCL_CHECK(error, error = q.enqueueMigrateMemObjects({output_matrix_1_ptr}, CL_MIGRATE_MEM_OBJECT_HOST));

        // Wait for kernel to finish
        OCL_CHECK(error, error = q.finish());

        unflatten_matrix(A_tiled, flattened_matrix_read_1, idx_mat_1);

    } else if (type_op == TTQRT) {
        // Flatten input matrices
        flatten_matrix(A_tiled, flattened_matrix_write_1, idx_mat_1);
        flatten_matrix(A_tiled, flattened_matrix_write_2, idx_mat_2);

        // Send input buffer to kernel
        OCL_CHECK(error, error = q.enqueueMigrateMemObjects({input_matrix_1_ptr, input_matrix_2_ptr}, 0));
        // OCL_CHECK(error, error = q.enqueueMigrateMemObjects({input_matrix_2_ptr}, 0));

        // Execute kernel
        OCL_CHECK(error, error = q.enqueueTask(qrd_kernel));

        // Wait for kernel to finish
        OCL_CHECK(error, error = q.finish());

        // Write output buffer data back to A_tiled
        OCL_CHECK(error, error = q.enqueueMigrateMemObjects({output_matrix_1_ptr, output_matrix_2_ptr}, CL_MIGRATE_MEM_OBJECT_HOST));
        // OCL_CHECK(error, error = q.enqueueMigrateMemObjects({output_matrix_2_ptr}, CL_MIGRATE_MEM_OBJECT_HOST));

        // Wait for kernel to finish
        OCL_CHECK(error, error = q.finish());

        unflatten_matrix(A_tiled, flattened_matrix_read_1, idx_mat_1);
        unflatten_matrix(A_tiled, flattened_matrix_read_2, idx_mat_2);
    }
}

bool init_matrix(data_t matrix[TAM][TAM], std::fstream *file) {
    if (!file->is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return false;
    } else {
        std::cout << "Opened file" << std::endl;
    }

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
