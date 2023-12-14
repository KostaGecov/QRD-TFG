#include <fstream>
#include <iostream>

#include "xcl2.hpp"

void init_matrix(data_t matrix[TAM][TAM], std::fstream *file) {
    if (!file->is_open()) {
        std::cerr << "Error opening file" << std::endl;
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
}

void init_matrix(float matrix[TAM][TAM], std::fstream *file) {
    if (!file->is_open()) {
        std::cerr << "Error opening file" << std::endl;
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
}

float error(data_t A[TAM][TAM], float out_gold[TAM][TAM]) {
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

int main(int argc, char **argv) {
    // Implementation remains the same until the OpenCL setup

    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }

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
     * stores all 32 Q matrices needed for tiled operations
     *
     */
    data_t Q_tiled[NUM_TILED][TAM_TILED][TAM];
    data_t Q[TAM][TAM];

    /**
     * to store the output data gold
     *
     */
    float out_gold[TAM][TAM];

    std::string binaryFile = argv[1];
    cl_int err;
    cl::Context context;
    cl::CommandQueue q;
    cl::Kernel qrd;
    // std::vector<int, aligned_allocator<int> > host_memory(elements, 42);
    // std::vector<int, aligned_allocator<int> > host_memory2(elements, 15);

    // size_t size_in_bytes = host_memory.size() * sizeof(int);

    // OPENCL HOST CODE AREA START
    auto start1 = chrono::high_resolution_clock::now();
    // get_xil_devices() is a utility API which will find the xilinx
    // platforms and will return list of devices connected to Xilinx platform
    auto devices = xcl::get_xil_devices();
    // read_binary_file() is a utility API which will load the binaryFile
    // and will return the pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    bool valid_device = false;
    for (unsigned int i = 0; i < devices.size(); i++) {
        // for (auto device : devices) {
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
            OCL_CHECK(err, qrd = cl::Kernel(program, "Kernel", &err));
            valid_device = true;
            break;  // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    auto end1 = chrono::high_resolution_clock::now();
    auto diff1 = end1 - start1;
    std::cout << "FPGA programming time: " << chrono::duration<double, milli>(diff1).count() << std::endl;

    std::fstream data_in("data_in.dat", std::ios::in);
    init_matrix(A, &data_in);

    std::fstream data_out_gold("data_out_gold.dat", std::ios::in);
    init_matrix(out_gold, &data_out_gold);

initialize_Q:
    for (uint16_t r = 0; r < TAM; r++) {
        for (uint16_t c = 0; c < TAM; c++) {
            if (r == c) {
                Q[r][c] = 1;
            } else {
                Q[r][c] = 0;
            }
        }
    }

divide_matrices_row_for:
    for (uint16_t r = 0; r < TAM; r++) {
        int tile_index = r / TAM_TILED;
        int tile_offset = r % TAM_TILED;

    divide_matrices_col_for:
        for (uint16_t c = 0; c < TAM; c++) {
            A_tiled[tile_index][tile_offset][c] = A[r][c];
            Q_tiled[tile_index][tile_offset][c] = Q[r][c];
        }
    }

    /**
     * @todo: continuar
     *
     */

    return 0;
}
