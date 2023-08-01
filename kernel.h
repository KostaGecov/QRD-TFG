#define TAM_INDEX 9
#define FIXED_POINT 32
#define FX_POINT_INT 16

#define TAM_TILED 8
#define TAM 256
#define N_ITER 20
#define NUM_OPERACIONES 65  //(63 + 2(offset))

#define GEQRT 0
#define TTQRT 1

#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_stream.h>

#include <iostream>

/**
 * @brief 24 bits fixed point data, 10 for integer value and 14 for decimals.
 * The more bits it has, more shifts can be performed later, so the approximation to 0 will be more precise.
 * For bigger matrices, we need bigger data formats to be able to calculate the right result.
 *
 */
typedef ap_fixed<FIXED_POINT, FX_POINT_INT, AP_RND> data_t;

/**
 * @brief data type used for indexes variables in for loops
 *
 */
typedef ap_uint<TAM_INDEX> index_t;

/**
 * @brief Scale factor to compensate rotations
 *
 */
const data_t SCALE_FACTOR = 0.6072529;

class Rotator {
   public:
    // after data is read from an hls::stream<>, it cannot be read again
    // Stream is FIFO type
    hls::stream<data_t, TAM> row_x_in, row_y_in;
    hls::stream<data_t, TAM> row_x_out, row_y_out;

    hls::stream<data_t, TAM> q_u_in, q_v_in;
    hls::stream<data_t, TAM> q_u_out, q_v_out;

    int row_x, row_y, col;  // posiciones de las filas a rotar. En teor�a las
                            // columnas son las mismas en ambas filas

   public:
    /**
     * @brief Construct a new Rotator object
     *
     * @param x
     * @param y
     * @param c
     */
    Rotator(int x, int y, int c);

    /**
     * @brief Perform a givens rotation to two input streams and stores the results in other two streams
     *
     * @param row_x_in
     * @param row_y_in
     * @param row_x_out
     * @param row_y_out
     * @param col_rotator From where the rotation should start
     */
    void givens_rotation(hls::stream<data_t, TAM>& row_x_in,
                         hls::stream<data_t, TAM>& row_y_in,
                         hls::stream<data_t, TAM>& row_x_out,
                         hls::stream<data_t, TAM>& row_y_out,
                         hls::stream<data_t, TAM>& q_u_in,
                         hls::stream<data_t, TAM>& q_v_in,
                         hls::stream<data_t, TAM>& q_u_out,
                         hls::stream<data_t, TAM>& q_v_out,
                         int col_rotator);
};

/**
 * @brief Read input rows using blocking write to streams
 *
 * @param A
 * @param row_in_1
 * @param row_in_2
 * @param row_in_3
 * @param row_in_4
 * @param row_in_5
 * @param row_in_6
 * @param row_in_7
 * @param row_in_8
 */
void read_input_rows(data_t Matrix[TAM][TAM],
                     hls::stream<data_t, TAM>& row_in_1,
                     hls::stream<data_t, TAM>& row_in_2,
                     hls::stream<data_t, TAM>& row_in_3,
                     hls::stream<data_t, TAM>& row_in_4,
                     hls::stream<data_t, TAM>& row_in_5,
                     hls::stream<data_t, TAM>& row_in_6,
                     hls::stream<data_t, TAM>& row_in_7,
                     hls::stream<data_t, TAM>& row_in_8);

extern "C" {
/**
 * @brief Performs the givens rotation of the tiles. It is the top function
 *
 * @param A_tiled_1
 * @param A_tiled_2 In case 'type_op' is equal to GEQRT, this parameter is not used
 * @param Q_tiled_1
 * @param Q_tiled_2
 * @param type_op It can be GEQRT or TTQRT
 * @param col_offset Offset used to avoid reading the positions that had already become 0
 */
void krnl_givens_rotation(data_t A_tiled_1[TAM_TILED][TAM], data_t A_tiled_2[TAM_TILED][TAM],
                          data_t Q_tiled_1[TAM_TILED][TAM], data_t Q_tiled_2[TAM_TILED][TAM],
                          index_t type_op, index_t col_offset);
}
