#ifndef UTILS_H
#define UTILS_H

/* Variables used. */
#ifndef MAT_SIZE
#define MAT_SIZE (64)
#endif

/* PMSIS includes */
#include "pmsis.h"

/**
 * Args struct
*/
struct cl_args_s
{
    uint32_t size;
    unsigned short *l1_in1;
    unsigned short *l1_in2;
    unsigned short *l1_out;
    unsigned short *l2_in1;
    unsigned short *l2_in2;
    unsigned short *l2_out;
    unsigned short *l2_conv;
    uint16_t nb_block;
};



/*Display matrix
* A : pointer to matrix
* size : size of matrix*/
void mat_display(unsigned short *A, int16_t size);

/**
 * Move a square block of a matrix from L2 to L1 and from L1 to L2
 * i, j : bloc coord in matrix . If matrix is devided into 9 bloc i is 0 1 or 2 and j 0 1 or 2
 * size : size of the block
 * dir : direction of tranfert
 *
 */
void move_block(u_int32_t i, u_int16_t j, uint16_t size, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc);

/**
 * Move a set of lines of a matrix to L1 and back to L2
 * nb_of_lines : Number of lines to tranfert
 * dir : direction
 * ext : adressof beginning of first line
 * loc : adress of local matrix
*/
void move_lines(u_int16_t nb_of_lines, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc);

/**
 * Move a set of cols of a matrix to L1 and back to L2
 * nb_of_lines : Number of cols to tranfert
 * dir : direction
 * ext : adressof beginning of first col
 * loc : adress of local matrix
*/
void move_columns(u_int16_t nb_of_cols, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc);

/*matrix product*/
void product(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2);

/* convolution of block padd with 0*/
void convolution(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2);

/**
 * Convolution of a 4*MAT_SIZE matirix, no padding
 * Function used to correct transition between blocks
*/
void convolution_correct_line(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2);

/**
 * Convolution of a MAT_SIZE*4 matirix, no padding
 * Function used to correct transition between blocks
*/
void convolution_correct_column(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2);

#endif