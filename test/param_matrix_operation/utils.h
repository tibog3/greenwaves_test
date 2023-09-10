#ifndef UTILS_H
#define UTILS_H

/* Variables used. */
#ifndef MAT_SIZE
#define MAT_SIZE (64)
#endif

/* PMSIS includes */
#include "pmsis.h"

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




void mat_display(unsigned short *A, int16_t size);

/***
 * i, j : bloc coord in matrix . If matrix is devided into 9 bloc i is 0 1 or 2 and j 0 1 or 2
 * size : size of the block
 * dir : direction of tranfert
 *
 */
void move_block(u_int32_t i, u_int16_t j, uint16_t size, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc);

void move_lines(u_int16_t nb_of_lines, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc);

void move_columns(u_int16_t nb_of_cols, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc);

/*matrix product*/
void product(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2);

/* convolution of block padd with 0*/
void convolution(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2);

void convolution_correct_line(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2);

void convolution_correct_column(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2);

#endif