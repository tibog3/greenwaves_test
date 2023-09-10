/* PMSIS includes */
#include "pmsis.h"

#include "utils.h"


void mat_display(unsigned short *A, int16_t size)
{
    int i, j, t = 0;

    printf("      ");
    for (j = 0; j < size; j++)
        printf("%i ", j);
    printf("\n");
    printf("    __");
    for (j = 0; j < size; j++)
        printf("____");
    printf("_\n");

    for (i = 0; i < size; i++)
    {
        printf("%i | ", i);
        for (j = 0; j < size; j++)
            printf("%2u ", A[t++]);
        printf("|\n");
    }
    printf("    --");
    for (j = 0; j < size; j++)
        printf("----");
    printf("-\n");
}

/***
 * i, j : bloc coord in matrix . If matrix is devided into 9 bloc i is 0 1 or 2 and j 0 1 or 2
 * size : size of the block
 * dir : direction of tranfert
 *
 */
void move_block(u_int32_t i, u_int16_t j, uint16_t size, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc)
{
    pi_cl_dma_copy_2d_t copy;
    copy.dir = dir;
    copy.merge = 0;
    copy.size = size * size * sizeof(unsigned short);
    copy.stride = MAT_SIZE * sizeof(unsigned short);
    copy.length = size * sizeof(unsigned short);
    copy.id = 0;
    copy.ext = (uint32_t)(ext + i * MAT_SIZE * size + j * size);
    copy.loc = (uint32_t)loc;

    pi_cl_dma_memcpy_2d(&copy);
    // printf("Transfer Matrix done.\n");
    pi_cl_dma_wait(&copy);
}

void move_lines(u_int16_t nb_of_lines, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc)
{
    pi_cl_dma_copy_t copy;
    copy.dir = dir;
    copy.merge = 0;
    copy.size = nb_of_lines * MAT_SIZE * sizeof(unsigned short);
    copy.id = 0;
    copy.ext = (uint32_t)(ext);
    copy.loc = (uint32_t)loc;

    pi_cl_dma_memcpy(&copy);
    printf("Send corrected line.\n");
    pi_cl_dma_wait(&copy);
}

void move_columns(u_int16_t nb_of_cols, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc)
{
    pi_cl_dma_copy_2d_t copy;
    copy.dir = PI_CL_DMA_DIR_EXT2LOC;
    copy.merge = 0;
    copy.size = MAT_SIZE * nb_of_cols * sizeof(unsigned short);
    copy.stride = MAT_SIZE * sizeof(unsigned short);
    copy.length = nb_of_cols * sizeof(unsigned short);
    copy.id = 0;
    copy.ext = (uint32_t)(ext);
    copy.loc = (uint32_t)loc;

    pi_cl_dma_memcpy_2d(&copy);
    printf("Get column to be corrected.\n");
    pi_cl_dma_wait(&copy);
}

/*matrix product*/
void product(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2)
{
    for (int k = 0; k < size_block; k++)
    {
        for (int i = start; i < end; i++)
        {
            for (int j = 0; j < size_block; j++)
            {
                l1_out[size_block * i + j] += l1_in1[size_block * i + k] * l1_in2[size_block * k + j];
            }
        }
    }
}

/* convolution of block padd with 0*/
void convolution(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2)
{
    /* Convolution step 1 : apply on each line*/
    for (uint16_t i = start; i < end; i++)
    {
        for (uint16_t j = 0; j < size_block; j++)
        {
            if (j > 0 && j < size_block - 1)
            {
                // normal
                l1_in1[size_block * i + j] = l1_out[size_block * i + j - 1] + 2 * l1_out[size_block * i + j] + l1_out[size_block * i + j + 1];
            }
            else if (j == 0)
            {
                // padd with 0 left
                l1_in1[size_block * i + j] = 2 * l1_out[size_block * i + j] + l1_out[size_block * i + j + 1];
            }
            else
            {
                // padd with 0 right
                l1_in1[size_block * i + j] = l1_out[size_block * i + j - 1] + 2 * l1_out[size_block * i + j];
            }
        }
    }
    /* Barrier synchronisation to wait for all cores. */
    pi_cl_team_barrier(0);
    /* Convolution step 2 : apply on each column*/
    for (uint16_t i = start; i < end; i++)
    {
        for (uint16_t j = 0; j < size_block; j++)
        {
            if (i > 0 && i < size_block - 1)
            {
                // normal
                l1_in2[size_block * i + j] = l1_in1[size_block * i + j + size_block] - l1_in1[size_block * i + j - size_block];
            }
            else if (i == 0)
            {
                // padd with 0 up
                l1_in2[size_block * i + j] = l1_in1[size_block * i + j + size_block];
            }
            else
            {
                // padd with 0 down
                l1_in2[size_block * i + j] = (-l1_in1[size_block * i + j - size_block]);
            }
        }
    }
}

void convolution_correct_line(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2)
{
    // convolution with no padding for a 4 x SIZE_MAT matrix
    for (uint16_t i = 0; i < 4; i++)
    {
        for (uint16_t j = start; j < end; j++)
        {
            if (j > 0 && j < MAT_SIZE - 1)
            {
                // normal
                l1_in1[i * MAT_SIZE + j] = l1_out[i * MAT_SIZE + j - 1] + 2 * l1_out[i * MAT_SIZE + j] + l1_out[i * MAT_SIZE + j];
            }
            else if (j == 0)
            {
                l1_in1[i * MAT_SIZE + j] = 2 * l1_out[i * MAT_SIZE + j] + l1_out[i * MAT_SIZE + j];
            }
            else
            {
                l1_in1[i * MAT_SIZE + j] = l1_out[i * MAT_SIZE + j - 1] + 2 * l1_out[i * MAT_SIZE + j];
            }
        }
    }
    pi_cl_team_barrier(0);

    for (uint16_t i = 0; i < 2; i++)
    {
        for (uint16_t j = start; j < end; j++)
        {
            l1_in2[MAT_SIZE * i + j] = l1_in1[MAT_SIZE * (i + 2) + j] - l1_in1[MAT_SIZE * i + j];
        }
    }
}

void convolution_correct_column(uint32_t size_block, uint32_t start, uint32_t end, unsigned short *l1_out, unsigned short *l1_in1, unsigned short *l1_in2)
{
    // convolution with no padding for a SIZE_MAT x 4 matrix
    for (uint16_t i = start; i < end; i++)
    {
        for (uint16_t j = 0; j < 2; j++)
        {
            l1_in1[i * 2 + j] = l1_out[i * 4 + j] + 2 * l1_out[i * 4 + j + 1] + l1_out[i * 4 + j + 2];
        }
    }
    pi_cl_team_barrier(0);

    for (uint16_t i = start; i < end; i++)
    {
        for (uint16_t j = 0; j < 2; j++)
        {
            if (i > 0 && i < MAT_SIZE - 1)
            {
                // normal
                l1_in2[2 * i + j] = l1_in1[2 * (i + 1) + j] - l1_in1[2 * (i - 1) + j];
            }
            else if (i == 0)
            {
                l1_in2[2 * i + j] = l1_in1[2 * i + j + 2];
            }
            else
            {
                l1_in2[2 * i + j] = (-l1_in1[2 * i + j - 2]);
            }
        }
    }
}
