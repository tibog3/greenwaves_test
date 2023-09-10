/* PMSIS includes */
#include "pmsis.h"

/* Variables used. */
#ifndef MAT_SIZE
#define MAT_SIZE (64)
#endif

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

PI_L2 static struct cl_args_s cl_arg;

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
    //printf("Transfer Matrix done.\n");
    pi_cl_dma_wait(&copy);
}

/* Task executed by cluster cores. */
void cluster_addition(void *arg)
{
    struct cl_args_s *op_args = (struct cl_args_s *)arg;
    unsigned short *l1_in1 = op_args->l1_in1;
    unsigned short *l1_in2 = op_args->l1_in2;
    unsigned short *l1_out = op_args->l1_out;
    unsigned short *l2_in1 = op_args->l2_in1;
    unsigned short *l2_in2 = op_args->l2_in2;
    unsigned short *l2_out = op_args->l2_out;
    unsigned short *l2_conv = op_args->l2_conv;
    uint32_t size_block = op_args->size;
    uint32_t coreid = pi_core_id(), start = 0, end = 0, start_glob = 0, end_glob = 0;
    uint16_t nb_block = op_args->nb_block;
    if (!coreid)
    printf("COmputing devided into %i blocks.\n", nb_block*nb_block);

    start = (coreid * (size_block / pi_cl_cluster_nb_pe_cores()));
    end = (start + (size_block / pi_cl_cluster_nb_pe_cores()));
    start_glob = (coreid * (size_block * size_block / pi_cl_cluster_nb_pe_cores()));
    end_glob = (start_glob + (size_block * size_block / pi_cl_cluster_nb_pe_cores()));

    // iterate over all output block
    for (uint16_t i_block = 0; i_block < nb_block; i_block++)
    {
        for (uint16_t j_block = 0; j_block < nb_block; j_block++)
        {
            pi_cl_team_barrier(0);
            // initialise block to 0
            for (uint32_t i = start_glob; i < end_glob; i++)
            {
                l1_out[i] = 0;
            }

            // ------------------------------------------------------ //
            /*-----------compute one matrix product-------------------*/
            // ------------------------------------------------------ //
            for (uint16_t k_block = 0; k_block < nb_block; k_block++)
            {

                /* Core 0 of cluster initiates DMA transfer from L2 to L1. */
                if (!coreid)
                {
                    printf("Get new blocks : requesting DMA transfer from l2_in to l1_buffer.\n");
                    move_block(i_block, k_block, size_block, PI_CL_DMA_DIR_EXT2LOC, l2_in1, l1_in1);
                    move_block(k_block, j_block, size_block, PI_CL_DMA_DIR_EXT2LOC, l2_in2, l1_in2);
                }

                /* Barrier synchronisation before starting to compute. */
                pi_cl_team_barrier(0);

                /* Operation addition: Each core computes on specific portion of buffer. */

                /* Operation product: Each core computes on specific portion of buffer. */
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

                /* Barrier synchronisation to wait for all cores. */
                pi_cl_team_barrier(0);
            }

            // ------------------------------------------------------ //
            /*-----------compute convolution on the block-------------*/
            // ------------------------------------------------------ //

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
            pi_cl_team_barrier(0);

            /*
            Core 0 of cluster initiates DMA transfer from L1 to L2.
            Put block result in L2 memory
            */
            if (!coreid)
            {
                printf("Send block result : requesting DMA transfer from l1_buffer to l2_out.\n");
                move_block(i_block, j_block, size_block, PI_CL_DMA_DIR_LOC2EXT, l2_out, l1_out);
                move_block(i_block, j_block, size_block, PI_CL_DMA_DIR_LOC2EXT, l2_conv, l1_in2);
                printf("Transfer output done.\n");
            }
        }
    }


    // ----------------------------------------------------------------//
    /*-----Complete convolution - recompute borders between blocks-----*/
    // --------------------------------------------------------------- //


    // fill borders inside blocks (k-1 colun) and (k-1) line
    start = (coreid * (MAT_SIZE / pi_cl_cluster_nb_pe_cores()));
    end = (start + (MAT_SIZE / pi_cl_cluster_nb_pe_cores()));
    for (uint16_t k = 1; k < nb_block; k++)
    {
        // get collumn
        if (!coreid)
        {
            pi_cl_dma_copy_2d_t copy;
            copy.dir = PI_CL_DMA_DIR_EXT2LOC;
            copy.merge = 0;
            copy.size = MAT_SIZE * 4 * sizeof(unsigned short);
            copy.stride = MAT_SIZE * sizeof(unsigned short);
            copy.length = 4 * sizeof(unsigned short);
            copy.id = 0;
            copy.ext = (uint32_t)(l2_out + k * size_block - 2);
            copy.loc = (uint32_t)l1_out;

            pi_cl_dma_memcpy_2d(&copy);
            printf("Get column to be corrected.\n");
            pi_cl_dma_wait(&copy);
        }
        // compute
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
        pi_cl_team_barrier(0);

        // send collumn
        if (!coreid)
        {
            pi_cl_dma_copy_2d_t copy;
            copy.dir = PI_CL_DMA_DIR_LOC2EXT;
            copy.merge = 0;
            copy.size = MAT_SIZE * 2 * sizeof(unsigned short);
            copy.stride = MAT_SIZE * sizeof(unsigned short);
            copy.length = 2 * sizeof(unsigned short);
            copy.id = 0;
            copy.ext = (uint32_t)(l2_conv + k * size_block - 1);
            copy.loc = (uint32_t)l1_in2;

            pi_cl_dma_memcpy_2d(&copy);
            printf("Send corrected column\n");
            pi_cl_dma_wait(&copy);
        }
        pi_cl_team_barrier(0);
    }

    for (uint16_t k = 1; k < nb_block; k++)
    {
        // get line

        if (!coreid)
        {
            pi_cl_dma_copy_t copy;
            copy.dir = PI_CL_DMA_DIR_EXT2LOC;
            copy.merge = 0;
            copy.size = 4 * MAT_SIZE * sizeof(unsigned short);
            copy.id = 0;
            copy.ext = (uint32_t)(l2_out + k * size_block * MAT_SIZE - 2 * MAT_SIZE);
            copy.loc = (uint32_t)l1_out;

            pi_cl_dma_memcpy(&copy);
            printf("Get line.\n");
            pi_cl_dma_wait(&copy);
        }

        // compute

        pi_cl_team_barrier(0);
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

        pi_cl_team_barrier(0);

        // send line
        if (!coreid)
        {
            pi_cl_dma_copy_t copy;
            copy.dir = PI_CL_DMA_DIR_LOC2EXT;
            copy.merge = 0;
            copy.size = 2 * MAT_SIZE * sizeof(unsigned short);
            copy.id = 0;
            copy.ext = (uint32_t)(l2_conv + k * size_block * MAT_SIZE - MAT_SIZE);
            copy.loc = (uint32_t)l1_in2;
            // mat_display(l1_in2);
            // mat_display(l1_in2, 16);

            pi_cl_dma_memcpy(&copy);
            printf("Send corrected line.\n");
            pi_cl_dma_wait(&copy);
        }
    }
}

/* Cluster main entry, executed by core 0. */
void master_entry(void *arg)
{
    printf("Cluster master core entry\n");
    /* Task dispatch to cluster cores. */
    pi_cl_team_fork(pi_cl_cluster_nb_pe_cores(), cluster_addition, arg);
    printf("Cluster master core exit\n");
}

void test_cluster_operation(void)
{
    printf("Entering main controller\n");
    struct pi_device cluster_dev;
    struct pi_cluster_conf conf;

    uint32_t nb_cl_pe_cores = pi_cl_cluster_nb_pe_cores();

    uint32_t block_size = 64;

    uint32_t m_size = (uint32_t)(MAT_SIZE * MAT_SIZE * sizeof(unsigned short));
    unsigned short *l2_in1 = pi_l2_malloc(m_size);
    unsigned short *l2_in2 = pi_l2_malloc(m_size);
    unsigned short *l2_out = pi_l2_malloc(m_size);
    unsigned short *l2_conv = pi_l2_malloc(m_size);

    if (l2_in1 == NULL || l2_in2 == NULL || l2_out == NULL)
    {
        printf("Matrix alloc in l2 failed !\n");
        pmsis_exit(-1);
    }

    /* Matrix Init. */
    for (uint32_t i = 0; i < (MAT_SIZE * MAT_SIZE); i++)
    {
        l2_in2[i] = 1;
        l2_in1[i] = 1;
        l2_out[i] = 0;
    }

    /* Init cluster configuration structure. */
    pi_cluster_conf_init(&conf);
    conf.id = 0; /* Set cluster ID. */
    /* Configure & open cluster. */
    pi_open_from_conf(&cluster_dev, &conf);
    if (pi_cluster_open(&cluster_dev))
    {
        printf("Cluster open failed !\n");
        pmsis_exit(-3);
    }

    uint32_t bytes_block_size = (uint32_t)(block_size * block_size * sizeof(unsigned short));
    unsigned short *l1_in1 = pi_cl_l1_malloc(&cluster_dev, bytes_block_size);
    unsigned short *l1_in2 = pi_cl_l1_malloc(&cluster_dev, bytes_block_size);
    unsigned short *l1_out = pi_cl_l1_malloc(&cluster_dev, bytes_block_size);
    if (l1_in1 == NULL || l1_in2 == NULL || l1_out == NULL)
    {
        printf("l1 alloc failed !\n");
        pi_cluster_close(&cluster_dev);
        pmsis_exit(-2);
    }

    /* Init arg struct. */
    cl_arg.size = block_size;
    cl_arg.l1_in1 = l1_in1;
    cl_arg.l1_in2 = l1_in2;
    cl_arg.l1_out = l1_out;
    cl_arg.l2_in1 = l2_in1;
    cl_arg.l2_in2 = l2_in2;
    cl_arg.l2_out = l2_out;
    cl_arg.l2_conv = l2_conv;
    cl_arg.nb_block = MAT_SIZE / 64;
    /* Prepare cluster task and send it to cluster. */
    struct pi_cluster_task *task = pi_l2_malloc(sizeof(struct pi_cluster_task));
    if (task == NULL)
    {
        printf("Cluster task alloc failed !\n");
        pi_cluster_close(&cluster_dev);
        pmsis_exit(-4);
    }

    pi_cluster_task(task, master_entry, &cl_arg);
    printf("Sending task.\n");
    pi_cluster_send_task_to_cl(&cluster_dev, task);

    // Free memory
    pi_l2_free(task, sizeof(struct pi_cluster_task));
    pi_cl_l1_free(&cluster_dev, l1_in1, bytes_block_size);
    pi_cl_l1_free(&cluster_dev, l1_in2, bytes_block_size);
    pi_cl_l1_free(&cluster_dev, l1_out, bytes_block_size);

    printf("Close cluster after end of computation.\n");
    pi_cluster_close(&cluster_dev);

    // mat_display(l2_out);
    mat_display(l2_conv, MAT_SIZE);

    // program end, free memory
    pi_l2_free(l2_out, m_size);
    pi_l2_free(l2_in1, m_size);
    pi_l2_free(l2_in2, m_size);

    printf("\nCluster operation done !\n");

    pmsis_exit(0);
}

/* Program Entry. */
int main(void)
{
    printf("\n\n\t *** PMSIS Cluster operation ***\n\n");
    return pmsis_kickoff((void *)test_cluster_operation);
}
