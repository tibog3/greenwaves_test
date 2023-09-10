/* PMSIS includes */
#include "pmsis.h"

#include "utils.h"



PI_L2 static struct cl_args_s cl_arg;


/* Task executed by cluster cores. */
void cluster_operation(void *arg)
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
        printf("COmputing devided into %i blocks.\n", nb_block * nb_block);

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

                /* Operation product: Each core computes on specific portion of buffer. */
                product(size_block, start, end, l1_out, l1_in1, l1_in2);

                /* Barrier synchronisation to wait for all cores. */
                pi_cl_team_barrier(0);
            }

            // ------------------------------------------------------ //
            /*-----------compute convolution on the block-------------*/
            // ------------------------------------------------------ //

            convolution(size_block, start, end, l1_out, l1_in1, l1_in2);

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

    // fill borders bitween matrix blocks (k-1 colun) and (k-1) line
    start = (coreid * (MAT_SIZE / pi_cl_cluster_nb_pe_cores()));
    end = (start + (MAT_SIZE / pi_cl_cluster_nb_pe_cores()));
    for (uint16_t k = 1; k < nb_block; k++)
    {
        // get collumn to be corrected
        if (!coreid)
            move_columns(4, PI_CL_DMA_DIR_EXT2LOC, l2_out + k * size_block - 2, l1_out);

        // compute
        pi_cl_team_barrier(0);
        convolution_correct_column(size_block, start, end, l1_out, l1_in1, l1_in2);
        pi_cl_team_barrier(0);

        // send collumn
        if (!coreid)
            move_lines(2, PI_CL_DMA_DIR_LOC2EXT, l2_conv + k * size_block - 1, l1_in2);

        pi_cl_team_barrier(0);
    }

    for (uint16_t k = 1; k < nb_block; k++)
    {
        // get line
        if (!coreid)
            move_lines(4, PI_CL_DMA_DIR_EXT2LOC, l2_out + k * size_block * MAT_SIZE - 2 * MAT_SIZE, l1_out);

        // compute
        pi_cl_team_barrier(0);
        convolution_correct_line(size_block, start, end, l1_out, l1_in1, l1_in2);
        pi_cl_team_barrier(0);

        // send line
        if (!coreid)
            move_lines(2, PI_CL_DMA_DIR_LOC2EXT, l2_conv + k * size_block * MAT_SIZE - MAT_SIZE, l1_in2);
    }
}

/* Cluster main entry, executed by core 0. */
void master_entry(void *arg)
{
    printf("Cluster master core entry\n");
    /* Task dispatch to cluster cores. */
    pi_cl_team_fork(pi_cl_cluster_nb_pe_cores(), cluster_operation, arg);
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
