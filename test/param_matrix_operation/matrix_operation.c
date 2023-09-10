/* PMSIS includes */
#include "pmsis.h"

/* Variables used. */
#define MAT_SIZE (128)

struct cl_args_s
{
    uint32_t size;
    unsigned short *l1_in1;
    unsigned short *l1_in2;
    unsigned short *l1_out;
    unsigned short *l2_in1;
    unsigned short *l2_in2;
    unsigned short *l2_out;
};

PI_L2 static struct cl_args_s cl_arg;

void mat_display(unsigned short *A)
{
    int i, j, t = 0;

    printf("      ");
    for (j = 0; j < MAT_SIZE; j++)
        printf("%i ", j);
    printf("\n");
    printf("    __");
    for (j = 0; j < MAT_SIZE; j++)
        printf("____");
    printf("_\n");

    for (i = 0; i < MAT_SIZE; i++)
    {
        printf("%i | ", i);
        for (j = 0; j < MAT_SIZE; j++)
            printf("%2u ", A[t++]);
        printf("|\n");
    }
    printf("    --");
    for (j = 0; j < MAT_SIZE; j++)
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
    printf("Transfer Matrix done.\n");
    pi_cl_dma_wait(&copy);
}


void move_block_out(u_int32_t i, u_int16_t j, uint16_t size, pi_cl_dma_dir_e dir, unsigned short *ext, unsigned short *loc)
{
    for(uint16_t k=0; k<size; k++){
        pi_cl_dma_copy_t copy;
        copy.dir = dir;
        copy.merge = 0;
        copy.size = size * sizeof(unsigned short);
        copy.id = 0;
        copy.ext = (uint32_t)(ext + (i*size+k) * MAT_SIZE + j * size);
        copy.loc = (uint32_t)(loc);

        pi_cl_dma_memcpy(&copy);
    
        pi_cl_dma_wait(&copy);
    }

    printf("Transfer Matrix done.\n");
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
    uint32_t size_block = op_args->size;
    uint32_t coreid = pi_core_id(), start = 0, end = 0, start_glob = 0, end_glob = 0;
    uint16_t nb_block = 2;
    
    // iterate over all output block
    for (uint16_t i_block = 0; i_block < nb_block; i_block++){
        for (uint16_t j_block = 0; j_block < nb_block; j_block++){
            pi_cl_team_barrier(0);
            start_glob = (coreid * (size_block * size_block / pi_cl_cluster_nb_pe_cores()));
            end_glob = (start_glob + (size_block * size_block / pi_cl_cluster_nb_pe_cores()));
            for (uint32_t i = start_glob; i < end_glob; i++){
                    l1_out[i] = 0;
                }

            // compute one matrix product
            for (uint16_t k_block = 0; k_block < nb_block; k_block++){

                /* Core 0 of cluster initiates DMA transfer from L2 to L1. */
                if (!coreid){
                    printf("Core %d requesting DMA transfer from l2_in to l1_buffer.\n", coreid);
                    move_block(i_block, k_block, size_block, PI_CL_DMA_DIR_EXT2LOC, l2_in1, l1_in1);
                    move_block(k_block, j_block, size_block, PI_CL_DMA_DIR_EXT2LOC, l2_in2, l1_in2);
                }
                
                
                

                start = (coreid * (size_block / pi_cl_cluster_nb_pe_cores()));
                end = (start + (size_block / pi_cl_cluster_nb_pe_cores()));

                /* Barrier synchronisation before starting to compute. */
                pi_cl_team_barrier(0);

                /* Operation addition: Each core computes on specific portion of buffer. */
                

                

                /* Operation product: Each core computes on specific portion of buffer. */
                for (int k = 0; k < size_block; k++){
                    for (int i = start; i < end; i++){
                        for (int j = 0; j < size_block; j++){
                            l1_out[size_block * i + j] += l1_in1[size_block * i + k] * l1_in2[size_block * k + j];
                        }
                    }
                }

                /* Barrier synchronisation to wait for all cores. */
                pi_cl_team_barrier(0);
                }

                /* Core 0 of cluster initiates DMA transfer from L1 to L2. */
            if (!coreid){
                printf("Core %d requesting DMA transfer from l1_buffer to l2_out.\n", coreid);
                move_block_out(1, 1, size_block, PI_CL_DMA_DIR_LOC2EXT, l2_out, l1_out);
                printf("Core %d : Transfer output done.\n", coreid);
                
            }
        }
    }
}



/* Cluster main entry, executed by core 0. */
void master_entry(void *arg){
    printf("Cluster master core entry\n");
    /* Task dispatch to cluster cores. */
    pi_cl_team_fork(pi_cl_cluster_nb_pe_cores(), cluster_addition, arg);
    printf("Cluster master core exit\n");
}

void test_cluster_operation(void){
    printf("Entering main controller\n");
    struct pi_device cluster_dev;
    struct pi_cluster_conf conf;

    uint32_t nb_cl_pe_cores = pi_cl_cluster_nb_pe_cores();
    
    uint32_t block_size = 64;

    uint32_t m_size = (uint32_t)(MAT_SIZE * MAT_SIZE * sizeof(unsigned short));
    unsigned short *l2_in1 = pi_l2_malloc(m_size);
    unsigned short *l2_in2 = pi_l2_malloc(m_size);
    unsigned short *l2_out = pi_l2_malloc(m_size);

    if (l2_in1 == NULL || l2_in2 == NULL || l2_out == NULL){
        printf("Matrix alloc in l2 failed !\n");
        pmsis_exit(-1);
    }

    /* Matrix Init. */
    for (uint32_t i = 0; i < (MAT_SIZE * MAT_SIZE); i++){
        l2_in2[i] = i%64;
        l2_in1[i] = 1;
        l2_out[i] = 0;
    }

    /* Init cluster configuration structure. */
    pi_cluster_conf_init(&conf);
    conf.id = 0; /* Set cluster ID. */
    /* Configure & open cluster. */
    pi_open_from_conf(&cluster_dev, &conf);
    if (pi_cluster_open(&cluster_dev)){
        printf("Cluster open failed !\n");
        pmsis_exit(-3);
    }

    
    uint32_t bytes_block_size = (uint32_t) (block_size * block_size * sizeof(unsigned short));
    unsigned short *l1_in1 = pi_cl_l1_malloc(&cluster_dev, bytes_block_size);
    unsigned short *l1_in2 = pi_cl_l1_malloc(&cluster_dev, bytes_block_size);
    unsigned short *l1_out = pi_cl_l1_malloc(&cluster_dev, bytes_block_size);
    if (l1_in1 == NULL || l1_in2 == NULL || l1_out == NULL){
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
    /* Prepare cluster task and send it to cluster. */
    struct pi_cluster_task *task = pi_l2_malloc(sizeof(struct pi_cluster_task));
    if (task == NULL){
        printf("Cluster task alloc failed !\n");
        pi_cluster_close(&cluster_dev);
        pmsis_exit(-4);
    }

    pi_cluster_task(task, master_entry, &cl_arg);
    printf("Sending task.\n");
    pi_cluster_send_task_to_cl(&cluster_dev, task);

    // pi_l2_free(task, sizeof(struct pi_cluster_task));

    pi_cl_l1_free(&cluster_dev, l1_in1, bytes_block_size);
    pi_cl_l1_free(&cluster_dev, l1_in2, bytes_block_size);
    pi_cl_l1_free(&cluster_dev, l1_out, bytes_block_size);

    printf("Close cluster after end of computation.\n");
    pi_cluster_close(&cluster_dev);

    mat_display(l2_out);

    pi_l2_free(l2_out, m_size);
    pi_l2_free(l2_in1, m_size);
    pi_l2_free(l2_in2, m_size);

    printf("\nCluster operation done !\n");

    pmsis_exit(0);
}

/* Program Entry. */
int main(void){
    printf("\n\n\t *** PMSIS Cluster operation ***\n\n");
    return pmsis_kickoff((void *)test_cluster_operation);
}
