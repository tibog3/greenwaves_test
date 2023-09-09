/* PMSIS includes */
#include "pmsis.h"

/* Variables used. */
#define MAT_SIZE ( 80 )

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
  int i,j,t=0;

  printf("      ");
  for(j=0;j<MAT_SIZE;j++)
    printf("%i ",j);
  printf("\n");
  printf("    __");
  for(j=0;j<MAT_SIZE;j++)
    printf("____");
  printf("_\n");

  for(i=0;i<MAT_SIZE;i++)
    {
      printf("%i | ",i);
      for(j=0;j<MAT_SIZE;j++)
      printf("%2u ",A[t++]);
      printf("|\n");
    }
  printf("    --");
  for(j=0;j<MAT_SIZE;j++)
    printf("----");
  printf("-\n");
}

/* Task executed by cluster cores. */
void cluster_addition(void *arg)
{
    struct cl_args_s *op_args = (struct cl_args_s *) arg;
    unsigned short *l1_in1 = op_args->l1_in1;
    unsigned short *l1_in2 = op_args->l1_in2;
    unsigned short *l1_out = op_args->l1_out;
    unsigned short *l2_in1 = op_args->l2_in1;
    unsigned short *l2_in2 = op_args->l2_in2;
    unsigned short *l2_out = op_args->l2_out;
    uint32_t size = op_args->size;
    uint32_t coreid = pi_core_id(), start = 0, end = 0;

    /* Core 0 of cluster initiates DMA transfer from L2 to L1. */
    if (!coreid)
    {
        printf("Core %d requesting DMA transfer from l2_in to l1_buffer.\n", coreid);
        pi_cl_dma_copy_t copy;
        copy.dir = PI_CL_DMA_DIR_EXT2LOC;
        copy.merge = 0;
        copy.size = size;
        copy.id = 0;
        copy.ext = (uint32_t) l2_in1;
        copy.loc = (uint32_t) l1_in1;

        pi_cl_dma_memcpy(&copy);
        pi_cl_dma_wait(&copy);
        printf("Core %d : Transfer Matrix 1 done.\n", coreid);


        copy.dir = PI_CL_DMA_DIR_EXT2LOC;
        copy.merge = 0;
        copy.size = size;
        copy.id = 0;
        copy.ext = (uint32_t) l2_in2;
        copy.loc = (uint32_t) l1_in2;

        pi_cl_dma_memcpy(&copy);
        pi_cl_dma_wait(&copy);
        printf("Core %d : Transfer Matrix 2 done.\n", coreid);
    }

    start = (coreid * (MAT_SIZE * MAT_SIZE / pi_cl_cluster_nb_pe_cores()));
    end = (start  + (MAT_SIZE* MAT_SIZE / pi_cl_cluster_nb_pe_cores()));

    /* Barrier synchronisation before starting to compute. */
    pi_cl_team_barrier(0);


    /* Operation addition: Each core computes on specific portion of buffer. */
    for (uint32_t i=start; i<end; i++)
    {
        l1_out[i] = l1_in2[i];
        l1_out[i] = 0;
    }

    start = (coreid * (MAT_SIZE / pi_cl_cluster_nb_pe_cores()));
    end = (start  + (MAT_SIZE / pi_cl_cluster_nb_pe_cores()));

    /* Operation multiplication: Each core computes on specific portion of buffer. */
    for(int k = 0; k<MAT_SIZE; k++){
        for(int i = start; i<end; i++){
            for(int j = 0; j<MAT_SIZE; j++){
                l1_out[MAT_SIZE*i+j] += l1_in1[MAT_SIZE*i+k] * l1_in2[MAT_SIZE*k+j]; 
            }
        }
    }

    uint32_t core_id = pi_core_id(), cluster_id = pi_cluster_id();
    printf("[%d %d] Hello World!\n", start, end);

    /* Barrier synchronisation to wait for all cores. */
    pi_cl_team_barrier(0);

    /* Core 0 of cluster initiates DMA transfer from L1 to L2. */
    if (!coreid)
    {
        printf("Core %d requesting DMA transfer from l1_buffer to l2_out.\n", coreid);
        pi_cl_dma_copy_t copy;
        copy.dir = PI_CL_DMA_DIR_LOC2EXT;
        copy.merge = 0;
        copy.size = size;
        copy.id = 0;
        copy.ext = (uint32_t) l2_out;
        copy.loc = (uint32_t) l1_out;

        pi_cl_dma_memcpy(&copy);
        pi_cl_dma_wait(&copy);
        printf("Core %d : Transfer output done.\n", coreid);
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
 
    uint32_t m_size = (uint32_t) (MAT_SIZE * MAT_SIZE * sizeof(unsigned short));
    unsigned short *l2_in1 = pi_l2_malloc(m_size);
    unsigned short *l2_in2 = pi_l2_malloc(m_size);
    unsigned short *l2_out = pi_l2_malloc(m_size);

    if (l2_in1 == NULL || l2_in2 == NULL || l2_out == NULL)
    {
        printf("Matrix alloc in l2 failed !\n");
        pmsis_exit(-1);
    }

    /* Matrix Init. */
    for (uint32_t i=0; i<(MAT_SIZE*MAT_SIZE); i++)
    {
        l2_in1[i] = 1;
        l2_in2[i] = i%8;
        l2_out[i] = 1;
    }
    

    /* Init cluster configuration structure. */
    pi_cluster_conf_init(&conf);
    conf.id = 0;                /* Set cluster ID. */
    /* Configure & open cluster. */
    pi_open_from_conf(&cluster_dev, &conf);
    if (pi_cluster_open(&cluster_dev))
    {
        printf("Cluster open failed !\n");
        pmsis_exit(-3);
    }

    unsigned short *l1_in1 = pi_cl_l1_malloc(&cluster_dev, m_size);
    unsigned short *l1_in2 = pi_cl_l1_malloc(&cluster_dev, m_size);
    unsigned short *l1_out = pi_cl_l1_malloc(&cluster_dev, m_size);
    if (l1_in1 == NULL || l1_in2 == NULL || l1_out == NULL)
    {
        printf("l1 alloc failed !\n");
        pi_cluster_close(&cluster_dev);
        pmsis_exit(-2);
    }

    /* Init arg struct. */
    cl_arg.size = m_size;
    cl_arg.l1_in1 = l1_in1;
    cl_arg.l1_in2 = l1_in2;
    cl_arg.l1_out = l1_out;
    cl_arg.l2_in1 = l2_in1;
    cl_arg.l2_in2 = l2_in2;
    cl_arg.l2_out = l2_out;
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


    //pi_l2_free(task, sizeof(struct pi_cluster_task));

    pi_cl_l1_free(&cluster_dev, l1_in1, m_size);
    pi_cl_l1_free(&cluster_dev, l1_in2, m_size);
    pi_cl_l1_free(&cluster_dev, l1_out, m_size);

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
int main(void)
{
    printf("\n\n\t *** PMSIS Cluster operation ***\n\n");
    return pmsis_kickoff((void *) test_cluster_operation);
}
