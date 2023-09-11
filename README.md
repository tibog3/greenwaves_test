# greenwaves_test
Test for Greenwaves Technologies

The work is located on the `test` branch in te test directory.


## Folder Descriptions:

- The `helloworld` directory allows printing "Hello, World!" using all available cores, including the main core and the 8 cluster cores. Use `./launch.sh` to compile and run the code on gvsoc. (question 1 and 2)

- The `matrix_operation` directory implements matrix addition and multiplication for matrices of size 64x64. It also applies convolution with the 3x3 mask [-1, -2, -1, 0, 0, 1, 2, 1] to the result of multiplication. This folder is designed for handling simple mountain manipulations and getting a grasp of workload distribution among the cores. (question 3&5)

- The `param_matrix_operation` directory contains similar code to the one above but adds the option to use a matrix size of 128. On the GAP 8 platform, the available L1 memory is insufficient for processing matrices of this size. Therefore, the workload needs to be divided, resulting in multiple memory accesses, making convolution and multiplication more complex. **I focused on the task of memory distribution by dividing the work into blocks of 64. My code can handle only matrices with sizes that are multiples of 64**, such as 64, 128, 192 etc. (Question 4 + adapt of 3&5)
Some improvements would be needed to handle all other matrix sizes efficiently.

- - While I didn't address the last question specifically, I took care of memory usage when building my algorithm, and all calculations are parallelized across the cluster cores. I tried to minimize data exchanges between L2 and L1 memory as much as possible. Additionally, I aimed for a more general approach, which is why I decided to process data in blocks.

  To handle other matrix sizes (not multiple of 64), the code would need to be adapted to work with variable-sized rectangular blocks. It would also require careful consideration of task division if the workload cannot be divided by 8 cores. This issue appears to be less complex to resolve than memory management. Therefore, I chose to focus on the essentials.
