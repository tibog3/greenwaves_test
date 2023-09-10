# greenwaves_test
Test for Greenwaves Technologies

The work is located on the `test` branch in te test directory.


## Folder Descriptions:

- The `helloworld` directory allows printing "Hello, World!" using all available cores, including the main core and the 8 cluster cores. Use `./launch.sh` to compile and run the code on gvsoc. (question 1 and 2)

- The `matrix_operation` directory implements matrix addition and multiplication for matrices of size 64x64. It also applies convolution with the mask [-1, -2, -1, 0, 0, 1, 2, 1] to the result of multiplication. This folder is designed for handling simple mountain manipulations and getting a grasp of workload distribution among the cores. (question 3-4-6)

- The `param_matrix_operation` directory contains similar code to the one above but adds the option to use a matrix size of 128. On the GAP 8 platform, the available L1 memory is insufficient for processing matrices of this size. Therefore, the workload needs to be divided, resulting in multiple memory accesses, making convolution and multiplication more complex. I focused on the task of memory distribution by dividing the work into blocks of 64. My code can handle matrices with sizes that are multiples of 64, such as 128, 192, 256, etc. Some improvements would be needed to handle all other matrix sizes efficiently. (Question 5 + addaptation of 3-4-6)