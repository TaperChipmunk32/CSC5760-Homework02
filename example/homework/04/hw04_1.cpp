#include "mpi.h"
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //1. Demonstrate the use of MPI Bcast and MPI Reduce to achieve the same result as MPI Allreduce.
    int data_1 = rank + 1;
    int result_1;

    MPI_Bcast(&data_1, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Reduce(&data_1, &result_1, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Result: " << result_1 << std::endl;
    }

    int data_2 = rank + 1;
    int result_2;
    
    MPI_Allreduce(&data_2, &result_2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Result: " << result_2 << std::endl;
    }

    MPI_Finalize();
    return 0;
}