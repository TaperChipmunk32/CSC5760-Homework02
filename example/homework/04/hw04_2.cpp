#include <iostream>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 2. Demonstrate the use of MPI Gather and MPI Bcast to achieve the same result as MPI Allgather.
    int* gatherData = new int[size];
    int* allGatherData = new int[size];
    int data = rank + 1;

    MPI_Gather(&data, 1, MPI_INT, gatherData, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(gatherData, size, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Allgather(&data, 1, MPI_INT, allGatherData, 1, MPI_INT, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Result using MPI Gather and MPI Bcast: ";
        for (int i = 0; i < size; i++) {
            std::cout << gatherData[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Result using MPI Allgather: ";
        for (int i = 0; i < size; i++) {
            std::cout << allGatherData[i] << " ";
        }
        std::cout << std::endl;
    }

    delete[] gatherData;
    delete[] allGatherData;

    MPI_Finalize();

    return 0;
}