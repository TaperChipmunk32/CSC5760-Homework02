#include <iostream>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 3. Demonstrate the use of MPI Alltoall to send a personalized communication between each process in MPI COMM WORLD.
    int* sendBuffer = new int[size];
    int* recvBuffer = new int[size];

    for (int i = 0; i < size; i++) {
        sendBuffer[i] = rank + i;
    }

    MPI_Alltoall(sendBuffer, 1, MPI_INT, recvBuffer, 1, MPI_INT, MPI_COMM_WORLD);

    std::cout << "MPI Alltoall result for process " << rank << ": ";
    for (int i = 0; i < size; i++) {
        std::cout << recvBuffer[i] << " ";
    }
    std::cout << std::endl;

    delete[] sendBuffer;
    delete[] recvBuffer;

    MPI_Finalize();

    return 0;
}