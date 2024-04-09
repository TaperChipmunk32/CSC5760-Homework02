#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Request request = MPI_REQUEST_NULL;

    int N = 5; // Number of times the message goes around the ring
    int message = 0;

    for (int i = 0; i < N; i++) {
        if (rank == 0) {
            // Process 0 sends the message to the next process
            MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            std::cout << "Process " << rank << " sent message: " << message << std::endl;
        } else {
            // Processes other than 0 receive the message from the previous process
            MPI_Recv(&message, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "Process " << rank << " received message: " << message << std::endl;

            // Process sends the message to the next process
            MPI_Send(&message, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD);
            std::cout << "Process " << rank << " sent message: " << message << std::endl;
        }

        // Process 0 receives the message from the last process
        if (rank == 0) {
            MPI_Recv(&message, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "Process " << rank << " received message: " << message << std::endl;
            message++;
        }
    }

    MPI_Finalize();
    return 0;
}
