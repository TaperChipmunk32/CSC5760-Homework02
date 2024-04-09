#include "mpi.h"
#include <iostream>

int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <Q>" << std::endl;
        return 1;
    }
    int Q = atoi(argv[1]);
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_size % Q != 0) {
        std::cerr << "World size must be divisible by Q" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // First split based on ranks divided by Q
    int color1 = world_rank / Q;
    MPI_Comm comm1;
    MPI_Comm_split(MPI_COMM_WORLD, color1, world_rank, &comm1);

    // Second split based on ranks mod Q
    int color2 = world_rank % Q;
    MPI_Comm comm2;
    MPI_Comm_split(MPI_COMM_WORLD, color2, world_rank, &comm2);

    // Get the size and rank within the new communicators
    int comm1_size, comm1_rank;
    MPI_Comm_size(comm1, &comm1_size);
    MPI_Comm_rank(comm1, &comm1_rank);

    int comm2_size, comm2_rank;
    MPI_Comm_size(comm2, &comm2_size);
    MPI_Comm_rank(comm2, &comm2_rank);

    // Print the results
    std::cout << "World Rank/Size: " << world_rank << "/" << world_size << ", ";
    std::cout << "Comm1 Rank/Size: " << comm1_rank << "/" << comm1_size << ", ";
    std::cout << "Comm2 Rank/Size: " << comm2_rank << "/" << comm2_size << std::endl;

    MPI_Finalize();
    return 0;
}