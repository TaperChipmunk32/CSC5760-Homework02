#include <Kokkos_Core.hpp>
#include "mpi.h"
#include <cstdio>

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
  
  // Get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // Create a View on the send side
  Kokkos::View<double*> send_view("send_view", size);
  
  // Initialize values on the send side
  if (rank == 0) {
    for (int i = 0; i < size; i++) {
      send_view(i) = i + 1;
    }
  }
  
  // Create a receive buffer on the recv side
  Kokkos::View<double*> recv_view("recv_view", size);
  
  // Send View values with MPI functions
  MPI_Allgather(&send_view(0), size, MPI_DOUBLE, &recv_view(0), size, MPI_DOUBLE, MPI_COMM_WORLD);
  
  // Output on recv side
  printf("Rank %d received the following values:\n", rank);
  for (int i = 0; i < size; i++) {
    printf("%f ", recv_view(i));
  }
  printf("\n");
  }
  Kokkos::finalize();
  MPI_Finalize();
}
