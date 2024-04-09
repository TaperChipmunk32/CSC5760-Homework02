#include <Kokkos_Core.hpp>
#include <cstdio>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
  // Make View of length n > 1000 and add values
  int n = 1001;
  Kokkos::View<int*> view("view", n);
  Kokkos::parallel_for("fill_view", n, KOKKOS_LAMBDA(const int i) {
    view(i) = i * i;
  });

  // create two additional views of same size and datatype
  Kokkos::View<int*> view2("view2", n);
  Kokkos::View<int*> view3("view3", n);

  // deep_copy
  Kokkos::Timer deepCopyTimer;
  Kokkos::deep_copy(view2, view);
  double deepCopyTime = deepCopyTimer.seconds();

  // user copy
  Kokkos::Timer userCopyTimer;
  for (int i = 0; i < n; i++) {
    view3(i) = view(i);
  }
  double userCopyTime = userCopyTimer.seconds();

  // Output times
  printf("Deep Copy Time: %f seconds\n", deepCopyTime);
  printf("User Copy Time: %f seconds\n", userCopyTime);

  }
  Kokkos::finalize();
}
