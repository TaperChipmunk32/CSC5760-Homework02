#include "mpi.h"
#include <iostream>

using namespace std;

#define GLOBAL_PRINT

class Domain
{
public:
  Domain(int _M, int _N, const char *_name="") : domain(new char[(_M+1)*(_N+1)]), M(_M), N(_N), name(_name)  {}
  virtual ~Domain() {delete[] domain;}
  char &operator()(int i, int j) {return domain[i*N+j];}
  char operator()(int i, int j)const {return domain[i*N+j];}

  int rows() const {return M;}
  int cols() const {return N;}

  const string & myname() const {return name;}

  char *rawptr() {return domain;}
  
protected:
  char *domain;
  int M;
  int N;

  string name;
};

void zero_domain(Domain &domain);
void print_domain(Domain &domain, int rank);
void update_domain(Domain &new_domain, Domain &old_domain, int world_rank, int world_size);
void game_of_life(int world_rank, int world_size, int N, int P, int Q, int iterations);

int main(int argc, char** argv) {

    int N, P, Q;
    int iterations;

    if(argc < 4)
    {
        cout << "usage: " << argv[0] << " N Q iterations" << endl;
        exit(0);
    }

    N = atoi(argv[1]); Q = atoi(argv[2]); iterations = atoi(argv[3]);

    int world_size, world_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int array[4];
    if(world_rank == 0)
    {
        N = atoi(argv[1]); Q = atoi(argv[2]); iterations = atoi(argv[3]);
        P = world_size / Q;
        array[0] = N;
        array[1] = P;
        array[2] = Q;
        array[3] = iterations;
        
    }
    MPI_Bcast(array, 4, MPI_INT, 0, MPI_COMM_WORLD);
    if(world_rank != 0)
    {
        N = array[0];
        P = array[1];
        Q = array[2];
        iterations = array[3];
    }

    game_of_life(world_rank, world_size, N, P, Q, iterations);

    MPI_Finalize();
    return 0;
}

void game_of_life(int world_rank, int world_size, int N, int P, int Q, int iterations)
{
    int nominal_rows = N / P;
    int nominal_cols = N / Q;
    int extra_rows = N % P;
    int extra_cols = N % Q;

    int p = (world_rank < extra_rows) ? (nominal_rows+1) : nominal_rows;
    int q = (world_rank % Q < extra_cols) ? (nominal_cols+1) : nominal_cols;

    // Create the domains
    Domain even_domain(p,q,"even Domain");
    Domain odd_domain(p,q,"odd Domain");

    zero_domain(even_domain);
    zero_domain(odd_domain);

    if (world_rank == 0)
    {
        // fill in even_domain with something meaningful (initial state)
        // this requires min size for default values to fit:
        if((q >= 8) && (p >= 10))
        {
            // blinker at top left, touching right...
            even_domain(0,(q-1)) = 1;
            even_domain(0,0)     = 1;
            even_domain(0,1)     = 1;

            // and a glider:
            even_domain(8,5)     = 1;
            even_domain(8,6)     = 1;
            even_domain(8,7)     = 1;
            even_domain(7,7)     = 1;
            even_domain(6,6)     = 1;
        }
        cout << "Initial State:" << 0 << endl;
        print_domain(even_domain, world_rank);
    }

    Domain *odd, *even; // pointer swap magic
    odd = &odd_domain;
    even = &even_domain;

    for(int i = 0; i < iterations; ++i)
    {
        update_domain(*odd, *even, world_size, world_rank);
        cout << "Iteration #" << i << endl; 
        print_domain(*odd, world_rank);

        // swap pointers:
        Domain *temp = odd;
        odd  = even;
        even = temp;
    }
}

void zero_domain(Domain &domain)
{
  for(int i = 0; i < domain.rows(); ++i)
    for(int j = 0; j < domain.cols(); ++j)
      domain(i,j) = 0;
}

void print_domain(Domain &domain, int rank)
{
  cout << rank << ": " << domain.myname() << ":" <<endl;
  // this is naive; it doesn't understand big domains at all 
  for(int i = 0; i < domain.rows(); ++i)
  {
    for(int j = 0; j < domain.cols(); ++j)
      cout << (domain(i,j) ? "*" : ".");
    cout << endl;
  }
}

// this was added to reduce code cloning...
inline char update_the_cell(char cell, int neighbor_count)
{
  char newcell;
  if(cell == 0) // dead now
    newcell = (neighbor_count == 3) ? 1 : 0;
  else // was live, what about now?
    newcell = ((neighbor_count == 2)||(neighbor_count == 3)) ? 1 : 0;
									    return newcell;
}

void update_domain(Domain &new_domain, Domain &old_domain, int world_rank, int world_size)
{
    MPI_Request request[8];
    
    int m = new_domain.rows();
    int n = new_domain.cols();

    // if we wanted to be much more efficient, we would do this
    // allocation once in an object, and not allocate and deallocate
    // each iteration... for a test program, this is OK.
    char *top_row = new char[n];
    char *bottom_row = new char[n];
    char *left_column = new char[m];
    char *right_column = new char[m];

    char *top_halo = new char[n];
    char *bottom_halo = new char[n];
    char *left_halo = new char[m];
    char *right_halo = new char[m];

    const int top_row_index = 0; 
    const int bottom_row_index = m-1;
    const int left_column_index = 0;
    const int right_column_index = n-1;

    const int TOP_HALO = 0, BOTTOM_HALO = 1;
    const int LEFT_HALO = 2, RIGHT_HALO = 3;
            // so order works with 1 process.

    // int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)

    // 1. post receives for the halosprocesses
    MPI_Irecv(top_halo, n, MPI_CHAR,    (world_rank+world_size-1)%world_size, TOP_HALO, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(bottom_halo, n, MPI_CHAR, (world_rank+1)%world_size, BOTTOM_HALO, MPI_COMM_WORLD, &request[3]);
    MPI_Irecv(left_halo, m, MPI_CHAR,    (world_rank+world_size-1)%world_size, LEFT_HALO, MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(right_halo, m, MPI_CHAR, (world_rank+1)%world_size, RIGHT_HALO, MPI_COMM_WORLD, &request[5]);

    // 2. gather+send my edges to adjacent process
    for(int j = 0; j < n; ++j)   // fill the top row
    {
        top_row[j] = old_domain(top_row_index,j);
    }
    MPI_Isend(top_row, n, MPI_CHAR, (world_rank-1+world_size)%world_size, BOTTOM_HALO, MPI_COMM_WORLD, &request[0]);

    for(int j = 0; j < n; ++j) // fill in the bottom row
    {
        bottom_row[j] = old_domain(bottom_row_index,j);
    }
    MPI_Isend(bottom_row, n, MPI_CHAR, (world_rank+1)%world_size, TOP_HALO, MPI_COMM_WORLD, &request[1]);

    for(int j = 0; j < m; ++j)   // fill the left column
    {
        left_column[j] = old_domain(j, left_column_index);
    }
    MPI_Isend(left_column, m, MPI_CHAR, (world_rank-1+world_size)%world_size, RIGHT_HALO, MPI_COMM_WORLD, &request[6]);

    for(int j = 0; j < m; ++j) // fill in the right column
    {
        right_column[j] = old_domain(j, right_column_index);
    }
    MPI_Isend(right_column, m, MPI_CHAR, (world_rank+1)%world_size, LEFT_HALO, MPI_COMM_WORLD, &request[7]);
    // complete all 8 transfers
    MPI_Waitall(8, request, MPI_STATUSES_IGNORE);
    
    // at this point, I have halos from my neighbor processes

    // work on the top row with the halo received from process above:
    for(int j = 0; j < n; ++j)
    {
        int neighbor_count = 0;

        for(int delta_i = 0; delta_i <= 1; ++delta_i)
        {
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
        if(delta_i == 0 && delta_j == 0) //skip self
            continue;

        if(old_domain((top_row_index+delta_i),
                (j+delta_j+old_domain.cols())%old_domain.cols()))
            ++neighbor_count;
        }
        }
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
            if(top_halo[(j+delta_j+old_domain.cols())%old_domain.cols()])
        ++neighbor_count;
        }
        
        new_domain(top_row_index,j)
        = update_the_cell(old_domain(top_row_index,j),neighbor_count);
    } // int j

    // work on the bottom row with the halo received from the process below.
    for(int j = 0; j < new_domain.cols(); ++j)
    {
        int neighbor_count = 0;

        for(int delta_i = -1; delta_i <= 0; ++delta_i)
        {
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
        if(delta_i == 0 && delta_j == 0) //skip self
            continue;

        if(old_domain((bottom_row_index+delta_i),
                (j+delta_j+old_domain.cols())%old_domain.cols()))
            ++neighbor_count;
        }
        }
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
            if(bottom_halo[(j+delta_j+old_domain.cols())%old_domain.cols()])
        ++neighbor_count;
        }
        
        new_domain(bottom_row_index,j)
        = update_the_cell(old_domain(bottom_row_index,j), neighbor_count);
    } // int j

    // work on the left column with the halo received from process to the left:
    for(int j = 0; j < n; ++j)
    {
        int neighbor_count = 0;

        for(int delta_i = 0; delta_i <= 1; ++delta_i)
        {
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
        if(delta_i == 0 && delta_j == 0) //skip self
            continue;

        if(old_domain((left_column_index+delta_i),
                (j+delta_j+old_domain.rows())%old_domain.rows()))
            ++neighbor_count;
        }
        }
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
            if(left_column[(j+delta_j+old_domain.rows())%old_domain.rows()])
        ++neighbor_count;
        }
        
        new_domain(j, left_column_index)
        = update_the_cell(old_domain(j, left_column_index),neighbor_count);
    } // int j

    // work on the right column with the halo received from the process to the right.
    for(int j = 0; j < new_domain.rows(); ++j)
    {
        int neighbor_count = 0;

        for(int delta_i = -1; delta_i <= 0; ++delta_i)
        {
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
        if(delta_i == 0 && delta_j == 0) //skip self
            continue;

        if(old_domain((right_column_index+delta_i),
                (j+delta_j+old_domain.rows())%old_domain.rows()))
            ++neighbor_count;
        }
        }
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
            if(right_halo[(j+delta_j+old_domain.rows())%old_domain.rows()])
        ++neighbor_count;
        }
        
        new_domain(right_column_index,j)
        = update_the_cell(old_domain(j, right_column_index), neighbor_count);
    } // int j

    // the interior of the domain updates as in sequential case:
    for(int i = (top_row_index+1); i < bottom_row_index ; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
        int neighbor_count = 0;
        for(int delta_i = -1; delta_i <= 1; ++delta_i)
        {
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
        if(delta_i == 0 && delta_j == 0) //skip self
            continue;

        if(old_domain((i+delta_i+old_domain.rows())%old_domain.rows(),
                (j+delta_j+old_domain.cols())%old_domain.cols()))
            ++neighbor_count;
        }
        }
        new_domain(i,j) = update_the_cell(old_domain(i,j), neighbor_count);
        } // int j
    } // int i

    for(int i = (left_column_index+1); i < right_column_index ; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
        int neighbor_count = 0;
        for(int delta_i = -1; delta_i <= 1; ++delta_i)
        {
        for(int delta_j = -1; delta_j <= 1; ++delta_j)
        {
        if(delta_i == 0 && delta_j == 0) //skip self
            continue;

        if(old_domain((i+delta_i+old_domain.cols())%old_domain.cols(),
                (j+delta_j+old_domain.rows())%old_domain.rows()))
            ++neighbor_count;
        }
        }
        new_domain(i,j) = update_the_cell(old_domain(i,j), neighbor_count);
        } // int j
    } // int i

    // remember, in a performant code, we would encapsulate the
    // dynamic memory allocation once level higher in the code...
    delete[] bottom_halo;
    delete[] top_halo;
    delete[] right_halo;
    delete[] left_halo;
    delete[] right_column;
    delete[] left_column;
}