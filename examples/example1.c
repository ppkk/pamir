#include <assert.h>

#include "pamir.h"

#define SCENARIO 5

void print_my_structures(p4est_wrapper_t* p4estw, double coordinates[][NUM_DIMS]);

/** The main function of the step1 example program.
 *
 * It creates a connectivity and forest, refines it, and writes a VTK file.
 */
int
main (int argc, char **argv)
{
  p4est_wrapper_t    *p4estw;
  int                 refinement_list[2], refine_list_len;
  int                 mpi_rank, mpi_num_proc, mpiret;
  int                 num_nodes, num_elements;
  int                *lnodes;

  // TODO: unable to pass dynamicaly alocated 2D array
  //double             coordinates[][NUM_DIMS];
  double             coordinates[100][NUM_DIMS];

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;
  sc_MPI_Comm_rank(mpicomm, &mpi_rank);
  sc_MPI_Comm_size(mpicomm, &mpi_num_proc);

  // so far all scenarios for 2 processors only
  assert(mpi_num_proc == 2);

  p4estw = p4estw_create(mpicomm);

  p4estw_refine_all(p4estw);
  p4estw_partition(p4estw);

  refinement_list[0] = refinement_list[1] = 0;

  if(SCENARIO == 1)
  {
    if(mpi_rank == 1)
    {
      refinement_list[0] = refinement_list[1] = 1;
      refine_list_len = 2;
    }
  }
  else if(SCENARIO == 2)
  {
    if(mpi_rank == 1)
    {
      refinement_list[0] = 1;
      refine_list_len = 1;
    }
  }
  else if((SCENARIO == 3) || (SCENARIO == 5))
  {
    if(mpi_rank == 0)
    {
      refinement_list[1] = 1;
      refine_list_len = 1;
    }
  }

  if(SCENARIO != 4)
    p4estw_refine_selected(p4estw, refine_list_len, refinement_list);

  if(SCENARIO == 5)
    p4estw_partition(p4estw);

  p4estw_get_dimensions(p4estw, &num_elements, &num_nodes);
  lnodes = malloc(num_elements * NUM_ELEM_NODES * sizeof(int));

  //coordinates = calloc(num_nodes * NUM_DIMS,  sizeof *coordinates);

  p4estw_get_data(p4estw, lnodes, coordinates);

  p4estw_destroy(p4estw);

  print_my_structures(p4estw, coordinates);

  /* Verify that allocations internal to p4est and sc do not leak memory.
   * This should be called if sc_init () has been called earlier. */
  sc_finalize ();

  /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
