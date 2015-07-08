#include <stdio.h>
#include <assert.h>

#ifndef P4_TO_P8
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>
#include <p4est_ghost.h>
#else
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p8est_iterate.h>
#include <p8est_lnodes.h>
#include <p8est_ghost.h>
#endif

#include "pamir.h"

/***************************************************************************************************/
/***************************************************************************************************/
// my structures describing the mesh
// all structures are extremely stupid, slow, large, etc... Just a "proof of concept"

#define NUM_ELEM_FACES  4
#define NUM_ELEM_CHILDREN  4
#define MAX_NODES 1000
#define MAX_ELEMENTS 1000
#define MAX_ADJACENCIES 3

typedef struct adjacency
{
  int ext_process;
  int ext_global_element;
  int ext_local_node;
}
adjacency_t;

typedef struct element
{
  int global_idx;
  int vertices[NUM_ELEM_NODES];
  int is_constrained[NUM_ELEM_NODES];
  int constraining_nodes[NUM_ELEM_NODES][2];
  int is_on_boundary;
  int num_adjacencies[NUM_ELEM_NODES];
  adjacency_t adjacencies[NUM_ELEM_NODES][MAX_ADJACENCIES];
}
element_t;


element_t elements[MAX_ELEMENTS];

// todo: replace by some kind of hash table
int hanging_indices[MAX_NODES][MAX_NODES];

// todo: neco
int n_only_hang, n_only_reg, n_both_hang_reg;
int *num_hanging_references, *num_regular_references;


/***************************************************************************************************/
// initialization, print, etc

void init_my_structures()
{
  int i, j, k;

  for(i = 0; i < MAX_NODES; i++)
  {
    for(j = 0; j < MAX_NODES; j++)
    {
      hanging_indices[i][j] = -1;
    }
  }

  for(i = 0; i < MAX_ELEMENTS; i++)
  {
    for(j = 0; j < NUM_ELEM_NODES; j++)
    {
      elements[i].vertices[j] = -1;
      elements[i].is_constrained[j] = 0;
      elements[i].is_on_boundary = 0;
      elements[i].num_adjacencies[j] = 0;
      for(k = 0; k < MAX_ADJACENCIES; k++)
      {
        elements[i].adjacencies[j][k].ext_global_element = -1999;
        elements[i].adjacencies[j][k].ext_local_node = -1999;
      }
    }
  }
  for(i = 0; i < MAX_NODES; i++)
    for(j = 0; j < MAX_NODES; j++)
      hanging_indices[i][j] = -1;
}

void print_my_structures(p4est_wrapper_t* p4estw, double coordinates[][NUM_DIMS])
{
  FILE *f;
  int node, elem, i_adj;
  char str[10000], filename[30];
  sprintf(str, "\n**********************************************************************************************\n");
  sprintf(str, "%s********                       DATA ON PARTITION %d                                  *********\n", str, p4estw->mpi_rank);
  sprintf(str, "%s**********************************************************************************************\n\n", str);
//  sprintf(str, "%s  vertex: ", str);
  for(node = 0; node < p4estw->num_full_nodes; node++)
  {
    sprintf(str, "%s  node %d -> (%3.2lf, %3.2lf)\n", str, node, coordinates[node][0], coordinates[node][1]);
  }
  if(p4estw->num_nodes > p4estw->num_full_nodes)
    sprintf(str, "%shanging nodes added (not in p4est):\n", str);
  for(node = p4estw->num_full_nodes; node < p4estw->num_nodes; node++)
  {
    sprintf(str, "%s  node %d -> (%3.2lf, %3.2lf)\n", str, node, coordinates[node][0], coordinates[node][1]);
  }
  sprintf(str, "%s\n", str);

  for(elem = 0; elem < p4estw->num_elements; elem++)
  {
    sprintf(str, "%s  elem %d (loc %d) : vertices [%d, %d, %d, %d] -> [(%3.2lf, %3.2lf), (%3.2lf, %3.2lf), (%3.2lf, %3.2lf), (%3.2lf, %3.2lf)]\n", str, elements[elem].global_idx, elem,
            elements[elem].vertices[0], elements[elem].vertices[1], elements[elem].vertices[2], elements[elem].vertices[3],
        coordinates[elements[elem].vertices[0]][0], coordinates[elements[elem].vertices[0]][1],
        coordinates[elements[elem].vertices[1]][0], coordinates[elements[elem].vertices[1]][1],
        coordinates[elements[elem].vertices[2]][0], coordinates[elements[elem].vertices[2]][1],
        coordinates[elements[elem].vertices[3]][0], coordinates[elements[elem].vertices[3]][1]);

    for(node = 0; node < NUM_ELEM_NODES; node++)
    {
      if(elements[elem].is_constrained[node])
      {
        sprintf(str, "%s      node %d constrained by %d and %d\n", str, elements[elem].vertices[node], elements[elem].constraining_nodes[node][0], elements[elem].constraining_nodes[node][1]);
      }
    }
    for(node = 0; node < NUM_ELEM_NODES; node++)
    {
      if(elements[elem].num_adjacencies[node] > 0)
      {
        sprintf(str, "%s      node %d on an interface to other processes\n", str, elements[elem].vertices[node]);
        for(i_adj = 0; i_adj < elements[elem].num_adjacencies[node]; i_adj++)
        {
          sprintf(str, "%s          connected to element %d through node with local id %d\n", str, elements[elem].adjacencies[node][i_adj].ext_global_element, elements[elem].adjacencies[node][i_adj].ext_local_node);
        }
      }
    }
  }

  sprintf(filename, "list_%d.txt", p4estw->mpi_rank);
  f = fopen(filename, "w");
  fprintf(f, "%s", str);
  fclose(f);
}



/***************************************************************************************************/
/***************************************************************************************************/

static void
init_element_data (p4est_t * p4est, p4est_topidx_t which_tree,
                              p4est_quadrant_t * q)
{
  /* the data associated with a quadrant is accessible by p.user_data */
  element_data_t       *data = (element_data_t *) q->p.user_data;

  data->global_idx = -12345;
  data->local_idx = -12345;
  data->refinement_flag = REF_FLAG_INVALID;
}

int get_global_offset(p4est_wrapper_t *p4estw, int local_num_elems)
{
  int send_buf[1], rec_buf[p4estw->mpi_num_proc];
  int i, offset;

  send_buf[0] = local_num_elems;
  sc_MPI_Allgather(send_buf, 1, sc_MPI_INT, rec_buf, 1, sc_MPI_INT, p4estw->p4est->mpicomm);

  offset = 0;
  for(i = 0; i < p4estw->mpi_rank; i++)
    offset += rec_buf[i];

  return offset;
}



static void
face_iter_callback (p4est_iter_face_info_t * info, void *user_data)
{
  int                 i, j;
  p4est_t            *p4est = info->p4est;
//  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
//  step3_data_t       *ghost_data = (step3_data_t *) user_data;
//  step3_data_t       *udata;
  double              q;
  double              h, facearea;
  int                 which_face[2] = {-88, -88}, is_hanging[2], quad_x[2] = {-55,-55}, quad_y[2] = {-55,-55};
  p4est_quadrant_t   *quad[2];
  p4est_iter_face_side_t *side[2];
  sc_array_t         *sides = &(info->sides);
  p4est_locidx_t      quadids[2];      /**< index in tree or ghost array */


  /* which of the quadrant's faces the interface touches */
  for(i = 0; i < sides->elem_count; i++)
  {
    side[i] = p4est_iter_fside_array_index_int (sides, i);
    which_face[i] = side[i]->face;
    is_hanging[i] = side[i]->is_hanging;
    if(is_hanging[i])
    {
    }
    else
    {
      quad[i] = side[i]->is.full.quad;
      quadids[i] = side[i]->is.full.quadid;
//      quad_x[i] = quad[i]->x;
//      quad_y[i] = quad[i]->y;
    }
  }
//  printf("    FACE  (MPI RANK %d) NUM SIDES %d \n", mpi_rank, sides->elem_count);
//  for(i = 0; i < sides->elem_count; i++)
//  {
//    printf("         ->   SIDE %d : QUADID: %d, loc face: %d, hang:  %d,     (%d, %d) \n", i,  quadids[i], which_face[i], is_hanging[i], quad_x[i], quad_y[i]);
//  }

}

static void
vertex_iter_callback (p4est_iter_corner_info_t * info, void *user_data)
{
  int                 i, j, i_loc, i_ghost, adj_id;
  p4est_t            *p4est = info->p4est;
//  step3_ctx_t        *ctx = (step3_ctx_t *) p4est->user_pointer;
  element_data_t       *ghost_data = (element_data_t *) user_data;
  element_data_t       *udata;
  p4est_iter_corner_side_t *side;
  sc_array_t         *sides = &(info->sides);
  p4est_locidx_t      local_quadids[4], ghost_quadids[4];
  int                 local_loc_ver[4], ghost_loc_ver[4];
  int                 local_global_id[4], ghost_global_id[4];
  int                 num_local_elems, num_ghost_elems;


  num_local_elems = num_ghost_elems = 0;
  for(i = 0; i < 4; i++)
    local_loc_ver[i] = ghost_loc_ver[i] = local_quadids[i] = ghost_quadids[i] = -99;

  for(i = 0; i < sides->elem_count; i++)
  {
    side = p4est_iter_cside_array_index_int (sides, i);

    if(side->is_ghost)
    {
      udata = &ghost_data[side->quadid];
      ghost_quadids[num_ghost_elems] = side->quadid;
      ghost_global_id[num_ghost_elems] = udata->global_idx;
      ghost_loc_ver[num_ghost_elems] = side->corner;
      num_ghost_elems++;
    }
    else
    {
      udata = (element_data_t *) side->quad->p.user_data;
      local_quadids[num_local_elems] = side->quadid;
      local_global_id[num_local_elems] = udata->global_idx;
      local_loc_ver[num_local_elems] = side->corner;
      num_local_elems++;
    }
  }

  assert(num_local_elems > 0);
  if(num_ghost_elems > 0)
  {
    for(i_loc = 0; i_loc < num_local_elems; i_loc++)
    {
      for(i_ghost = 0; i_ghost < num_ghost_elems; i_ghost++)
      {
        adj_id = elements[local_quadids[i_loc]].num_adjacencies[local_loc_ver[i_loc]];
        elements[local_quadids[i_loc]].adjacencies[local_loc_ver[i_loc]][adj_id].ext_global_element = ghost_global_id[i_ghost];
        elements[local_quadids[i_loc]].adjacencies[local_loc_ver[i_loc]][adj_id].ext_local_node = ghost_loc_ver[i_ghost];

        elements[local_quadids[i_loc]].num_adjacencies[local_loc_ver[i_loc]] = adj_id + 1;
      }
    }
  }

//  printf("    VERTEX  (MPI RANK %d) NUM SIDES %d \n", mpi_rank, sides->elem_count);
//  for(i = 0; i < sides->elem_count; i++)
//  {
//    printf("         ->   SIDE %d : QUADID: %d, loc face: %d, hang:  %d,     (%d, %d) \n", i,  quadids[i], which_face[i], is_hanging[i], quad_x[i], quad_y[i]);
//  }

}


static int
refine_all_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  return 1;
}


// TODO: is there another way to get information to callback than global data?
// TODO: or I could base the decision on global element idx, which is available in quadrant
int* refine_array;
int total_to_be_refined;

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  int will_refine;
  int refine_flag;
  element_data_t     *elem_data;

  elem_data = (element_data_t* ) quadrant->p.user_data;
  //printf("refine     %d, %d, %d, %d\n", quadrant->x, quadrant->y, elem_data->local_idx, elem_data->refinement_flag);

  will_refine = 0;
//  refine_flag = refine_array[elem_data->local_idx];
  refine_flag = elem_data->refinement_flag;

  if( refine_flag == REF_FLAG_REFINE)
  {
    will_refine = 1;
    total_to_be_refined = total_to_be_refined + 1;
  }
  else if((refine_flag == REF_FLAG_NOTHING) || (refine_flag == REF_FLAG_COARSEN))
  {
  }
  else
  {
   // assert(0);
  }

  return  will_refine;
}


static int
coarsen_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrants[])
{
  int will_coarsen;
  int refine_flag;
  element_data_t     *elem_data;
  int idx;

  will_coarsen = 1;
  //printf("\ncoarsen:\n");
  for(idx = 0; idx < NUM_ELEM_CHILDREN; idx++)
  {
    elem_data = (element_data_t* ) ((quadrants[idx])->p.user_data);
    //printf("     %d, %d, %d, %d\n", quadrants[idx]->x, quadrants[idx]->y, elem_data->local_idx, elem_data->refinement_flag);
    //refine_flag = refine_array[elem_data->local_idx];
    refine_flag = elem_data->refinement_flag;

    if(refine_flag == REF_FLAG_COARSEN)
    {
    }
    else if((refine_flag == REF_FLAG_NOTHING) || (refine_flag == REF_FLAG_REFINE))
    {
      will_coarsen = 0;
    }
    else
    {
      assert(0);
    }

  }

  if(will_coarsen)
    printf("\nCOARSENED!!!!\n\n");

  return will_coarsen;
}


void destroy_additional_data(p4est_wrapper_t* p4estw)
{
  if(p4estw->ghost != NULL)
  {
    p4est_ghost_destroy (p4estw->ghost);
    p4estw->ghost = NULL;
  }

  if(p4estw->ghost_data != NULL)
  {
    P4EST_FREE (p4estw->ghost_data);
    p4estw->ghost_data = NULL;
  }

  if(p4estw->lnodes != NULL)
  {
    p4est_lnodes_destroy (p4estw->lnodes);
    p4estw->lnodes = NULL;
  }
}

void create_ghost_data(p4est_wrapper_t* p4estw)
{
  if((p4estw->ghost == NULL) || (p4estw->ghost_data == NULL))
  {
    assert(p4estw->ghost == NULL);
    assert(p4estw->ghost_data == NULL);

    /* Create the ghost layer to learn about parallel neighbors. */
    p4estw->ghost = p4est_ghost_new (p4estw->p4est, P4EST_CONNECT_FULL);
    /* create space for storing the ghost data */
    p4estw->ghost_data = P4EST_ALLOC (element_data_t, p4estw->ghost->ghosts.elem_count);
    /* synchronize the ghost data */
    p4est_ghost_exchange_data (p4estw->p4est, p4estw->ghost, p4estw->ghost_data);
  }
}



/********************************************************************************************/
//           inteface to Fempar:
/********************************************************************************************/

p4est_wrapper_t* p4estw_create(sc_MPI_Comm mpicomm)
{
  p4est_connectivity_t *conn;
  p4est_wrapper_t *p4estw = malloc(sizeof(p4est_wrapper_t));
  //sc_MPI_Comm mpicomm = MPI_Comm_f2c(fortran_mpicomm);

  printf("*********** creating p4est\n");
  conn = p4est_connectivity_new_unitsquare ();
  p4estw->p4est = p4est_new (mpicomm, conn, sizeof(element_data_t), init_element_data, NULL);
  p4estw->ghost = NULL;
  p4estw->ghost_data = NULL;
  p4estw->lnodes = NULL;
  p4estw->info_found = 0;
  sc_MPI_Comm_rank(p4estw->p4est->mpicomm, &p4estw->mpi_rank);
  sc_MPI_Comm_size(p4estw->p4est->mpicomm, &p4estw->mpi_num_proc);

  return p4estw;
}

int p4estw_refine_all(p4est_wrapper_t* p4estw)
{
  destroy_additional_data(p4estw);

  int recursive = 0;

  printf("*********** refine all \n");
  p4est_refine (p4estw->p4est, recursive, refine_all_fn, init_element_data);
  p4est_balance (p4estw->p4est, P4EST_CONNECT_FACE, init_element_data);
  return 0;
}

int p4estw_refine_selected(p4est_wrapper_t* p4estw, int what_refine_size, int* what_refine)
{
  p4est_topidx_t      tt;       /* Connectivity variables have this type. */
  p4est_locidx_t      k, q, Q;  /* Process-local counters have this type. */
  p4est_tree_t       *tree;     /* Pointer to one octree */
  p4est_quadrant_t   *quad;
  sc_array_t         *tquadrants;       /* Quadrant array for one tree */
  element_data_t     *elem_data;

  for (tt = p4estw->p4est->first_local_tree, k = 0;
       tt <= p4estw->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4estw->p4est->trees, tt);   /* Current tree */
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t) tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k) {
      /* This is now a loop over all local elements.
       * Users might aggregate the above code into a more compact iterator. */
      quad = p4est_quadrant_array_index (tquadrants, q);
      elem_data = (element_data_t* ) quad->p.user_data;
//      printf("(%d, %d)\n", elem_data->local_idx, what_refine_size);
      assert((elem_data->local_idx >= 0) && (elem_data->local_idx < what_refine_size));
      elem_data->refinement_flag = what_refine[elem_data->local_idx];
    }
  }



  destroy_additional_data(p4estw);

  int recursive = 0;

  refine_array = what_refine;
  total_to_be_refined = 0;

  p4est_coarsen (p4estw->p4est, recursive, coarsen_fn, init_element_data);
  p4est_refine (p4estw->p4est, recursive, refine_fn, init_element_data);
  p4est_balance (p4estw->p4est, P4EST_CONNECT_FACE, init_element_data);

  printf("P4est: refined %d\n", total_to_be_refined);

  return 0;
}

int p4estw_partition(p4est_wrapper_t* p4estw)
{
  destroy_additional_data(p4estw);

  int partforcoarsen = 0;
  p4est_partition (p4estw->p4est, partforcoarsen, NULL);
  return 0;
}

int p4estw_get_num_elements(p4est_wrapper_t *p4estw, int *num_elements)
{
  *num_elements = p4estw->p4est->local_num_quadrants;
  return 0;
}

int p4estw_get_dimensions(p4est_wrapper_t* p4estw, int* num_elements, int* num_nodes)
{
  int                 anyhang, hanging_face[P4EST_CHILDREN], hanging_edge[12];
  p4est_topidx_t      tt;       /* Connectivity variables have this type. */
  p4est_locidx_t      k, q, Q;  /* Process-local counters have this type. */
  p4est_tree_t       *tree;     /* Pointer to one octree */
  p4est_quadrant_t   *quad, node_quadrant;
  sc_array_t         *tquadrants;       /* Quadrant array for one tree */
  int node, face, constraining0, constraining1, constrained, local_constrained_idx,  dim;
  double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
  element_data_t     *udata; // user data
  int num_hanging_nodes = 0;
  int localnode, elem, idx = 0;

  create_ghost_data(p4estw);
  assert(p4estw->lnodes == NULL);
  /* Create a node numbering for continuous linear finite elements. */
  p4estw->lnodes = p4est_lnodes_new (p4estw->p4est, p4estw->ghost, 1);

  init_my_structures();

  idx = 0;
  for(elem = 0; elem < p4estw->lnodes->num_local_elements; elem++)
  {
    for(localnode = 0; localnode < p4estw->lnodes->vnodes; localnode++)
    {
      elements[elem].vertices[localnode] = p4estw->lnodes->element_nodes[idx];
      idx++;
    }
  }


  num_hanging_references = malloc(p4estw->lnodes->num_local_nodes * sizeof(int));
  num_regular_references = malloc(p4estw->lnodes->num_local_nodes * sizeof(int));
  memset(num_hanging_references, 0, p4estw->lnodes->num_local_nodes * sizeof(int));
  memset(num_regular_references, 0, p4estw->lnodes->num_local_nodes * sizeof(int));

  for (tt = p4estw->p4est->first_local_tree, k = 0;
       tt <= p4estw->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4estw->p4est->trees, tt);   /* Current tree */
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t) tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k) {
      /* This is now a loop over all local elements.
       * Users might aggregate the above code into a more compact iterator. */
      quad = p4est_quadrant_array_index (tquadrants, q);

      /* We need to determine whether any node on this element is hanging. */
      anyhang = p4est_lnodes_decode (p4estw->lnodes->face_code[q], hanging_face
#ifdef P4_TO_P8
                                     , hanging_edge
#endif
                                     );

#ifdef P4_TO_P8
      assert(0);  // TODO: distinguis middle of the edge and middle of the face
#endif

      for(face = 0; face < NUM_ELEM_FACES; face++)
      {
        constraining0 = elements[k].vertices[p4est_face_corners[face][0]];
        constraining1 = elements[k].vertices[p4est_face_corners[face][1]];
        if((!anyhang) || (hanging_face[face] == -1))
        {
          num_regular_references[constraining0]++;
          num_regular_references[constraining1]++;
          printf("        *** ceddss  elem %d, face %d     (%d, %d), (reg, reg)    on rank %d\n", k, face, constraining0, constraining1, p4estw->mpi_rank);

        }
        else if(hanging_face[face] == 0)
        {
          num_regular_references[constraining0]++;
          num_hanging_references[constraining1]++;
          printf("        *** ceddss  elem %d, face %d     (%d, %d), (reg, hang)    on rank %d\n", k, face, constraining0, constraining1, p4estw->mpi_rank);
        }
        else if(hanging_face[face] == 1)
        {
          num_hanging_references[constraining0]++;
          num_regular_references[constraining1]++;
          printf("        *** ceddss  elem %d, face %d     (%d, %d), (hang, reg)    on rank %d\n", k, face, constraining0, constraining1, p4estw->mpi_rank);
        }
        else
        {
          printf(" ffff   %d\n", hanging_face[face]);
          assert(0);
        }
      }
    }
  }

  char str[2000];
  sprintf(str, "-----------------rank %d : ", p4estw->mpi_rank);
  for(idx = 0; idx < p4estw->lnodes->num_local_nodes; idx++)
  {
    sprintf(str, "%s(%d, %d), ", str, num_regular_references[idx], num_hanging_references[idx]);
  }
  printf("%s\n", str);

  n_only_hang = n_only_reg = n_both_hang_reg = 0;
  for(idx = 0; idx < p4estw->lnodes->num_local_nodes; idx++)
  {
    if(num_regular_references[idx] > 0)
    {
      if(num_hanging_references[idx] > 0)
        n_both_hang_reg++;
      else
        n_only_reg++;
    }
    else
    {
      if(num_hanging_references[idx] > 0)
        n_only_hang++;
      else
        assert(0);
    }
  }
  p4estw->num_full_nodes = n_only_hang + n_only_reg;
  *num_nodes = p4estw->num_nodes = n_only_hang + n_only_reg + n_both_hang_reg;
  *num_elements = p4estw->num_elements = p4estw->lnodes->num_local_elements;
  printf("++++++++++++++++++++rank %d : only hang %d, only reg  %d, both %d\n", p4estw->mpi_rank, n_only_hang, n_only_reg, n_both_hang_reg);

  return 0;
}

int p4estw_get_data(p4est_wrapper_t* p4estw, int* list_nodes, double coordinates[][NUM_DIMS])
{
  int                 anyhang, hanging_face[P4EST_CHILDREN], hanging_edge[12];
  p4est_topidx_t      tt;       /* Connectivity variables have this type. */
  p4est_locidx_t      k, q, Q;  /* Process-local counters have this type. */
  p4est_tree_t       *tree;     /* Pointer to one octree */
  p4est_quadrant_t   *quad, node_quadrant;
  sc_array_t         *tquadrants;       /* Quadrant array for one tree */
  int node, face, constraining0, constraining1, constrained, local_constrained_idx,  dim;
  double              vxyz[3];  /* We embed the 2D vertices into 3D space. */
  element_data_t     *udata; // user data
  int localnode, elem, idx = 0;
  int correction, num_nodes;

  num_nodes = p4estw->lnodes->num_local_nodes;

  for (tt = p4estw->p4est->first_local_tree, k = 0;
       tt <= p4estw->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4estw->p4est->trees, tt);   /* Current tree */
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t) tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k) {
      /* This is now a loop over all local elements.
       * Users might aggregate the above code into a more compact iterator. */
      quad = p4est_quadrant_array_index (tquadrants, q);

      /* We need to determine whether any node on this element is hanging. */
      anyhang = p4est_lnodes_decode (p4estw->lnodes->face_code[q], hanging_face
#ifdef P4_TO_P8
                                     , hanging_edge
#endif
                                     );
      if (anyhang)
      {
         printf("hanging %d -> [%d, %d, %d, %d]\n", q, hanging_face[0], hanging_face[1], hanging_face[2], hanging_face[3]);
         for(face = 0; face < NUM_ELEM_FACES; face++)
         {
           if(hanging_face[face] > -1)
           {
             constraining0 = elements[k].vertices[p4est_face_corners[face][0]];
             constraining1 = elements[k].vertices[p4est_face_corners[face][1]];
             if(hanging_indices[constraining0][constraining1] == -1)
             {
               assert(num_regular_references[constraining0] + num_regular_references[constraining1] > 0);
               if(num_regular_references[constraining0] == 0)
               {
                 constrained = constraining0;
               }
               else if(num_regular_references[constraining1] == 0)
               {
                 constrained = constraining1;
               }
               else
               {
                 constrained = num_nodes;
                 num_nodes++;
               }
               hanging_indices[constraining0][constraining1] =
                                  hanging_indices[constraining1][constraining0] = constrained;
             }
             else
             {
               constrained = hanging_indices[constraining0][constraining1];
             }

             if(hanging_face[face] == 0)
               local_constrained_idx = p4est_face_corners[face][1];
             else if(hanging_face[face] == 1)
               local_constrained_idx = p4est_face_corners[face][0];
             else
               assert(0);

             elements[k].vertices[local_constrained_idx] = constrained;
             elements[k].is_constrained[local_constrained_idx] = 1;
             elements[k].constraining_nodes[local_constrained_idx][0] = constraining0;
             elements[k].constraining_nodes[local_constrained_idx][1] = constraining1;

           }
         }
      }
    }
  }

  p4estw->global_elem_idx_offset = get_global_offset(p4estw, k);

  // second loop through elements
  for (tt = p4estw->p4est->first_local_tree, k = 0;
       tt <= p4estw->p4est->last_local_tree; ++tt) {
    tree = p4est_tree_array_index (p4estw->p4est->trees, tt);   /* Current tree */
    tquadrants = &tree->quadrants;
    Q = (p4est_locidx_t) tquadrants->elem_count;
    for (q = 0; q < Q; ++q, ++k) {
      /* This is now a loop over all local elements.
       * Users might aggregate the above code into a more compact iterator. */
      quad = p4est_quadrant_array_index (tquadrants, q);

      // assign global index to the user data
      udata = (element_data_t *) quad->p.user_data;
      udata->global_idx = k + p4estw->global_elem_idx_offset;
      elements[k].global_idx = k + p4estw->global_elem_idx_offset;

      for(node = 0; node < NUM_ELEM_NODES; node++)
      {
        p4est_quadrant_corner_node (quad, node, &node_quadrant);

        p4est_qcoord_to_vertex (p4estw->p4est->connectivity, tt, node_quadrant.x, node_quadrant.y,
      #ifdef P4_TO_P8
                                node_quadrant.z,
      #endif
                                vxyz);
        for(dim = 0; dim < NUM_DIMS; dim++)
        {
          coordinates[elements[k].vertices[node]][dim] = vxyz[dim];
        }
      }

    }
  }

  /* synchronize the ghost data */
  p4est_ghost_exchange_data (p4estw->p4est, p4estw->ghost, p4estw->ghost_data);


  /// 3. iterate to find additional info
  ///
  p4est_iterate (p4estw->p4est,  p4estw->ghost, (void *) p4estw->ghost_data,
                 NULL, //volume_iter_callback,    /* callback function that interpolate from the cell center to the cell corners, defined above */
                 face_iter_callback,          /* there is no callback for the faces between quadrants */
#ifdef P4_TO_P8
                 NULL,          /* there is no callback for the edges between quadrants */
#endif
                 vertex_iter_callback);         /* there is no callback for the corners between quadrants */




  idx = 0;
  for(elem = 0; elem < p4estw->lnodes->num_local_elements; elem++)
  {
    for(localnode = 0; localnode < p4estw->lnodes->vnodes; localnode++)
    {
      // fempar and p4est have different local nodes numbering
      // TODO : do it in a more systematic way
      if(NUM_DIMS == 2)
      {
        if(localnode == 2)
          correction = 1;
        else if(localnode == 3)
          correction = -1;
        else
          correction = 0;
      }
      else
        assert(0);

      // transform from 0-based to 1-based numbering
      //list_nodes[idx] = p4estw->lnodes->element_nodes[idx+correction] + 1;
      list_nodes[idx] = elements[elem].vertices[localnode + correction] + 1;
      idx++;
    }
  }

  return 0;
}


int p4estw_destroy(p4est_wrapper_t* p4estw)
{
  destroy_additional_data(p4estw);

  p4est_t            *p4est = p4estw->p4est;
  p4est_connectivity_t *conn = p4est->connectivity;
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  return 0;
}
