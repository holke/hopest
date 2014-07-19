// P4est Bindings
 
#include <sc_io.h>
#include <p8est_connectivity.h>
#include <p8est_bits.h>
#include <p8est_vtk.h>
#include <p8est_mesh.h>
// 3D mode
#include <p4est_to_p8est.h>


// extern"C"{
void test( double *areal, int *bint)
{
  P4EST_GLOBAL_PRODUCTIONF ("test areal= %f bint= %i \n",  *areal,  *bint);
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  if(which_tree == 0)  // | for naca which_tree == 147 | which_tree == 148 ) 
  {
    return 1;
  }
  else
  {
    return 0;
  }
}


void p4est_connectivity_treevertex (p4est_topidx_t * num_vertices_in,
                               p4est_topidx_t * num_trees_in        ,
                               double         *vertices             ,
                               p4est_topidx_t * tree_to_vertex      ,
                               int *p4est_out )
                               //p4est_t *p4est )
{
  p4est_topidx_t      tree;
  int                 face;
  p4est_topidx_t      num_vertices = *num_vertices_in;
  p4est_topidx_t      num_trees = *num_trees_in;

  p4est_connectivity_t *conn = NULL;

  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;



  /* Initialize MPI; see sc_mpi.h.
   * If configure --enable-mpi is given these are true MPI calls.
   * Else these are dummy functions that simulate a single-processor run. */
  mpicomm = sc_MPI_COMM_WORLD;

  /* These functions are optional.  If called they store the MPI rank as a
   * static variable so subsequent global p4est log messages are only issued
   * from processor zero.  Here we turn off most of the logging; see sc.h. */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_PRODUCTION);



  P4EST_GLOBAL_PRODUCTIONF ("creating connectivity from tree to vertex only...%d %d \n",num_vertices,num_trees);



  conn = p4est_connectivity_new (num_vertices, num_trees,
                                 0, 0,
                                 0, 0);
  conn->vertices = vertices;
  conn->tree_to_vertex = tree_to_vertex;

  /*
   * Fill tree_to_tree and tree_to_face to make sure we have a valid
   * connectivity.
   */
  for (tree = 0; tree < conn->num_trees; ++tree) {
    for (face = 0; face < P4EST_FACES; ++face) {
      conn->tree_to_tree[P4EST_FACES * tree + face] = tree;
      conn->tree_to_face[P4EST_FACES * tree + face] = face;
    }
  }

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));

  p4est_connectivity_complete (conn);

  P4EST_ASSERT (p4est_connectivity_is_valid (conn));


  P4EST_GLOBAL_PRODUCTIONF
    ("New connectivity with %lld trees and %lld vertices\n",
     (long long) conn->num_trees, (long long) conn->num_vertices);


  //return conn;

  /* Create a forest that is not refined; it consists of the root octant. */
  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

  *p4est_out = p4est;
  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG p4est %d  p4est_out %d \n", p4est,*p4est_out);

}

void p4est_refine_mesh ( int        *p4est_in,
                         int        *mesh_out ) 
{
  p4est_t            *p4est;
  p4est_mesh_t       *mesh;
  p4est_ghost_t      *ghost;
  p4est_connect_type_t mesh_btype;
  int                 level;
  int                 balance;
  static int          refine_level = 3;

  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG p4est_in %d  \n", *p4est_in);
  p4est = *p4est_in;

  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG:  %d \n",0);
  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: New connectivity with %lld trees and %lld vertices\n",
     (long long) p4est->connectivity->num_trees, (long long) p4est->connectivity->num_vertices);
  /* Refine the forest iteratively, load balancing at each iteration.
   * This is important when starting with an unrefined forest */
  for (level = 0; level < refine_level; ++level) {
    p4est_refine (p4est, 0, refine_fn, NULL);
    /* Refinement has lead to up to 8x more elements; redistribute them. */
    p4est_partition (p4est, 0, NULL);
  }

  /* If we call the 2:1 balance we ensure that neighbors do not differ in size
   * by more than a factor of 2.  This can optionally include diagonal
   * neighbors across edges or corners as well; see p4est.h.
   *
   * Note that this balance step is not strictly necessary since we are using
   * uniform refinement but may be required for other types of refinement.
   */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_afterrefine");
  balance = 1;
  if (balance) {
    p4est_balance (p4est, P4EST_CONNECT_FACE, NULL);
    p4est_partition (p4est, 0, NULL);
  }

  /* Write the forest to disk for visualization, one file per processor. */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_afterbalance");

  /* create ghost layer and mesh */
  // ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FACE);
  // mesh = p4est_mesh_new (p4est, ghost, mesh_btype);
  // *mesh_out=mesh;

 
  
}
//}
