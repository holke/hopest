// P4est Bindings
 
#include <sc_io.h>
#include <p8est_connectivity.h>
#include <p8est_bits.h>
#include <p8est_vtk.h>
#include <p8est_mesh.h>
#include <p8est.h>
// 3D mode
#include <p4est_to_p8est.h>


static int
refine_one (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant) 
{
  if(which_tree == 0)  
  //if( which_tree == 147 | which_tree == 148 ) 
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

static int
refine_all (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  return 1;
}


void p4est_connectivity_treevertex (p4est_topidx_t num_vertices,
                                    p4est_topidx_t num_trees,
                                    double         *vertices,
                                    p4est_topidx_t *tree_to_vertex,
                                    p4est_t        **p4est_out )
{
  p4est_t              *p4est;
  p4est_topidx_t        tree;
  int                   face,i;
  p4est_connectivity_t *conn = NULL;
  sc_MPI_Comm           mpicomm;


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
  for (i=0; i < 3*num_vertices; ++i){
    conn->vertices[i] = vertices[i];
  }
  for (i=0; i < 8*num_trees; ++i){
    conn->tree_to_vertex[i] = tree_to_vertex[i];
  }

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

  /* Create a forest that is not refined; it consists of the root octant. */
  p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);
  *p4est_out=*(&p4est);

}

void p4est_refine_mesh ( p4est_t        *p4est,
                         int             refine_level,
                         int             refine_elem,
                         p4est_mesh_t   **mesh_out)
{
  p4est_mesh_t       *mesh;
  p4est_ghost_t      *ghost;
  p4est_connect_type_t mesh_btype;
  int                 level;
  int                 balance;

  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: refine_level  %d \n",refine_level);
  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: refine_elem  %d \n",refine_elem);
  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: New connectivity with %lld trees and %lld vertices\n",
     (long long) p4est->connectivity->num_trees, (long long) p4est->connectivity->num_vertices);
  /* Refine the forest iteratively, load balancing at each iteration.
   * This is important when starting with an unrefined forest */
  for (level = 0; level < refine_level; ++level) {
    if(refine_elem < 0 )
    {
      p4est_refine (p4est, 0, refine_all, NULL);
    }
    else
    {
      p4est_refine (p4est, 0, refine_one, NULL);
    }
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
  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: before first vtk %i  \n",p4est);

  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_afterrefine");
  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: after first vtk %d  \n",0);

  balance = 1;
  if (balance) {
    p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
    p4est_partition (p4est, 0, NULL);
  }
  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: before vtk %d  \n",0);

  /* Write the forest to disk for visualization, one file per processor. */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_afterbalance");

  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: before ghosts %d  \n",0);

  /* create ghost layer and mesh */
  ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  mesh = p4est_mesh_new (p4est, ghost, mesh_btype);
  //return mesh as pointer adress;
  *mesh_out=(p4est_mesh_t *)mesh;

  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: REFINE FINISHED %d  \n",0);
}

void p4est_get_mesh_info ( p4est_t        *p4est,
                           p4est_mesh_t   *mesh,
                           int            *global_num_quadrants,
                           int            *num_half_faces )
{
  *global_num_quadrants = p4est->global_num_quadrants;
  *num_half_faces = mesh->quad_to_half->elem_count;      // big face with 4 small neighbours
  SC_CHECK_ABORTF (mesh->local_num_quadrants == p4est->global_num_quadrants,
                   "Global quads %d and local quads %d mismatch ! ",
                   p4est->global_num_quadrants,mesh->local_num_quadrants );

}

void p4est_get_quadrants ( p4est_t       *p4est,
                           p4est_mesh_t   *mesh,
                           int            global_num_quadrants,
                           int            num_half_faces,
                           double         *intsize,
                           p4est_topidx_t **quad_to_tree,
                           p4est_locidx_t **quad_to_quad,
                           int8_t         **quad_to_face, 
                           p4est_locidx_t **quad_to_half, 
                           p4est_qcoord_t *quadcoords,
                           int8_t         *quadlevel ) 
{
  int num_trees = p4est->connectivity->num_trees;
  int iquad,iquadloc,iface,i;
  p8est_tree_t       *tree;
  p8est_quadrant_t   *q;
  sc_array_t         *quadrants;
  p4est_locidx_t     *halfentries;

  *intsize = (double) 1.0/P4EST_ROOT_LEN;

  P4EST_ASSERT (global_num_quadrants == p4est->local_num_quadrants);
  for (iquad = 0; iquad < mesh->local_num_quadrants; iquad++) {
    tree = p8est_tree_array_index (p4est->trees,mesh->quad_to_tree[iquad]);
    quadrants = &(tree->quadrants);
    iquadloc = iquad - tree->quadrants_offset;
    q = p8est_quadrant_array_index(quadrants, iquadloc);
    quadlevel [iquad    ] = q->level;
    quadcoords[iquad*3  ] = q->x;
    quadcoords[iquad*3+1] = q->y;
    quadcoords[iquad*3+2] = q->z;
  }

  *quad_to_tree=mesh->quad_to_tree;
  *quad_to_quad=mesh->quad_to_quad;
  *quad_to_face=mesh->quad_to_face;

  if(num_half_faces>0) *quad_to_half=mesh->quad_to_half->array;
}

void p4est_save_all ( char    filename[],
                      p4est_t *p4est)
{
  p4est_save(filename,p4est,0);
}
