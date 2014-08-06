// P4est Bindings
 
#include <sc_io.h>
#include <p8est_connectivity.h>
#include <p8est_bits.h>
#include <p8est_vtk.h>
#include <p8est_mesh.h>
#include <p8est.h>
// 3D mode
#include <p4est_to_p8est.h>


/* Theses functions are required for performing the mesh refinement in p4est.
   The dummy function refine_hopest is used a refinement function for p4est.
   The actual refinement functions are in HOPEST, called through refine_f.
*/

int (*refine_f) (p4est_qcoord_t,p4est_qcoord_t,p4est_qcoord_t,p4est_topidx_t,int8_t,int);

static int
refine_hopest (p4est_t * p4est, p4est_topidx_t which_tree,
               p4est_quadrant_t * q)
{
  // Call HOPEST refinemtn routines
  return refine_f(q->x,q->y,q->z,which_tree+1,q->level,p8est_quadrant_child_id(q));
}


void p4_refine_mesh(p4est_t  *p4est,
                    int     (*myrefine_f)
                            (p4est_qcoord_t,p4est_qcoord_t,p4est_qcoord_t,p4est_topidx_t,int8_t,int),
                    int     refine_level,
                    p4est_mesh_t  **mesh_out )
{
  int                 level;
  int                 balance;

  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: refine_level  %d \n",refine_level);
  P4EST_GLOBAL_PRODUCTIONF
    ("DEBUG: New connectivity with %lld trees and %lld vertices\n",
     (long long) p4est->connectivity->num_trees, (long long) p4est->connectivity->num_vertices);

  // Set refinement function pointer called by refine_hopest to function provided by HOPEST
  refine_f=myrefine_f;
  
  /* Refine the forest iteratively, load balancing at each iteration.
   * This is important when starting with an unrefined forest */
  for (level = 0; level < refine_level; ++level) {
    p4est_refine (p4est, 0, &refine_hopest, NULL);
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
}

