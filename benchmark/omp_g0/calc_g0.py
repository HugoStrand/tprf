
# ----------------------------------------------------------------------

import sys
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.gf import MeshImFreq
from pytriqs.archive import HDFArchive
from pytriqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from pyed.OperatorUtils import quartic_tensor_from_operator
from pyed.OperatorUtils import quartic_permutation_symmetrize

# ----------------------------------------------------------------------

from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice import solve_rpa_PH
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    # ------------------------------------------------------------------
    # -- Discretizations
    
    n_k = (32, 32, 1)
    n_w = 400
    
    # ------------------------------------------------------------------
    # -- tight binding parameters

    beta = 10.0
    mu = 0.0
    t = 1.0
    
    h_loc = np.zeros((2, 2))        
    T = - t * np.eye(2)

    # ------------------------------------------------------------------
    # -- tight binding
    
    print '--> tight binding model'
    t_r = TBLattice(
        units = [(1, 0, 0), (0, 1, 0)],
        hopping = {
            # nearest neighbour hopping -t
            ( 0, 0): h_loc,
            ( 0,+1): T,
            ( 0,-1): T,
            (+1, 0): T,
            (-1, 0): T,
            },
        orbital_positions = [(0,0,0)]*2,
        orbital_names = ['up_0', 'do_0'],
        )

    print '--> e_k'
    e_k = t_r.on_mesh_brillouin_zone(n_k)    

    print '--> g0_wk'
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_w)

    import time

    t = time.time()
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)
    t_tot = time.time() - t
    
    nc = int(sys.argv[1])
    print 'nc =', nc
    print 'lattice_dyson_g0_wk time: ', t_tot, ' | time*cores ', t_tot * nc, 's'
