################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 by The Simons Foundation
# Author: H. U.R. Strand
#
# TPRF is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TPRF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TPRF. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import triqs_cthyb

from common import *


if mpi.is_master_node():
    with HDFArchive('data_B_0.000000.h5', 'r') as a:
        p = a['ps'].objects[-1]
else: p = None
p = mpi.bcast(p)

# -- Sample G2

p.init.n_iw = 128
p.init.delta_interface = True

p.solve = ParameterCollection(
    h_int = p.solve.h_int,
    h_loc0 = p.solve.h_loc0,
    n_cycles = int(1e8),
    #n_cycles = int(1e5),
    measure_G_l = False,
    measure_G_tau = True,
    measure_G2_iw_ph = True,
    measure_G2_blocks = set([('up','up'), ('up','do')]),
    measure_G2_n_bosonic = 1,
    measure_G2_n_fermionic = 20,
    )

cthyb = triqs_cthyb.Solver(**p.init.dict())

cthyb.Delta_tau['up'][0, 0] << p.delta_tau[0, 0]
cthyb.Delta_tau['do'][0, 0] << p.delta_tau[1, 1]

cthyb.solve(**p.solve.dict())

p.G_tau_cthyb = cthyb.G_tau.copy()
p.G2_iw_ph = cthyb.G2_iw_ph.copy()

# -- Compute DMFT impurity vertex

from triqs_tprf.linalg import inverse_PH
from triqs_tprf.chi_from_gg2 import chi0_from_gg2_PH

p.chi_m = p.G2_iw_ph[('up','up')] - p.G2_iw_ph[('up','do')]

g_w = make_gf_imfreq(p.G_w['up'], n_iw=p.init.n_iw)
p.chi0_m = chi0_from_gg2_PH(g_w, p.chi_m)
p.gamma_m = inverse_PH(p.chi0_m) - inverse_PH(p.chi_m)

del p.solve.measure_G2_blocks
if mpi.is_master_node():
    with HDFArchive('data_g2.h5', 'w') as a: a['p'] = p
