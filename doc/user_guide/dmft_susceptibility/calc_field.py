################################################################################
#
# TPRF: Two-Particle Response Function (TPRF) Toolbox for TRIQS
#
# Copyright (C) 2019 by The Simons Foundation
# Copyright (C) 2024 by Hugo U. R. Strand
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


""" Solve self-consistent DMFT for the Hubbard model on the square lattice
in a sweep over applied external magnetic fields B. """


from common import *


p = ParameterCollection(
    t=1.,
    B=0.,
    U=10.,
    mu=5.,
    n_k=16,
    n_iter=10,
    G_tol=1e-4,
    dlr_w_max=15.,
    dlr_eps=1e-12,
    n_iw=128,
    )

p.init = ParameterCollection(
    beta=1.,
    n_tau=127 * 4,
    gf_struct=[('up', 1), ('do', 1)],
    )

p.solve = ParameterCollection(
    length_cycle=5,
    n_warmup_cycles=int(1e7),
    n_cycles = int(1e8),
    move_double_insert_segment=False,
    move_double_remove_segment=False,
    measure_densities=True,
    measure_F_tau=True,
    )

p = setup_dmft_calculation(p)

for B in [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.4, 0.6, 0.8, 1.0]:
    p.B = B
    ps = solve_self_consistent_dmft(p)
    p = ps[-1]
    if mpi.is_master_node():
        with HDFArchive('data_B_{:f}.h5'.format(p.B), 'w') as a:
            a['ps'] = ParameterCollections(ps)
