
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

import copy

import numpy as np

import triqs.utility.mpi as mpi

from h5 import HDFArchive

from triqs.gf import MeshDLRImFreq
from triqs.gf import Gf, inverse, iOmega_n
from triqs.gf import make_gf_dlr_imfreq, make_gf_dlr_imtime
from triqs.gf import fit_gf_dlr, make_gf_imfreq, make_gf_imtime

from triqs_tprf.lattice import lattice_dyson_g_w
from triqs_tprf.ParameterCollection import ParameterCollection
from triqs_tprf.ParameterCollection import ParameterCollections
from triqs_tprf.utilities import BlockGf_data

from triqs_ctseg import Solver


def setup_dmft_calculation(p):

    p = copy.deepcopy(p)
    p.iter = 0

    # -- Local Hubbard interaction
    from triqs.operators import n, Operator
    
    p.h_loc0 = -0.5*p.U*(n('up', 0) + n('do', 0))
    #p.h_loc0 = Operator()
    p.solve.h_int = p.U*n('up', 0)*n('do', 0)

    # -- 2D square lattice with nearest neighbour hopping t

    from triqs_tprf.tight_binding import TBLattice
    
    T = -p.t * np.eye(2)
    p.H_loc = p.B * np.diag([+1, 1])
    H = TBLattice(
        units=[(1, 0, 0), (0, 1, 0)],
        orbital_positions=[(0,0,0)]*2,
        orbital_names=['up', 'do'],
        hopping={(0, +1) : T, (0, -1) : T, (+1, 0) : T, (-1, 0) : T,
                 #(0, 0): p.H_loc
                 })
    
    kmesh = H.get_kmesh(n_k=p.n_k)
    p.e_k = H.fourier(kmesh)

    # -- Initial zero guess for the self-energy

    p.wmesh = MeshDLRImFreq(p.init.beta, 'Fermion', p.dlr_w_max, p.dlr_eps)
    p.sigma_w = Gf(mesh=p.wmesh, target_shape=[2, 2])
    p.sigma_w.zero()
    p.sigma_w[:] += p.mu * np.eye(2)
    
    return p


def solve_self_consistent_dmft(p):

    ps = []
    for dmft_iter in range(p.n_iter):
        p = dmft_self_consistent_step(p)
        ps.append(p)
        mpi.report('-'*72)
        mpi.report(f'--> DMFT Iteration {p.iter}: dG = {p.dG:2.2E}, dF = {p.dF:2.2E}, dM = {p.dM:2.2E}')
        mpi.report('-'*72)
        if p.dG < p.G_tol: break

    mpi.report('='*72)
    if dmft_iter >= p.n_iter - 1:
        mpi.report(f'--> Warning: DMFT Not converged! dG = {p.dG:2.2E}, dF = {p.dF:2.2E}')
    else:
        mpi.report(f'--> DMFT Converged in {p.iter} iterations: dG = {p.dG:2.2E}, dF = {p.dF:2.2E}')
    mpi.report('='*72)
    
    return ps


def dmft_self_consistent_step(p):

    p = copy.deepcopy(p)
    p.iter += 1

    # -- Local lattice Green's function
    
    p.g_loc_w = lattice_dyson_g_w(p.mu, p.e_k, p.sigma_w - p.H_loc) # DEBUG
    #p.g_loc_w = lattice_dyson_g_w(p.mu, p.e_k, p.sigma_w - p.H_loc)
    p.rho_loc = p.g_loc_w.density()
    mpi.report(f'rho_loc =\n{p.rho_loc}')

    # -- Setup impurity solver

    from triqs.operators import n
    p.solve.h_loc0 =  p.h_loc0 - p.B * (n('up', 0) - n('do', 0))
    
    S = Solver(**p.init.dict())

    # -- Hybridization function
    
    p.delta_w = p.g_loc_w.copy()
    p.delta_w << iOmega_n - inverse(p.g_loc_w) - p.sigma_w + p.H_loc + p.mu * np.eye(2)
    
    p.delta_tau = make_gf_imtime(p.delta_w, n_tau=p.init.n_tau)
    p.delta_tau.data[:] = p.delta_tau.data.real                

    S.Delta_tau['up'][0, 0] << p.delta_tau[0, 0]
    S.Delta_tau['do'][0, 0] << p.delta_tau[1, 1]

    # -- Solve impurity problem
    
    S.solve(**p.solve.dict())

    # -- Get sampled propagators
    
    p.G_tau_raw = S.results.G_tau.copy()
    p.F_tau_raw = S.results.F_tau.copy()
    
    p.G_dlr = fit_gf_dlr(p.G_tau_raw, p.dlr_w_max, p.dlr_eps)
    p.F_dlr = fit_gf_dlr(p.F_tau_raw, p.dlr_w_max, p.dlr_eps)

    # -- Check convergence
    
    G_tau = make_gf_dlr_imtime(p.G_dlr)
    F_tau = make_gf_dlr_imtime(p.F_dlr)

    p.dG = np.max(np.abs(BlockGf_data(G_tau - p.G_tau))) if hasattr(p,'G_tau') else float('nan')
    p.dF = np.max(np.abs(BlockGf_data(F_tau - p.F_tau))) if hasattr(p,'F_tau') else float('nan')

    p.G_tau = G_tau
    p.F_tau = F_tau

    # -- Compute local self-energy from impurity propagators (\Sigma = F / G)
    
    p.G_w = make_gf_dlr_imfreq(p.G_dlr)
    p.F_w = make_gf_dlr_imfreq(p.F_dlr)
    
    p.g_imp_w = p.g_loc_w.copy()
    p.g_imp_w[0, 0] << p.G_w['up'][0, 0]
    p.g_imp_w[1, 1] << p.G_w['do'][0, 0]

    p.f_w = p.g_loc_w.copy()
    p.f_w[0, 0] << p.F_w['up'][0, 0]
    p.f_w[1, 1] << p.F_w['do'][0, 0]

    p.sigma_w << p.f_w * inverse(p.g_imp_w)

    # -- local observables
    
    p.N_up = S.results.densities['up'][0]
    p.N_do = S.results.densities['do'][0]
    M_old = p.M if hasattr(p, 'M') else float('nan')
    p.M = 0.5*(p.N_up - p.N_do).real
    p.dM = np.abs(p.M - M_old)

    mpi.report(f'N_up = {p.N_up}, N_do = {p.N_do}, M = {p.M}')
    mpi.report(f'densities = {S.results.densities}')
    
    return p
