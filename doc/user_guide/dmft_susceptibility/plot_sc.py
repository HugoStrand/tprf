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

from common import *
from triqs.plot.mpl_interface import oplot, oplotr, oploti, plt

with HDFArchive('data_sc.h5', 'r') as a: ps = a['ps']
p = ps.objects[-1]
    
plt.figure(figsize=(3.25*2, 8))
subp = [3, 2, 1]

plt.subplot(*subp); subp[-1] += 1
plt.plot(ps.iter, ps.dG, 's-')
plt.ylabel(r'$\max | \Delta G(\tau) |$')
plt.xlabel('Iteration')
plt.semilogy([], [])

plt.subplot(*subp); subp[-1] += 1
plt.plot(ps.iter, ps.N_up, 'v-', label=r'$N_\uparrow$')
plt.plot(ps.iter, ps.N_do, '^-', label=r'$N_\downarrow$')
plt.ylabel(r'$N_{\uparrow/\downarrow}$')
plt.xlabel('Iteration')
plt.ylim([0, 1])
plt.legend(loc='best')

plt.subplot(*subp); subp[-1] += 1
oplotr(p.G_tau_raw['up'], alpha=0.75, label=r'Binned $\uparrow$')
oplotr(p.G_tau_raw['do'], alpha=0.75, label=r'Binned $\downarrow$')
oplotr(p.G_tau['up'], label=r'DLR fit $\uparrow$')
oplotr(p.G_tau['do'], label=r'DLR fit $\downarrow$')
plt.legend(loc='best', fontsize=8)
plt.ylabel(r'$G(\tau)$')

plt.subplot(*subp); subp[-1] += 1
oplotr(p.F_tau_raw['up'], alpha=0.75, label=r'Binned $\uparrow$')
oplotr(p.F_tau_raw['do'], alpha=0.75, label=r'Binned $\downarrow$')
oplotr(p.F_tau['up'], label=r'DLR fit $\uparrow$')
oplotr(p.F_tau['do'], label=r'DLR fit $\downarrow$')
plt.legend(loc='best', fontsize=8)
plt.ylabel(r'$F(\tau)$')

plt.subplot(*subp); subp[-1] += 1
oplot(p.delta_w[0, 0])
oplot(p.delta_w[1, 1])
plt.ylabel(r'$\Delta(i \omega_n)$')
plt.xlim([-100, 100])

plt.subplot(*subp); subp[-1] += 1
oplot(p.sigma_w[0,0])
oplot(p.sigma_w[1,1])
plt.ylabel(r'$\Sigma(i\omega_n)$')
plt.xlim([-100, 100])
plt.ylim([-10, 10])

plt.tight_layout()
plt.savefig('figure_sc.svg')
plt.show()
