_dmft_susceptibility:

DMFT lattice susceptibility
===========================

In this guide we will compute the uniform magnetic susceptibility of the single band Hubbard model on the square lattice with nearest neighbour hopping using dynamical mean-field theory (DMFT).

First we compute the uniform magnetic susceptibility by appyling an external magnetic field :math:`B` and measuring the resulting magnetization :math:`M`. The susceptibility can then be obtained by extrapolating the first derivative at zero field, i.e.

.. math::
   \chi = \left. \frac{dM}{dB} \right|_{B \rightarrow 0}

To do this we setup a small framework for doing self consistent DMFT calculations.

First we set up the lattice and the initial (zero) guess for the self-energy :math:`\Sigma`.

.. literalinclude:: common.py
   :lines: 1-10		
