# Generated automatically using the command :
# c++2py ../../c++/triqs_tprf/lattice.hpp --members_read_only -N triqs_tprf -a triqs_tprf -m lattice -o lattice -C pytriqs --moduledoc="Lattice functionality" --cxxflags="-std=c++17"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "lattice", doc = "Lattice functionality", app_name = "triqs_tprf")

# Imports
import pytriqs.gf
import pytriqs.lattice

# Add here all includes
module.add_include("triqs_tprf/lattice.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/tuple.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_tprf;
""")

module.add_enum("Channel_t", ['Channel_t::PP', 'Channel_t::PH', 'Channel_t::PH_bar'], "triqs_tprf", """Two-particle channel enum class, PP (particle-particle), PH (particle-hole), PH_bar (particle-hole-bar)""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::lattice_dyson_g0_wk (double mu, triqs_tprf::e_k_cvt e_k, gf_mesh<triqs::gfs::imfreq> mesh)", doc = """Construct a non-interacting Matsubara frequency lattice Green\'s function :math:`G^{(0)}_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`\n\n  Computes\n\n  .. math::\n     G^{(0)}_{a\\bar{b}}(i\\omega_n, \\mathbf{k}) = \\left[ \n         (i\\omega_n + \\mu ) \\cdot \\mathbf{1}  - \\epsilon(\\mathbf{k})\n\t \\right]^{-1}_{a\\bar{b}},\n\n  using a discretized dispersion :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`, chemical potential :math:`\\mu`,\n  and a Matsubara frequency Green\'s function mesh.\n\n  :param mu: chemical potential :math:`\\mu`\n  :param e_k: discretized lattice dispersion :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`\n  :param mesh: imaginary frequency mesh\n  @return Matsubara frequency lattice Green\'s function :math:`G^{(0)}_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_w_cvt sigma_w)", doc = """Construct an interacting Matsubara frequency lattice Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`\n   \n Computes\n\n .. math::\n    G_{a\\bar{b}}(i\\omega_n, \\mathbf{k}) = \\left[ \n        (i\\omega_n + \\mu ) \\cdot \\mathbf{1}  - \\epsilon(\\mathbf{k}) - \\Sigma(i\\omega_n)\n\t\\right]^{-1}_{a\\bar{b}},\n\n using a discretized dispersion :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`, \n chemical potential :math:`\\mu`, and a momentum independent Matsubara frequency \n self energy :math:`\\Sigma_{\\bar{a}b}(i\\omega_n)`.\n\n :param mu: chemical potential :math:`\\mu`\n :param e_k: discretized lattice dispersion :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`\n :param sigma_w: imaginary frequency self-energy :math:`\\Sigma_{\\bar{a}b}(i\\omega_n)`\n @return Matsubara frequency lattice Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::lattice_dyson_g_wk (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_wk_cvt sigma_wk)", doc = """Construct an interacting Matsubara frequency lattice Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`\n   \n Computes\n\n .. math::\n    G_{a\\bar{b}}(i\\omega_n, \\mathbf{k}) = \\left[ \n        (i\\omega_n + \\mu ) \\cdot \\mathbf{1}  - \\epsilon(\\mathbf{k}) - \\Sigma(i\\omega_n, \\mathbf{k})\n\t\\right]^{-1}_{a\\bar{b}},\n\n using a discretized dispersion :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`, \n chemical potential :math:`\\mu`, and a momentum independent Matsubara frequency \n self energy :math:`\\Sigma_{\\bar{a}b}(i\\omega_n, \\mathbf{k})`.\n\n :param mu: chemical potential :math:`\\mu`\n :param e_k: discretized lattice dispersion :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`\n :param sigma_wk: imaginary frequency self-energy :math:`\\Sigma_{\\bar{a}b}(i\\omega_n, \\mathbf{k})`\n @return Matsubara frequency lattice Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`""")

module.add_function ("triqs_tprf::g_w_t triqs_tprf::lattice_dyson_g_w (double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_w_cvt sigma_w)", doc = """Construct an interacting Matsubara frequency local (:math:`\\mathbf{r}=\\mathbf{0}`) lattice Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n)`\n   \n Computes\n\n .. math::\n    G_{a\\bar{b}}(i\\omega_n) = \\frac{1}{N_k} \\sum_\\mathbf{k} \\left[ \n        (i\\omega_n + \\mu ) \\cdot \\mathbf{1}  - \\epsilon(\\mathbf{k}) - \\Sigma(i\\omega_n)\n\t\\right]^{-1}_{a\\bar{b}},\n\n using a discretized dispersion :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`, \n chemical potential :math:`\\mu`, and a momentum independent Matsubara frequency \n self energy :math:`\\Sigma_{\\bar{a}b}(i\\omega_n)`.\n\n :param mu: chemical potential :math:`\\mu`\n :param e_k: discretized lattice dispersion :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`\n :param sigma_w: imaginary frequency self-energy :math:`\\Sigma_{\\bar{a}b}(i\\omega_n)`\n @return Matsubara frequency lattice Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`""")

module.add_function ("triqs_tprf::g_wr_t triqs_tprf::fourier_wk_to_wr (triqs_tprf::g_wk_cvt g_wk)", doc = """Inverse fast fourier transform of imaginary frequency Green\'s function from k-space to real space\n\n    Computes: :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{r}) = \\mathcal{F}^{-1} \\left\\{ G_{a\\bar{b}}(i\\omega_n, \\mathbf{k}) \\right\\}`\n\n    :param g_wk: k-space imaginary frequency Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`\n    @return real-space imaginary frequency Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{r})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::fourier_wr_to_wk (triqs_tprf::g_wr_cvt g_wr)", doc = """Fast fourier transform of imaginary frequency Green\'s function from real-space to k-space\n\n    Computes: :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{k}) = \\mathcal{F} \\left\\{ G_{a\\bar{b}}(i\\omega_n, \\mathbf{r}) \\right\\}`\n\n    :param g_wr: real-space imaginary frequency Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{r})`\n    @return k-space imaginary frequency Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{k})`""")

module.add_function ("triqs_tprf::g_tr_t triqs_tprf::fourier_wr_to_tr (triqs_tprf::g_wr_cvt g_wr, int nt = -1)", doc = """Fast fourier transform of real-space Green\'s function from Matsubara frequency to imaginary time\n\n    Computes: :math:`G_{a\\bar{b}}(\\tau, \\mathbf{r}) = \\mathcal{F} \\left\\{ G_{a\\bar{b}}(i\\omega_n, \\mathbf{r}) \\right\\}`\n\n    :param g_wr: real-space imaginary frequency Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{r})`\n    @return real-space imaginary time Green\'s function :math:`G_{a\\bar{b}}(\\tau, \\mathbf{r})`""")

module.add_function ("triqs_tprf::g_wr_t triqs_tprf::fourier_tr_to_wr (triqs_tprf::g_tr_cvt g_tr, int nw = -1)", doc = """Fast fourier transform of real-space Green\'s function from imaginary time to Matsubara frequency\n\n    Computes: :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{r}) = \\mathcal{F} \\left\\{ G_{a\\bar{b}}(\\tau, \\mathbf{r}) \\right\\}`\n\n    :param g_tr: real-space imaginary time Green\'s function :math:`G_{a\\bar{b}}(\\tau, \\mathbf{r})`\n    @return real-space Matsubara frequency Green\'s function :math:`G_{a\\bar{b}}(i\\omega_n, \\mathbf{r})`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::lindhard_chi00_wk (triqs_tprf::e_k_cvt e_k, int nw, double beta, double mu)", doc = """Generalized Lindhard susceptibility in the particle-hole channel :math:`\\chi^{(00)}_{\\bar{a}b\\bar{c}d}(i\\omega_n, \\mathbf{q})`.\n   \n   Analytic calculation of the generalized (non-interacting) Lindhard susceptibility \n   in the particle-hole channel. The analytic expression is obtained using residue calculus \n   to explicitly evaluate the matsubara sum of the fourier transformed imaginary time\n   bubble product of two non-interacting single-particle Green\'s functions.\n\n   .. math::\n      G^{(0)}_{a\\bar{b}}(\\mathbf{k}, i\\omega_n) =\n      \\left[ i\\omega_n \\cdot \\mathbf{1} - \\epsilon(\\mathbf{k}) \\right]^{-1} .\n\n   The analytic evaluation of the bubble diagram gives\n\n   .. math::\n        \\chi^{(00)}_{\\bar{a}b\\bar{c}d}(i\\omega_n, \\mathbf{q}) \\equiv \n        \\mathcal{F} \\left\\{\n\t  - G^{(0)}_{d\\bar{a}}(\\tau, \\mathbf{r}) G^{(0)}_{b\\bar{c}}(-\\tau, -\\mathbf{r})\n\t\\right\\}\n        = \n\t- \\frac{1}{N_k} \\sum_{\\nu} \\sum_{\\mathbf{k}}\n          G^{(0)}_{d\\bar{a}}(\\nu, \\mathbf{k}) \n\t  G^{(0)}_{b\\bar{c}}(\\nu + \\omega, \\mathbf{k} + \\mathbf{q})\n\t\\\\ =\n\t- \\frac{1}{N_k} \\sum_{\\nu} \\sum_{\\mathbf{k}} \n\t  \\left( \\sum_{i}\n          U^\\dagger_{di}(\\mathbf{k}) \\frac{1}{i\\nu - \\epsilon_{\\mathbf{k}, i}} U_{i\\bar{a}}(\\mathbf{k}) \n\t  \\right)\n\t  \\left( \\sum_j\n\t  U^\\dagger_{bj}(\\mathbf{k} + \\mathbf{q})\n\t  \\frac{1}{i\\nu + i\\omega - \\epsilon_{\\mathbf{k} + \\mathbf{q}, j}}\n\t  U_{j\\bar{c}}(\\mathbf{k} + \\mathbf{q}) \n\t  \\right)\n\t\\\\ =  \n\t\\frac{1}{N_k} \\sum_{\\mathbf{k}} \\sum_{ij} \n\t  \\left(\n\t    [1-\\delta(\\epsilon_{\\mathbf{k},i} - \\epsilon_{\\mathbf{k}+\\mathbf{q}, j})]\n\t    \\frac{ f(\\epsilon_{\\mathbf{k}, i}) - f(\\epsilon_{\\mathbf{k}+\\mathbf{q}, j}) }\n\t         {i\\omega_n + \\epsilon_{\\mathbf{k} + \\mathbf{q}, j} - \\epsilon_{\\mathbf{k}, i}}\n\t    +\n  \t    \\delta(\\epsilon_{\\mathbf{k},i} - \\epsilon_{\\mathbf{k}+\\mathbf{q}, j})\n\t    \\frac{\\beta}{4 \\cosh^2 (\\beta \\epsilon_{\\mathbf{k}, i} / 2) }\n\t  \\right)\n\t  \\\\ \\times\n\t  U_{i\\bar{a}}(\\mathbf{k}) U^\\dagger_{di}(\\mathbf{k}) \n\t  U_{j\\bar{c}}(\\mathbf{k} + \\mathbf{q}) U^\\dagger_{bj}(\\mathbf{k} + \\mathbf{q})\n\n   where the :math:`U(\\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued \n   dispersion relation :math:`\\epsilon_{\\bar{a}b}(\\mathbf{k})`, i.e.\n\n   .. math::\n      \\sum_{\\bar{a}b} U_{i\\bar{a}}(\\mathbf{k}) \\epsilon_{\\bar{a}b}(\\mathbf{k}) U^\\dagger_{bj} (\\mathbf{k})\n      = \\delta_{ij} \\epsilon_{\\mathbf{k}, i}\n\n   .. note::\n      The analytic formula is sub-optimal in terms of performance for higher temperatures. The evaluation\n      scales as :math:`\\mathcal{O}(N_k^2)` which is worse than computing the bubble explicitly in imaginary \n      time, with scaling :math:`\\mathcal{O}(N_k N_\\tau \\log(N_k N_\\tau)` for :math:`N_k \\gg N_\\tau`.\n\n   .. note::\n      Care must be taken when evaluating the fermionic Matsubara frequency sum of the\n      product of two simple poles. By extending the sum to an integral over the complex \n      plane the standard expression for the Lindhard response is obtained when the \n      poles are non-degenerate. The degenerate case produces an additional frequency independent\n      contribution (the last term on the last row).""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::solve_rpa_PH (triqs_tprf::chi_wk_vt chi0, array_view<std::complex<double>,4> U)", doc = """Random Phase Approximation (RPA) in the particle-hole channel\n   \n     Computes the equation\n\n     .. math::\n         \\chi(\\bar{a}b\\bar{c}d) = \\big(\n         \\mathbb{1} \n         - \\chi^{(0)}(\\bar{a}b\\bar{B}A) U(A\\bar{B}D\\bar{C})\n         \\big)^{-1} \\chi^{(0)}(\\bar{C}D\\bar{c}d)\\,.\n     \n     :param chi0: bare particle-hole bubble :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\mathbf{k}, i\\omega_n)`\n     :param U: RPA static vertex as obtained from triqs_tprf.rpa_tensor.get_rpa_tensor :math:`U_{a\\bar{b}c\\bar{d}}`\n     @return RPA suceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\mathbf{k}, i\\omega_n)`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::retarded_screened_interaction_Wr_wk (triqs_tprf::chi_wk_cvt PI_wk, triqs_tprf::chi_k_cvt V_k)", doc = """Screened interaction :math:`W(i\\omega_n, \\mathbf{k})` calculator \n    for static momentum dependent interaction :math:`V(\\mathbf{k})`\n\n    .. math::\n        W_{abcd}(i\\omega_n, \\mathbf{k}) = \n          V_{abcd}(\\mathbf{k}) +\n\t  \\sum_{efgh} V_{abef}(\\mathbf{k}) \\cdot\n          \\Pi_{fegh}(i\\omega_n, \\mathbf{k}) \\cdot\n          W_{hgcd}(i\\omega_n, \\mathbf{k})\n\n    Instead of returning :math:`W` we return the retarded part\n    :math:`W^{(r)}` (with zero high-frequency offset)\n    \n    .. math::\n        W^{(r)}_{abcd}(i\\omega_n, \\mathbf{k}) = \n            W_{abcd}(i\\omega_n, \\mathbf{k}) - V_{abcd}(\\mathbf{k})\n\n    :param PI_wk: polarization bubble :math:`\\Pi_{abcd}(i\\omega_n, \\mathbf{k})`\n    :param V_k: static interaction :math:`V_{abcd}(\\mathbf{k})`\n    @return retarded screened interaction :math:`W^{(r)}_{abcd}(i\\omega_n, \\mathbf{k})`""")

module.add_function ("triqs_tprf::g_wk_t triqs_tprf::gw_self_energy (triqs_tprf::chi_wk_cvt Wr_wk, triqs_tprf::g_wk_cvt g_wk)", doc = """GW self energy :math:`\\Sigma(i\\omega_n, \\mathbf{k})` calculator \n\n    Fourier transforms the screened interaction and the single-particle\n    Green\'s function to imagiary time and real space.\n\n    .. math::\n        G_{ab}(\\tau, \\mathbf{r}) = \\mathcal{F}^{-1}\n          \\left\\{ G_{ab}(i\\omega_n, \\mathbf{k}) \\right\\}\n\n    .. math::\n        W^{(r)}_{abcd}(\\tau, \\mathbf{r}) = \\mathcal{F}^{-1}\n          \\left\\{ W^{(r)}_{abcd}(i\\omega_n, \\mathbf{k}) \\right\\}\n\n    computes the GW self-energy as the product\n\n    .. math::\n        \\Sigma_{ab}(\\tau, \\mathbf{r}) =\n          \\sum_{cd} W^{(r)}_{abcd}(\\tau, \\mathbf{r}) G_{cd}(\\tau, \\mathbf{r})\n\n    and transforms back to frequency and momentum\n\n    .. math::\n        \\Sigma_{ab}(i\\omega_n, \\mathbf{k}) =\n          \\mathcal{F} \\left\\{ \\Sigma_{ab}(\\tau, \\mathbf{r}) \\right\\}\n\n    :param V_k: static bare interaction :math:`V_{abcd}(\\mathbf{k})`\n    :param Wr_wk: retarded screened interaction :math:`W^{(r)}_{abcd}(i\\omega_n, \\mathbf{k})`\n    :param g_wk: single particle Green\'s function :math:`G_{ab}(i\\omega_n, \\mathbf{k})`\n    @return GW self-energy :math:`\\Sigma_{ab}(i\\omega_n, \\mathbf{k})`""")

module.add_function ("triqs_tprf::g_tr_t triqs_tprf::gw_sigma_tr (triqs_tprf::chi_tr_cvt Wr_tr, triqs_tprf::g_tr_cvt g_tr)", doc = """GW self energy :math:`\\Sigma(\\tau, \\mathbf{r})` calculator \n\n    Computes the GW self-energy as the product\n\n    .. math::\n        \\Sigma_{ab}(\\tau, \\mathbf{r}) =\n          \\sum_{cd} W^{(r)}_{abcd}(\\tau, \\mathbf{r}) G_{cd}(\\tau, \\mathbf{r})\n\n    :param Wr_tr: retarded screened interaction :math:`W^{(r)}_{abcd}(\\tau, \\mathbf{r})`\n    :param g_tr: single particle Green\'s function :math:`G_{ab}(\\tau, \\mathbf{r})`\n    @return GW self-energy :math:`\\Sigma_{ab}(\\tau, \\mathbf{r})`""")

module.add_function ("triqs_tprf::gk_iw_t triqs_tprf::eliashberg_product (triqs_tprf::chi_wk_vt Gamma_pp, triqs_tprf::gk_iw_vt g_wk, triqs_tprf::gk_iw_vt delta_wk)", doc = """Linearized Eliashberg product\n\n     Computes the product\n\n     .. math::\n         \\Delta^{(out)}_{\\bar{a}\\bar{b}}(\\mathbf{k},i\\nu) =  -\\frac{1}{N_k \\beta}\\sum_{\\mathbf{k}\'} \\sum_{i\\nu\'}\n\t \\Gamma_{A\\bar{a}B\\bar{b}}(\\mathbf{k}-\\mathbf{k}\', i\\nu - i\\nu\')\n\t \\\\ \\times\n\t G_{A\\bar{c}}(\\mathbf{k}\', i\\nu\')\n\t \\Delta_{\\bar{c}\\bar{d}}(\\mathbf{k}\', i\\nu\')\n\t G_{B\\bar{d}}(-\\mathbf{k}\', -i\\nu\')\n\n     :param chi_pp: particle-particle vertex :math:`\\Gamma^{(pp)}_{a\\bar{b}c\\bar{d}}(\\mathbf{k}, i\\nu_n)`\n     :param g_kw: single particle Green\'s function :math:`G_{a\\bar{b}}(\\mathbf{k}, i\\nu_n)`\n     :param delta_kw: pairing self-energy :math:`\\Delta_{\\bar{a}\\bar{b}}(\\mathbf{k}, i\\nu_n)`\n     @return Gives the result of the product :math:`\\Delta^{(out)} \\sim \\Gamma^{(pp)}GG \\Delta`""")

module.add_function ("triqs_tprf::gk_iw_t triqs_tprf::eliashberg_product_fft (triqs_tprf::chi_tr_vt Gamma_pp_dyn_tr, triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::gk_iw_vt g_wk, triqs_tprf::gk_iw_vt delta_wk)", doc = """Linearized Eliashberg product via FFT\n\n     Computes the product\n\n     .. math::\n         \\Delta^{(out)}_{\\bar{a}\\bar{b}}(\\mathbf{k},i\\nu) =  -\\frac{1}{N_k \\beta}\\sum_{\\mathbf{k}\'} \\sum_{i\\nu\'}\n\t \\Gamma_{A\\bar{a}B\\bar{b}}(\\mathbf{k}-\\mathbf{k}\', i\\nu - i\\nu\')\n\t \\\\ \\times\n\t G_{A\\bar{c}}(\\mathbf{k}\', i\\nu\')\n\t \\Delta_{\\bar{c}\\bar{d}}(\\mathbf{k}\', i\\nu\')\n\t G_{B\\bar{d}}(-\\mathbf{k}\', -i\\nu\')\\,,\n\n     by taking advantage of the convolution theorem.\n\n     We therefore first calculate\n\n     .. math::\n        \\Delta^{(out)}_{\\bar{a}\\bar{b}}(\\mathbf{r}, \\tau) = \n\t -\\Gamma_{A\\bar{a}B\\bar{b}}(\\mathbf{r}, \\tau) F_{AB}(\\mathbf{r}, \\tau) \\,,\n\n     where \n\n     .. math::\n        F_{AB}(\\mathbf{r}, \\tau)  = \n        \\mathcal{F}\\big(G_{A\\bar{c}}(\\mathbf{k}\', i\\nu\')\n\t \\Delta_{\\bar{c}\\bar{d}}(\\mathbf{k}\', i\\nu\')\n\t G_{B\\bar{d}}(-\\mathbf{k}\', -i\\nu\')\\big)\\,.\n\n     Then we Fourier transform \n\n     .. math::\n          \\Delta^{(out)}_{\\bar{a}\\bar{b}}(\\mathbf{k},i\\nu) = \n          \\mathcal{F}\\big(\\Delta^{(out)}_{\\bar{a}\\bar{b}}(\\mathbf{r}, \\tau)\\big)\\,,\n\n    to get the same result, but with far less computational effort.\n\n     :param chi_rt: dynamic part of the particle-particle vertex :math:`\\Gamma^{(pp)}_{a\\bar{b}c\\bar{d}}(\\mathbf{r}, \\tau)`\n     :param chi_r: constant part of the particle-particle vertex :math:`\\Gamma^{(pp)}_{a\\bar{b}c\\bar{d}}(\\mathbf{r})`\n     :param g_kw: single particle Green\'s function :math:`G_{a\\bar{b}}(\\mathbf{k}, i\\nu_n)`\n     :param delta_kw: pairing self-energy :math:`\\Delta_{\\bar{a}\\bar{b}}(\\mathbf{k}, i\\nu_n)`\n     @return Gives the result of the product :math:`\\Delta^{(out)} \\sim \\Gamma^{(pp)}GG \\Delta`""")

module.add_function ("triqs_tprf::gk_iw_t triqs_tprf::eliashberg_g_delta_g_product (triqs_tprf::gk_iw_vt g_wk, triqs_tprf::gk_iw_vt delta_wk)", doc = """""")

module.add_function ("std::tuple<chi_wk_vt,chi_k_vt> triqs_tprf::split_into_dynamic_wk_and_constant_k (triqs_tprf::chi_wk_vt Gamma_pp)", doc = """""")

module.add_function ("std::tuple<chi_tr_vt,chi_r_vt> triqs_tprf::dynamic_and_constant_to_tr (triqs_tprf::chi_wk_vt Gamma_pp_dyn_wk, triqs_tprf::chi_k_vt Gamma_pp_const_k)", doc = """""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::gamma_PP_singlet (triqs_tprf::chi_wk_vt chi_c, triqs_tprf::chi_wk_vt chi_s, array_view<std::complex<double>,4> U_c, array_view<std::complex<double>,4> U_s)", doc = """Gamma particle-particle singlet\n\n     Computes the particle-particle vertex for singlet pairing in the RPA limit\n\n    .. math::\n        \\Gamma^{(\\mathrm{singlet})}(a\\bar{b}c\\bar{d}) =\n        \\frac{3}{2} U^{(\\mathrm{s})}(a\\bar{b}A\\bar{B}) \\chi^{(\\mathrm{s})}(\\bar{B}A\\bar{C}D) \n        U^{(\\mathrm{s})}(D\\bar{C}c\\bar{d}) \\\\\n        -\\frac{1}{2} U^{(\\mathrm{c})}(a\\bar{b}A\\bar{B}) \\chi^{(\\mathrm{c})}(\\bar{B}A\\bar{C}D) \n        U^{(\\mathrm{c})}(D\\bar{C}c\\bar{d}) \\\\\n       + \\frac{1}{2}\\big(U^{(\\mathrm{s})}(a\\bar{b}c\\bar{d})+\n        U^{(\\mathrm{c})}(a\\bar{b}c\\bar{d})\\big)\n\n     :param chi_c: charge susceptibility  :math:`\\chi^{(\\mathrm{c})}_{\\bar{a}b\\bar{c}d}(\\mathbf{k}, i\\omega_n)`\n     :param chi_s: spin susceptibility  :math:`\\chi^{(\\mathrm{s})}_{\\bar{a}b\\bar{c}d}(\\mathbf{k}, i\\omega_n)`\n     :param U_c: charge interaction  :math:`U^{(\\mathrm{c})}_{a\\bar{b}c\\bar{d}}`\n     :param U_s: spin interaction  :math:`U^{(\\mathrm{s})}_{a\\bar{b}c\\bar{d}}`\n     @return :math:`\\Gamma^{(\\mathrm{singlet})}_{a\\bar{b}c\\bar{d}}(\\mathbf{k}, i\\omega_n)`""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::gamma_PP_triplet (triqs_tprf::chi_wk_vt chi_c, triqs_tprf::chi_wk_vt chi_s, array_view<std::complex<double>,4> U_c, array_view<std::complex<double>,4> U_s)", doc = """Gamma particle-particle triplet\n\n     Computes the particle-particle vertex for triplet pairing in the RPA limit\n\n    .. math::\n        \\Gamma^{(\\mathrm{triplet})}(a\\bar{b}c\\bar{d}) =\n        -\\frac{1}{2} U^{(\\mathrm{s})}(a\\bar{b}A\\bar{B}) \\chi^{(\\mathrm{s})}(\\bar{B}A\\bar{C}D) \n        U^{(\\mathrm{s})}(D\\bar{C}c\\bar{d}) \\\\\n        -\\frac{1}{2} U^{(\\mathrm{c})}(a\\bar{b}A\\bar{B}) \\chi^{(\\mathrm{c})}(\\bar{B}A\\bar{C}D) \n        U^{(\\mathrm{c})}(D\\bar{C}c\\bar{d}) \\\\\n       + \\frac{1}{2}\\big(U^{(\\mathrm{s})}(a\\bar{b}c\\bar{d})+\n        U^{(\\mathrm{c})}(a\\bar{b}c\\bar{d})\\big)\n\n     :param chi_c: charge susceptibility  :math:`\\chi^{(\\mathrm{c})}_{\\bar{a}b\\bar{c}d}(\\mathbf{k}, i\\omega_n)`\n     :param chi_s: spin susceptibility  :math:`\\chi^{(\\mathrm{s})}_{\\bar{a}b\\bar{c}d}(\\mathbf{k}, i\\omega_n)`\n     :param U_c: charge interaction  :math:`U^{(\\mathrm{c})}_{a\\bar{b}c\\bar{d}}`\n     :param U_s: spin interaction  :math:`U^{(\\mathrm{s})}_{a\\bar{b}c\\bar{d}}`\n     @return :math:`\\Gamma^{(\\mathrm{triplet})}_{a\\bar{b}c\\bar{d}}(\\mathbf{k}, i\\omega_n)`""")

module.add_function ("array<std::complex<double>,6> triqs_tprf::cluster_mesh_fourier_interpolation (array<double,2> k_vecs, triqs_tprf::chi_wr_cvt chi)", doc = """""")

module.add_function ("triqs_tprf::chi_tr_t triqs_tprf::chi0_tr_from_grt_PH (triqs_tprf::g_tr_cvt g_tr)", doc = """Generalized susceptibility imaginary time bubble in the particle-hole channel :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r})`\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r}) =\n     - G_{d\\bar{a}}(\\tau, \\mathbf{r}) G_{b\\bar{c}}(-\\tau, -\\mathbf{r})\n\n  :param g_tr: Imaginary time Green\'s function in real-space, :math:`G_{a\\bar{b}}(\\tau, \\mathbf{r})`.\n  @return Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r})` in imaginary time and real-space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi0_w0r_from_grt_PH (triqs_tprf::g_tr_cvt g_tr)", doc = """Generalized susceptibility zero imaginary frequency bubble in the particle-hole channel :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega=0, \\mathbf{r})`\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\mathbf{r}) =\n     - \\int_0^\\beta d\\tau \\,\n     G_{d\\bar{a}}(\\tau, \\mathbf{r}) G_{b\\bar{c}}(-\\tau, -\\mathbf{r})\n\n  :param g_tr: Imaginary time Green\'s function in real-space, :math:`G_{a\\bar{b}}(\\tau, \\mathbf{r})`.\n  @return Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\mathbf{r})` in real-space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi_w0r_from_chi_tr (triqs_tprf::chi_tr_cvt chi_tr)", doc = """Static susceptibility calculation :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega=0, \\mathbf{r})`\n   \n  Explicit calculation of the static, zero frequency response, by 2nd order trapetzoidal \n  integration in imaginary time, i.e.\n  \n  .. math::\n     \\chi_{\\bar{a}b\\bar{c}d}(\\omega=0, \\mathbf{r}) =\n         \\int_0^\\beta d\\tau \\, \\chi_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r}) \n\n  :param chi_tr: Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r})` \n                in imaginary time and real space.\n  @return Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega=0, \\mathbf{r})` \n          at zero Matsubara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi_wr_from_chi_tr (triqs_tprf::chi_tr_cvt chi_tr, int nw)", doc = """Parallell Fourier transform from  :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r})` to :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r})`\n\n  Computes\n\n  .. math::\n     \\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r}) =\n         \\mathcal{F}_{\\tau \\rightarrow \\omega} \\left\\{\n         \\chi_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r}) \n         \\right\\}\n\n  :param chi_tr: Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r})` \n                in imaginary time and real space.\n  @return Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r})` \n          in Matsubara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_tr_t triqs_tprf::chi_tr_from_chi_wr (triqs_tprf::chi_wr_cvt chi_wr, int ntau = -1)", doc = """Fourier transform from :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r})` to :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r})`\n\n  Computes\n\n  .. math::\n         \\chi_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r}) =\n         \\mathcal{F}_{\\omega \\rightarrow \\tau} \\left\\{\n\t \\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r}) =\n         \\right\\}\n\n  :param chi_tr: Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\tau, \\mathbf{r})` \n                in imaginary time and real space.\n  @return Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r})` \n          in Matsubara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::chi_wk_from_chi_wr (triqs_tprf::chi_wr_cvt chi_wr)", doc = """Parallell Fourier transform from :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r})` to :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})`\n\n  Computes\n\n  .. math::\n     \\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k}) =\n         \\mathcal{F}_{\\mathbf{r} \\rightarrow \\mathbf{k}} \\left\\{\n         \\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r}) \n         \\right\\}\n\n  :param chi_wr: Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r})` \n                in Matsubara frequency and real space.\n  @return Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})` \n          in Matsubara frequency and momentum space.""")

module.add_function ("triqs_tprf::chi_wr_t triqs_tprf::chi_wr_from_chi_wk (triqs_tprf::chi_wk_cvt chi_wk)", doc = """Parallell Fourier transform from :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})` to :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r})`\n\n  Computes\n\n  .. math::\n     \\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r}) =\n         \\mathcal{F}_{\\mathbf{k} \\rightarrow \\mathbf{r}} \\left\\{\n         \\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k}) \n         \\right\\}\n\n  :param chi_wr: Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})` \n                in imaginary time and momentum space.\n  @return Generalized susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{r})` \n          in Matsubara frequency and real space.""")

module.add_function ("chi_t_t::zero_t triqs_tprf::chi_trapz_tau (triqs_tprf::chi_t_cvt chi_t)", doc = """""")

module.add_function ("triqs_tprf::chi_wnr_t triqs_tprf::chi0r_from_gr_PH (int nw, int nn, triqs_tprf::g_wr_cvt g_nr)", doc = """Generalized susceptibility bubble in the particle-hole channel :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})`.\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r}) =\n     - \\beta G_{d\\bar{a}}(\\nu, \\mathbf{r}) \\cdot G_{b\\bar{c}}(\\nu + \\omega, -\\mathbf{r})\n\n  :param nw: Number of bosonic Matsubara freqiencies.\n  :param nn: Number of fermionic Matsubara freqiencies.\n  :param g_tr: Imaginary time Green\'s function in real-space, :math:`G_{a\\bar{b}}(\\nu, \\mathbf{r})`.\n  @return Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wnr_t triqs_tprf::chi0r_from_gr_PH_nompi (int nw, int nn, triqs_tprf::g_wr_cvt g_nr)", doc = """Generalized susceptibility bubble in the particle-hole channel :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})` without MPI parallellization.\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r}) =\n     - \\beta G_{d\\bar{a}}(\\nu, \\mathbf{r}) \\cdot G_{b\\bar{c}}(\\nu + \\omega, -\\mathbf{r})\n\n  :param nw: Number of bosonic Matsubara freqiencies.\n  :param nn: Number of fermionic Matsubara freqiencies.\n  :param g_tr: Imaginary time Green\'s function in real-space, :math:`G_{a\\bar{b}}(\\nu, \\mathbf{r})`.\n  @return Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wnk_t triqs_tprf::chi0q_from_g_wk_PH (int nw, int nn, triqs_tprf::g_wk_cvt g_wk)", doc = """Generalized susceptibility bubble in the particle-hole channel :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{q})` with convolution in k-space.\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{q}) =\n     - \\frac{\\beta}{N_k} \\sum_\\mathbf{k} \n     G_{d\\bar{a}}(\\nu, \\mathbf{k}) \\cdot G_{b\\bar{c}}(\\nu + \\omega, \\mathbf{k} - \\mathbf{q})\n\n  :param nw: Number of bosonic Matsubara freqiencies.\n  :param nn: Number of fermionic Matsubara freqiencies.\n  :param g_tr: Imaginary time Green\'s function in real-space, :math:`G_{a\\bar{b}}(\\nu, \\mathbf{r})`.\n  @return Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.""")

module.add_function ("triqs_tprf::chi_wnr_t triqs_tprf::chi0r_from_chi0q (triqs_tprf::chi_wnk_cvt chi_wnk)", doc = """Fourier transform of the generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{q})` in momentum-space to :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})` in real-space.\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r}) =\n     \\mathcal{F}_{\\mathbf{q} \\rightarrow \\mathbf{r}} \\left\\{\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{q})\n     \\right\\}\n\n  :param chi_wnk: Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.\n  @return Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.""")

module.add_function ("triqs_tprf::chi_wnk_t triqs_tprf::chi0q_from_chi0r (triqs_tprf::chi_wnr_cvt chi_wnr)", doc = """Fourier transform of the generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})` in real space to :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{q})` in momentum space.\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{q}) =\n     \\mathcal{F}_{\\mathbf{r} \\rightarrow \\mathbf{q}} \\left\\{\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})\n     \\right\\}\n\n  :param chi_wnr: Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.\n  @return Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::chi0q_sum_nu (triqs_tprf::chi_wnk_cvt chi_wnk)", doc = """Sum over fermionic frequency in the generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{k})`. (NB! without tail corrections)\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k}) =\n     \\frac{1}{\\beta^2} \\sum_{\\nu=\\nu_{min}}^\\nu_{max} \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{k})\n\n  :param chi_wnk: Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.\n  @return Susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})` in one bosonic Matsubara frequency and momentum space.""")

module.add_function ("triqs_tprf::chi_wk_t triqs_tprf::chi0q_sum_nu_tail_corr_PH (triqs_tprf::chi_wnk_cvt chi_wnk)", doc = """Sum over fermionic frequency in the generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{k})` using higher order tail corrections when summing to infinity.\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k}) =\n     \\frac{1}{\\beta^2} \\sum_{\\nu=-\\infty}^\\infty \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{k})\n\n  :param chi_wnk: Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.\n  @return Susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})` in one bosonic Matsubara frequency and momentum space.""")

module.add_function ("triqs_tprf::chi_w_t triqs_tprf::chi0q_sum_nu_q (triqs_tprf::chi_wnk_cvt chi_wnk)", doc = """Sum over fermionic frequency and momentum in the generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{k})`. (NB! without tail corrections)\n\n  Computes\n\n  .. math::\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k}) =\n     \\frac{1}{N_k} \\sum_\\matbf{k} \\frac{1}{\\beta^2} \\sum_{\\nu=\\nu_{min}}^\\nu_{max}\n     \\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\mathbf{k})\n\n  :param chi_wnk: Generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.\n  @return Susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega)` in one bosonic Matsubara frequency.""")

module.add_function ("triqs_tprf::chi_kwnn_t triqs_tprf::chiq_from_chi0q_and_gamma_PH (triqs_tprf::chi_wnk_cvt chi0_wnk, triqs_tprf::chi_wnn_cvt gamma_ph_wnn)", doc = """Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\nu\', \\mathbf{k})`.\n\n  Computes\n\n  .. math::\n     \\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\nu\', \\mathbf{k}) =\n     \\chi^{(0)} \\left[ 1 - \\Gamma^{(PH)} \\chi^{(0)} \\right]^{-1}\n\n  :param chi0_wnk: Generalized lattice bubble susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})`.\n  :param gamma_ph_wnn: Local particle-hole vertex function :math:`\\Gamma^{(PH)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\nu\')`.\n  @return Generalized lattice susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\nu\', \\mathbf{k})`.""")

module.add_function ("triqs_tprf::chi_kw_t triqs_tprf::chiq_sum_nu_from_chi0q_and_gamma_PH (triqs_tprf::chi_wnk_cvt chi0_wnk, triqs_tprf::chi_wnn_cvt gamma_ph_wnn)", doc = """Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})`.\n\n  Computes\n\n  .. math::\n     \\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k}) =\n     \\chi^{(0)} \\left[ 1 - \\Gamma^{(PH)} \\chi^{(0)} \\right]^{-1}\n\n  :param chi0_wnk: Generalized lattice bubble susceptibility :math:`\\chi^{(0)}_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})`.\n  :param gamma_ph_wnn: Local particle-hole vertex function :math:`\\Gamma^{(PH)}_{\\bar{a}b\\bar{c}d}(\\omega, \\nu, \\nu\')`.\n  @return Generalized lattice susceptibility :math:`\\chi_{\\bar{a}b\\bar{c}d}(\\omega, \\mathbf{k})`.""")

module.add_function ("gf<cartesian_product<triqs::lattice::brillouin_zone,triqs::gfs::imfreq>,tensor_valued<4>> triqs_tprf::chiq_sum_nu_from_g_wk_and_gamma_PH (triqs_tprf::gk_iw_t g_wk, triqs_tprf::g2_iw_vt gamma_ph_wnn, int tail_corr_nwf = -1)", doc = """""")

module.add_function ("gf<cartesian_product<triqs::lattice::brillouin_zone,triqs::gfs::imfreq>,tensor_valued<4>> triqs_tprf::chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH (double mu, triqs_tprf::ek_vt e_k, triqs_tprf::g_iw_vt sigma_w, triqs_tprf::g2_iw_vt gamma_ph_wnn, int tail_corr_nwf = -1)", doc = """""")

module.add_function ("gf<cartesian_product<triqs::lattice::brillouin_zone,triqs::gfs::imfreq>,tensor_valued<4>> triqs_tprf::chiq_sum_nu (triqs_tprf::chiq_t chiq)", doc = """""")

module.add_function ("gf<triqs::gfs::imfreq,tensor_valued<4>> triqs_tprf::chiq_sum_nu_q (triqs_tprf::chiq_t chiq)", doc = """""")



module.generate_code()