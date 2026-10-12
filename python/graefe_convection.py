"""Free-convection boundary-layer conductance of inclined elliptical leaves.

Implementation of Graefe, Grosch & Bitterlich (2022): "The boundary layer
conductance of inclined elliptical leaves under free convection", Agricultural
and Forest Meteorology 317, 108884, https://doi.org/10.1016/j.agrformet.2022.108884
(CC BY 4.0).  Equation numbers below refer to that paper.

Model summary
-------------
* geometry (Section 2.2): area S = pi/4 L W, dissection index DI (Eq. 8),
  characteristic length of horizontal leaves Lh = S/P (Eq. 7), mean chord length
  LCo(phi) (Eq. 9 / A5) along the steepest-ascent direction (angle phi to the short
  axis, Eq. 10/11, Appendix B) and Li = (0.04511 n + 1.006) LCo, n = 3 * 0.25;
* flux-proportional Nusselt numbers NuF (m-1) for the upper (U) and lower (L)
  leaf side (Eq. 3a-3g, Table 3 parameters p1..p6); for a cooled leaf
  (dT < 0) the right-hand sides of U and L are exchanged;
* boundary-layer conductance gb = kappa * NuF (Eq. 1) and sensible heat flux
  H = gb rho cp (TL - TA) per side (Eq. 14);
* optional buoyancy correction of the leaf-air temperature difference for vapour
  pressure differences (Eq. 18);
* extension (not part of the paper): laminar forced convection with the
  streamline ("flowline") length along the wind direction projected onto the leaf
  plane, obtained from the exact solution of Eq. 4 with n = 0.5, and mixed
  convection by cubic blending of the free and forced Nusselt numbers.

All functions are vectorised over numpy arrays of temperatures; the geometry is
scalar.  Temperatures in Kelvin, lengths in metres, angles in degrees.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

__version__ = "1.1.0"

# Table 3 (Graefe et al. 2022)
P1, P2, P3, P4, P5, P6 = 6.0, 1.42, 0.513, 0.4, 0.202, 0.492
PR = 0.71          # Prandtl number of air (Table 1)
G = 9.81           # m s-2
CP_AIR = 1004.0    # J kg-1 K-1
R_DRY = 287.05     # J kg-1 K-1
N_RA = 0.25        # Rayleigh exponent of Eq. 3b/3c
N_LI = 3.0 * N_RA  # exponent n in Eq. 4 / Eq. 9


# ----------------------------------------------------------------- air
def air_properties(T, p=101325.0):
    """Dry-air density rho (kg m-3), kinematic viscosity nu (m2 s-1), thermal
    diffusivity kappa (m2 s-1) and conductivity k (W m-1 K-1) at temperature T (K)."""
    T = np.asarray(T, dtype=float)
    rho = np.asarray(p, dtype=float) / (R_DRY * T)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)          # Sutherland
    k = 0.02624 * (T / 300.0) ** 0.8646
    nu = mu / rho
    kappa = k / (rho * CP_AIR)
    return rho, nu, kappa, k


# ------------------------------------------------------------ geometry
@dataclass(frozen=True)
class LeafGeometry:
    L: float
    W: float
    alpha_x_deg: float
    alpha_y_deg: float
    S: float          # one-sided area, m2
    AR: float         # aspect ratio L/W
    DI: float         # dissection index (Eq. 8)
    P: float          # perimeter implied by DI, m
    Lh: float         # characteristic length, horizontal leaf (Eq. 7)
    theta: float      # leaf normal angle to the vertical, rad (Eq. 11)
    phi: float        # streamline angle to the short axis, rad (Eq. 10)
    LCo: float        # mean chord length along the streamline, m (Eq. 9)
    Li: float         # characteristic length, inclined leaf (Eq. 9)


def dissection_index(L, W):
    """Eq. 8 (max. error < 0.04 % for AR 1..4)."""
    return 1.112 * ((L + W) / math.sqrt(L * W)) ** 0.783 - 0.913


def mean_chord_length(L, W, phi):
    """LCo(phi), Eq. 9 / A5, written with sin/cos to stay finite at phi = pi/2."""
    s, c = math.sin(phi), math.cos(phi)
    return math.pi * L * W / (4.0 * math.sqrt((s * W) ** 2 + (c * L) ** 2))


def streamline_angle(alpha_x, alpha_y):
    """Angle phi between the steepest-ascent vector s and the rotated short axis w
    (Appendix B).  The leaf (long axis along y, short axis along x) is rotated first
    about y by alpha_y, then about x by alpha_x (radians).  Returns phi in [0, pi/2];
    for a horizontal leaf the direction is undefined and pi/2 is returned."""
    cy, sy = math.cos(alpha_y), math.sin(alpha_y)
    cx, sx = math.cos(alpha_x), math.sin(alpha_x)
    Ry = np.array([[cy, 0.0, -sy], [0.0, 1.0, 0.0], [sy, 0.0, cy]])     # Eq. B2
    Rx = np.array([[1.0, 0.0, 0.0], [0.0, cx, -sx], [0.0, sx, cx]])
    R = Rx @ Ry
    n = R @ np.array([0.0, 0.0, 1.0])                                   # leaf normal (Eq. B1)
    w = R @ np.array([1.0, 0.0, 0.0])                                   # short axis
    z = np.array([0.0, 0.0, 1.0])
    s = z - np.dot(z, n) * n                                            # projection of the vertical onto the plane
    norm = np.linalg.norm(s)
    if norm < 1e-12:
        return 0.5 * math.pi
    return math.acos(min(1.0, abs(float(np.dot(s, w)) / norm)))


def characteristic_length(L, W, phi, n):
    """Exact flux-weighted characteristic length of Eq. 4 for an ellipse.

    With chords parallel to the streamline direction (angle phi to the short axis)
    the chord length varies as D(t) = D0 sqrt(1 - (t/t_max)^2), D0 = (4/pi) LCo(phi),
    so that Eq. 4 / A6 has the closed form

        Lc^(n-1) = D0^(n-1) * (2/pi) * B(1/2, n/2 + 1)

    with the Beta function B.  n = 3*0.25 reproduces Eq. 9 (free convection) to
    within 0.2 %; n = 0.5 gives the length for laminar forced convection
    (Nu ~ Re^0.5, i.e. h ~ L^-0.5).
    """
    if abs(n - 1.0) < 1e-12:
        return mean_chord_length(L, W, phi)
    D0 = 4.0 / math.pi * mean_chord_length(L, W, phi)
    beta = math.gamma(0.5) * math.gamma(0.5 * n + 1.0) / math.gamma(0.5 * n + 1.5)
    return D0 * (2.0 / math.pi * beta) ** (1.0 / (n - 1.0))


def _leaf_frame_vectors(alpha_x, alpha_y):
    """Leaf normal n and short-axis vector w after the rotations Rx(alpha_x) Ry(alpha_y)
    (Appendix B), expressed in the unrotated leaf frame (x = short axis, y = long
    axis, z = up)."""
    cy, sy = math.cos(alpha_y), math.sin(alpha_y)
    cx, sx = math.cos(alpha_x), math.sin(alpha_x)
    Ry = np.array([[cy, 0.0, -sy], [0.0, 1.0, 0.0], [sy, 0.0, cy]])
    Rx = np.array([[1.0, 0.0, 0.0], [0.0, cx, -sx], [0.0, sx, cx]])
    R = Rx @ Ry
    return R @ np.array([0.0, 0.0, 1.0]), R @ np.array([1.0, 0.0, 0.0])


def wind_streamline_angle(alpha_x_deg, alpha_y_deg, leaf_azimuth_deg, wind_from_deg):
    """Angle phi_w (rad, 0..pi/2) between the wind streamlines over the leaf and the
    short leaf axis, plus the fraction of the wind speed lying in the leaf plane.

    The horizontal wind vector (blowing towards azimuth wind_from + 180) is
    projected onto the leaf plane.  ``leaf_azimuth_deg`` is the compass direction
    the upper leaf surface faces (180 = south), i.e. the horizontal direction of
    the leaf normal; the long axis then runs towards leaf_azimuth + 180 (up-slope).
    Returns (phi_w, in_plane_fraction)."""
    n, w = _leaf_frame_vectors(math.radians(alpha_x_deg), math.radians(alpha_y_deg))
    A = math.radians(leaf_azimuth_deg)
    beta = math.radians(wind_from_deg + 180.0)                  # direction the wind blows to
    d_world = np.array([math.sin(beta), math.cos(beta)])        # (east, north)
    y_leaf = np.array([-math.sin(A), -math.cos(A)])             # long axis (horizontal, before tilting)
    x_leaf = np.array([-math.cos(A), math.sin(A)])              # short axis = y_leaf x z
    d = np.array([float(d_world @ x_leaf), float(d_world @ y_leaf), 0.0])
    d_p = d - float(d @ n) * n
    norm = float(np.linalg.norm(d_p))
    if norm < 1e-9:
        return 0.5 * math.pi, 0.0
    return math.acos(min(1.0, abs(float(d_p @ w)) / norm)), norm


def leaf_geometry(L, W, alpha_x_deg, alpha_y_deg=0.0) -> LeafGeometry:
    """Geometric descriptors of an elliptical leaf with length L >= width W (m),
    rotated by alpha_x (inclination from horizontal about the short axis) and
    alpha_y (pre-rotation about the long axis), degrees."""
    if L < W:
        raise ValueError("L must be the longer axis (L >= W)")
    ax, ay = math.radians(alpha_x_deg), math.radians(alpha_y_deg)
    S = 0.25 * math.pi * L * W
    DI = dissection_index(L, W)
    P = 2.0 * DI * math.sqrt(math.pi * S)
    Lh = S / P                                                          # = 0.5 sqrt(S) / (DI sqrt(pi))
    theta = math.acos(max(-1.0, min(1.0, math.cos(ax) * math.cos(ay))))  # Eq. 11
    phi = streamline_angle(ax, ay)
    LCo = mean_chord_length(L, W, phi)
    Li = (0.04511 * N_LI + 1.006) * LCo
    return LeafGeometry(L, W, alpha_x_deg, alpha_y_deg, S, L / W, DI, P, Lh, theta, phi, LCo, Li)


# ------------------------------------------------------------ Nusselt
def conduction_limit(geom: LeafGeometry) -> float:
    """Nuc, Eq. 3f (m-1)."""
    return 5.3 * math.sin(0.0338 * geom.AR + 0.607) / math.sqrt(geom.S)


def rayleigh(dT, TA, nu):
    """Modified Rayleigh number without length scale, Eq. 3g (m-3)."""
    return G * np.abs(dT) * PR / (np.asarray(TA, dtype=float) * np.asarray(nu, dtype=float) ** 2)


def virtual_temperature_difference(TL, TA, eL, eA, p=101325.0):
    """Eq. 18: leaf-air temperature difference including the vapour density effect."""
    TL = np.asarray(TL, dtype=float); TA = np.asarray(TA, dtype=float)
    return (TL - TA) + 0.38 * (np.asarray(eL, dtype=float) * TL - np.asarray(eA, dtype=float) * TA) / p


def nusselt_free(dT, TA, nu, geom: LeafGeometry):
    """Flux-proportional Nusselt numbers (NuF_U, NuF_L) in m-1, Eq. 3a-3g.

    dT: (virtual) leaf minus air temperature difference (K), negative for a cooled
    leaf; TA: air temperature (K); nu: kinematic viscosity at film temperature.
    Also returns Nuc and Ra."""
    dT = np.asarray(dT, dtype=float)
    Ra = rayleigh(dT, TA, nu)
    st, ct = math.sin(geom.theta), math.cos(geom.theta)
    ra25 = Ra ** 0.25
    nu_i_u = P3 * geom.Li ** (-0.25) * st ** P4 * ra25            # Eq. 3b
    nu_i_l = P3 * geom.Li ** (-0.25) * st ** 0.25 * ra25          # Eq. 3c
    nu_h_u = P5 * geom.Lh ** (-0.25) * ct ** 0.25 * ra25          # Eq. 3d
    nu_h_l = P6 * geom.Lh ** (-0.4) * ct ** 0.2 * Ra ** 0.2       # Eq. 3e
    nuc = conduction_limit(geom)
    cooled = dT < 0.0
    iu = np.where(cooled, nu_i_l, nu_i_u); il = np.where(cooled, nu_i_u, nu_i_l)
    hu = np.where(cooled, nu_h_l, nu_h_u); hl = np.where(cooled, nu_h_u, nu_h_l)

    def blend(ni, nh):                                               # Eq. 3a
        return ((ni ** P1 + nh ** P1) ** (P2 / P1) + nuc ** P2) ** (1.0 / P2)

    return blend(iu, hu), blend(il, hl), nuc, Ra


C_FREE_TURB = 0.14      # turbulent free convection Nu_L = 0.14 Ra_L^(1/3) (upper side / vertical, Ra > ~1e7)
M_FREE_TURB = 6.0       # blending exponent laminar/turbulent free convection


def nusselt_free_turbulent(dT, TA, nu, geom: LeafGeometry, c=C_FREE_TURB):
    """Turbulent free-convection flux-proportional Nusselt number (m-1) for the side
    facing the buoyant flow (upper side of a heated or lower side of a cooled leaf,
    and vertical leaves): NuF = c Ra^(1/3), independent of the leaf size.  The
    stably stratified side (lower side of a heated leaf) keeps its laminar value.
    Weighted by sin(theta) + cos(theta) contributions as in Eq. 3 (c applies to both)."""
    Ra = rayleigh(dT, TA, nu)
    return c * Ra ** (1.0 / 3.0)


def free_convection_transition(TL, TA, L, W, alpha_x_deg, alpha_y_deg=0.0, p=101325.0, eL=None, eA=None,
                               c_turb=C_FREE_TURB, m=M_FREE_TURB):
    """free_convection() with a laminar-turbulent transition of the free convection:
    on the buoyantly unstable side NuF = (NuF_lam^m + NuF_turb^m)^(1/m) with
    NuF_turb = c_turb Ra^(1/3); the stable side stays laminar.  Extension, not part of
    Graefe et al. (2022)."""
    r = free_convection(TL, TA, L, W, alpha_x_deg, alpha_y_deg, p, eL, eA)
    g = r["geometry"]
    turb = nusselt_free_turbulent(r["dT"], TA, r["nu"], g, c_turb)
    heated = np.asarray(r["dT"]) >= 0
    vertical = abs(math.sin(g.theta) - 1.0) < 1e-9
    unstable_u = heated | vertical           # upper side buoyantly unstable (heated) or vertical
    unstable_l = (~heated) | vertical
    nf_u = np.where(unstable_u, (r["NuF_U"] ** m + turb ** m) ** (1.0 / m), r["NuF_U"])
    nf_l = np.where(unstable_l, (r["NuF_L"] ** m + turb ** m) ** (1.0 / m), r["NuF_L"])
    r["NuF_lam_U"], r["NuF_lam_L"] = r["NuF_U"], r["NuF_L"]
    r["NuF_turb"] = turb
    # turbulent share of the free convection per side (power-mean partition)
    r["free_turb_share_U"] = np.where(unstable_u, turb ** m / (r["NuF_U"] ** m + turb ** m), 0.0)
    r["free_turb_share_L"] = np.where(unstable_l, turb ** m / (r["NuF_L"] ** m + turb ** m), 0.0)
    r["NuF_U"], r["NuF_L"] = nf_u, nf_l
    r["gb_U"], r["gb_L"] = r["kappa"] * nf_u, r["kappa"] * nf_l
    TLa = np.asarray(TL, dtype=float); TAa = np.asarray(TA, dtype=float)
    r["H_U"] = r["rho"] * CP_AIR * r["gb_U"] * (TLa - TAa)
    r["H_L"] = r["rho"] * CP_AIR * r["gb_L"] * (TLa - TAa)
    return r


def free_convection(TL, TA, L, W, alpha_x_deg, alpha_y_deg=0.0, p=101325.0, eL=None, eA=None):
    """Boundary-layer conductances and sensible heat fluxes of both leaf sides.

    Parameters
    ----------
    TL, TA : leaf and air temperature (K), scalars or broadcastable arrays
    L, W : leaf length and width (m); alpha_x_deg, alpha_y_deg : rotations (deg)
    p : air pressure (Pa); eL, eA : leaf-surface and air vapour pressure (Pa),
        both required to apply the buoyancy correction of Eq. 18.

    Returns a dict with gb_U, gb_L (m s-1), H_U, H_L (W m-2 of one-sided leaf
    area, positive from leaf to air), NuF_U, NuF_L, Nuc (m-1), Ra (m-3), dT (K),
    TF, rho, nu, kappa and the geometry.
    """
    TL = np.asarray(TL, dtype=float); TA = np.asarray(TA, dtype=float)
    geom = leaf_geometry(L, W, alpha_x_deg, alpha_y_deg)
    TF = 0.5 * (TL + TA)
    rho, nu, kappa, _ = air_properties(TF, p)
    if eL is not None and eA is not None:
        dT = virtual_temperature_difference(TL, TA, eL, eA, p)
    else:
        dT = TL - TA
    dT = np.where(np.abs(dT) < 1e-9, 1e-9, dT)
    nf_u, nf_l, nuc, Ra = nusselt_free(dT, TA, nu, geom)
    gb_u, gb_l = kappa * nf_u, kappa * nf_l                          # Eq. 1
    H_u = rho * CP_AIR * gb_u * (TL - TA)                            # Eq. 14
    H_l = rho * CP_AIR * gb_l * (TL - TA)
    return {"gb_U": gb_u, "gb_L": gb_l, "H_U": H_u, "H_L": H_l, "NuF_U": nf_u, "NuF_L": nf_l,
            "Nuc": nuc, "Ra": Ra, "dT": dT, "TF": TF, "rho": rho, "nu": nu, "kappa": kappa, "geometry": geom}


# --------------------------------------------------- forced / mixed convection
N_FORCED = 0.5          # Nu ~ Re^0.5  ->  h ~ L^-0.5, exponent n of Eq. 4 for the forced streamline length
C_FORCED = 0.664        # laminar flat plate (Pohlhausen); Campbell & Norman (1998) use 0.6
C_TURB_LOCAL = 0.0296   # local turbulent flat plate Nu_x = 0.0296 Re_x^0.8 Pr^0.6 (Lienhard 2020)
RE_L_SMOOTH = 5.0e5     # transition onset for a disturbance-free stream
RE_L_MIN = 2.0e3        # lower bound of the onset Reynolds number at very high turbulence


def transition_onset_reynolds(TI_pct):
    """Onset of laminar-turbulent transition Re_l as a function of the free-stream
    turbulence intensity TI (%): Mayle (1991) Re_theta,t = 400 TI^-5/8 converted with the
    Blasius momentum thickness Re_theta = 0.664 Re_x^0.5, i.e. Re_l = 3.63e5 TI^-1.25,
    bounded by RE_L_SMOOTH (TI -> 0) and RE_L_MIN."""
    TI = np.maximum(np.asarray(TI_pct, dtype=float), 1e-3)
    return np.clip(3.63e5 * TI ** (-1.25), RE_L_MIN, RE_L_SMOOTH)


def forced_nusselt_parts(Re, Pr=PR, TI_pct=1.0, bc="isothermal", unheated_fraction=0.0, Re_l=None):
    """The three additive contributions (laminar, transitional, turbulent) of the
    Lienhard (2020) average Nusselt number; their sum is forced_nusselt_lienhard()."""
    Re = np.asarray(Re, dtype=float)
    a_l = 0.332
    hfa = 0.680 / 0.664 if bc == "uniform_flux" else 1.0
    f_sl = (1.0 - unheated_fraction ** 0.75) ** (-1.0 / 3.0) if unheated_fraction > 0 else 1.0
    if Re_l is None:
        Re_l = transition_onset_reynolds(TI_pct)
    Re_l = np.asarray(Re_l, dtype=float)
    c = 0.9922 * np.log10(Re_l) - 3.013
    Re_u = ((Re_l ** (0.5 - c) * a_l / C_TURB_LOCAL) * Pr ** (1.0 / 3.0 - 0.6)) ** (1.0 / (0.8 - c))
    Re_u = np.maximum(Re_u, Re_l * 1.0001)
    lam = 2.0 * a_l * f_sl * hfa * Pr ** (1.0 / 3.0) * np.minimum(Re, Re_l) ** 0.5
    nu_lam_local_l = a_l * f_sl * hfa * Re_l ** 0.5 * Pr ** (1.0 / 3.0)
    re_mid = np.clip(Re, Re_l, Re_u)
    trans = np.where(Re > Re_l, nu_lam_local_l / c * ((re_mid / Re_l) ** c - 1.0), 0.0)
    turb = np.where(Re > Re_u, 0.037 * Pr ** 0.6 * np.maximum(Re ** 0.8 - Re_u ** 0.8, 0.0), 0.0)
    return lam, trans, turb


def forced_nusselt_lienhard(Re, Pr=PR, TI_pct=1.0, bc="isothermal", unheated_fraction=0.0, Re_l=None):
    """Average Nusselt number of one side of a flat plate including the laminar-turbulent
    transition after Lienhard (2020, J. Heat Transfer 142, 061805).

    Local laws: laminar Nu_x = 0.332 Re_x^0.5 Pr^(1/3) (uniform heat flux: average
    0.680 instead of 0.664 with the mean temperature difference), transitional Nu_x = Nu_x(Re_l) (Re_x/Re_l)^c with
    c = 0.9922 log10(Re_l) - 3.013, turbulent Nu_x = 0.0296 Re_x^0.8 Pr^0.6.  The end
    of transition Re_u follows from continuity of the local Nusselt number; the
    average over the plate is the exact integral of the three pieces.  An unheated
    starting length multiplies the laminar part by (1 - (x0/L)^0.75)^(-1/3).

    Returns the average Nu based on the plate length (dimensionless, not flux-proportional)."""
    Re = np.asarray(Re, dtype=float)
    a_l = 0.332
    # uniform heat flux: average Nu based on the mean temperature difference is 0.680 Re^0.5 Pr^(1/3)
    # (Lienhard 2020), i.e. 2.4 % above the isothermal plate; the local law keeps a_l = 0.332
    hfa = 0.680 / 0.664 if bc == "uniform_flux" else 1.0
    f_sl = (1.0 - unheated_fraction ** 0.75) ** (-1.0 / 3.0) if unheated_fraction > 0 else 1.0
    if Re_l is None:
        Re_l = transition_onset_reynolds(TI_pct)
    Re_l = np.asarray(Re_l, dtype=float)
    c = 0.9922 * np.log10(Re_l) - 3.013
    # Re_u: a_l Re_l^0.5 Pr^(1/3) (Re_u/Re_l)^c = 0.0296 Re_u^0.8 Pr^0.6
    Re_u = ((Re_l ** (0.5 - c) * a_l / C_TURB_LOCAL) * Pr ** (1.0 / 3.0 - 0.6)) ** (1.0 / (0.8 - c))
    Re_u = np.maximum(Re_u, Re_l * 1.0001)
    lam_avg = 2.0 * a_l * f_sl * hfa * Pr ** (1.0 / 3.0)    # 0.664 (isothermal) or 0.680 (uniform flux)
    nu_lam_local_l = a_l * f_sl * hfa * Re_l ** 0.5 * Pr ** (1.0 / 3.0)   # local Nu at onset
    re_mid = np.clip(Re, Re_l, Re_u)
    trans = nu_lam_local_l / c * ((re_mid / Re_l) ** c - 1.0)
    turb = 0.037 * Pr ** 0.6 * np.maximum(Re ** 0.8 - Re_u ** 0.8, 0.0)
    lam = lam_avg * np.minimum(Re, Re_l) ** 0.5
    nu = lam + np.where(Re > Re_l, trans, 0.0) + np.where(Re > Re_u, turb, 0.0)
    return nu


def forced_nusselt_turbulent(Re, Pr=PR):
    """Average Nusselt number of a plate that is turbulent from the leading edge
    (0.037 Re^0.8 Pr^0.6): blunt leading edges with Re_t = u t / nu >= 2000 (Kottke,
    Blenke & Schmidt 1977) trip the boundary layer immediately."""
    return 0.037 * np.asarray(Re, dtype=float) ** 0.8 * Pr ** 0.6


RE_THICKNESS_BLUNT = 2000.0


def forced_regime(Re, TI_pct=1.0, Re_l=None):
    """'laminar' | 'transitional' | 'turbulent' classification used in the reports."""
    Re = np.asarray(Re, dtype=float)
    Re_l = transition_onset_reynolds(TI_pct) if Re_l is None else np.asarray(Re_l, float)
    c = 0.9922 * np.log10(Re_l) - 3.013
    Re_u = ((Re_l ** (0.5 - c) * 0.332 / C_TURB_LOCAL) * PR ** (1.0 / 3.0 - 0.6)) ** (1.0 / (0.8 - c))
    return np.where(Re <= Re_l, "laminar", np.where(Re <= Re_u, "transitional", "turbulent"))


def forced_length(L, W, alpha_x_deg, alpha_y_deg, leaf_azimuth_deg, wind_from_deg):
    """Characteristic length for laminar forced convection along the wind streamlines
    (flowline concept of Section 2.2 applied with n = 0.5).  Returns (Lf, phi_w,
    in_plane_fraction).  For wind normal to the leaf (no in-plane component) the
    orientation-averaged chord length pi S / P (Eq. 5) is used."""
    phi_w, frac = wind_streamline_angle(alpha_x_deg, alpha_y_deg, leaf_azimuth_deg, wind_from_deg)
    if frac < 1e-6:
        g = leaf_geometry(L, W, alpha_x_deg, alpha_y_deg)
        D0 = 4.0 * g.S / g.P                                     # (4/pi) * LCo with LCo = pi S / P (Eq. 5)
        beta = math.gamma(0.5) * math.gamma(0.5 * N_FORCED + 1.0) / math.gamma(0.5 * N_FORCED + 1.5)
        return D0 * (2.0 / math.pi * beta) ** (1.0 / (N_FORCED - 1.0)), phi_w, frac
    return characteristic_length(L, W, phi_w, N_FORCED), phi_w, frac


# ------------------------------------------------ angle of attack (parsimonious leaf model)
# Fitted 2026-10-12 to 815 thin-plate/leaf measurements (docs/validation.md): Chen 1988a/b,
# Haseba 1973/1984, Sparrow 1977/1979, Parkhurst 1968, Powell 1940, Schuepp 1993,
# Ramirez 2002, Kobus 2006; rRMSE 24 %, bias -4 %.
ANGLE_PARAMS = {"f_e": 1.161, "a_W": 0.226, "a_R": 0.574, "beta_s_deg": 31.2, "c_w": 0.833}


def angle_of_attack_factors(beta_deg, params=ANGLE_PARAMS):
    """Multipliers of the parallel-flow Nusselt number for the windward (gW) and the
    leeward (gB) side of a thin plate at the angle of attack beta (0 = parallel,
    90 = normal to the flow):

        gW = 1 + a_W sin(beta)                                  accelerated outer flow
        gB = w (1 + a_R sin(2 beta)) + (1 - w) c_w,  w = exp(-(beta/beta_s)^2)
             reattaching separation bubble (enhancement) -> wake regime (c_w < 1)

    Returns (gW, gB); both are 1 for parallel flow."""
    b = np.radians(np.asarray(beta_deg, dtype=float))
    gW = 1.0 + params["a_W"] * np.sin(b)
    w = np.exp(-(np.degrees(b) / params["beta_s_deg"]) ** 2)
    gB = w * (1.0 + params["a_R"] * np.sin(2.0 * b)) + (1.0 - w) * params["c_w"]
    return gW, gB


def forced_nusselt_leaf(Re, beta_deg=0.0, Pr=PR, TI_pct=1.0, bc="isothermal", transition=True, params=ANGLE_PARAMS):
    """Average Nusselt numbers (dimensionless, on the streamline length) of the windward
    and leeward side of a thin plate or leaf: f_e * Nu_plate(Re) * (gW, gB), where
    Nu_plate is the Lienhard (2020) laminar-transitional-turbulent average (or the
    laminar law when transition=False) and f_e the finite-plate enhancement."""
    Re = np.asarray(Re, dtype=float)
    base = forced_nusselt_lienhard(Re, Pr, TI_pct, bc) if transition else 0.664 * (0.680 / 0.664 if bc == "uniform_flux" else 1.0) * Re ** 0.5 * Pr ** (1.0 / 3.0)
    gW, gB = angle_of_attack_factors(beta_deg, params)
    return params["f_e"] * base * gW, params["f_e"] * base * gB


def wind_attack_angle(alpha_x_deg, alpha_y_deg, leaf_azimuth_deg, wind_from_deg):
    """Angle of attack beta (deg, 0 = parallel) of the horizontal wind on the leaf plane
    and the side it hits first ('U' upper or 'L' lower).  Same frame as
    wind_streamline_angle()."""
    n, _ = _leaf_frame_vectors(math.radians(alpha_x_deg), math.radians(alpha_y_deg))
    A = math.radians(leaf_azimuth_deg)
    beta = math.radians(wind_from_deg + 180.0)
    d_world = np.array([math.sin(beta), math.cos(beta)])
    y_leaf = np.array([-math.sin(A), -math.cos(A)])
    x_leaf = np.array([-math.cos(A), math.sin(A)])
    d = np.array([float(d_world @ x_leaf), float(d_world @ y_leaf), 0.0])
    dn = float(d @ n)
    return math.degrees(math.asin(min(1.0, abs(dn)))), ("U" if dn < 0 else "L")


def nusselt_forced(u, nu, Lf, c=C_FORCED, transition=False, TI_pct=1.0, bc="isothermal"):
    """Forced convection per leaf side, flux-proportional (m-1).

    transition=False: laminar flat plate NuF = c Re^0.5 Pr^(1/3) / Lf, Re = u Lf / nu.
    transition=True : Lienhard (2020) laminar-transitional-turbulent average with the
    onset Reynolds number from the turbulence intensity TI (%), divided by Lf."""
    u = np.asarray(u, dtype=float)
    Re = np.maximum(u * Lf / np.asarray(nu, dtype=float), 0.0)
    if transition:
        return forced_nusselt_lienhard(Re, PR, TI_pct, bc) / Lf
    return c * np.sqrt(Re) * PR ** (1.0 / 3.0) / Lf


def mixed_convection(TL, TA, u, L, W, alpha_x_deg, alpha_y_deg=0.0, leaf_azimuth_deg=180.0, wind_from_deg=None,
                     p=101325.0, eL=None, eA=None, c_forced=C_FORCED, transition=False, TI_pct=1.0, bc="isothermal",
                     angle_model=False, free_transition=False):
    """Mixed free/forced convection of both leaf sides.

    Free convection follows Graefe et al. (2022); forced convection is the laminar
    flat-plate relation with the streamline length Lf along the wind direction
    projected onto the leaf plane (wind_from_deg = meteorological wind direction;
    None = orientation-averaged chord length).  Both are combined per side with the
    cubic blending NuF = (NuF_free^3 + NuF_forced^3)^(1/3) (no distinction of
    assisting/opposing flows).  With ``transition=True`` the forced part includes the
    laminar-turbulent transition (Lienhard 2020) with the onset Reynolds number from
    the turbulence intensity ``TI_pct``.  With ``angle_model=True`` the forced part is
    multiplied by the finite-plate enhancement f_e and the windward/leeward angle-of-
    attack factors (angle_of_attack_factors, fitted to the thin-plate literature data);
    the angle of attack and the windward side follow from the wind direction (beta = 0
    when the direction is unknown).  Returns the dict of free_convection() extended by
    NuF_forced, gb_forced, Lf, phi_w (rad), NuF_mix_U/L, gb_U/gb_L (mixed), H_U/H_L
    (mixed), forced_share and gb_free_U/L."""
    if free_transition:
        free = free_convection_transition(TL, TA, L, W, alpha_x_deg, alpha_y_deg, p, eL, eA)
    else:
        free = free_convection(TL, TA, L, W, alpha_x_deg, alpha_y_deg, p, eL, eA)
        free["NuF_lam_U"], free["NuF_lam_L"] = free["NuF_U"], free["NuF_L"]
        free["free_turb_share_U"] = free["free_turb_share_L"] = 0.0
    if wind_from_deg is None:
        # orientation-averaged chord length (Eq. 5) when the wind direction is unknown
        g = free["geometry"]
        D0 = 4.0 * g.S / g.P
        beta = math.gamma(0.5) * math.gamma(0.5 * N_FORCED + 1.0) / math.gamma(0.5 * N_FORCED + 1.5)
        Lf = D0 * (2.0 / math.pi * beta) ** (1.0 / (N_FORCED - 1.0))
        phi_w, frac = float("nan"), float("nan")
    else:
        Lf, phi_w, frac = forced_length(L, W, alpha_x_deg, alpha_y_deg, leaf_azimuth_deg, wind_from_deg)
    nf_forced = nusselt_forced(u, free["nu"], Lf, c_forced, transition=transition, TI_pct=TI_pct, bc=bc)
    out = dict(free)
    out["Re_f"] = np.asarray(u, dtype=float) * Lf / free["nu"]
    out["forced_regime"] = forced_regime(out["Re_f"], TI_pct) if transition else "laminar"
    out["gb_free_U"], out["gb_free_L"] = free["gb_U"], free["gb_L"]
    out["H_free_U"], out["H_free_L"] = free["H_U"], free["H_L"]
    # angle of attack: finite-plate enhancement f_e and windward/leeward factors (thin-plate fit)
    nf_forced_u, nf_forced_l = nf_forced, nf_forced
    beta_deg, windward = 0.0, "U"
    if angle_model:
        if wind_from_deg is not None:
            beta_deg, windward = wind_attack_angle(alpha_x_deg, alpha_y_deg, leaf_azimuth_deg, wind_from_deg)
        gW, gB = angle_of_attack_factors(beta_deg)
        fe = ANGLE_PARAMS["f_e"]
        if windward == "U":
            nf_forced_u, nf_forced_l = fe * nf_forced * gW, fe * nf_forced * gB
        else:
            nf_forced_u, nf_forced_l = fe * nf_forced * gB, fe * nf_forced * gW
    out["attack_deg"], out["windward_side"] = beta_deg, windward
    out["NuF_forced"] = nf_forced
    out["NuF_forced_U"], out["NuF_forced_L"] = nf_forced_u, nf_forced_l
    out["gb_forced"] = free["kappa"] * nf_forced
    out["gb_forced_U"], out["gb_forced_L"] = free["kappa"] * nf_forced_u, free["kappa"] * nf_forced_l
    out["Lf"], out["phi_w"], out["wind_in_plane_fraction"] = Lf, phi_w, frac
    mix_u = (free["NuF_U"] ** 3 + nf_forced_u ** 3) ** (1.0 / 3.0)
    mix_l = (free["NuF_L"] ** 3 + nf_forced_l ** 3) ** (1.0 / 3.0)
    out["NuF_mix_U"], out["NuF_mix_L"] = mix_u, mix_l
    out["gb_U"], out["gb_L"] = free["kappa"] * mix_u, free["kappa"] * mix_l
    TLa = np.asarray(TL, dtype=float); TAa = np.asarray(TA, dtype=float)
    out["H_U"] = free["rho"] * CP_AIR * out["gb_U"] * (TLa - TAa)
    out["H_L"] = free["rho"] * CP_AIR * out["gb_L"] * (TLa - TAa)
    share_u = nf_forced_u ** 3 / (free["NuF_U"] ** 3 + nf_forced_u ** 3)
    share_l = nf_forced_l ** 3 / (free["NuF_L"] ** 3 + nf_forced_l ** 3)
    out["forced_share"] = 0.5 * (share_u + share_l)
    return out


# ----------------------------------------------------- scalars: heat, H2O, CO2
# diffusivity exponent k of each transport mode: g_s = g_heat (D_s / kappa)^k, k = 1 - (Pr exponent)
K_COND, K_FREE_LAM, K_FREE_TURB, K_FORCED_LAM, K_FORCED_TURB = 1.0, 0.75, 2.0 / 3.0, 2.0 / 3.0, 0.4
D_H2O_0, D_CO2_0 = 2.12e-5, 1.39e-5    # m2/s at 273.15 K, 101325 Pa (Massman 1998)


def scalar_diffusivities(T, p=101325.0):
    """Molecular diffusivities of water vapour and CO2 in air (m2 s-1) at T (K), p (Pa)."""
    T = np.asarray(T, dtype=float)
    f = (T / 273.15) ** 1.75 * 101325.0 / np.asarray(p, dtype=float)
    return D_H2O_0 * f, D_CO2_0 * f


def scalar_exponent(nuc, nuf_free_lam, free_turb_share, forced_turb_share, forced_share):
    """Effective diffusivity exponent of one leaf side from the shares of the transport
    modes (generalisation of Eq. 17 of Graefe et al. 2022):
    laminar free convection k = 0.75 interpolated to 1 by the conduction share
    Nuc/NuF_lam, turbulent free convection 2/3, forced laminar/transitional 2/3,
    forced turbulent 0.4; weighted by the turbulent shares and the forced share."""
    s_c = np.clip(np.asarray(nuc) / np.maximum(np.asarray(nuf_free_lam), np.asarray(nuc)), 0.0, 1.0)
    k_free_lam = K_FREE_LAM + (K_COND - K_FREE_LAM) * s_c                  # Eq. 17
    k_free = (1.0 - free_turb_share) * k_free_lam + free_turb_share * K_FREE_TURB
    k_forced = (1.0 - forced_turb_share) * K_FORCED_LAM + forced_turb_share * K_FORCED_TURB
    return (1.0 - forced_share) * k_free + forced_share * k_forced


def leaf_conductances(TL, TA, u, L, W, alpha_x_deg, alpha_y_deg=0.0, leaf_azimuth_deg=180.0, wind_from_deg=None,
                      p=101325.0, eL=None, eA=None, TI_pct=1.0, bc="isothermal", transition=True,
                      free_transition=True, angle_model=True):
    """General convection model: boundary-layer conductances (m s-1) of both leaf sides
    for heat, water vapour and CO2 with free (laminar/turbulent), forced
    (laminar/transitional/turbulent) convection, angle of attack and conduction limit.

    Returns a dict with, per side X in (U, L): gb_X (heat), gv_X (H2O), gc_X (CO2),
    k_X (effective diffusivity exponent), shares cond_share_X, free_share_X,
    free_turb_share_X, forced_turb_share_X, plus the fields of mixed_convection()."""
    out = mixed_convection(TL, TA, u, L, W, alpha_x_deg, alpha_y_deg, leaf_azimuth_deg, wind_from_deg, p, eL, eA,
                           transition=transition, TI_pct=TI_pct, bc=bc, angle_model=angle_model, free_transition=free_transition)
    nuc = out["Nuc"]
    # forced: turbulent share of the Lienhard average (additive parts)
    if transition:
        lam, trans, turb = forced_nusselt_parts(out["Re_f"], PR, TI_pct, bc)
        tot = lam + trans + turb
        f_turb = np.where(tot > 0, turb / np.maximum(tot, 1e-300), 0.0)
    else:
        f_turb = np.zeros_like(np.asarray(out["Re_f"], dtype=float))
    TF = out["TF"]
    D_v, D_c = scalar_diffusivities(TF, p)
    for side in ("U", "L"):
        nf_free_lam = out.get(f"NuF_lam_{side}", out[f"NuF_free_{side}"] if f"NuF_free_{side}" in out else out[f"gb_free_{side}"] / out["kappa"])
        free_turb_share = out.get(f"free_turb_share_{side}", 0.0)
        nf_free = out[f"gb_free_{side}"] / out["kappa"]
        nf_forced = out[f"NuF_forced_{side}"]
        forced_share = nf_forced ** 3 / (nf_free ** 3 + nf_forced ** 3)
        k = scalar_exponent(nuc, nf_free_lam, free_turb_share, f_turb, forced_share)
        out[f"k_{side}"] = k
        out[f"gv_{side}"] = out[f"gb_{side}"] * (D_v / out["kappa"]) ** k
        out[f"gc_{side}"] = out[f"gb_{side}"] * (D_c / out["kappa"]) ** k
        out[f"cond_share_{side}"] = np.clip(nuc / np.maximum(nf_free_lam, nuc), 0.0, 1.0)
        out[f"free_share_{side}"] = 1.0 - forced_share
        out[f"free_turb_share_{side}"] = free_turb_share
        out[f"forced_turb_share_{side}"] = f_turb
    out["D_H2O"], out["D_CO2"] = D_v, D_c
    return out


def monteith_gb(TL, TA):
    """Simple empirical relation of Monteith (1973), Eq. 2 (m s-1, both sides)."""
    return np.abs(np.asarray(TL, dtype=float) - np.asarray(TA, dtype=float)) ** 0.33 / 400.0


if __name__ == "__main__":
    # Example: 8 x 4 cm leaf, 60 deg inclined, 3 K warmer than 20 degC air
    res = free_convection(296.15, 293.15, 0.08, 0.04, 60.0)
    g = res["geometry"]
    print(f"Lh={g.Lh*1e3:.2f} mm  Li={g.Li*1e3:.2f} mm  phi={math.degrees(g.phi):.1f} deg  theta={math.degrees(g.theta):.1f} deg")
    print(f"gb_U={res['gb_U']*1e3:.3f} mm/s  gb_L={res['gb_L']*1e3:.3f} mm/s  H_U={res['H_U']:.2f}  H_L={res['H_L']:.2f} W/m2")
    mix = mixed_convection(296.15, 293.15, 1.0, 0.08, 0.04, 60.0, 0.0, 180.0, wind_from_deg=270.0)
    print(f"wind 1 m/s from W: Lf={mix['Lf']*1e3:.2f} mm  phi_w={math.degrees(mix['phi_w']):.1f} deg  "
          f"gb_U={mix['gb_U']*1e3:.3f} mm/s  forced share={mix['forced_share']:.2f}")
    lc = leaf_conductances(296.15, 293.15, 1.0, 0.08, 0.04, 60.0, wind_from_deg=270.0, TI_pct=20.0)
    for s in ("U", "L"):
        print(f"side {s}: heat {lc[f'gb_{s}']*1e3:.2f}  H2O {lc[f'gv_{s}']*1e3:.2f}  CO2 {lc[f'gc_{s}']*1e3:.2f} mm/s  k={float(lc[f'k_{s}']):.3f}")
