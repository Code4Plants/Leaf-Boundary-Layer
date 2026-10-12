"""pytest suite for graefe_convection.py (run: python -m pytest python -q)."""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from graefe_convection import (C_FORCED, P3, PR, air_properties, characteristic_length, forced_length, free_convection,
                               leaf_geometry, mean_chord_length, mixed_convection, monteith_gb, nusselt_forced,
                               nusselt_free, rayleigh, wind_streamline_angle)

REF = Path(__file__).resolve().parent.parent / "reference" / "reference_values.json"


def test_geometry_circle():
    g = leaf_geometry(0.084, 0.084, 0.0)
    assert abs(g.DI - 1.0) < 0.005                      # circle: DI = 1
    assert abs(g.Lh - 0.084 / 4.0) < 1e-4               # S/P = D/4
    assert abs(g.S - math.pi * 0.042 ** 2) < 1e-12


def test_geometry_ellipse_axes_and_li():
    g = leaf_geometry(0.08, 0.04, 60.0)
    assert g.AR == 2.0
    assert abs(math.degrees(g.phi) - 90.0) < 1e-9       # alpha_y = 0: flow along the long axis
    assert abs(g.LCo - math.pi * 0.08 / 4.0) < 1e-12
    assert abs(g.Li - 0.816 * 0.08) < 0.0005            # paper: Li = 0.816 L
    assert abs(math.degrees(g.theta) - 60.0) < 1e-9
    g2 = leaf_geometry(0.08, 0.04, 60.0, 45.0)
    assert 0.0 < math.degrees(g2.phi) < 90.0            # pre-rotated leaf: oblique streamlines
    assert g2.Li < g.Li
    assert abs(math.cos(g2.theta) - math.cos(math.radians(60)) * math.cos(math.radians(45))) < 1e-12
    with pytest.raises(ValueError):
        leaf_geometry(0.04, 0.08, 0.0)


def test_air_properties_20C():
    rho, nu, kappa, k = air_properties(293.15)
    assert abs(rho - 1.204) < 0.005
    assert abs(nu - 1.51e-5) < 0.05e-5
    assert abs(kappa - 2.12e-5) < 0.1e-5
    assert abs(k - 0.0257) < 0.0005


def test_vertical_plate_matches_p3():
    """For a vertical leaf (theta = 90 deg) the horizontal term vanishes, so
    NuF * Li = p3 Ra_Li^0.25 apart from the conduction-limit blending."""
    g = leaf_geometry(0.08, 0.04, 90.0)
    TA = 293.15; dT = 12.0
    _, nu, _, _ = air_properties(TA + dT / 2)
    nf_u, nf_l, nuc, Ra = nusselt_free(dT, TA, nu, g)
    assert abs(nf_u - nf_l) < 1e-12                     # both sides identical when vertical
    nu_lam = P3 * g.Li ** (-0.25) * Ra ** 0.25
    assert nu_lam < nf_u < 1.3 * nu_lam                 # conduction limit adds a little
    Ra_L = Ra * g.Li ** 3
    assert abs(nu_lam * g.Li - P3 * Ra_L ** 0.25) < 1e-9


def test_cooled_leaf_swaps_sides():
    g = leaf_geometry(0.08, 0.04, 0.0)
    TA = 293.15
    _, nu, _, _ = air_properties(TA)
    hu, hl, _, _ = nusselt_free(6.0, TA, nu, g)
    cu, cl, _, _ = nusselt_free(-6.0, TA, nu, g)
    assert abs(hu - cl) < 1e-12 and abs(hl - cu) < 1e-12
    assert hu < hl                                       # heated horizontal leaf: p5 < p6 (upper side weaker)


def test_fig6_magnitudes_and_convergence():
    """Reproduces the behaviour of Fig. 6b: 8 x 4 cm leaf, dT = 3 and 12 K."""
    TA = 293.15
    for dT in (3.0, 12.0):
        prev_gap = None
        for ax in (0.0, 30.0, 60.0, 90.0):
            r = free_convection(TA + dT, TA, 0.08, 0.04, ax)
            assert 2.0 < r["H_U"] < 120.0 and 2.0 < r["H_L"] < 120.0
            gap = abs(r["H_U"] - r["H_L"]) / max(r["H_U"], r["H_L"])
            if prev_gap is not None:
                assert gap <= prev_gap + 1e-9            # sides converge with inclination
            prev_gap = gap
        assert prev_gap < 1e-9                           # identical at 90 deg
    r3 = free_convection(TA + 3.0, TA, 0.08, 0.04, 60.0)
    r12 = free_convection(TA + 12.0, TA, 0.08, 0.04, 60.0)
    assert r12["H_U"] > 3.0 * r3["H_U"]                  # H ~ dT^1.25


def test_vectorised_and_eq18():
    TA = 293.15
    dT = np.array([-5.0, -1.0, 1.0, 5.0])
    r = free_convection(TA + dT, TA, 0.08, 0.04, 60.0)
    assert r["gb_U"].shape == (4,) and np.all(r["gb_U"] > 0)
    assert np.all(np.sign(r["H_U"]) == np.sign(dT))
    # wet, transpiring leaf cooler than the air can still be buoyant (Eq. 18)
    r_wet = free_convection(TA - 0.5, TA, 0.08, 0.04, 60.0, eL=2300.0, eA=800.0)
    assert r_wet["dT"] > 0.0
    assert abs(monteith_gb(TA + 3.0, TA) - 3.0 ** 0.33 / 400.0) < 1e-12


def test_exact_characteristic_length_matches_eq9():
    """Closed form of Eq. 4 with n = 0.75 reproduces the linear fit of Eq. 9 (paper:
    Li = 0.816 L for flow along the long axis) and the rectangle limit."""
    for L, W in ((0.08, 0.04), (0.12, 0.04), (0.084, 0.084)):
        for phi in (0.0, math.pi / 4, math.pi / 2):
            exact = characteristic_length(L, W, phi, 0.75)
            fit = (0.04511 * 0.75 + 1.006) * mean_chord_length(L, W, phi)
            assert abs(exact / fit - 1.0) < 0.003
    assert abs(characteristic_length(0.08, 0.04, math.pi / 2, 0.75) - 0.816 * 0.08) < 1e-4
    # forced convection (n = 0.5): Lf = 1.0286 LCo
    assert abs(characteristic_length(0.08, 0.04, math.pi / 2, 0.5) / mean_chord_length(0.08, 0.04, math.pi / 2) - 1.0286) < 0.002
    # n = 1 degenerates to the mean chord length
    assert abs(characteristic_length(0.08, 0.04, 0.3, 1.0) - mean_chord_length(0.08, 0.04, 0.3)) < 1e-15


def test_wind_streamline_angle():
    # horizontal leaf facing south: long axis runs N-S -> north wind flows along the long axis (phi = 90 deg)
    phi, frac = wind_streamline_angle(0.0, 0.0, 180.0, 0.0)
    assert abs(math.degrees(phi) - 90.0) < 1e-9 and abs(frac - 1.0) < 1e-12
    phi, frac = wind_streamline_angle(0.0, 0.0, 180.0, 90.0)     # east wind: along the short axis
    assert abs(math.degrees(phi)) < 1e-9
    phi, _ = wind_streamline_angle(0.0, 0.0, 180.0, 45.0)
    assert abs(math.degrees(phi) - 45.0) < 1e-9
    # vertical leaf facing south: north/south wind hits the leaf normal -> no in-plane component
    phi, frac = wind_streamline_angle(90.0, 0.0, 180.0, 180.0)
    assert frac < 1e-9
    # inclined leaf, wind along the slope keeps phi = 90 but reduces the in-plane fraction to cos(60)
    phi, frac = wind_streamline_angle(60.0, 0.0, 180.0, 180.0)
    assert abs(math.degrees(phi) - 90.0) < 1e-9 and abs(frac - math.cos(math.radians(60))) < 1e-9
    Lf, phi_w, frac = forced_length(0.08, 0.04, 60.0, 0.0, 180.0, 0.0)
    assert abs(Lf - characteristic_length(0.08, 0.04, math.pi / 2, 0.5)) < 1e-15
    Lf_n, _, frac_n = forced_length(0.08, 0.04, 90.0, 0.0, 180.0, 180.0)
    assert frac_n < 1e-9 and 0.02 < Lf_n < 0.08     # orientation-averaged fallback


def test_forced_and_mixed():
    TA = 293.15
    _, nu, _, _ = air_properties(TA)
    Lf = characteristic_length(0.08, 0.04, math.pi / 2, 0.5)
    nf = nusselt_forced(1.0, nu, Lf)
    Re = 1.0 * Lf / nu
    assert abs(nf * Lf - C_FORCED * Re ** 0.5 * PR ** (1 / 3)) < 1e-9
    assert nusselt_forced(0.0, nu, Lf) == 0.0
    still = mixed_convection(TA + 3.0, TA, 0.0, 0.08, 0.04, 60.0, wind_from_deg=270.0)
    free = free_convection(TA + 3.0, TA, 0.08, 0.04, 60.0)
    assert abs(still["gb_U"] - free["gb_U"]) < 1e-15 and still["forced_share"] == 0.0
    windy = mixed_convection(TA + 3.0, TA, 3.0, 0.08, 0.04, 60.0, wind_from_deg=270.0)
    assert windy["gb_U"] > free["gb_U"] and windy["forced_share"] > 0.9
    assert windy["gb_U"] ** 3 >= windy["gb_forced"] ** 3 + free["gb_U"] ** 3 - 1e-18
    # wind along the long axis gives a longer streamline and thus a smaller forced conductance than across it
    along = mixed_convection(TA + 3.0, TA, 2.0, 0.08, 0.04, 0.0, wind_from_deg=0.0)
    across = mixed_convection(TA + 3.0, TA, 2.0, 0.08, 0.04, 0.0, wind_from_deg=90.0)
    assert along["Lf"] > across["Lf"] and along["gb_forced"] < across["gb_forced"]
    # unknown wind direction -> orientation-averaged length between the two
    unknown = mixed_convection(TA + 3.0, TA, 2.0, 0.08, 0.04, 0.0)
    assert across["Lf"] < unknown["Lf"] < along["Lf"]
    # vectorised wind speeds
    vec = mixed_convection(TA + 3.0, TA, np.array([0.1, 1.0, 5.0]), 0.08, 0.04, 60.0, wind_from_deg=200.0)
    assert np.all(np.diff(vec["gb_U"]) > 0)


def test_lienhard_transition():
    from graefe_convection import forced_nusselt_lienhard, forced_regime, transition_onset_reynolds, free_convection_transition
    # onset: smooth stream capped at 5e5, 5 % turbulence ~ 4.8e4, floor 2e3
    assert transition_onset_reynolds(0.1) == 5e5
    assert abs(transition_onset_reynolds(5.0) / 4.8e4 - 1) < 0.05
    assert transition_onset_reynolds(100.0) == 2e3
    Re = np.array([1e3, 1e4, 1e5, 1e6, 3e6])
    nu = forced_nusselt_lienhard(Re, 0.71, TI_pct=0.1)       # Re_l = 5e5
    lam = 0.664 * Re ** 0.5 * 0.71 ** (1 / 3)
    assert np.allclose(nu[:3], lam[:3])                       # laminar below onset
    assert nu[3] > lam[3] and nu[4] > lam[4]                  # enhanced beyond
    turb_only = 0.037 * 0.71 ** 0.6 * 3e6 ** 0.8
    assert 0.7 * turb_only < nu[4] < turb_only                # below the all-turbulent plate (laminar start)
    assert list(forced_regime(Re, 0.1))[-1] == "turbulent" and list(forced_regime(Re, 0.1))[0] == "laminar"
    # blunt leading edge: fully turbulent average
    from graefe_convection import forced_nusselt_turbulent
    assert abs(float(forced_nusselt_turbulent(3e6)) - turb_only) < 1e-9
    # continuity in Re and monotonic increase
    grid = np.logspace(3, 7, 400)
    val = forced_nusselt_lienhard(grid, 0.71, TI_pct=2.0)
    assert np.all(np.diff(val) > 0)
    assert np.max(np.abs(np.diff(np.log(val)))) < 0.05
    # higher turbulence -> earlier transition -> larger Nu at 1e5
    assert forced_nusselt_lienhard(1e5, 0.71, 10.0) > forced_nusselt_lienhard(1e5, 0.71, 1.0)
    assert forced_nusselt_lienhard(4e5, 0.71, 1.0) > forced_nusselt_lienhard(4e5, 0.71, 0.1)
    # uniform heat flux larger than isothermal in the laminar range
    assert forced_nusselt_lienhard(1e4, 0.71, 1.0, "uniform_flux") > forced_nusselt_lienhard(1e4, 0.71, 1.0)
    # mixed convection with transition: regime reported and conductance grows with u
    m = mixed_convection(293.15 + 3, 293.15, np.array([0.5, 5.0, 20.0]), 0.08, 0.04, 60.0, wind_from_deg=270.0, transition=True, TI_pct=10.0)
    assert np.all(np.diff(m["gb_U"]) > 0) and m["forced_regime"][0] == "laminar"
    # free-convection transition: identical at low Ra, larger at very high Ra on the unstable side
    small = free_convection_transition(293.15 + 1, 293.15, 0.08, 0.04, 0.0)
    base = free_convection(293.15 + 1, 293.15, 0.08, 0.04, 0.0)
    assert abs(small["gb_U"] / base["gb_U"] - 1) < 0.02
    big = free_convection_transition(293.15 + 40, 293.15, 1.5, 1.5, 0.0)
    base_big = free_convection(293.15 + 40, 293.15, 1.5, 1.5, 0.0)
    assert big["gb_U"] > 1.3 * base_big["gb_U"] and abs(big["gb_L"] - base_big["gb_L"]) < 1e-12


def test_angle_of_attack_model():
    from graefe_convection import ANGLE_PARAMS, angle_of_attack_factors, forced_nusselt_leaf, wind_attack_angle
    gW, gB = angle_of_attack_factors(0.0)
    assert abs(gW - 1.0) < 1e-12 and abs(gB - 1.0) < 1e-12
    gW90, gB90 = angle_of_attack_factors(90.0)
    assert abs(gW90 - (1 + ANGLE_PARAMS["a_W"])) < 1e-12 and abs(gB90 - ANGLE_PARAMS["c_w"]) < 0.01
    b = np.arange(0, 91, 5.0)
    _, gBs = angle_of_attack_factors(b)
    assert gBs.max() > 1.15 and 10 <= b[np.argmax(gBs)] <= 35     # reattachment peak at small angles
    assert gBs[-1] < 1.0                                           # wake regime below parallel flow
    nuW, nuB = forced_nusselt_leaf(1e4, 20.0, TI_pct=1.0)
    base = 0.664 * 1e4 ** 0.5 * 0.71 ** (1 / 3)
    assert abs(nuW / base - ANGLE_PARAMS["f_e"] * angle_of_attack_factors(20.0)[0]) < 1e-9
    assert nuB > nuW                                               # reattaching leeward side at 20 deg
    # wind geometry: horizontal leaf -> beta = 0 for any direction; vertical south-facing leaf, north wind hits the lower side
    assert abs(wind_attack_angle(0.0, 0.0, 180.0, 90.0)[0]) < 1e-9
    beta, side = wind_attack_angle(90.0, 0.0, 180.0, 180.0)        # wind from the south blows onto the upper (south-facing) side
    assert abs(beta - 90.0) < 1e-9 and side == "U"
    beta, side = wind_attack_angle(90.0, 0.0, 180.0, 0.0)
    assert abs(beta - 90.0) < 1e-9 and side == "L"
    beta, side = wind_attack_angle(60.0, 0.0, 180.0, 180.0)
    assert abs(beta - 60.0) < 1e-9 and side == "U"
    m = mixed_convection(296.15, 293.15, 3.0, 0.08, 0.04, 60.0, wind_from_deg=180.0, transition=True, TI_pct=10.0, angle_model=True)
    assert m["windward_side"] == "U" and abs(m["attack_deg"] - 60.0) < 1e-9
    assert m["NuF_forced_U"] > m["NuF_forced_L"]                   # windward side stronger at 60 deg
    m0 = mixed_convection(296.15, 293.15, 3.0, 0.08, 0.04, 60.0, wind_from_deg=180.0, transition=True, TI_pct=10.0)
    assert m0["NuF_forced_U"] == m0["NuF_forced_L"]                # without the angle model both sides equal


def test_reference_values_angle():
    data = json.loads(REF.read_text(encoding="utf-8"))
    for c in data["angle_cases"]:
        m = mixed_convection(296.15, 293.15, c["u"], 0.08, 0.04, c["alpha_x"], c["alpha_y"], c["leaf_azimuth"],
                             wind_from_deg=c["wind_from"], transition=True, TI_pct=c["TI"], angle_model=True)
        for key in ("attack_deg", "Re_f", "NuF_forced_U", "NuF_forced_L", "gb_U", "gb_L", "H_U", "H_L"):
            assert abs(float(m[key]) - c[key]) <= 1e-9 * max(1.0, abs(c[key])), key
        assert m["windward_side"] == c["windward_side"]


def test_scalar_conductances():
    from graefe_convection import (K_COND, K_FORCED_LAM, K_FORCED_TURB, K_FREE_LAM, leaf_conductances,
                                   scalar_diffusivities, scalar_exponent)
    D_v, D_c = scalar_diffusivities(293.15)
    assert abs(D_v / D_c - 2.12 / 1.39) < 1e-12 and 2.3e-5 < D_v < 2.5e-5
    # pure modes
    assert abs(scalar_exponent(0.0, 1.0, 0.0, 0.0, 0.0) - K_FREE_LAM) < 1e-12
    assert abs(scalar_exponent(1.0, 1.0, 0.0, 0.0, 0.0) - K_COND) < 1e-12
    assert abs(scalar_exponent(0.0, 1.0, 0.0, 0.0, 1.0) - K_FORCED_LAM) < 1e-12
    assert abs(scalar_exponent(0.0, 1.0, 0.0, 1.0, 1.0) - K_FORCED_TURB) < 1e-12
    TA = 293.15
    still = leaf_conductances(TA + 3, TA, 0.0, 0.08, 0.04, 60.0, wind_from_deg=180.0)
    assert K_FREE_LAM < float(still["k_U"]) < K_COND and float(still["cond_share_U"]) > 0.2
    assert float(still["gv_U"]) > float(still["gb_U"]) > float(still["gc_U"])       # D_H2O > kappa > D_CO2
    assert abs(float(still["gv_U"]) / float(still["gb_U"]) - (still["D_H2O"] / still["kappa"]) ** float(still["k_U"])) < 1e-12
    windy = leaf_conductances(TA + 3, TA, 3.0, 0.08, 0.04, 60.0, wind_from_deg=180.0, TI_pct=1.0)
    assert abs(float(windy["k_U"]) - K_FORCED_LAM) < 0.01 and float(windy["free_share_U"]) < 0.05
    # large leaf, strong turbulent wind: turbulent share pushes k towards 0.4
    big = leaf_conductances(TA + 3, TA, 15.0, 0.6, 0.3, 0.0, wind_from_deg=90.0, TI_pct=40.0)
    assert float(big["forced_turb_share_U"]) > 0.3 and float(big["k_U"]) < K_FORCED_LAM
    # hot large leaf in still air: free turbulent share > 0 lowers k below the laminar value
    hot = leaf_conductances(TA + 15, TA, 0.0, 0.6, 0.3, 0.0, wind_from_deg=90.0)
    assert float(hot["free_turb_share_U"]) > 0.3 and float(hot["k_U"]) < float(still["k_U"])
    # vectorised wind speeds, sides independent
    vec = leaf_conductances(TA + 3, TA, np.array([0.0, 1.0, 5.0]), 0.08, 0.04, 60.0, wind_from_deg=180.0)
    assert vec["gv_U"].shape == (3,) and np.all(np.diff(vec["gc_L"]) > 0)


def test_reference_values_scalar():
    from graefe_convection import leaf_conductances
    data = json.loads(REF.read_text(encoding="utf-8"))
    for c in data["scalar_cases"]:
        r = leaf_conductances(293.15 + c["dT"], 293.15, c["u"], c["L"], c["W"], c["alpha_x"], 0.0, 180.0,
                              wind_from_deg=c["wind_from"], TI_pct=c["TI"])
        for key in ("gb_U", "gb_L", "gv_U", "gv_L", "gc_U", "gc_L", "k_U", "k_L", "cond_share_U", "free_share_U",
                    "free_turb_share_U", "forced_turb_share_U"):
            assert abs(float(r[key]) - c[key]) <= 1e-9 * max(1.0, abs(c[key])), key


def test_reference_values_mixed():
    data = json.loads(REF.read_text(encoding="utf-8"))
    for c in data["mixed_cases"]:
        m = mixed_convection(c["TL"], c["TA"], c["u"], c["L"], c["W"], c["alpha_x"], c["alpha_y"], c["leaf_azimuth"],
                             wind_from_deg=c["wind_from"])
        for key in ("Lf", "NuF_forced", "NuF_mix_U", "NuF_mix_L", "gb_U", "gb_L", "H_U", "H_L", "forced_share"):
            assert abs(float(m[key]) - c[key]) <= 1e-9 * max(1.0, abs(c[key])), key


def test_reference_values():
    data = json.loads(REF.read_text(encoding="utf-8"))
    for c in data["cases"]:
        r = free_convection(c["TL"], c["TA"], c["L"], c["W"], c["alpha_x"], c["alpha_y"],
                            eL=c.get("eL"), eA=c.get("eA"))
        for key in ("gb_U", "gb_L", "H_U", "H_L", "NuF_U", "NuF_L", "Nuc", "Ra"):
            assert abs(float(r[key]) - c[key]) <= 1e-9 * max(1.0, abs(c[key])), key
        assert abs(r["geometry"].Li - c["Li"]) < 1e-12
