"""Comparison of the laminar free-convection core with Ryms, Tesch & Lewandowski (2025).

Ryms, M., Tesch, K., Lewandowski, W.M. (2025). Experimental and numerical study of free
convective heat transfer in air performed simultaneously from the upper and lower surfaces
of a horizontal plate heated on both sides. International Journal of Heat and Mass Transfer
245, 126973, https://doi.org/10.1016/j.ijheatmasstransfer.2025.126973

A horizontal 150 x 75 mm plate heated on both sides ("Type II") is the situation of a
horizontal leaf; the paper gives Nusselt-Rayleigh correlations per side from experiments
(Eqs. 29/32), from CFD of the same two-sided plate (Eqs. 41/42) and of a one-sided plate
(Eqs. 39/40), all in the form Nu_b = C Ra_b^0.2 with the plate width b = 75 mm as length.

Only leaf-air temperature differences up to DT_MAX = 10 K are considered (the leaf range);
the points of the paper with larger differences (VI/7 at 10.7 / 11.1 K, IX/14 at 57 K) are
left out.

Run from the repository root:  python python/validate_ryms2025.py
Prints Markdown tables; the numbers are quoted in docs/validation.md.
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from graefe_convection import conduction_limit, free_convection, leaf_geometry  # noqa: E402

DT_MAX = 10.0                                           # K, upper limit of the temperature differences considered
B, LEN = 0.075, 0.150                                   # plate width (characteristic length) and length, m
S_RECT = B * LEN
L_ELL = math.sqrt(4.0 * S_RECT * (LEN / B) / math.pi)   # equal-area ellipse with the same aspect ratio
W_ELL = L_ELL * B / LEN

# Nu_b = C Ra_b^0.2 (n = 1/5 imposed), Table 7 of the paper
C_PAPER = {"Exp. II (29/32)": (0.804, 0.744), "CFD Type II (41/42)": (0.601, 0.812),
           "CFD Type I (39/40)": (0.941, 0.808), "MLC.I (53/54)": (1.093, 0.660)}

# Table 1: series/point, Tw (degC), Tinf (degC), Nu (Eq. 17), Ra (Eq. 18); filtered to Tw - Tinf <= DT_MAX
TABLE1_ALL = {"U": [("III/5", 23.16, 20.07, 7.633, 1.333e5), ("VI/7", 30.62, 19.90, 9.583, 4.365e5),
                    ("IX/14", 78.28, 21.27, 15.446, 1.614e6), ("XIV/10", 26.73, 19.83, 9.028, 2.985e5)],
              "L": [("III/5", 23.45, 20.07, 5.688, 1.456e5), ("VI/7", 30.99, 19.90, 8.742, 4.501e5),
                    ("IX/14", 78.20, 21.27, 14.330, 1.613e6), ("XIV/10", 26.60, 19.83, 7.476, 2.843e5)]}
TABLE1 = {side: [pt for pt in pts if pt[1] - pt[2] <= DT_MAX] for side, pts in TABLE1_ALL.items()}


def model_c(dT, TA=293.15):
    """C_U, C_L = Nu_b / Ra_b^0.2 and Ra_b of the laminar model for the equal-area ellipse."""
    r = free_convection(TA + dT, TA, L_ELL, W_ELL, 0.0)
    Ra_b = float(r["Ra"]) * B ** 3
    return float(r["NuF_U"]) * B / Ra_b ** 0.2, float(r["NuF_L"]) * B / Ra_b ** 0.2, Ra_b


def main():
    g = leaf_geometry(L_ELL, W_ELL, 0.0)
    print(f"Plate 150 x 75 mm -> equal-area ellipse L = {L_ELL * 1e3:.1f} mm, W = {W_ELL * 1e3:.1f} mm, "
          f"S = {g.S * 1e4:.2f} cm2, Lh = {g.Lh * 1e3:.2f} mm (rectangle A/P = {S_RECT / (2 * (B + LEN)) * 1e3:.1f} mm), "
          f"Nu_c b = {conduction_limit(g) * B:.2f}; temperature differences up to {DT_MAX:g} K")
    print("\nCoefficients C = Nu_b / Ra_b^0.2 (n = 1/5), Ra_b = Ra b^3, b = 75 mm\n")
    print("| Source | C_U | C_L | U/L |\n| --- | --- | --- | --- |")
    for name, (cu, cl) in C_PAPER.items():
        print(f"| {name} | {cu:.3f} | {cl:.3f} | {cu / cl:.2f} |")
    for dT in (3.0, 5.0, 8.0, 10.0):
        cu, cl, Ra_b = model_c(dT)
        print(f"| model, laminar (Graefe 2022), dT = {dT:g} K (Ra_b = {Ra_b:.1e}) | {cu:.3f} | {cl:.3f} | {cu / cl:.2f} |")
    print(f"\nTable 1 points with dT <= {DT_MAX:g} K (Nu_model / Nu_measured, laminar core)\n")
    print("| Point | Side | dT (K) | Ra_b (paper) | Ra_b (model) | Nu meas. | Nu model | model/meas |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    ratios = {"U": [], "L": []}
    for side, pts in TABLE1.items():
        for pt, Tw, Ti, Nu, Ra in pts:
            TL, TA = Tw + 273.15, Ti + 273.15
            r = free_convection(TL, TA, L_ELL, W_ELL, 0.0)
            lam = float(r[f"NuF_{side}"]) * B
            ratios[side].append(lam / Nu)
            print(f"| {pt} | {side} | {Tw - Ti:.1f} | {Ra:.3e} | {float(r['Ra']) * B ** 3:.3e} | {Nu:.3f} | {lam:.3f} | "
                  f"{lam / Nu:.2f} |")
    for side in ("U", "L"):
        print(f"geometric mean {side}: {math.exp(np.mean(np.log(ratios[side]))):.2f}")
    # Rayleigh numbers of leaves against the turbulent onset (~1e7 horizontal with L_h, ~1e9 vertical with L_i)
    print("\nRayleigh numbers of heated leaves at 20 degC (laminar range; no turbulent branch in the model)\n")
    print("| Leaf | dT (K) | Ra (m-3) | Ra_Lh (horizontal) | Ra_Li (vertical) |")
    print("| --- | --- | --- | --- | --- |")
    TA = 293.15
    for L, W, lab in ((0.08, 0.04, "8 x 4 cm"), (0.20, 0.10, "20 x 10 cm"), (0.30, 0.15, "30 x 15 cm"), (2.0, 0.5, "200 x 50 cm")):
        gl = leaf_geometry(L, W, 0.0)
        for dT in (3.0, 10.0):
            Ra = float(free_convection(TA + dT, TA, L, W, 0.0)["Ra"])
            print(f"| {lab} | {dT:g} | {Ra:.1e} | {Ra * gl.Lh ** 3:.1e} | {Ra * gl.Li ** 3:.1e} |")


if __name__ == "__main__":
    main()
