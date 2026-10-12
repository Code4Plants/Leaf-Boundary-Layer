"""Generate reference values (JSON) used by the Python and MATLAB cross-check tests."""
import json
import math
from pathlib import Path

from graefe_convection import free_convection, leaf_conductances, mixed_convection

cases = []
for L, W in ((0.08, 0.04), (0.084, 0.084), (0.12, 0.04)):
    for ax in (0.0, 30.0, 60.0, 90.0):
        for ay in (0.0, 45.0):
            for dT in (-12.0, -3.0, 1.0, 3.0, 12.0):
                TA = 293.15
                r = free_convection(TA + dT, TA, L, W, ax, ay)
                g = r["geometry"]
                cases.append({"L": L, "W": W, "alpha_x": ax, "alpha_y": ay, "TA": TA, "TL": TA + dT,
                              "Lh": g.Lh, "Li": g.Li, "phi_deg": math.degrees(g.phi), "theta_deg": math.degrees(g.theta),
                              "DI": g.DI, "Nuc": float(r["Nuc"]), "Ra": float(r["Ra"]),
                              "NuF_U": float(r["NuF_U"]), "NuF_L": float(r["NuF_L"]),
                              "gb_U": float(r["gb_U"]), "gb_L": float(r["gb_L"]),
                              "H_U": float(r["H_U"]), "H_L": float(r["H_L"])})
# one case with the vapour correction (Eq. 18): wet leaf at 25 degC, air 20 degC / 50 % RH
r = free_convection(298.15, 293.15, 0.08, 0.04, 60.0, 0.0, 101325.0, eL=3169.0, eA=1169.0)
cases.append({"L": 0.08, "W": 0.04, "alpha_x": 60.0, "alpha_y": 0.0, "TA": 293.15, "TL": 298.15, "eL": 3169.0, "eA": 1169.0,
              "dT": float(r["dT"]), "gb_U": float(r["gb_U"]), "gb_L": float(r["gb_L"]), "H_U": float(r["H_U"]), "H_L": float(r["H_L"]),
              "NuF_U": float(r["NuF_U"]), "NuF_L": float(r["NuF_L"]), "Nuc": float(r["Nuc"]), "Ra": float(r["Ra"]),
              "Lh": r["geometry"].Lh, "Li": r["geometry"].Li, "phi_deg": math.degrees(r["geometry"].phi),
              "theta_deg": math.degrees(r["geometry"].theta), "DI": r["geometry"].DI})

# mixed convection cases: wind speed, direction and leaf azimuth
mixed = []
for L, W in ((0.08, 0.04), (0.12, 0.04)):
    for ax, ay in ((0.0, 0.0), (60.0, 0.0), (60.0, 45.0), (90.0, 0.0)):
        for az in (180.0, 90.0):
            for wd in (0.0, 90.0, 180.0, 225.0, None):
                for u in (0.3, 2.0):
                    for dT in (-3.0, 3.0):
                        TA = 293.15
                        m = mixed_convection(TA + dT, TA, u, L, W, ax, ay, az, wind_from_deg=wd)
                        mixed.append({"L": L, "W": W, "alpha_x": ax, "alpha_y": ay, "leaf_azimuth": az,
                                      "wind_from": wd, "u": u, "TA": TA, "TL": TA + dT,
                                      "Lf": float(m["Lf"]), "phi_w_deg": None if wd is None else math.degrees(m["phi_w"]),
                                      "NuF_forced": float(m["NuF_forced"]), "NuF_mix_U": float(m["NuF_mix_U"]),
                                      "NuF_mix_L": float(m["NuF_mix_L"]), "gb_U": float(m["gb_U"]), "gb_L": float(m["gb_L"]),
                                      "H_U": float(m["H_U"]), "H_L": float(m["H_L"]), "forced_share": float(m["forced_share"])})
# angle-of-attack / transition cases: wind direction, turbulence intensity, angle model on
angle = []
for ax, ay, az in ((0.0, 0.0, 180.0), (60.0, 0.0, 180.0), (60.0, 45.0, 90.0), (90.0, 0.0, 180.0)):
    for wd in (0.0, 45.0, 90.0, 180.0, 270.0):
        for u, ti in ((0.5, 1.0), (3.0, 10.0), (12.0, 30.0)):
            m = mixed_convection(296.15, 293.15, u, 0.08, 0.04, ax, ay, az, wind_from_deg=wd, transition=True, TI_pct=ti, angle_model=True)
            angle.append({"alpha_x": ax, "alpha_y": ay, "leaf_azimuth": az, "wind_from": wd, "u": u, "TI": ti,
                          "attack_deg": float(m["attack_deg"]), "windward_side": m["windward_side"], "Re_f": float(m["Re_f"]),
                          "forced_regime": str(m["forced_regime"]), "NuF_forced_U": float(m["NuF_forced_U"]),
                          "NuF_forced_L": float(m["NuF_forced_L"]), "gb_U": float(m["gb_U"]), "gb_L": float(m["gb_L"]),
                          "H_U": float(m["H_U"]), "H_L": float(m["H_L"])})
# scalar conductance cases (general model: free transition, forced transition, angle model)
scalar = []
for dT in (-6.0, 3.0, 12.0):
    for u, ti in ((0.0, 1.0), (0.4, 5.0), (2.5, 20.0), (12.0, 40.0)):
        for ax, wd in ((0.0, 90.0), (60.0, 180.0), (90.0, 0.0)):
            for L, W in ((0.08, 0.04), (0.6, 0.3)):
                r = leaf_conductances(293.15 + dT, 293.15, u, L, W, ax, 0.0, 180.0, wind_from_deg=wd, TI_pct=ti)
                scalar.append({"dT": dT, "u": u, "TI": ti, "alpha_x": ax, "wind_from": wd, "L": L, "W": W,
                               **{k: float(r[k]) for k in ("gb_U", "gb_L", "gv_U", "gv_L", "gc_U", "gc_L", "k_U", "k_L",
                                                            "cond_share_U", "free_share_U", "free_turb_share_U", "forced_turb_share_U")}})
out = Path(__file__).resolve().parent.parent / "reference" / "reference_values.json"
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps({"source": "python/graefe_convection.py", "cases": cases, "mixed_cases": mixed, "angle_cases": angle,
                           "scalar_cases": scalar}, indent=1), encoding="utf-8")
print(f"{len(cases)} free + {len(mixed)} mixed + {len(angle)} angle + {len(scalar)} scalar cases written to {out}")
