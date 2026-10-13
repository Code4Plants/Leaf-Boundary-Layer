"""Figures and numerical tables of the model documentation (docs/manual).

Run from the repository root:  python docs/manual/make_figures.py
Writes PDF figures to docs/manual/fig/ and LaTeX macros/tables to docs/manual/generated/.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "python"))
import graefe_convection as gc  # noqa: E402

FIG = HERE / "fig"; FIG.mkdir(exist_ok=True)
GEN = HERE / "generated"; GEN.mkdir(exist_ok=True)

plt.rcParams.update({"font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True, "axes.unicode_minus": False, "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 1.4, "savefig.bbox": "tight"})
C1, C2, C3, C4 = "#1f4e79", "#c0504d", "#4f8a3c", "#8064a2"
TA = 293.15
macros = {}


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- Fig. geometry sketch
fig, ax = plt.subplots(figsize=(3.3, 3.0))
L, W = 2.0, 1.0
t = np.linspace(0, 2 * np.pi, 400)
ax.plot(0.5 * W * np.cos(t), 0.5 * L * np.sin(t), color="k")
ax.annotate("", xy=(0.75, 0), xytext=(-0.75, 0), arrowprops=dict(arrowstyle="->", color="0.4"))
ax.annotate("", xy=(0, 1.2), xytext=(0, -1.2), arrowprops=dict(arrowstyle="->", color="0.4"))
ax.text(0.78, -0.08, r"$w$ (short axis, $W$)", fontsize=10, va="top")
ax.text(0.05, 1.18, r"$y$ (long axis, $L$)", fontsize=10)
phi = math.radians(55)
d = np.array([math.cos(phi), math.sin(phi)]); nrm = np.array([-d[1], d[0]])
for off in np.linspace(-0.55, 0.55, 7):
    p0 = off * nrm
    # chord of the ellipse along d through p0
    a, b = 0.5 * W, 0.5 * L
    A = (d[0] / a) ** 2 + (d[1] / b) ** 2
    B = 2 * (p0[0] * d[0] / a ** 2 + p0[1] * d[1] / b ** 2)
    Cq = (p0[0] / a) ** 2 + (p0[1] / b) ** 2 - 1
    disc = B * B - 4 * A * Cq
    if disc <= 0:
        continue
    s1, s2 = (-B - math.sqrt(disc)) / (2 * A), (-B + math.sqrt(disc)) / (2 * A)
    ax.plot(*zip(p0 + s1 * d, p0 + s2 * d), color=C1, lw=0.8)
    ax.annotate("", xy=p0 + s2 * d, xytext=p0 + (s2 - 0.12) * d, arrowprops=dict(arrowstyle="->", color=C1, lw=0.8))
arc = np.linspace(0, phi, 30)
ax.plot(0.35 * np.cos(arc), 0.35 * np.sin(arc), color=C2, lw=1)
ax.text(0.38, 0.14, r"$\varphi$", color=C2)
ax.text(-1.05, -1.42, "streamlines (chords) at angle $\\varphi$\nto the short axis", fontsize=10, color=C1)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.6); ax.set_ylim(-1.45, 1.35)
save(fig, "geometry")

# ---------------------------------------------------------------- Fig. characteristic lengths
phis = np.linspace(0, 90, 91)
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.6), gridspec_kw={"wspace": 0.38})
worst = 0.0
for AR, col in ((1.0, C1), (2.0, C2), (3.0, C3), (4.0, C4)):
    Lx, Wx = 0.08, 0.08 / AR
    lco = np.array([gc.mean_chord_length(Lx, Wx, math.radians(p)) for p in phis])
    exact = np.array([gc.characteristic_length(Lx, Wx, math.radians(p), 0.75) for p in phis])
    worst = max(worst, np.max(np.abs(exact / ((0.04511 * 0.75 + 1.006) * lco) - 1)))
    axs[0].plot(phis, lco / Lx, color=col, label=f"AR = {AR:g}")
axs[0].set_xlabel(r"streamline angle $\varphi$ (deg)"); axs[0].set_ylabel(r"$L_{Co}(\varphi)\,/\,L$")
axs[0].set_xlim(0, 90); axs[0].set_xticks([0, 30, 60, 90])
axs[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.92)); axs[0].set_title("(a) mean chord length, Eq. (9)", fontsize=9, loc="left")
# (b) L_c / L_Co depends only on n (self-similar chord distribution of the ellipse)
ns = np.linspace(0.3, 0.99, 200)
ratio_n = np.array([gc.characteristic_length(0.08, 0.04, 0.3, n) / gc.mean_chord_length(0.08, 0.04, 0.3) for n in ns])
axs[1].plot(ns, ratio_n, color=C1, label="exact, Eq. (11)")
axs[1].plot(ns, 0.04511 * ns + 1.006, color=C2, ls="--", label="linear fit, Eq. (10)")
for n_mark, lab in ((0.5, "forced"), (0.75, "free")):
    r_mark = gc.characteristic_length(0.08, 0.04, 0.3, n_mark) / gc.mean_chord_length(0.08, 0.04, 0.3)
    axs[1].plot([n_mark], [r_mark], "o", color="k", ms=3.5)
    axs[1].annotate(f"$n$ = {n_mark:g} ({lab})", (n_mark, r_mark), xytext=(6, -12), textcoords="offset points", fontsize=8)
axs[1].set_xlabel(r"exponent $n$"); axs[1].set_ylabel(r"$L_c\,/\,L_{Co}$")
axs[1].legend(frameon=False, loc="upper left"); axs[1].set_title("(b) flux-weighted length, all $\\varphi$ and AR", fontsize=9, loc="left")
save(fig, "charlength")
macros["MaxDevEqNine"] = f"{worst * 100:.2f}"
lco = gc.mean_chord_length(0.08, 0.04, math.pi / 2)
macros["LfOverLCo"] = f"{gc.characteristic_length(0.08, 0.04, math.pi / 2, 0.5) / lco:.4f}"
macros["LiOverLCo"] = f"{gc.characteristic_length(0.08, 0.04, math.pi / 2, 0.75) / lco:.4f}"
macros["LiFitOverLCo"] = f"{0.04511 * 0.75 + 1.006:.4f}"

# ---------------------------------------------------------------- Fig. free convection (cf. Fig. 6b of the paper)
axs_deg = np.arange(0, 91, 2.0)
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.5))
for dT, ls in ((3.0, "-"), (12.0, "--")):
    HU = [float(gc.free_convection(TA + dT, TA, 0.08, 0.04, a)["H_U"]) for a in axs_deg]
    HL = [float(gc.free_convection(TA + dT, TA, 0.08, 0.04, a)["H_L"]) for a in axs_deg]
    axs[0].plot(axs_deg, HU, ls, color=C1, label=f"upper, $\\Delta T$ = {dT:g} K")
    axs[0].plot(axs_deg, HL, ls, color=C2, label=f"lower, $\\Delta T$ = {dT:g} K")
axs[0].set_xlabel(r"inclination $\alpha_x$ (deg)"); axs[0].set_ylabel(r"$H$ (W m$^{-2}$)")
axs[0].legend(frameon=False, ncol=1); axs[0].set_title("(a) sensible heat flux, 8 cm $\\times$ 4 cm", fontsize=9, loc="left")
dTs = np.linspace(0.1, 15, 150)
for a, col in ((0.0, C1), (45.0, C3), (90.0, C2)):
    r = gc.free_convection(TA + dTs, TA, 0.08, 0.04, a)
    axs[1].plot(dTs, r["gb_U"] * 1e3, color=col, label=f"upper, $\\alpha_x$ = {a:g}$^\\circ$")
    axs[1].plot(dTs, r["gb_L"] * 1e3, ":", color=col, label=f"lower, $\\alpha_x$ = {a:g}$^\\circ$")
axs[1].set_xlabel(r"$T_L - T_A$ (K)"); axs[1].set_ylabel(r"$g_b$ (mm s$^{-1}$)")
axs[1].legend(frameon=False, fontsize=7, ncol=2); axs[1].set_title("(b) conductance per side", fontsize=9, loc="left")
save(fig, "free")

# ---------------------------------------------------------------- Fig. Lienhard transition
Re = np.logspace(3, 6.5, 400)
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.5))
lam = 0.664 * Re ** 0.5 * gc.PR ** (1 / 3)
axs[0].loglog(Re, lam, color="0.5", lw=1, label="laminar 0.664 Re$^{1/2}$Pr$^{1/3}$")
axs[0].loglog(Re, gc.forced_nusselt_turbulent(Re), color="0.5", lw=1, ls=":", label="turbulent from edge")
for ti, col in ((0.1, C1), (2.0, C3), (10.0, C2), (30.0, C4)):
    axs[0].loglog(Re, gc.forced_nusselt_lienhard(Re, gc.PR, ti), color=col, label=f"TI = {ti:g} %")
axs[0].set_xlabel(r"$Re = u L_f/\nu$"); axs[0].set_ylabel(r"$\overline{Nu}$")
axs[0].legend(frameon=False, fontsize=7); axs[0].set_title("(a) average flat-plate Nusselt number", fontsize=9, loc="left")
TI = np.logspace(-1, 2, 200)
axs[1].loglog(TI, gc.transition_onset_reynolds(TI), color=C1)
axs[1].set_xlabel("turbulence intensity TI (%)"); axs[1].set_ylabel(r"onset $Re_l$")
axs[1].set_title("(b) transition onset (Mayle 1991)", fontsize=9, loc="left")
save(fig, "lienhard")
for ti in (0.1, 1.0, 5.0, 10.0, 20.0, 50.0):
    v = float(gc.transition_onset_reynolds(ti)); e = int(math.floor(math.log10(v)))
    macros["ReL" + {0.1: "a", 1.0: "b", 5.0: "c", 10.0: "d", 20.0: "e", 50.0: "f"}[ti]] = f"\\ensuremath{{{v / 10 ** e:.3g}\\times10^{{{e}}}}}"

# ---------------------------------------------------------------- Fig. angle of attack
b = np.linspace(0, 90, 181)
gW, gB = gc.angle_of_attack_factors(b)
fe = gc.ANGLE_PARAMS["f_e"]
fig, ax = plt.subplots(figsize=(3.4, 2.5))
ax.plot(b, fe * gW, color=C2, label=r"windward $f_e\,g_W$")
ax.plot(b, fe * gB, color=C1, label=r"leeward $f_e\,g_B$")
ax.plot(b, fe * 0.5 * (gW + gB), color="0.4", ls="--", label="two-sided mean")
ax.axhline(1, color="0.6", lw=0.6)
ax.set_xlabel(r"angle of attack $\beta$ (deg)"); ax.set_ylabel(r"$\overline{Nu}\,/\,\overline{Nu}_{\rm plate}$")
ax.legend(frameon=False)
save(fig, "angle")
i = int(np.argmax(gB)); macros["GBmax"] = f"{gB[i]:.3f}"; macros["GBmaxAt"] = f"{b[i]:.1f}"

# ---------------------------------------------------------------- Fig. general model: conductances vs wind speed
u = np.concatenate([[0.0], np.logspace(-2, 1, 120)])
fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.5))
r = gc.leaf_conductances(TA + 3, TA, u, 0.08, 0.04, 60.0, wind_from_deg=270.0, TI_pct=20.0)
for key, col, lab in (("gb", C1, "heat"), ("gv", C3, "H$_2$O"), ("gc", C2, "CO$_2$")):
    axs[0].plot(u[1:], r[f"{key}_U"][1:] * 1e3, color=col, label=f"{lab}, upper")
    axs[0].plot(u[1:], r[f"{key}_L"][1:] * 1e3, ":", color=col, label=f"{lab}, lower")
axs[0].set_xscale("log"); axs[0].set_xlabel(r"wind speed $u$ (m s$^{-1}$)"); axs[0].set_ylabel(r"$g$ (mm s$^{-1}$)")
axs[0].legend(frameon=False, fontsize=7, ncol=2); axs[0].set_title("(a) conductances, 8 cm $\\times$ 4 cm, 60$^\\circ$", fontsize=9, loc="left")
axs[1].plot(u[1:], r["k_U"][1:], color=C1, label="$k$, upper")
axs[1].plot(u[1:], r["k_L"][1:], ":", color=C1, label="$k$, lower")
axs[1].plot(u[1:], r["free_share_U"][1:], color=C3, label="free share, upper")
axs[1].plot(u[1:], float(r["cond_share_U"]) * np.ones_like(u[1:]), color="0.5", ls="--", label="conduction share (laminar free)")
axs[1].set_xscale("log"); axs[1].set_xlabel(r"wind speed $u$ (m s$^{-1}$)"); axs[1].set_ylabel("exponent, share")
axs[1].set_ylim(0, 1.02); axs[1].legend(frameon=False, fontsize=7)
axs[1].set_title("(b) effective exponent and mode shares", fontsize=9, loc="left")
save(fig, "general")

# ---------------------------------------------------------------- worked example table
rows = []
f = gc.free_convection(TA + 3, TA, 0.08, 0.04, 60.0)
g = f["geometry"]
ex = {"S": g.S * 1e4, "DI": g.DI, "Lh": g.Lh * 1e3, "LCo": g.LCo * 1e3, "Li": g.Li * 1e3, "phi": math.degrees(g.phi),
      "theta": math.degrees(g.theta), "Ra": float(f["Ra"]), "Nuc": float(f["Nuc"]), "NuFU": float(f["NuF_U"]),
      "NuFL": float(f["NuF_L"]), "gbU": float(f["gb_U"]) * 1e3, "gbL": float(f["gb_L"]) * 1e3,
      "HU": float(f["H_U"]), "HL": float(f["H_L"]), "nu": float(f["nu"]) * 1e5, "kappa": float(f["kappa"]) * 1e5,
      "rho": float(f["rho"])}
m = gc.mixed_convection(TA + 3, TA, 1.0, 0.08, 0.04, 60.0, 0.0, 180.0, wind_from_deg=270.0)
ex.update({"Lf": m["Lf"] * 1e3, "phiw": math.degrees(m["phi_w"]), "Ref": f"{float(m['Re_f']):.0f}",
           "gbmixU": float(m["gb_U"]) * 1e3, "gbmixL": float(m["gb_L"]) * 1e3, "fshare": float(m["forced_share"])})
lc0 = gc.leaf_conductances(TA + 3, TA, 0.0, 0.08, 0.04, 60.0, wind_from_deg=270.0, TI_pct=20.0)
lc3 = gc.leaf_conductances(TA + 3, TA, 3.0, 0.08, 0.04, 60.0, wind_from_deg=270.0, TI_pct=20.0)
for tag, lc in (("Still", lc0), ("Wind", lc3)):
    for s in ("U", "L"):
        ex[f"gb{tag}{s}"] = float(lc[f"gb_{s}"]) * 1e3
        ex[f"gv{tag}{s}"] = float(lc[f"gv_{s}"]) * 1e3
        ex[f"gc{tag}{s}"] = float(lc[f"gc_{s}"]) * 1e3
        ex[f"k{tag}{s}"] = float(lc[f"k_{s}"])
    ex[f"att{tag}"] = float(lc["attack_deg"]); ex[f"ww{tag}"] = lc["windward_side"]
ex["regimeWind"] = str(lc3["forced_regime"])
ex["RefWind"] = f"{float(lc3['Re_f']):.0f}"


def fmt(v, nd=3):
    if isinstance(v, str):
        return v
    if abs(v) >= 1e4 or (abs(v) < 1e-2 and v != 0):
        mant, exp = f"{v:.{nd}e}".split("e")
        return f"{mant}\\times10^{{{int(exp)}}}"
    return f"{v:.{nd}f}"


macro_lines = [f"\\newcommand{{\\Ex{k}}}{{{fmt(v)}}}" for k, v in ex.items()]
macro_lines += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in macros.items()]

# ---------------------------------------------------------------- reference case statistics
ref = json.loads((ROOT / "reference" / "reference_values.json").read_text(encoding="utf-8"))
counts = {k: len(ref[k]) for k in ("cases", "mixed_cases", "angle_cases", "scalar_cases")}
macro_lines += [f"\\newcommand{{\\NRef{k.replace('_', '').replace('cases', 'Cases')}}}{{{v}}}" for k, v in counts.items()]
macro_lines.append(f"\\newcommand{{\\NRefTotal}}{{{sum(counts.values())}}}")
regimes = {}
for c in ref["angle_cases"]:
    regimes[c["forced_regime"]] = regimes.get(c["forced_regime"], 0) + 1
macro_lines.append("\\newcommand{\\AngleRegimes}{" + ", ".join(f"{v} {k}" for k, v in sorted(regimes.items())) + "}")

# Python recomputation of the reference values (regression check)
dev = {}
d = 0.0
for c in ref["cases"]:
    r = gc.free_convection(c["TL"], c["TA"], c["L"], c["W"], c["alpha_x"], c["alpha_y"], eL=c.get("eL"), eA=c.get("eA"))
    for k in ("gb_U", "gb_L", "H_U", "H_L", "NuF_U", "NuF_L", "Nuc", "Ra"):
        d = max(d, abs(float(r[k]) - c[k]) / max(1.0, abs(c[k])))
dev["free"] = d
d = 0.0
for c in ref["mixed_cases"]:
    r = gc.mixed_convection(c["TL"], c["TA"], c["u"], c["L"], c["W"], c["alpha_x"], c["alpha_y"], c["leaf_azimuth"], wind_from_deg=c["wind_from"])
    for k in ("Lf", "NuF_forced", "NuF_mix_U", "NuF_mix_L", "gb_U", "gb_L", "H_U", "H_L", "forced_share"):
        d = max(d, abs(float(r[k]) - c[k]) / max(1.0, abs(c[k])))
dev["mixed"] = d
d = 0.0
for c in ref["angle_cases"]:
    r = gc.mixed_convection(296.15, 293.15, c["u"], 0.08, 0.04, c["alpha_x"], c["alpha_y"], c["leaf_azimuth"],
                            wind_from_deg=c["wind_from"], transition=True, TI_pct=c["TI"], angle_model=True)
    for k in ("attack_deg", "Re_f", "NuF_forced_U", "NuF_forced_L", "gb_U", "gb_L", "H_U", "H_L"):
        d = max(d, abs(float(r[k]) - c[k]) / max(1.0, abs(c[k])))
dev["angle"] = d
d = 0.0
for c in ref["scalar_cases"]:
    r = gc.leaf_conductances(293.15 + c["dT"], 293.15, c["u"], c["L"], c["W"], c["alpha_x"], 0.0, 180.0,
                             wind_from_deg=c["wind_from"], TI_pct=c["TI"])
    for k in ("gb_U", "gb_L", "gv_U", "gv_L", "gc_U", "gc_L", "k_U", "k_L", "cond_share_U", "free_share_U",
              "free_turb_share_U", "forced_turb_share_U"):
        d = max(d, abs(float(r[k]) - c[k]) / max(1.0, abs(c[k])))
dev["scalar"] = d
for k, v in dev.items():
    macro_lines.append(f"\\newcommand{{\\PyDev{k.capitalize()}}}{{{fmt(v, 1) if v > 0 else '0'}}}")

# MATLAB cross-check (written by matlab_crosscheck.m, optional)
mfile = HERE / "generated" / "matlab_crosscheck.json"
if mfile.exists():
    mres = json.loads(mfile.read_text())
    for k in ("free", "mixed", "angle", "scalar"):
        v = mres[k]["maxRelDev"]
        macro_lines.append(f"\\newcommand{{\\MlDev{k.capitalize()}}}{{{fmt(v, 1) if v > 0 else '0'}}}")
        macro_lines.append(f"\\newcommand{{\\MlVals{k.capitalize()}}}{{{mres[k]['values']}}}")
    macro_lines.append(f"\\newcommand{{\\MatlabVersion}}{{{mres['matlabVersion']}}}")

macro_lines.append(f"\\newcommand{{\\PkgVersion}}{{{gc.__version__}}}")
macro_lines.append(f"\\newcommand{{\\NumpyVersion}}{{{np.__version__}}}")
macro_lines.append(f"\\newcommand{{\\PythonVersion}}{{{sys.version.split()[0]}}}")
(GEN / "values.tex").write_text("\n".join(macro_lines) + "\n", encoding="utf-8")
print("\n".join(macro_lines))
