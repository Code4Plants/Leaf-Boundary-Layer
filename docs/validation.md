# Validation of the forced-convection extensions

The free-convection core (Graefe, Grosch & Bitterlich 2022) was calibrated and validated
in the paper against CFD simulations and the measurements of Hassani & Hollands (1989).
The extensions of this repository (transition, angle of attack, finite-plate enhancement)
were fitted and validated on 2026-10-12 against a data set of 1169 digitised literature
measurements of mean Nusselt numbers of finite flat plates, disks and leaf models in
forced convection (fit set: 815 rows of thin plates with Re·t/L < 2000 and wind-tunnel
blockage ≤ 0.12; the thicker bodies of aus der Wiesche & Helcig 2016 and Motwani et al.
1985 are outside the fit).

Sources (full references in `references.md` / `references.bib`): Chen, Ibbetson & Milford
(1988a, b), Haseba (1973), Haseba & Ito (1984), Sparrow & Tien (1977), Sparrow, Ramsey &
Mass (1979), Parkhurst et al. (1968), Powell (1940), Schuepp (1993), Ramirez, Murray &
Fitzpatrick (2002), Kobus & Shumway (2006); blunt bodies from Motwani, Gaitonde & Sukhatme
(1985) and aus der Wiesche & Helcig (2016); Re 1.3e3–1.6e5, angle of attack 0–90°,
windward, leeward and two-sided measurements, turbulence intensities 0.2–22 %, heat and
mass (naphthalene, evaporation, electrochemical) transfer converted to Nu at Pr = 0.71.
Free-convection checks: Hassani & Hollands (1989), Kobus & Wedekind (2001, 2002), Ryms, Tesch &
Lewandowski (2025, horizontal plate heated on both sides; see the last section).
Model ingredients: Lienhard (2020), Mayle (1991), Kottke, Blenke & Schmidt (1977),
Yovanovich & Jafarpur (1993), Massman (1998).

## Model

```
Nu_0  = f_e · Nu_Lienhard(Re, TI, bc)                       Lienhard (2020) transition, onset Re_l(TI) after Mayle (1991)
windward: Nu_W = Nu_0 (1 + a_W sin β)
leeward:  Nu_B = Nu_0 [w (1 + a_R sin 2β) + (1 − w) c_w],    w = exp(−(β/β_s)²)
two-sided measurements: mean of both sides; conduction limit Nu_c blended (Yovanovich & Jafarpur 1993)
```

| Parameter | Value | Std. error | Meaning |
| --- | --- | --- | --- |
| f_e | 1.161 | 0.016 | finite-plate enhancement (edges, three-dimensional flow) |
| a_W | 0.226 | 0.029 | windward enhancement by the accelerated outer flow |
| a_R | 0.574 | 0.095 | leeward enhancement by the reattaching separation bubble |
| β_s | 31.2° | 2.6° | transition from the reattaching to the wake regime |
| c_w | 0.833 | 0.021 | wake-regime level relative to parallel flow |

Fit: least squares on log(Nu_pred/Nu_obs); rRMSE 23.9 % (plain laminar plate without
enhancement and angle terms: 37.3 %), bias −3.8 % (−21.9 %), log-sd 0.19. Parallel flow
only (143 rows): bias −1.9 % (−15.3 % without f_e), rRMSE 11.4 % (20.2 %).

| Source | Face | n | rRMSE % | Bias % |
| --- | --- | --- | --- | --- |
| Chen 1988a | both | 42 | 7.8 | −2.7 |
| Chen 1988b | both / windward / leeward | 99 / 39 / 33 | 16.9 / 15.1 / 18.8 | −6.5 / +14.3 / −1.6 |
| Chen et al. 1988a, b (second blocks) | both | 40 / 36 | 11.0 / 15.8 | −6.4 / −10.0 |
| Haseba 1973 | windward / leeward | 70 / 78 | 11.3 / 13.2 | −10.6 / −5.4 |
| Haseba 1984 | windward / leeward | 51 / 34 | 26.4 / 30.8 | −5.7 / +2.4 |
| Sparrow 1977 | windward | 34 | 7.1 | +5.3 |
| Sparrow 1979 | windward | 33 | 22.0 | −10.4 |
| Parkhurst 1968 | both | 36 | 15.6 | +9.5 |
| Powell 1940 | windward / leeward | 12 / 52 | 32.6 / 35.3 | +28.2 / −11.7 |
| Schuepp 1993 | leeward | 72 | 30.0 | −1.9 |
| Ramirez 2002 (thin subset) | windward / leeward | 8 / 7 | 17.4 / 34.5 | −14.4 / −5.3 |
| Kobus 2006 | both | 17 | 36.3 | +19.8 |

The rows listed above sum to 793 of the 815 fit rows; the remaining 22 rows are not
itemised per source (the digitised data set is not part of this repository).

Geometric mean of predicted/observed by angle-of-attack bin (windward / leeward / both):
0–15° 1.04 / 1.08 / 0.98; 15–30° 1.04 / 1.02 / 0.85; 30–45° 0.96 / 0.99 / 1.09;
45–60° 0.97 / 0.97 / 1.07; 60–75° 0.98 / 0.96 / 1.06; 75–90° 0.90 / 1.04 / 0.98.

## Notes

- Blunt bodies (Re·t/L ≥ 2000) are not described by the model: the assumption "turbulent
  from the leading edge" overpredicts thick disks by up to a factor of two.
- The free-convection data of Hassani & Hollands (1989, disks and squares, Ra 18–7.6e7) and
  Kobus & Wedekind (2002, inclined 2 cm disks) are reproduced by the laminar model within
  −20 % / −10 % (edges carry no transfer in the model). The CFD correlations of Ryms et
  al. (2025) for a horizontal plate heated on both sides are met by the laminar core
  within 4 % (see below).
- No measured mixed-convection data (simultaneous Re and Gr) were available; the mixed
  model is checked for its limits only. Kobus & Wedekind (1996, vertical plate, assisting
  and opposing flow, Ri 0.01–1000) is the candidate data set.

## Free convection: horizontal plate heated on both sides (Ryms et al. 2025)

Ryms, Tesch & Lewandowski (2025) measured free convection in air simultaneously from the upper
and lower surface of a horizontal 150 × 75 × 2.1 mm plate heated on both sides ("Type II"),
Ra_b = 1e5–1.6e6 with the width b = 75 mm as length scale (laminar range), and computed the
same plate (CFD, no turbulence model) with both sides heated (Type II) and with one side heated
(Type I). A plate heated on both sides is the situation of a horizontal leaf, so the data are
an independent check of the two-sided laminar core, whose upper-side coefficient p5 = 0.202 is
far below the single-sided literature value. Only leaf-air temperature differences up to 10 K
(the leaf range) are considered: ΔT = 3, 5, 8 and 10 K (Ra_b 1.3e5–4.2e5) and the two tabulated
points per side with ΔT ≤ 10 K (III/5 and XIV/10); the points VI/7 (10.7 / 11.1 K) and IX/14
(57 K) are left out. The model is evaluated for the equal-area ellipse of aspect ratio 2
(169 × 85 mm, L_h = 27.4 mm; rectangle A/P = 25 mm) and converted with Nu_b = NuF·b,
Ra_b = Ra·b³ (`python/validate_ryms2025.py`).

Coefficients C = Nu_b / Ra_b^0.2 (n = 1/5 imposed, Table 7 of the paper):

| Source | C_U | C_L | U/L |
| --- | --- | --- | --- |
| Experiment Type II (Eqs. 29/32) | 0.804 | 0.744 | 1.08 |
| CFD Type II, two-sided (Eqs. 41/42) | 0.601 | 0.812 | 0.74 |
| CFD Type I, one-sided (Eqs. 39/40) | 0.941 | 0.808 | 1.16 |
| Literature mean, one-sided plates (Eqs. 53/54) | 1.093 | 0.660 | 1.66 |
| Model, laminar core, ΔT = 3 / 5 / 8 / 10 K (Ra_b 1.3e5–4.2e5) | 0.577 / 0.574 / 0.574 / 0.574 | 0.828 / 0.816 / 0.806 / 0.802 | 0.70–0.72 |

The laminar core reproduces the CFD of the two-sided plate within −4 to −5 % (upper) and
−1 to +2 % (lower) and the central finding of the paper that the upper side of a plate heated on
both sides is weaker than the lower side (the rising flow from the lower side inhibits the upper
side, CFD −36 % against the one-sided plate). Against the measurements the upper side is 20 %
too low (Table 1 points 0.80 and 0.79; 29 % below the fitted correlation), the lower side
8–11 % high on the fitted correlation and 34–56 % high on the two Table 1 points (1.56, 1.34; the
measured lower side scatters strongly, fitted exponent 0.33). The authors themselves find their
upper-side measurement 34 % above their CFD and cannot resolve the difference (the radiative
correction amounts to two thirds of the heater power).

## Free convection: laminar range of leaves (no turbulent branch)

Turbulent free convection sets in at Ra_L ≈ 1e7 for upward-facing horizontal plates (L = A/P)
and at Ra_L ≈ 1e9 for vertical plates. With ΔT ≤ 10 K at 20 °C the Rayleigh number per volume is
Ra ≤ 1e9 m⁻³, so leaves stay far below the onset (heated leaf, AR = 2, L_h ≈ 0.16 L, L_i = 0.816 L):

| Leaf | ΔT (K) | Ra_Lh (horizontal) | Ra_Li (vertical) |
| --- | --- | --- | --- |
| 8 × 4 cm | 3 | 6.7e2 | 8.6e4 |
| 8 × 4 cm | 10 | 2.2e3 | 2.8e5 |
| 20 × 10 cm | 3 | 1.1e4 | 1.3e6 |
| 30 × 15 cm | 10 | 1.1e5 | 1.5e7 |
| 200 × 50 cm | 10 | 6.1e6 | 4.3e9 |

Only metre-sized leaves at ΔT ≈ 10 K approach the onset. The laminar core therefore covers the
leaf range, and version 1.2.0 removed the optional turbulent branch of versions 1.1.x
(`free_convection_transition`, 0.14 Ra^(1/3), power mean m = 6, default on in the general
model). That branch had no Rayleigh onset: its single-sided coefficient exceeded the two-sided
laminar upper-side term 0.202 Ra_Lh^0.25 from Ra_Lh ≈ 80 on (≈ 1e3–1e4 including the conduction
term) instead of from ≈ 5e6 (0.54 vs 0.15 for one-sided plates), raised the upper-side
conductance of a horizontal 8 × 4 cm leaf by 4 % at ΔT = 3 K and 11 % at 10 K and of a
20 × 10 cm leaf by 30 % at 3 K, and gave a slope of 1/3 in the range of Ryms et al. (2025),
which is laminar (CFD slope 0.17–0.20, experiment 0.27).

Regression tests `test_ryms2025_two_sided_plate` (Python) and `ryms2025TwoSidedPlate` (MATLAB)
hold the laminar core within 10 % of the two-sided CFD correlations and NuF_U < NuF_L.
