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
Free-convection checks: Hassani & Hollands (1989), Kobus & Wedekind (2001, 2002).
Model ingredients: Lienhard (2020), Mayle (1991), Kottke, Blenke & Schmidt (1977),
Yovanovich & Jafarpur (1993), Massman (1998).

## Model

```
Nu_0  = f_e · Nu_Lienhard(Re, TI, bc)                       Lienhard (2020) transition, onset Re_l(TI) after Mayle (1991)
windward: Nu_W = Nu_0 (1 + a_W sin β)
leeward:  Nu_B = Nu_0 [w (1 + a_R sin 2β) + (1 − w) c_w],    w = exp(−(β/β_s)²)
two-sided measurements: mean of both sides; conduction limit Nu_c blended (Yovanovich & Jafarpur 1994)
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

Geometric mean of predicted/observed by angle-of-attack bin (windward / leeward / both):
0–15° 1.04 / 1.08 / 0.98; 15–30° 1.04 / 1.02 / 0.85; 30–45° 0.96 / 0.99 / 1.09;
45–60° 0.97 / 0.97 / 1.07; 60–75° 0.98 / 0.96 / 1.06; 75–90° 0.90 / 1.04 / 0.98.

## Notes

- Blunt bodies (Re·t/L ≥ 2000) are not described by the model: the assumption "turbulent
  from the leading edge" overpredicts thick disks by up to a factor of two.
- The free-convection data of Hassani & Hollands (1989, disks and squares, Ra 18–7.6e7) and
  Kobus & Wedekind (2002, inclined 2 cm disks) are reproduced by the laminar model within
  −20 % / −10 % (edges carry no transfer in the model); the Ra^(1/3) branch reduces the
  Hassani deviation to −14 %.
- No measured mixed-convection data (simultaneous Re and Gr) were available; the mixed
  model is checked for its limits only. Kobus & Wedekind (1996, vertical plate, assisting
  and opposing flow, Ri 0.01–1000) is the candidate data set.
