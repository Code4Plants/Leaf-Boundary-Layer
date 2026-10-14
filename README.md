# leaf-boundary-layer

Python and MATLAB implementations of the boundary-layer conductance of leaves under
free, forced and mixed convection, per leaf side and for heat, water vapour and CO2.
The free-convection core is the model of

> Graefe, J., Grosch, R., Bitterlich, M. (2022). The boundary layer conductance of
> inclined elliptical leaves under free convection. *Agricultural and Forest
> Meteorology* 317, 108884. https://doi.org/10.1016/j.agrformet.2022.108884 (CC BY 4.0)

extended (sections marked *Extension*) by the flowline length for laminar forced
convection, the laminar-turbulent transition (Lienhard 2020), a parsimonious
angle-of-attack model fitted to 815 literature measurements, mixed convection and the
mode-dependent conversion to scalar conductances.

Both implementations share the same formulation and are cross-checked against a
common set of 573 reference cases (`reference/reference_values.json`, relative
tolerance 1e-9). Install: `pip install .` (module `graefe_convection`) or add `matlab/`
to the MATLAB path. Tests: `python -m pytest python -q`, `runTests` in MATLAB.

Quick start (general model):

```python
from graefe_convection import leaf_conductances
r = leaf_conductances(TL=296.15, TA=293.15, u=1.5, L=0.08, W=0.04, alpha_x_deg=60,
                      leaf_azimuth_deg=180, wind_from_deg=270, TI_pct=20)
r["gb_U"], r["gv_U"], r["gc_U"], r["k_U"]      # heat, H2O, CO2 conductance (m/s) and exponent, upper side
```

```matlab
r = leafConductances(296.15, 293.15, 1.5, 0.08, 0.04, 60, 0, 180, 270, [], [], [], 20);
```

## Model

For a flat elliptical leaf of length *L* and width *W* (area *S* = π/4 *L W*,
aspect ratio *AR* = *L*/*W*) rotated about its long axis by *α_y* and about its
short axis by *α_x* (inclination from the horizontal), the leaf-normal angle to the
vertical is cos *θ* = cos *α_x* cos *α_y* (Eq. 11). The characteristic lengths are

* horizontal leaf: *L_h* = *S*/*P* = 0.5 √*S* / (*DI* √π), with the dissection index
  *DI* ≈ 1.112 [(*L*+*W*)/√(*L W*)]^0.783 − 0.913 (Eq. 7, 8);
* inclined leaf: *L_i* = (0.04511 *n* + 1.006) *L_Co*(*φ*), *n* = 0.75, where
  *L_Co*(*φ*) = π *L W* √(tan²*φ* + 1) / (4 √(tan²*φ* *W*² + *L*²)) is the mean chord
  length along the steepest-ascent direction, which subtends the angle *φ* with the
  short axis (Eq. 9, 10; Appendix A/B). For *α_y* = 0 the flow runs along the long
  axis and *L_i* = 0.816 *L*.

Flux-proportional Nusselt numbers (m⁻¹) for the upper (U) and lower (L) side of a
heated leaf (Eq. 3a–3g, parameters of Table 3: p1 = 6, p2 = 1.42, p3 = 0.513,
p4 = 0.4, p5 = 0.202, p6 = 0.492, Pr = 0.71):

```
Ra      = g |ΔT| Pr / (T_A ν²)                      (per m³, ν at film temperature)
Nu_i,U  = p3 L_i^-0.25 sin(θ)^p4  Ra^0.25
Nu_i,L  = p3 L_i^-0.25 sin(θ)^0.25 Ra^0.25
Nu_h,U  = p5 L_h^-0.25 cos(θ)^0.25 Ra^0.25
Nu_h,L  = p6 L_h^-0.4  cos(θ)^0.2  Ra^0.2
Nu_c    = 5.3 sin(0.0338 AR + 0.607) / √S
NuF     = { [Nu_i^p1 + Nu_h^p1]^(p2/p1) + Nu_c^p2 }^(1/p2)
```

For a cooled leaf (ΔT < 0) the right-hand sides of U and L are exchanged. The
boundary-layer conductance is *g_b* = κ NuF (Eq. 1) and the sensible heat flux per
side *H* = *g_b* ρ c_p (T_L − T_A) (Eq. 14). Optionally the buoyancy term uses the
virtual difference ΔT' = ΔT + 0.38 (e_L T_L − e_A T_A)/p (Eq. 18).

Air properties (ρ, ν via Sutherland, κ = k/(ρ c_p) with k = 0.02624 (T/300)^0.8646)
are evaluated at the film temperature (T_L + T_A)/2.

### Exact characteristic length (Eq. 4)

For an ellipse the chords parallel to a streamline direction vary as
*D*(t) = *D*₀ √(1 − (t/t_max)²) with *D*₀ = (4/π) *L_Co*(φ), so the flux-weighted
length of Eq. 4 has the closed form

```
Lc^(n-1) = D0^(n-1) · (2/π) · B(1/2, n/2 + 1)        (B = Beta function)
```

`characteristic_length(L, W, phi, n)` / `characteristicLength.m` implement it.
For *n* = 0.75 it reproduces the linear fit of Eq. 9 within 0.03 % (L_c/L_Co depends only on *n*, not on φ or the aspect ratio).

### Extension: laminar forced and mixed convection (flowline concept)

Not part of the paper, added as a feature. Forced convection over a leaf side is
the laminar flat-plate relation NuF = c Re^0.5 Pr^(1/3) / *L_f*, Re = *u L_f*/ν,
c = 0.664 (Pohlhausen; Campbell & Norman use 0.6), where the streamline length
*L_f* follows the flowline concept: the horizontal wind vector (meteorological
direction `wind_from`) is projected onto the leaf plane, its angle φ_w to the short
axis gives *L_f* = characteristic_length(L, W, φ_w, n = 0.5) because *h* ∝ *L*^−0.5.
The leaf azimuth is the compass direction the upper surface faces (180 = south);
the long axis runs up-slope towards azimuth + 180. With unknown wind direction
(or wind normal to the leaf) the orientation-averaged chord length π *S*/*P*
(Eq. 5) is used. Free and forced Nusselt numbers are blended per side with
NuF = (NuF_free³ + NuF_forced³)^(1/3); assisting/opposing flows are not
distinguished.

```python
m = mixed_convection(TL, TA, u=1.5, L=0.08, W=0.04, alpha_x_deg=60, leaf_azimuth_deg=180, wind_from_deg=270)
m["gb_U"], m["gb_L"], m["Lf"], m["forced_share"]
```

```matlab
m = graefeMixedConvection(TL, TA, 1.5, 0.08, 0.04, 60, 0, 180, 270);
```

### Extension: laminar-turbulent transition

Forced convection: `forced_nusselt_lienhard(Re, Pr, TI_pct, bc)` /
`forcedNusseltLienhard.m` implement the average flat-plate Nusselt number of Lienhard
(2020, *J. Heat Transfer* 142:061805): laminar Nu_x = a_l Re_x^0.5 Pr^(1/3)
(a_l = 0.332 isothermal, 0.453 uniform heat flux), transitional
Nu_x = Nu_x(Re_l)(Re_x/Re_l)^c with c = 0.9922 log10 Re_l − 3.013, turbulent
Nu_x = 0.0296 Re_x^0.8 Pr^0.6, end of transition Re_u from continuity, exact integral
over the plate. The onset Re_l follows the free-stream turbulence intensity TI (%)
after Mayle (1991): Re_l = 3.63e5 TI^−1.25, bounded to [2e3, 5e5]
(`transition_onset_reynolds`). Blunt leading edges with Re·t/L ≥ 2000 (Kottke et al.
1977) are turbulent from the edge: `forced_nusselt_turbulent`. In
`mixed_convection(..., transition=True, TI_pct=..)` / `graefeMixedConvection(..., true, TIpct)`
the forced part uses this model and `forced_regime` reports laminar / transitional /
turbulent.

### Extension: angle of attack and finite-plate enhancement (parsimonious leaf model)

Five parameters fitted to 815 thin-plate / leaf measurements (Chen 1988, Haseba 1973/84,
Sparrow 1977/79, Parkhurst 1968, Powell 1940, Schuepp 1993, Ramirez 2002, Kobus 2006;
Re 1e3–1.6e5, angle of attack 0–90°, both faces; rRMSE 24 %, bias −4 %;
see `docs/validation.md`):

```
Nu_0  = f_e · Nu_Lienhard(Re, TI, bc)                          f_e = 1.161 (edges, 3-D flow)
windward: Nu_W = Nu_0 (1 + a_W sin β)                          a_W = 0.226 (accelerated outer flow)
leeward:  Nu_B = Nu_0 [w (1 + a_R sin 2β) + (1 − w) c_w],       a_R = 0.574, c_w = 0.833
          w = exp(−(β/β_s)²), β_s = 31.2°                      (reattaching bubble → wake)
```

β is the angle of attack (0 = parallel, 90 = normal). `angle_of_attack_factors`,
`forced_nusselt_leaf`, `wind_attack_angle` / `angleOfAttackFactors.m`,
`forcedNusseltLeaf.m`, `windAttackAngle.m`; in `mixed_convection(..., angle_model=True)` /
`graefeMixedConvection(..., angleModel)` the windward side and β follow from the wind
direction and the leaf orientation. Blunt bodies (Re·t/L ≥ 2000) are outside the fit.

### General model: heat, water vapour and CO2 conductances per leaf side

`leaf_conductances(...)` / `leafConductances.m` combine everything (free laminar and
turbulent, forced laminar / transitional / turbulent with the Lienhard transition,
angle of attack, conduction limit) and return per side X ∈ {U, L}: `gb_X` (heat),
`gv_X` (H2O), `gc_X` (CO2), the effective diffusivity exponent `k_X` and the mode shares
(`cond_share_X`, `free_share_X`, `free_turb_share_X`, `forced_turb_share_X`).

The scalar conductances follow g_s = g_heat (D_s/κ)^k. Each transport mode has its own
exponent k = 1 − (Pr exponent of its law): conduction 1, laminar free convection 0.75
(Ra^1/4), turbulent free convection 2/3 (Ra^1/3), forced laminar / transitional 2/3
(Pr^1/3), forced turbulent 0.4 (Pr^0.6). The effective k of a side is the share-weighted
mean, generalising Eq. 17 of the paper (laminar free value raised towards 1 by the
conduction share Nu_c/NuF): k_free = (1 − s_t)[0.75 + 0.25 Nu_c/NuF_lam] + s_t·2/3,
k_forced = (1 − f_t)·2/3 + f_t·0.4, k = (1 − s_f) k_free + s_f k_forced with the
power-mean shares s_t (free turbulent), f_t (forced turbulent, additive Lienhard parts)
and s_f (forced, cubic blending). Diffusivities D_H2O = 2.12e-5, D_CO2 = 1.39e-5 m²/s at
0 °C scaled with (T/273.15)^1.75 (p0/p). Example (8×4 cm leaf, 60°, +3 K, 20 °C):
still air k = 0.85, g_heat 3.5 / g_H2O 3.8 / g_CO2 2.7 mm/s; 3 m/s wind k = 0.67.

```python
r = leaf_conductances(TL, TA, u, L, W, alpha_x_deg, wind_from_deg=270, TI_pct=20)
r["gb_U"], r["gv_U"], r["gc_U"], r["k_U"]
```

```matlab
r = leafConductances(TL, TA, u, L, W, alphaX, 0, 180, 270, [], [], [], 20);
```

Free convection: `free_convection_transition()` adds a turbulent branch
NuF = 0.14 Ra^(1/3) on the buoyantly unstable side (and both sides of a vertical
leaf), blended with the laminar Graefe model as (NuF_lam^6 + NuF_turb^6)^(1/6); the
stably stratified side stays laminar. Both extensions are validated against the
literature data set described in `docs/validation.md`.

## Python

```python
from graefe_convection import free_convection, leaf_geometry
r = free_convection(TL=296.15, TA=293.15, L=0.08, W=0.04, alpha_x_deg=60)
r["gb_U"], r["gb_L"], r["H_U"], r["H_L"]      # m/s and W/m2 per side
```

Requires numpy (tests: pytest). Run `python -m pytest python -q`.
`python/make_reference.py` regenerates the reference values.

## MATLAB

```matlab
addpath matlab
r = graefeFreeConvection(296.15, 293.15, 0.08, 0.04, 60);   % struct with gb_U, gb_L, H_U, H_L, ...
g = leafGeometry(0.08, 0.04, 60);                            % Lh, Li, phi, theta, ...
runTests                                                     % matlab.unittest suite incl. Python cross-check
exampleFig6                                                  % H over inclination (cf. Fig. 6b)
```

Inputs are elementwise over arrays of leaf/air temperatures; the geometry is scalar.

## Files

```
python/graefe_convection.py        model (single file, numpy): free, forced, mixed convection, transition, angle of attack, scalar conductances
python/test_graefe_convection.py   pytest suite
python/make_reference.py           writes reference/reference_values.json
matlab/graefeFreeConvection.m      free convection (paper)
matlab/graefeMixedConvection.m     forced + mixed convection extension
matlab/leafConductances.m          general model: heat / H2O / CO2 per side with mode-dependent exponents
matlab/forcedNusseltLienhard.m, forcedNusseltLeaf.m, angleOfAttackFactors.m, windAttackAngle.m, freeConvectionTransition.m
matlab/leafGeometry.m, airProperties.m, characteristicLength.m, windStreamlineAngle.m
matlab/tests/testGraefeFreeConvection.m, matlab/runTests.m, matlab/exampleFig6.m
reference/reference_values.json    cross-check cases (free + mixed)
```

## References

`docs/references.md` lists all sources of the model and of the validation data with
verified DOIs (`docs/references.bib` for BibTeX); `docs/validation.md` describes the fit.

## License

GPL-3.0 (see LICENSE). Please cite the paper above when using the model
(`CITATION.cff`).
