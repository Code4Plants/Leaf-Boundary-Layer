# Changelog

## Unreleased

- License changed from GPL-3.0-or-later to LGPL-3.0-or-later (`COPYING.LESSER`, with the
  GNU GPL v3 text in `COPYING`).

## 1.1.0 (2026-10-12)

- General model `leaf_conductances` / `leafConductances.m`: free (laminar, turbulent
  Ra^1/3 branch), forced (laminar, transitional, turbulent after Lienhard 2020), angle of
  attack, conduction limit; conductances per side for heat, H2O and CO2 with the effective
  diffusivity exponent as the share-weighted mean of the mode exponents (generalised Eq. 17).
- Parsimonious angle-of-attack model (f_e = 1.161, a_W = 0.226, a_R = 0.574, beta_s = 31.2 deg,
  c_w = 0.833) fitted to 815 thin-plate measurements (rRMSE 24 %, bias -4 %).
- Laminar-turbulent transition with turbulence-dependent onset (Mayle 1991); uniform-heat-flux
  average coefficient 0.680.
- Flowline length for laminar forced convection (exact Eq. 4, n = 0.5) along the wind
  direction projected onto the leaf plane; mixed convection by cubic blending.
- 573 Python/MATLAB cross-check cases, pyproject for `pip install`.

## 1.0.0 (2026-10-12)

- Free convection of inclined elliptical leaves after Graefe, Grosch & Bitterlich (2022),
  Python and MATLAB, with tests.
