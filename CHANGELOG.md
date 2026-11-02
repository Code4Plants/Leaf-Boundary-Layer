# Changelog

## 1.2.0 (2026-11-02)

- Turbulent branch of free convection removed. Leaves with ΔT ≤ 10 K and L ≤ 30 cm stay two
  orders of magnitude below the turbulent onset (Ra_Lh ≤ ~1e5 against ~1e7 horizontal,
  Ra_Li ≤ ~1.5e7 against ~1e9 vertical), so the laminar model of Graefe et al. (2022) covers
  the leaf range. The former 0.14 Ra^(1/3) branch had no Rayleigh onset and, being on by
  default in the general model, raised the upper-side conductance of horizontal leaves within
  the laminar range (8 × 4 cm, 3 K: +4 %; 20 × 10 cm, 3 K: +30 %), contradicting the two-sided
  plate data of Ryms et al. (2025). In still air the general model now equals the laminar
  free-convection model.
- API: `free_convection_transition()` / `freeConvectionTransition.m` removed; parameter
  `free_transition` / `freeTransition` of `mixed_convection`, `leaf_conductances` /
  `graefeMixedConvection.m`, `leafConductances.m` removed (`angleModel` is now the 16th
  argument of `leafConductances.m`); output fields `free_turb_share_X`, `NuF_lam_X`,
  `NuF_turb` and the constants `C_FREE_TURB`, `M_FREE_TURB`, `K_FREE_TURB` removed;
  `scalar_exponent(nuc, nuf_free, forced_turb_share, forced_share)` has four parameters.
- Reference values (`scalar_cases`) regenerated with the laminar free convection.
- Validation of the laminar free-convection core against Ryms, Tesch & Lewandowski (2025,
  horizontal 150 × 75 mm plate heated on both sides): the two-sided CFD correlations are met
  within 4 %, the upper side is weaker than the lower side as observed; new regression tests
  and `python/validate_ryms2025.py` (`docs/validation.md`).
- License changed from GPL-3.0-or-later to LGPL-2.1-or-later (license text in `COPYING`).

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
