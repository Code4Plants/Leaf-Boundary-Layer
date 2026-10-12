function r = freeConvectionTransition(TL, TA, L, W, alphaXdeg, alphaYdeg, p, eL, eA, cTurb, m)
%FREECONVECTIONTRANSITION graefeFreeConvection with a laminar-turbulent branch.
%   On the buoyantly unstable side (upper side of a heated leaf, lower side of a
%   cooled leaf, both sides of a vertical leaf) NuF = (NuF_lam^m + NuF_turb^m)^(1/m)
%   with NuF_turb = cTurb Ra^(1/3) (default 0.14, m = 6); the stable side stays laminar.
%   Extension, not part of Graefe et al. (2022).  Adds NuF_lam_U/L, NuF_turb and the
%   turbulent shares free_turb_share_U/L to the output struct.
if nargin < 6 || isempty(alphaYdeg), alphaYdeg = 0; end
if nargin < 7 || isempty(p), p = 101325.0; end
if nargin < 9, eL = []; eA = []; end
if nargin < 10 || isempty(cTurb), cTurb = 0.14; end
if nargin < 11 || isempty(m), m = 6; end
r = graefeFreeConvection(TL, TA, L, W, alphaXdeg, alphaYdeg, p, eL, eA);
turb = cTurb .* r.Ra.^(1/3);
heated = r.dT >= 0;
vertical = abs(sin(r.geometry.theta) - 1) < 1e-9;
unstU = heated | vertical; unstL = (~heated) | vertical;
r.NuF_lam_U = r.NuF_U; r.NuF_lam_L = r.NuF_L; r.NuF_turb = turb;
r.free_turb_share_U = unstU .* turb.^m ./ (r.NuF_U.^m + turb.^m);
r.free_turb_share_L = unstL .* turb.^m ./ (r.NuF_L.^m + turb.^m);
r.NuF_U = unstU .* (r.NuF_U.^m + turb.^m).^(1/m) + (~unstU) .* r.NuF_U;
r.NuF_L = unstL .* (r.NuF_L.^m + turb.^m).^(1/m) + (~unstL) .* r.NuF_L;
r.gb_U = r.kappa .* r.NuF_U; r.gb_L = r.kappa .* r.NuF_L;
cp = 1004.0;
r.H_U = r.rho .* cp .* r.gb_U .* (TL - TA);
r.H_L = r.rho .* cp .* r.gb_L .* (TL - TA);
end
