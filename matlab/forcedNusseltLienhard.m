function [Nu, regime, ReL, ReU, lam, trans, turb] = forcedNusseltLienhard(Re, Pr, TIpct, bc, unheatedFraction, ReL)
%FORCEDNUSSELTLIENHARD Average flat-plate Nusselt number with laminar-turbulent transition.
%   [Nu, regime, ReL, ReU] = forcedNusseltLienhard(Re, Pr, TIpct, bc, unheatedFraction, ReL)
%
%   Lienhard (2020, J. Heat Transfer 142:061805): local laws laminar
%   Nu_x = a_l Re_x^0.5 Pr^(1/3) (a_l = 0.332 isothermal, 0.453 uniform heat flux),
%   transitional Nu_x = Nu_x(ReL) (Re_x/ReL)^c with c = 0.9922 log10(ReL) - 3.013,
%   turbulent Nu_x = 0.0296 Re_x^0.8 Pr^0.6; ReU from continuity of Nu_x; the
%   average is the exact integral of the three pieces.  Transition onset ReL from
%   the turbulence intensity (Mayle 1991): ReL = 3.63e5 TI^-1.25, bounded to
%   [2e3, 5e5], unless ReL is given.  Elementwise over Re.
%   bc: 'isothermal' (default) or 'uniform_flux'; unheatedFraction = x0/L.
if nargin < 2 || isempty(Pr), Pr = 0.71; end
if nargin < 3 || isempty(TIpct), TIpct = 1.0; end
if nargin < 4 || isempty(bc), bc = 'isothermal'; end
if nargin < 5 || isempty(unheatedFraction), unheatedFraction = 0; end
if nargin < 6 || isempty(ReL)
    TI = max(TIpct, 1e-3);
    ReL = min(max(3.63e5 .* TI.^(-1.25), 2e3), 5e5);
end
aT = 0.0296; aL = 0.332;
% uniform heat flux: average Nu with the mean temperature difference is 0.680 (Lienhard 2020)
if strcmp(bc, 'uniform_flux'), hfa = 0.680 / 0.664; else, hfa = 1; end
if unheatedFraction > 0, fSL = (1 - unheatedFraction^0.75)^(-1/3); else, fSL = 1; end
c = 0.9922 .* log10(ReL) - 3.013;
ReU = ((ReL.^(0.5 - c) .* aL ./ aT) .* Pr^(1/3 - 0.6)).^(1 ./ (0.8 - c));
ReU = max(ReU, ReL * 1.0001);
lam = 2 * aL * fSL * hfa * Pr^(1/3) .* min(Re, ReL).^0.5;
nuLocalL = aL * fSL * hfa .* ReL.^0.5 .* Pr^(1/3);
reMid = min(max(Re, ReL), ReU);
trans = (Re > ReL) .* nuLocalL ./ c .* ((reMid ./ ReL).^c - 1);
turb = (Re > ReU) .* 0.037 * Pr^0.6 .* max(Re.^0.8 - ReU.^0.8, 0);
Nu = lam + trans + turb;   % the three additive parts are also returned
regime = repmat({'laminar'}, size(Re));
regime(Re > ReL) = {'transitional'};
regime(Re > ReU) = {'turbulent'};
end
