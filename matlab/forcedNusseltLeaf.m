function [NuW, NuB] = forcedNusseltLeaf(Re, betaDeg, Pr, TIpct, bc, transition)
%FORCEDNUSSELTLEAF Windward and leeward Nusselt numbers of a thin plate / leaf.
%   [NuW, NuB] = forcedNusseltLeaf(Re, betaDeg, Pr, TIpct, bc, transition)
%   = f_e * Nu_plate(Re) * [gW, gB] with Nu_plate the Lienhard (2020) average
%   (transition = true, default) or the laminar law, and the factors of
%   angleOfAttackFactors.  Dimensionless, based on the streamline length.
if nargin < 2 || isempty(betaDeg), betaDeg = 0; end
if nargin < 3 || isempty(Pr), Pr = 0.71; end
if nargin < 4 || isempty(TIpct), TIpct = 1.0; end
if nargin < 5 || isempty(bc), bc = 'isothermal'; end
if nargin < 6 || isempty(transition), transition = true; end
if transition
    base = forcedNusseltLienhard(Re, Pr, TIpct, bc);
else
    if strcmp(bc, 'uniform_flux'), cl = 0.680; else, cl = 0.664; end
    base = cl .* Re.^0.5 .* Pr^(1/3);
end
[gW, gB, prm] = angleOfAttackFactors(betaDeg);
NuW = prm.f_e .* base .* gW;
NuB = prm.f_e .* base .* gB;
end
