function [gW, gB, prm] = angleOfAttackFactors(betaDeg, prm)
%ANGLEOFATTACKFACTORS Windward / leeward multipliers of the parallel-flow Nusselt number.
%   [gW, gB] = angleOfAttackFactors(betaDeg) for the angle of attack beta (deg,
%   0 = parallel flow, 90 = normal to the flow) of a thin plate or leaf:
%       gW = 1 + a_W sin(beta)                                   accelerated outer flow
%       gB = w (1 + a_R sin(2 beta)) + (1 - w) c_w, w = exp(-(beta/beta_s)^2)
%            reattaching separation bubble -> wake regime (c_w < 1)
%   The finite-plate enhancement f_e (prm.f_e) multiplies the base Nusselt number
%   (see forcedNusseltLeaf).  Parameters fitted 2026-10-12 to 815 thin-plate
%   measurements (see docs/validation.md).
if nargin < 2 || isempty(prm)
    prm = struct('f_e', 1.161, 'a_W', 0.226, 'a_R', 0.574, 'beta_s_deg', 31.2, 'c_w', 0.833);
end
b = deg2rad(betaDeg);
gW = 1 + prm.a_W .* sin(b);
w = exp(-(betaDeg ./ prm.beta_s_deg).^2);
gB = w .* (1 + prm.a_R .* sin(2 * b)) + (1 - w) .* prm.c_w;
end
