function r = graefeFreeConvection(TL, TA, L, W, alphaXdeg, alphaYdeg, p, eL, eA)
%GRAEFEFREECONVECTION Free-convection conductance of an inclined elliptical leaf.
%   r = graefeFreeConvection(TL, TA, L, W, alphaXdeg, alphaYdeg, p, eL, eA)
%
%   Implements Graefe, Grosch & Bitterlich (2022), "The boundary layer conductance
%   of inclined elliptical leaves under free convection", Agric. For. Meteorol.
%   317:108884, doi:10.1016/j.agrformet.2022.108884 (Eq. 1, 3a-3g, 14, 18; Table 3).
%
%   Inputs (temperatures in K, elementwise over arrays of TL/TA):
%     TL, TA      leaf and air temperature
%     L, W        leaf length and width (m), L >= W
%     alphaXdeg   inclination from horizontal about the short axis (deg)
%     alphaYdeg   pre-rotation about the long axis (deg), default 0
%     p           air pressure (Pa), default 101325
%     eL, eA      leaf-surface and air vapour pressure (Pa); if both given the
%                 buoyancy correction of Eq. 18 is applied
%
%   Output struct r: gb_U, gb_L (m s-1), H_U, H_L (W m-2 one-sided area, positive
%   leaf -> air), NuF_U, NuF_L, Nuc (m-1), Ra (m-3), dT (K), TF, rho, nu, kappa,
%   geometry (see leafGeometry).
if nargin < 6 || isempty(alphaYdeg), alphaYdeg = 0; end
if nargin < 7 || isempty(p), p = 101325.0; end
if nargin < 9, eL = []; eA = []; end
% Table 3
p1 = 6.0; p2 = 1.42; p3 = 0.513; p4 = 0.4; p5 = 0.202; p6 = 0.492;
Pr = 0.71; g0 = 9.81; cp = 1004.0;

geom = leafGeometry(L, W, alphaXdeg, alphaYdeg);
TF = 0.5 * (TL + TA);
[rho, nu, kappa] = airProperties(TF, p);
if ~isempty(eL) && ~isempty(eA)
    dT = (TL - TA) + 0.38 * (eL .* TL - eA .* TA) ./ p;      % Eq. 18
else
    dT = TL - TA;
end
dT(abs(dT) < 1e-9) = 1e-9;
Ra = g0 .* abs(dT) .* Pr ./ (TA .* nu.^2);                   % Eq. 3g
st = sin(geom.theta); ct = cos(geom.theta);
ra25 = Ra.^0.25;
nuIU = p3 * geom.Li^(-0.25) * st^p4 .* ra25;                 % Eq. 3b
nuIL = p3 * geom.Li^(-0.25) * st^0.25 .* ra25;               % Eq. 3c
nuHU = p5 * geom.Lh^(-0.25) * ct^0.25 .* ra25;               % Eq. 3d
nuHL = p6 * geom.Lh^(-0.4) * ct^0.2 .* Ra.^0.2;              % Eq. 3e
nuc = 5.3 * sin(0.0338 * geom.AR + 0.607) / sqrt(geom.S);    % Eq. 3f
cooled = dT < 0;                                             % swap sides for cooled leaves
iu = nuIU; iu(cooled) = nuIL(cooled); il = nuIL; il(cooled) = nuIU(cooled);
hu = nuHU; hu(cooled) = nuHL(cooled); hl = nuHL; hl(cooled) = nuHU(cooled);
blend = @(ni, nh) ((ni.^p1 + nh.^p1).^(p2 / p1) + nuc^p2).^(1 / p2);   % Eq. 3a
NuF_U = blend(iu, hu);
NuF_L = blend(il, hl);
gb_U = kappa .* NuF_U;                                       % Eq. 1
gb_L = kappa .* NuF_L;
H_U = rho .* cp .* gb_U .* (TL - TA);                        % Eq. 14
H_L = rho .* cp .* gb_L .* (TL - TA);
r = struct('gb_U', gb_U, 'gb_L', gb_L, 'H_U', H_U, 'H_L', H_L, 'NuF_U', NuF_U, ...
    'NuF_L', NuF_L, 'Nuc', nuc, 'Ra', Ra, 'dT', dT, 'TF', TF, 'rho', rho, 'nu', nu, ...
    'kappa', kappa, 'geometry', geom);
end
