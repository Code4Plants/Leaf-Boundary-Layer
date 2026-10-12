function r = graefeMixedConvection(TL, TA, u, L, W, alphaXdeg, alphaYdeg, leafAzimuthDeg, windFromDeg, p, eL, eA, cForced, transition, TIpct, bc, angleModel, freeTransition)
%GRAEFEMIXEDCONVECTION Mixed free/forced convection of an inclined elliptical leaf.
%   r = graefeMixedConvection(TL, TA, u, L, W, alphaX, alphaY, leafAzimuth, windFrom, p, eL, eA, cForced)
%
%   Free convection: Graefe, Grosch & Bitterlich (2022) via graefeFreeConvection.
%   Forced convection (extension, not part of the paper): laminar flat plate
%   NuF = cForced Re^0.5 Pr^(1/3) / Lf with Re = u Lf / nu, where the streamline
%   length Lf follows the flowline concept of Section 2.2 with the exact Eq. 4
%   solution for n = 0.5 along the wind direction projected onto the leaf plane
%   (see windStreamlineAngle, characteristicLength).  Both are blended per side:
%   NuF = (NuF_free^3 + NuF_forced^3)^(1/3).  Assisting/opposing flow effects are
%   neglected.
%
%   Inputs: u = wind speed at the leaf (m s-1, elementwise with TL/TA);
%   leafAzimuth = compass direction the upper surface faces (default 180);
%   windFrom = meteorological wind direction (deg); [] = unknown direction ->
%   orientation-averaged chord length (Eq. 5).  p, eL, eA as in graefeFreeConvection;
%   cForced default 0.664 (Pohlhausen; Campbell & Norman use 0.6).
%
%   Output struct: fields of graefeFreeConvection (with gb_U/gb_L/H_U/H_L replaced
%   by the mixed values) plus gb_free_U/L, H_free_U/L, NuF_forced, gb_forced, Lf,
%   phi_w (rad, NaN if windFrom empty), wind_in_plane_fraction, NuF_mix_U/L,
%   forced_share.
if nargin < 7 || isempty(alphaYdeg), alphaYdeg = 0; end
if nargin < 8 || isempty(leafAzimuthDeg), leafAzimuthDeg = 180; end
if nargin < 9, windFromDeg = []; end
if nargin < 10 || isempty(p), p = 101325.0; end
if nargin < 12, eL = []; eA = []; end
if nargin < 13 || isempty(cForced), cForced = 0.664; end
if nargin < 14 || isempty(transition), transition = false; end
if nargin < 15 || isempty(TIpct), TIpct = 1.0; end
if nargin < 16 || isempty(bc), bc = 'isothermal'; end
if nargin < 17 || isempty(angleModel), angleModel = false; end
if nargin < 18 || isempty(freeTransition), freeTransition = false; end
%   freeTransition = true uses freeConvectionTransition (Ra^(1/3) branch) for the free part.
%   transition = true uses the Lienhard (2020) laminar-transitional-turbulent average
%   (forcedNusseltLienhard) with the onset from the turbulence intensity TIpct.
%   angleModel = true multiplies the forced part by the finite-plate enhancement f_e
%   and the windward / leeward angle-of-attack factors (angleOfAttackFactors); the
%   angle of attack and the windward side follow from the wind direction
%   (windAttackAngle; beta = 0 when windFromDeg is empty).
nForced = 0.5; Pr = 0.71; cp = 1004.0;

if freeTransition
    free = freeConvectionTransition(TL, TA, L, W, alphaXdeg, alphaYdeg, p, eL, eA);
else
    free = graefeFreeConvection(TL, TA, L, W, alphaXdeg, alphaYdeg, p, eL, eA);
    free.NuF_lam_U = free.NuF_U; free.NuF_lam_L = free.NuF_L;
    free.free_turb_share_U = zeros(size(free.NuF_U)); free.free_turb_share_L = zeros(size(free.NuF_L));
end
g = free.geometry;
if isempty(windFromDeg)
    D0 = 4 * g.S / g.P;                                        % (4/pi) LCo, LCo = pi S / P (Eq. 5)
    B = gamma(0.5) * gamma(0.5 * nForced + 1) / gamma(0.5 * nForced + 1.5);
    Lf = D0 * (2 / pi * B)^(1 / (nForced - 1));
    phiW = NaN; frac = NaN;
else
    [phiW, frac] = windStreamlineAngle(alphaXdeg, alphaYdeg, leafAzimuthDeg, windFromDeg);
    if frac < 1e-6
        D0 = 4 * g.S / g.P;
        B = gamma(0.5) * gamma(0.5 * nForced + 1) / gamma(0.5 * nForced + 1.5);
        Lf = D0 * (2 / pi * B)^(1 / (nForced - 1));
    else
        Lf = characteristicLength(L, W, phiW, nForced);
    end
end
Re = max(u .* Lf ./ free.nu, 0);
if transition
    [NuAvg, regime] = forcedNusseltLienhard(Re, Pr, TIpct, bc);
    NuF_forced = NuAvg ./ Lf;
else
    NuF_forced = cForced .* sqrt(Re) .* Pr^(1 / 3) ./ Lf;
    regime = repmat({'laminar'}, size(Re));
end
NuF_forced_U = NuF_forced; NuF_forced_L = NuF_forced;
betaDeg = 0; windward = 'U';
if angleModel
    if ~isempty(windFromDeg)
        [betaDeg, windward] = windAttackAngle(alphaXdeg, alphaYdeg, leafAzimuthDeg, windFromDeg);
    end
    [gW, gB, prm] = angleOfAttackFactors(betaDeg);
    if windward == 'U'
        NuF_forced_U = prm.f_e .* NuF_forced .* gW; NuF_forced_L = prm.f_e .* NuF_forced .* gB;
    else
        NuF_forced_U = prm.f_e .* NuF_forced .* gB; NuF_forced_L = prm.f_e .* NuF_forced .* gW;
    end
end
NuF_mix_U = (free.NuF_U.^3 + NuF_forced_U.^3).^(1 / 3);
NuF_mix_L = (free.NuF_L.^3 + NuF_forced_L.^3).^(1 / 3);
r = free;
r.gb_free_U = free.gb_U; r.gb_free_L = free.gb_L;
r.H_free_U = free.H_U; r.H_free_L = free.H_L;
r.NuF_forced = NuF_forced;
r.NuF_forced_U = NuF_forced_U; r.NuF_forced_L = NuF_forced_L;
r.attack_deg = betaDeg; r.windward_side = windward;
r.Re_f = Re; r.forced_regime = regime;
r.gb_forced = free.kappa .* NuF_forced;
r.gb_forced_U = free.kappa .* NuF_forced_U; r.gb_forced_L = free.kappa .* NuF_forced_L;
r.Lf = Lf; r.phi_w = phiW; r.wind_in_plane_fraction = frac;
r.NuF_mix_U = NuF_mix_U; r.NuF_mix_L = NuF_mix_L;
r.gb_U = free.kappa .* NuF_mix_U;
r.gb_L = free.kappa .* NuF_mix_L;
r.H_U = free.rho .* cp .* r.gb_U .* (TL - TA);
r.H_L = free.rho .* cp .* r.gb_L .* (TL - TA);
shareU = NuF_forced_U.^3 ./ (free.NuF_U.^3 + NuF_forced_U.^3);
shareL = NuF_forced_L.^3 ./ (free.NuF_L.^3 + NuF_forced_L.^3);
r.forced_share = 0.5 * (shareU + shareL);
end
