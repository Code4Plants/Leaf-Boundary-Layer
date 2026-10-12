function r = leafConductances(TL, TA, u, L, W, alphaXdeg, alphaYdeg, leafAzimuthDeg, windFromDeg, p, eL, eA, TIpct, bc, transition, freeTransition, angleModel)
%LEAFCONDUCTANCES General convection model: heat, H2O and CO2 conductances per leaf side.
%   r = leafConductances(TL, TA, u, L, W, alphaX, alphaY, leafAzimuth, windFrom, p, eL, eA, TIpct, bc, transition, freeTransition, angleModel)
%   Free (laminar/turbulent) + forced (laminar/transitional/turbulent) convection with
%   angle of attack and conduction limit (graefeMixedConvection).  The scalar
%   conductances follow g_s = g_heat (D_s/kappa)^k with the effective diffusivity
%   exponent k of each side as the share-weighted mean of the transport modes
%   (generalisation of Eq. 17 of Graefe et al. 2022): conduction 1, laminar free 0.75,
%   turbulent free 2/3, forced laminar/transitional 2/3, forced turbulent 0.4.
%   Output adds per side X in {U, L}: gv_X (H2O), gc_X (CO2), k_X, cond_share_X,
%   free_share_X, free_turb_share_X, forced_turb_share_X; and D_H2O, D_CO2.
if nargin < 7 || isempty(alphaYdeg), alphaYdeg = 0; end
if nargin < 8 || isempty(leafAzimuthDeg), leafAzimuthDeg = 180; end
if nargin < 9, windFromDeg = []; end
if nargin < 10 || isempty(p), p = 101325.0; end
if nargin < 12, eL = []; eA = []; end
if nargin < 13 || isempty(TIpct), TIpct = 1.0; end
if nargin < 14 || isempty(bc), bc = 'isothermal'; end
if nargin < 15 || isempty(transition), transition = true; end
if nargin < 16 || isempty(freeTransition), freeTransition = true; end
if nargin < 17 || isempty(angleModel), angleModel = true; end
kCond = 1.0; kFreeLam = 0.75; kFreeTurb = 2/3; kForcedLam = 2/3; kForcedTurb = 0.4;
r = graefeMixedConvection(TL, TA, u, L, W, alphaXdeg, alphaYdeg, leafAzimuthDeg, windFromDeg, p, eL, eA, [], transition, TIpct, bc, angleModel, freeTransition);
if transition
    [~, ~, ~, ~, lam, trans, turb] = forcedNusseltLienhard(r.Re_f, 0.71, TIpct, bc);
    tot = lam + trans + turb;
    fTurb = zeros(size(tot)); nz = tot > 0; fTurb(nz) = turb(nz) ./ tot(nz);
else
    fTurb = zeros(size(r.Re_f));
end
fD = (r.TF ./ 273.15).^1.75 .* 101325 ./ p;
r.D_H2O = 2.12e-5 .* fD; r.D_CO2 = 1.39e-5 .* fD;
for s = {'U', 'L'}
    X = s{1};
    nfFreeLam = r.(['NuF_lam_' X]);
    nfFree = r.(['gb_free_' X]) ./ r.kappa;
    nfForced = r.(['NuF_forced_' X]);
    forcedShare = nfForced.^3 ./ (nfFree.^3 + nfForced.^3);
    freeTurbShare = r.(['free_turb_share_' X]);
    sC = min(max(r.Nuc ./ max(nfFreeLam, r.Nuc), 0), 1);
    kFreeL = kFreeLam + (kCond - kFreeLam) .* sC;
    kFree = (1 - freeTurbShare) .* kFreeL + freeTurbShare .* kFreeTurb;
    kForced = (1 - fTurb) .* kForcedLam + fTurb .* kForcedTurb;
    k = (1 - forcedShare) .* kFree + forcedShare .* kForced;
    r.(['k_' X]) = k;
    r.(['gv_' X]) = r.(['gb_' X]) .* (r.D_H2O ./ r.kappa).^k;
    r.(['gc_' X]) = r.(['gb_' X]) .* (r.D_CO2 ./ r.kappa).^k;
    r.(['cond_share_' X]) = sC;
    r.(['free_share_' X]) = 1 - forcedShare;
    r.(['forced_turb_share_' X]) = fTurb;
end
end
