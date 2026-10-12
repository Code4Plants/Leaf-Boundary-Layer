function [rho, nu, kappa, k] = airProperties(T, p)
%AIRPROPERTIES Dry-air properties at temperature T (K) and pressure p (Pa).
%   [rho, nu, kappa, k] = airProperties(T, p) returns density (kg m-3),
%   kinematic viscosity (m2 s-1, Sutherland), thermal diffusivity (m2 s-1)
%   and thermal conductivity (W m-1 K-1).  Elementwise over arrays.
%   Identical formulation to python/graefe_convection.py.
if nargin < 2 || isempty(p), p = 101325.0; end
Rd = 287.05; cp = 1004.0;
rho = p ./ (Rd .* T);
mu = 1.458e-6 .* T.^1.5 ./ (T + 110.4);
k = 0.02624 .* (T ./ 300.0).^0.8646;
nu = mu ./ rho;
kappa = k ./ (rho .* cp);
end
