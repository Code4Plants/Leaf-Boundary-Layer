function g = leafGeometry(L, W, alphaXdeg, alphaYdeg)
%LEAFGEOMETRY Geometric descriptors of an inclined elliptical leaf.
%   g = leafGeometry(L, W, alphaXdeg, alphaYdeg) for a flat ellipse of length
%   L >= width W (m), rotated first about its long axis by alphaY and then about
%   its short axis by alphaX (inclination from horizontal), degrees.
%   Graefe, Grosch & Bitterlich (2022) AFM 317:108884, Section 2.2 / Appendix B.
%
%   Fields: S (area), AR, DI (Eq. 8), P (perimeter), Lh (Eq. 7), theta (leaf
%   normal angle to the vertical, Eq. 11), phi (streamline angle to the short
%   axis, Eq. 10), LCo (mean chord length, Eq. 9), Li (Eq. 9).
if nargin < 4 || isempty(alphaYdeg), alphaYdeg = 0; end
if L < W, error('leafGeometry:axes', 'L must be the longer axis (L >= W)'); end
nLi = 3 * 0.25;                               % n = 3 nRa
ax = deg2rad(alphaXdeg); ay = deg2rad(alphaYdeg);
S = 0.25 * pi * L * W;
DI = 1.112 * ((L + W) / sqrt(L * W))^0.783 - 0.913;   % Eq. 8
P = 2 * DI * sqrt(pi * S);
Lh = S / P;                                    % Eq. 7
theta = acos(max(-1, min(1, cos(ax) * cos(ay))));     % Eq. 11
phi = streamlineAngle(ax, ay);
LCo = meanChordLength(L, W, phi);
Li = (0.04511 * nLi + 1.006) * LCo;           % Eq. 9
g = struct('L', L, 'W', W, 'alphaXdeg', alphaXdeg, 'alphaYdeg', alphaYdeg, ...
    'S', S, 'AR', L / W, 'DI', DI, 'P', P, 'Lh', Lh, 'theta', theta, 'phi', phi, ...
    'LCo', LCo, 'Li', Li);
end

function LCo = meanChordLength(L, W, phi)
% Eq. 9 / A5 written with sin/cos (finite at phi = pi/2: LCo = pi L / 4)
s = sin(phi); c = cos(phi);
LCo = pi * L * W / (4 * sqrt((s * W)^2 + (c * L)^2));
end

function phi = streamlineAngle(ax, ay)
% Angle between the steepest-ascent vector s (projection of the vertical onto the
% leaf plane) and the rotated short axis w, Appendix B.  Returns phi in [0, pi/2];
% pi/2 for a horizontal leaf (direction undefined, Li not used).
cy = cos(ay); sy = sin(ay); cx = cos(ax); sx = sin(ax);
Ry = [cy 0 -sy; 0 1 0; sy 0 cy];              % Eq. B2
Rx = [1 0 0; 0 cx -sx; 0 sx cx];
R = Rx * Ry;
n = R * [0; 0; 1];
w = R * [1; 0; 0];
z = [0; 0; 1];
s = z - dot(z, n) * n;
nrm = norm(s);
if nrm < 1e-12
    phi = pi / 2;
else
    phi = acos(min(1, abs(dot(s, w)) / nrm));
end
end
