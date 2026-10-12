function Lc = characteristicLength(L, W, phi, n)
%CHARACTERISTICLENGTH Exact flux-weighted characteristic length of Eq. 4 (ellipse).
%   Lc = characteristicLength(L, W, phi, n) for chords parallel to the streamline
%   direction at angle phi (rad) to the short axis.  With D(t) = D0 sqrt(1-(t/tmax)^2)
%   and D0 = (4/pi) LCo(phi), Eq. 4 / A6 has the closed form
%       Lc^(n-1) = D0^(n-1) * (2/pi) * B(1/2, n/2 + 1)
%   (B = Beta function).  n = 0.75 reproduces Eq. 9 (free convection) within
%   0.2 %; n = 0.5 is the length for laminar forced convection (Nu ~ Re^0.5).
LCo = meanChordLength(L, W, phi);
if abs(n - 1) < 1e-12
    Lc = LCo;
    return
end
D0 = 4 / pi * LCo;
B = gamma(0.5) * gamma(0.5 * n + 1) / gamma(0.5 * n + 1.5);
Lc = D0 * (2 / pi * B)^(1 / (n - 1));
end

function LCo = meanChordLength(L, W, phi)
s = sin(phi); c = cos(phi);
LCo = pi * L * W / (4 * sqrt((s * W)^2 + (c * L)^2));
end
