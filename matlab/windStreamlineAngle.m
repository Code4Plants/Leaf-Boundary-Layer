function [phiW, frac] = windStreamlineAngle(alphaXdeg, alphaYdeg, leafAzimuthDeg, windFromDeg)
%WINDSTREAMLINEANGLE Angle of the wind streamlines over the leaf to the short axis.
%   [phiW, frac] = windStreamlineAngle(alphaX, alphaY, leafAzimuth, windFrom)
%   projects the horizontal wind vector (blowing towards windFrom + 180 deg) onto
%   the leaf plane.  leafAzimuth is the compass direction the upper leaf surface
%   faces (180 = south), i.e. the horizontal direction of the leaf normal; the long
%   axis runs towards leafAzimuth + 180 (up-slope).  phiW in [0, pi/2]; frac is the
%   fraction of the wind speed lying in the leaf plane (0 = wind normal to the leaf).
ax = deg2rad(alphaXdeg); ay = deg2rad(alphaYdeg);
cy = cos(ay); sy = sin(ay); cx = cos(ax); sx = sin(ax);
Ry = [cy 0 -sy; 0 1 0; sy 0 cy];
Rx = [1 0 0; 0 cx -sx; 0 sx cx];
R = Rx * Ry;
n = R * [0; 0; 1];
w = R * [1; 0; 0];
A = deg2rad(leafAzimuthDeg);
beta = deg2rad(windFromDeg + 180);
dWorld = [sin(beta); cos(beta)];              % (east, north)
yLeaf = [-sin(A); -cos(A)];                   % long axis (horizontal, before tilting)
xLeaf = [-cos(A); sin(A)];                    % short axis = yLeaf x z
d = [dot(dWorld, xLeaf); dot(dWorld, yLeaf); 0];
dP = d - dot(d, n) * n;
frac = norm(dP);
if frac < 1e-9
    phiW = pi / 2;
    frac = 0;
else
    phiW = acos(min(1, abs(dot(dP, w)) / frac));
end
end
