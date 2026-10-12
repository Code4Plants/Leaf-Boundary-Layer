function [betaDeg, windward] = windAttackAngle(alphaXdeg, alphaYdeg, leafAzimuthDeg, windFromDeg)
%WINDATTACKANGLE Angle of attack of the horizontal wind on the leaf plane and the side hit first.
%   [betaDeg, windward] = windAttackAngle(alphaX, alphaY, leafAzimuth, windFrom)
%   betaDeg: 0 = parallel flow, 90 = normal; windward: 'U' (upper side) or 'L'.
%   Same frame conventions as windStreamlineAngle.
ax = deg2rad(alphaXdeg); ay = deg2rad(alphaYdeg);
cy = cos(ay); sy = sin(ay); cx = cos(ax); sx = sin(ax);
R = [1 0 0; 0 cx -sx; 0 sx cx] * [cy 0 -sy; 0 1 0; sy 0 cy];
n = R * [0; 0; 1];
A = deg2rad(leafAzimuthDeg);
beta = deg2rad(windFromDeg + 180);
dWorld = [sin(beta); cos(beta)];
yLeaf = [-sin(A); -cos(A)];
xLeaf = [-cos(A); sin(A)];
d = [dot(dWorld, xLeaf); dot(dWorld, yLeaf); 0];
dn = dot(d, n);
betaDeg = rad2deg(asin(min(1, abs(dn))));
if dn < 0, windward = 'U'; else, windward = 'L'; end
end
