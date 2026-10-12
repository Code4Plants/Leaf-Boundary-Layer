%EXAMPLEFIG6 Sensible heat flux of an 8 x 4 cm leaf over the inclination (cf. Fig. 6b).
TA = 293.15;
ax = 0:10:90;
figure; hold on
for dT = [3 12]
    HU = zeros(size(ax)); HL = HU;
    for i = 1:numel(ax)
        r = graefeFreeConvection(TA + dT, TA, 0.08, 0.04, ax(i));
        HU(i) = r.H_U; HL(i) = r.H_L;
    end
    plot(ax, HU, '-o', 'DisplayName', sprintf('upper, \\DeltaT = %g K', dT));
    plot(ax, HL, '--s', 'DisplayName', sprintf('lower, \\DeltaT = %g K', dT));
end
xlabel('rotation about the x axis \alpha_x (deg)'); ylabel('H (W m^{-2})');
title('Free convection, elliptical leaf 8 x 4 cm (Graefe et al. 2022)'); legend('Location', 'best'); grid on
