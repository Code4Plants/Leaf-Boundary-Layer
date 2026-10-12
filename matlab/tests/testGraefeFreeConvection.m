classdef testGraefeFreeConvection < matlab.unittest.TestCase
    %TESTGRAEFEFREECONVECTION Unit tests (run: runTests from the matlab folder).

    methods (TestClassSetup)
        function addPaths(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here, '..'));
            tc.addTeardown(@() rmpath(fullfile(here, '..')));
        end
    end

    methods (Test)
        function geometryCircle(tc)
            g = leafGeometry(0.084, 0.084, 0);
            tc.verifyEqual(g.DI, 1.0, 'AbsTol', 0.005);
            tc.verifyEqual(g.Lh, 0.084 / 4, 'AbsTol', 1e-4);
            tc.verifyEqual(g.S, pi * 0.042^2, 'AbsTol', 1e-12);
        end

        function geometryEllipse(tc)
            g = leafGeometry(0.08, 0.04, 60);
            tc.verifyEqual(g.AR, 2.0);
            tc.verifyEqual(rad2deg(g.phi), 90, 'AbsTol', 1e-9);
            tc.verifyEqual(g.LCo, pi * 0.08 / 4, 'AbsTol', 1e-12);
            tc.verifyEqual(g.Li, 0.816 * 0.08, 'AbsTol', 5e-4);
            tc.verifyEqual(rad2deg(g.theta), 60, 'AbsTol', 1e-9);
            g2 = leafGeometry(0.08, 0.04, 60, 45);
            tc.verifyGreaterThan(rad2deg(g2.phi), 0);
            tc.verifyLessThan(rad2deg(g2.phi), 90);
            tc.verifyLessThan(g2.Li, g.Li);
            tc.verifyEqual(cos(g2.theta), cosd(60) * cosd(45), 'AbsTol', 1e-12);
            tc.verifyError(@() leafGeometry(0.04, 0.08, 0), 'leafGeometry:axes');
        end

        function airAt20C(tc)
            [rho, nu, kappa, k] = airProperties(293.15);
            tc.verifyEqual(rho, 1.204, 'AbsTol', 0.005);
            tc.verifyEqual(nu, 1.51e-5, 'AbsTol', 0.05e-5);
            tc.verifyEqual(kappa, 2.12e-5, 'AbsTol', 0.1e-5);
            tc.verifyEqual(k, 0.0257, 'AbsTol', 5e-4);
        end

        function verticalPlate(tc)
            r = graefeFreeConvection(293.15 + 12, 293.15, 0.08, 0.04, 90);
            tc.verifyEqual(r.NuF_U, r.NuF_L, 'AbsTol', 1e-12);
            nuLam = 0.513 * r.geometry.Li^(-0.25) * r.Ra^0.25;
            tc.verifyGreaterThan(r.NuF_U, nuLam);
            tc.verifyLessThan(r.NuF_U, 1.3 * nuLam);
        end

        function cooledLeafSwapsSides(tc)
            % heated and cooled leaf with the same |dT|: sides are exchanged (the film
            % temperatures differ by 6 K, hence ~1 % difference in nu and Ra)
            h = graefeFreeConvection(293.15 + 6, 293.15, 0.08, 0.04, 0);
            c = graefeFreeConvection(293.15 - 6, 293.15, 0.08, 0.04, 0);
            tc.verifyEqual(h.NuF_U, c.NuF_L, 'RelTol', 0.02);
            tc.verifyEqual(h.NuF_L, c.NuF_U, 'RelTol', 0.02);
            tc.verifyLessThan(h.NuF_U, h.NuF_L);
            % heated horizontal leaf: the upper side is the weaker one (p5 < p6), the cooled leaf the reverse
            tc.verifyGreaterThan(c.NuF_U, c.NuF_L);
        end

        function fig6Behaviour(tc)
            TA = 293.15;
            for dT = [3 12]
                prevGap = Inf;
                for ax = [0 30 60 90]
                    r = graefeFreeConvection(TA + dT, TA, 0.08, 0.04, ax);
                    tc.verifyGreaterThan([r.H_U r.H_L], 2);
                    tc.verifyLessThan([r.H_U r.H_L], 120);
                    gap = abs(r.H_U - r.H_L) / max(r.H_U, r.H_L);
                    tc.verifyLessThanOrEqual(gap, prevGap + 1e-9);
                    prevGap = gap;
                end
                tc.verifyLessThan(prevGap, 1e-9);
            end
        end

        function vectorisedAndEq18(tc)
            TA = 293.15; dT = [-5 -1 1 5];
            r = graefeFreeConvection(TA + dT, TA, 0.08, 0.04, 60);
            tc.verifySize(r.gb_U, [1 4]);
            tc.verifyEqual(sign(r.H_U), sign(dT));
            rw = graefeFreeConvection(TA - 0.5, TA, 0.08, 0.04, 60, 0, 101325, 2300, 800);
            tc.verifyGreaterThan(rw.dT, 0);
        end

        function exactCharacteristicLength(tc)
            for LW = {[0.08 0.04], [0.12 0.04], [0.084 0.084]}
                L = LW{1}(1); W = LW{1}(2);
                for phi = [0 pi/4 pi/2]
                    LCo = pi * L * W / (4 * sqrt((sin(phi) * W)^2 + (cos(phi) * L)^2));
                    exact = characteristicLength(L, W, phi, 0.75);
                    tc.verifyEqual(exact / ((0.04511 * 0.75 + 1.006) * LCo), 1, 'AbsTol', 0.003);
                end
            end
            tc.verifyEqual(characteristicLength(0.08, 0.04, pi/2, 0.75), 0.816 * 0.08, 'AbsTol', 1e-4);
            LCo = pi * 0.08 / 4;
            tc.verifyEqual(characteristicLength(0.08, 0.04, pi/2, 0.5) / LCo, 1.0286, 'AbsTol', 0.002);
            tc.verifyEqual(characteristicLength(0.08, 0.04, 0.3, 1.0), ...
                pi * 0.08 * 0.04 / (4 * sqrt((sin(0.3) * 0.04)^2 + (cos(0.3) * 0.08)^2)), 'AbsTol', 1e-15);
        end

        function windStreamlines(tc)
            [phi, frac] = windStreamlineAngle(0, 0, 180, 0);        % north wind along the long axis
            tc.verifyEqual(rad2deg(phi), 90, 'AbsTol', 1e-9); tc.verifyEqual(frac, 1, 'AbsTol', 1e-12);
            [phi, ~] = windStreamlineAngle(0, 0, 180, 90);          % east wind along the short axis
            tc.verifyEqual(rad2deg(phi), 0, 'AbsTol', 1e-9);
            [phi, ~] = windStreamlineAngle(0, 0, 180, 45);
            tc.verifyEqual(rad2deg(phi), 45, 'AbsTol', 1e-9);
            [~, frac] = windStreamlineAngle(90, 0, 180, 180);       % wind normal to a vertical leaf
            tc.verifyLessThan(frac, 1e-9);
            [phi, frac] = windStreamlineAngle(60, 0, 180, 180);
            tc.verifyEqual(rad2deg(phi), 90, 'AbsTol', 1e-9); tc.verifyEqual(frac, cosd(60), 'AbsTol', 1e-9);
        end

        function mixedConvection(tc)
            TA = 293.15;
            free = graefeFreeConvection(TA + 3, TA, 0.08, 0.04, 60);
            still = graefeMixedConvection(TA + 3, TA, 0, 0.08, 0.04, 60, 0, 180, 270);
            tc.verifyEqual(still.gb_U, free.gb_U, 'AbsTol', 1e-15);
            tc.verifyEqual(still.forced_share, 0);
            windy = graefeMixedConvection(TA + 3, TA, 3, 0.08, 0.04, 60, 0, 180, 270);
            tc.verifyGreaterThan(windy.gb_U, free.gb_U);
            tc.verifyGreaterThan(windy.forced_share, 0.9);
            along = graefeMixedConvection(TA + 3, TA, 2, 0.08, 0.04, 0, 0, 180, 0);
            across = graefeMixedConvection(TA + 3, TA, 2, 0.08, 0.04, 0, 0, 180, 90);
            unknown = graefeMixedConvection(TA + 3, TA, 2, 0.08, 0.04, 0, 0, 180, []);
            tc.verifyGreaterThan(along.Lf, across.Lf);
            tc.verifyLessThan(along.gb_forced, across.gb_forced);
            tc.verifyGreaterThan(unknown.Lf, across.Lf); tc.verifyLessThan(unknown.Lf, along.Lf);
            vec = graefeMixedConvection(TA + 3, TA, [0.1 1 5], 0.08, 0.04, 60, 0, 180, 200);
            tc.verifyGreaterThan(diff(vec.gb_U), 0);
        end

        function lienhardTransition(tc)
            Re = [1e3 1e4 1e5 1e6 3e6];
            [nu, regime, ReL] = forcedNusseltLienhard(Re, 0.71, 0.1);
            tc.verifyEqual(ReL, 5e5);
            lam = 0.664 * Re.^0.5 * 0.71^(1/3);
            tc.verifyEqual(nu(1:3), lam(1:3), 'RelTol', 1e-12);
            tc.verifyGreaterThan(nu(4:5), lam(4:5));
            tc.verifyEqual(regime{end}, 'turbulent');
            grid = logspace(3, 7, 400);
            v = forcedNusseltLienhard(grid, 0.71, 2.0);
            tc.verifyGreaterThan(diff(v), 0);
            tc.verifyLessThan(max(abs(diff(log(v)))), 0.05);
            tc.verifyGreaterThan(forcedNusseltLienhard(1e5, 0.71, 10), forcedNusseltLienhard(1e5, 0.71, 1));
            tc.verifyGreaterThan(forcedNusseltLienhard(1e4, 0.71, 1, 'uniform_flux'), forcedNusseltLienhard(1e4, 0.71, 1));
            % cross-check with the Python implementation (graefe_convection.forced_nusselt_lienhard)
            tc.verifyEqual(forcedNusseltLienhard(1e3, 0.71, 0.1), 18.73214577951894, 'RelTol', 1e-9);
            tc.verifyEqual(forcedNusseltLienhard(1e6, 0.71, 0.1), 834.3053288803666, 'RelTol', 1e-9);
            tc.verifyEqual(forcedNusseltLienhard(1e5, 0.71, 10.0), 282.81377236628475, 'RelTol', 1e-9);
            m = graefeMixedConvection(293.15 + 3, 293.15, [0.5 5 20], 0.08, 0.04, 60, 0, 180, 270, [], [], [], [], true, 10);
            tc.verifyGreaterThan(diff(m.gb_U), 0);
            tc.verifyEqual(m.forced_regime{1}, 'laminar');
        end

        function angleOfAttackModel(tc)
            [gW, gB] = angleOfAttackFactors(0);
            tc.verifyEqual([gW gB], [1 1], 'AbsTol', 1e-12);
            b = 0:5:90;
            [~, gBs] = angleOfAttackFactors(b);
            tc.verifyGreaterThan(max(gBs), 1.15);
            [~, i] = max(gBs); tc.verifyGreaterThanOrEqual(b(i), 10); tc.verifyLessThanOrEqual(b(i), 35);
            tc.verifyLessThan(gBs(end), 1);
            [nuW, nuB] = forcedNusseltLeaf(1e4, 20, 0.71, 1.0);
            tc.verifyGreaterThan(nuB, nuW);
            [beta, side] = windAttackAngle(90, 0, 180, 180);
            tc.verifyEqual(beta, 90, 'AbsTol', 1e-9); tc.verifyEqual(side, 'U');
            [beta, side] = windAttackAngle(60, 0, 180, 180);
            tc.verifyEqual(beta, 60, 'AbsTol', 1e-9); tc.verifyEqual(side, 'U');
            tc.verifyEqual(windAttackAngle(0, 0, 180, 90), 0, 'AbsTol', 1e-9);
        end

        function matchesPythonReferenceAngle(tc)
            here = fileparts(mfilename('fullpath'));
            ref = jsondecode(fileread(fullfile(here, '..', '..', 'reference', 'reference_values.json')));
            cases = ref.angle_cases;
            if isstruct(cases), cases = num2cell(cases); end
            for i = 1:numel(cases)
                c = cases{i};
                m = graefeMixedConvection(296.15, 293.15, c.u, 0.08, 0.04, c.alpha_x, c.alpha_y, c.leaf_azimuth, c.wind_from, [], [], [], [], true, c.TI, [], true);
                for key = {'attack_deg', 'Re_f', 'NuF_forced_U', 'NuF_forced_L', 'gb_U', 'gb_L', 'H_U', 'H_L'}
                    tc.verifyEqual(m.(key{1}), c.(key{1}), 'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf('angle case %d field %s', i, key{1}));
                end
                tc.verifyEqual(m.windward_side, c.windward_side);
            end
        end

        function scalarConductances(tc)
            TA = 293.15;
            still = leafConductances(TA + 3, TA, 0, 0.08, 0.04, 60, 0, 180, 180);
            tc.verifyGreaterThan(still.k_U, 0.75); tc.verifyLessThan(still.k_U, 1);
            tc.verifyGreaterThan(still.gv_U, still.gb_U); tc.verifyGreaterThan(still.gb_U, still.gc_U);
            windy = leafConductances(TA + 3, TA, 3, 0.08, 0.04, 60, 0, 180, 180, [], [], [], 1.0);
            tc.verifyEqual(windy.k_U, 2/3, 'AbsTol', 0.01);
            big = leafConductances(TA + 3, TA, 15, 0.6, 0.3, 0, 0, 180, 90, [], [], [], 40);
            tc.verifyGreaterThan(big.forced_turb_share_U, 0.3); tc.verifyLessThan(big.k_U, 2/3);
            hot = leafConductances(TA + 15, TA, 0, 0.6, 0.3, 0, 0, 180, 90);
            tc.verifyGreaterThan(hot.free_turb_share_U, 0.3);
        end

        function matchesPythonReferenceScalar(tc)
            here = fileparts(mfilename('fullpath'));
            ref = jsondecode(fileread(fullfile(here, '..', '..', 'reference', 'reference_values.json')));
            cases = ref.scalar_cases;
            if isstruct(cases), cases = num2cell(cases); end
            for i = 1:numel(cases)
                c = cases{i};
                r = leafConductances(293.15 + c.dT, 293.15, c.u, c.L, c.W, c.alpha_x, 0, 180, c.wind_from, [], [], [], c.TI);
                for key = {'gb_U', 'gb_L', 'gv_U', 'gv_L', 'gc_U', 'gc_L', 'k_U', 'k_L', 'cond_share_U', 'free_share_U', 'free_turb_share_U', 'forced_turb_share_U'}
                    tc.verifyEqual(r.(key{1}), c.(key{1}), 'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf('scalar case %d field %s', i, key{1}));
                end
            end
        end

        function matchesPythonReferenceMixed(tc)
            here = fileparts(mfilename('fullpath'));
            ref = jsondecode(fileread(fullfile(here, '..', '..', 'reference', 'reference_values.json')));
            cases = ref.mixed_cases;
            if isstruct(cases), cases = num2cell(cases); end
            for i = 1:numel(cases)
                c = cases{i};
                wd = c.wind_from;
                if isempty(wd) || (isnumeric(wd) && any(isnan(wd))), wd = []; end
                m = graefeMixedConvection(c.TL, c.TA, c.u, c.L, c.W, c.alpha_x, c.alpha_y, c.leaf_azimuth, wd);
                for key = {'Lf', 'NuF_forced', 'NuF_mix_U', 'NuF_mix_L', 'gb_U', 'gb_L', 'H_U', 'H_L', 'forced_share'}
                    tc.verifyEqual(m.(key{1}), c.(key{1}), 'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                        sprintf('mixed case %d field %s', i, key{1}));
                end
            end
        end

        function matchesPythonReference(tc)
            here = fileparts(mfilename('fullpath'));
            ref = jsondecode(fileread(fullfile(here, '..', '..', 'reference', 'reference_values.json')));
            cases = ref.cases;
            if isstruct(cases), cases = num2cell(cases); end
            for i = 1:numel(cases)
                c = cases{i};
                if isfield(c, 'eL') && ~isempty(c.eL)
                    r = graefeFreeConvection(c.TL, c.TA, c.L, c.W, c.alpha_x, c.alpha_y, 101325, c.eL, c.eA);
                else
                    r = graefeFreeConvection(c.TL, c.TA, c.L, c.W, c.alpha_x, c.alpha_y);
                end
                for key = {'gb_U', 'gb_L', 'H_U', 'H_L', 'NuF_U', 'NuF_L', 'Nuc', 'Ra'}
                    tc.verifyEqual(r.(key{1}), c.(key{1}), 'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                        sprintf('case %d field %s', i, key{1}));
                end
                tc.verifyEqual(r.geometry.Li, c.Li, 'AbsTol', 1e-12);
                tc.verifyEqual(rad2deg(r.geometry.phi), c.phi_deg, 'AbsTol', 1e-9);
            end
        end
    end
end
