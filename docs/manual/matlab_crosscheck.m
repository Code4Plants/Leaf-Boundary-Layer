function matlab_crosscheck(outFile)
%MATLAB_CROSSCHECK Maximum relative deviation of the MATLAB implementation from the
%   Python reference values (reference/reference_values.json), written as JSON to outFile.
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(fullfile(root, 'matlab'));
ref = jsondecode(fileread(fullfile(root, 'reference', 'reference_values.json')));
res = struct();

cases = asCell(ref.cases); dev = 0; nval = 0;
for i = 1:numel(cases)
    c = cases{i};
    if isfield(c, 'eL') && ~isempty(c.eL)
        r = graefeFreeConvection(c.TL, c.TA, c.L, c.W, c.alpha_x, c.alpha_y, 101325, c.eL, c.eA);
    else
        r = graefeFreeConvection(c.TL, c.TA, c.L, c.W, c.alpha_x, c.alpha_y);
    end
    for key = {'gb_U', 'gb_L', 'H_U', 'H_L', 'NuF_U', 'NuF_L', 'Nuc', 'Ra'}
        dev = max(dev, relDev(r.(key{1}), c.(key{1}))); nval = nval + 1;
    end
end
res.free = struct('n', numel(cases), 'values', nval, 'maxRelDev', dev);

cases = asCell(ref.mixed_cases); dev = 0; nval = 0;
for i = 1:numel(cases)
    c = cases{i}; wd = c.wind_from;
    if isempty(wd) || (isnumeric(wd) && any(isnan(wd))), wd = []; end
    m = graefeMixedConvection(c.TL, c.TA, c.u, c.L, c.W, c.alpha_x, c.alpha_y, c.leaf_azimuth, wd);
    for key = {'Lf', 'NuF_forced', 'NuF_mix_U', 'NuF_mix_L', 'gb_U', 'gb_L', 'H_U', 'H_L', 'forced_share'}
        dev = max(dev, relDev(m.(key{1}), c.(key{1}))); nval = nval + 1;
    end
end
res.mixed = struct('n', numel(cases), 'values', nval, 'maxRelDev', dev);

cases = asCell(ref.angle_cases); dev = 0; nval = 0;
for i = 1:numel(cases)
    c = cases{i};
    m = graefeMixedConvection(296.15, 293.15, c.u, 0.08, 0.04, c.alpha_x, c.alpha_y, c.leaf_azimuth, c.wind_from, [], [], [], [], true, c.TI, [], true);
    for key = {'attack_deg', 'Re_f', 'NuF_forced_U', 'NuF_forced_L', 'gb_U', 'gb_L', 'H_U', 'H_L'}
        dev = max(dev, relDev(m.(key{1}), c.(key{1}))); nval = nval + 1;
    end
end
res.angle = struct('n', numel(cases), 'values', nval, 'maxRelDev', dev);

cases = asCell(ref.scalar_cases); dev = 0; nval = 0;
for i = 1:numel(cases)
    c = cases{i};
    r = leafConductances(293.15 + c.dT, 293.15, c.u, c.L, c.W, c.alpha_x, 0, 180, c.wind_from, [], [], [], c.TI);
    for key = {'gb_U', 'gb_L', 'gv_U', 'gv_L', 'gc_U', 'gc_L', 'k_U', 'k_L', 'cond_share_U', 'free_share_U', 'free_turb_share_U', 'forced_turb_share_U'}
        dev = max(dev, relDev(r.(key{1}), c.(key{1}))); nval = nval + 1;
    end
end
res.scalar = struct('n', numel(cases), 'values', nval, 'maxRelDev', dev);
res.matlabVersion = version;
fid = fopen(outFile, 'w'); fprintf(fid, '%s', jsonencode(res)); fclose(fid);
end

function c = asCell(x)
if isstruct(x), c = num2cell(x); else, c = x; end
end

function d = relDev(a, b)
d = abs(double(a) - double(b)) / max(1, abs(double(b)));
end
