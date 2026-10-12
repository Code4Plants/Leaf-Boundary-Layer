function results = runTests()
%RUNTESTS Run the MATLAB unit tests of the free-convection model.
here = fileparts(mfilename('fullpath'));
addpath(here);
results = runtests(fullfile(here, 'tests'));
disp(table(results));
if any([results.Failed])
    error('runTests:failed', '%d test(s) failed', sum([results.Failed]));
end
end
