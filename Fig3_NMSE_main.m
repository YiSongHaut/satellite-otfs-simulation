function summary = Fig3_NMSE_main(nMC,maxIter,neighborhood,outputDir)
% Fig. 3 six-panel NMSE reproduction entry point.
% Diagnostic run:      Fig3_NMSE_main(5,2500,'paperCircular')
% Full run:            Fig3_NMSE_main(100,2500,'paperCircular')
% DS=1 us and 16 diffuse rays/tap are explicit assumptions because the paper does not report them.
here = fileparts(mfilename('fullpath'));
addpath(here);
cd(here);
if nargin < 1 || isempty(nMC), nMC = 100; end
if nargin < 2 || isempty(maxIter), maxIter = 2500; end
if nargin < 3 || isempty(neighborhood), neighborhood = 'paperCircular'; end
if nargin < 4 || isempty(outputDir), outputDir = fullfile(here,'otfs_fig3_converged'); end
overrides = struct;
overrides.snrDb = 4:4:24;
overrides.nMC = nMC;
overrides.maxIter = maxIter;
overrides.minIter = min(10,maxIter);
overrides.tol = 1e-5;
overrides.delaySpread = 1e-6;
overrides.raysPerTap = 16;
overrides.neighborhood = neighborhood;
overrides.outputDir = outputDir;
overrides.makePlots = true;
overrides.saveOutputs = true;
summary = run_otfs_nmse('compare',overrides);
end
