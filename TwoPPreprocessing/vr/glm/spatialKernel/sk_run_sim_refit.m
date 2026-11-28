function sk_run_sim_refit(iSess, lags, nWorkers)
% Refit simulated cells without a position term (sk_simulate_null) with the refit-1
% pipeline, for the GLM1 spatial cells of one session, at several timing lags (samples
% at 60 Hz). One file per lag in <outRoot>/simrefit/, skipped if present (resumable).
% Load the kernels with sk_kernelset('sim<lag>_vsip') etc.
C = sk_config();
if nargin < 1, iSess = 9; end                       % M25126_20260313
if nargin < 2, lags = [0 6 -6]; end
if nargin < 3, nWorkers = 8; end
T = readtable(C.tablePath); od = fullfile(C.outRoot, 'simrefit'); if ~exist(od, 'dir'), mkdir(od); end
if isempty(gcp('nocreate')), parpool(nWorkers); end
tag = sprintf('%s_%d', C.sessions{iSess,1}, C.sessions{iSess,2});
cells = sort(T.cellNum(T.sessionid == C.sessions{iSess,3} & T.spatialcells == 1))';
for lag = lags
    f = fullfile(od, sprintf('sim%s_%s.csv', lagName(lag), tag));
    if exist(f, 'file'), fprintf('%s exists, skipping\n', f); continue; end
    t0 = tic; [Y, phi] = sk_simulate_null(iSess, lag);
    [Rf, Ex] = sk_refit_session(iSess, cells, struct('Y', Y));
    Rf.phiSim = phi(cells)'; Ex.lagSamples = lag;
    save(strrep(f, '.csv', '.mat'), 'Ex'); writetable(Rf, f);
    fprintf('lag %+d: %d cells in %.1f h (%s)\n', lag, numel(cells), toc(t0)/3600, datestr(now));
end
end

function s = lagName(lag)
if lag == 0, s = '0'; elseif lag > 0, s = sprintf('lag%d', lag); else, s = sprintf('lead%d', -lag); end
end
