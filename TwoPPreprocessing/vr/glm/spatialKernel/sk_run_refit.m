function sk_run_refit(sessions)
% Refit Gaussian/Poisson, +/- position, +/- vision x speed for the spatial
% cells of each session (sk_refit_session). One CSV per session in
% <outRoot>\refit\, skipped if already present, so the run can be resumed;
% then all sessions are combined into refit_spatialcells.csv (CellIdx order).
C = sk_config();
if nargin < 1, sessions = 1:size(C.sessions,1); end
T = readtable(C.tablePath);
od = fullfile(C.outRoot, 'refit'); if ~exist(od,'dir'), mkdir(od); end
if isempty(gcp('nocreate')), parpool(16); end
for iS = sessions
    tag = sprintf('%s_%d', C.sessions{iS,1}, C.sessions{iS,2});
    f = fullfile(od, ['refit_' tag '.csv']);
    if exist(f,'file'), fprintf('%s exists, skipping\n', tag); continue; end
    t = T(T.sessionid == C.sessions{iS,3} & T.spatialcells == 1, :);
    t0 = tic;
    Rf = sk_refit_session(iS, sort(t.cellNum)');
    writetable(Rf, f);
    fprintf('%s: %d cells in %.1f h (%s)\n', tag, height(Rf), toc(t0)/3600, datestr(now));
end
% combine
d = dir(fullfile(od, 'refit_M*.csv'));
R = vertcat(arrayfun(@(x) readtable(fullfile(od, x.name)), d, 'UniformOutput', false));
R = vertcat(R{:});
[~, ia] = ismember([R.sessionid R.cellNum], [T.sessionid T.cellNum], 'rows');
R.CellIdx = T.CellIdx(ia); R.LLHrel_table = T.LLHrel(ia);
R.LLHi_vis_table = T.LLHi_vis(ia); R.LLHi_pos_table = T.LLHi_pos(ia);
R = movevars(sortrows(R, 'CellIdx'), 'CellIdx', 'Before', 1);
writetable(R, fullfile(C.outRoot, 'refit_spatialcells.csv'));
fprintf('combined %d cells from %d sessions\n', height(R), numel(d));
end
