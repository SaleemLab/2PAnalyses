function sk_run_refit2(sessions)
% Run sk_refit2_session (Poisson; ridge and elastic net; omit fitted after
% position; 1.5 s onset/offset temporal kernels) on the spatial cells of each
% session. One CSV + MAT per session in <outRoot>\refit2\, skipped if the CSV
% exists, so runs can be resumed; then combined in CellIdx order into
% refit2_spatialcells.csv and refit2_kernels.mat.
C = sk_config();
if nargin < 1, sessions = 1:size(C.sessions,1); end
T = readtable(C.tablePath);
od = fullfile(C.outRoot, 'refit2'); if ~exist(od,'dir'), mkdir(od); end
if isempty(gcp('nocreate')), parpool(16); end
for iS = sessions
    tag = sprintf('%s_%d', C.sessions{iS,1}, C.sessions{iS,2});
    f = fullfile(od, ['refit2_' tag '.csv']);
    if exist(f,'file'), fprintf('%s exists, skipping\n', tag); continue; end
    t = T(T.sessionid == C.sessions{iS,3} & T.spatialcells == 1, :);
    t0 = tic;
    [Rf, Ex] = sk_refit2_session(iS, sort(t.cellNum)');
    save(strrep(f, '.csv', '.mat'), 'Ex');
    writetable(Rf, f);                                   % written last: marks the session done
    fprintf('%s: %d cells in %.1f h (%s)\n', tag, height(Rf), toc(t0)/3600, datestr(now));
end

% combine everything finished so far, in CellIdx order
d = dir(fullfile(od, 'refit2_M*.csv'));
R = arrayfun(@(x) readtable(fullfile(od, x.name)), d, 'UniformOutput', false);
R = vertcat(R{:});
[~, ia] = ismember([R.sessionid R.cellNum], [T.sessionid T.cellNum], 'rows');
R.CellIdx = T.CellIdx(ia); R.LLHrel_table = T.LLHrel(ia);
R.LLHi_vis_table = T.LLHi_vis(ia); R.LLHi_pos_table = T.LLHi_pos(ia);
R = movevars(sortrows(R, 'CellIdx'), 'CellIdx', 'Before', 1);
writetable(R, fullfile(C.outRoot, 'refit2_spatialcells.csv'));

E = arrayfun(@(x) load(fullfile(od, strrep(x.name,'.csv','.mat'))), d);
sidc = arrayfun(@(e) repmat(e.Ex.sessionid, numel(e.Ex.cellNum), 1), E, 'UniformOutput', false);
key = [vertcat(sidc{:}), cell2mat(arrayfun(@(e) e.Ex.cellNum, E, 'UniformOutput', false))];
[~, ord] = ismember([R.sessionid R.cellNum], key, 'rows');
K = struct('CellIdx', R.CellIdx, 'sessionid', R.sessionid, 'cellNum', R.cellNum, ...
    'binCentres', E(1).Ex.binCentres, 'models', {E(1).Ex.models}, 'pens', {E(1).Ex.pens}, ...
    'onOffBasis', E(1).Ex.onOffBasis, 'onOffLagS', E(1).Ex.onOffLagS);
for fn = fieldnames(E(1).Ex.kernels)'
    A = cell2mat(arrayfun(@(e) e.Ex.kernels.(fn{1}), E, 'UniformOutput', false)); K.(['pos_' fn{1}]) = A(ord,:);
end
for fn = fieldnames(E(1).Ex.onoff)'
    A = cell2mat(arrayfun(@(e) e.Ex.onoff.(fn{1}), E, 'UniformOutput', false)); K.(['onoff_' fn{1}]) = A(ord,:);
end
F = cell2mat(arrayfun(@(e) reshape(e.Ex.foldLL, size(e.Ex.foldLL,1), []), E, 'UniformOutput', false));
K.foldLL = reshape(F(ord,:), [height(R), size(E(1).Ex.foldLL, 2:4)]);
save(fullfile(C.outRoot, 'refit2_kernels.mat'), '-struct', 'K');
fprintf('combined %d cells from %d sessions\n', height(R), numel(d));
end
