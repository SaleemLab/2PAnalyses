function P = sk_cell_components(rebuild)
% Per-cell quantities along the track that do not depend on which position-kernel
% fit is being studied, for every cell in V1_GLMData (6473 cells, 12 sessions).
% Built once from the *_binned.mat files and V1_GLMData.mat, then cached in
% <outRoot>/cell_components.mat. Rows follow CellIdx of V1_GLMData_table.csv.
%
%   mapObs     measured activity per position bin, base-condition trials (activity / sample)
%   mapObsA/B  the same from odd / even base trials (for map reliability)
%   relMap     correlation of mapObsA and mapObsB
%   mapFull    GLM1 full prediction, base trials
%   mapVS      GLM1 prediction without position (vision + speed + on/off)
%   mapPos     GLM1 prediction from the position component
%   visBG, visLM, visEOC   GLM1 background, landmark and end-of-corridor contributions, base condition
%   speedMap   session mean running speed per bin, base trials
%   speedCorr  correlation of activity with speed over trial x bin samples (all trials)
%   speedTune  activity in 8 speed bins (session speed octiles); speedTuneX = bin centres
% The GLM1 vision/speed terms are the only ones saved per cell; the refits store
% position kernels only, so these serve as the reference for every kernel set.
if nargin < 1, rebuild = false; end
C = sk_config(); f = fullfile(C.outRoot, 'cell_components.mat');
if exist(f, 'file') && ~rebuild, P = load(f); return; end

T = readtable(C.tablePath); n = height(T); nb = C.nBins;
z = nan(n, nb, 'single');
P = struct('CellIdx', T.CellIdx, 'x', C.binCentres, 'mapObs', z, 'mapObsA', z, 'mapObsB', z, 'relMap', nan(n,1), ...
    'mapFull', z, 'mapVS', z, 'mapPos', z, 'speedMap', z, 'speedCorr', nan(n,1), ...
    'speedTune', nan(n, 8, 'single'), 'speedTuneX', nan(n, 8, 'single'));
g = 'D:\DATA\Sonali\GLM\ForAS\V1_GLMData.mat';
P.visBG = single(h5read(g, '/Kernels/BGBase')); P.visLM = single(h5read(g, '/Kernels/LandmarkBase'));
P.visEOC = single(h5read(g, '/Kernels/EOCBase'));
assert(size(P.visBG, 1) == n, 'visual bases do not match the table');

for s = 1:size(C.sessions, 1)
    tag = sprintf('%s_%d', C.sessions{s,1}, C.sessions{s,2});
    B = load(fullfile(C.outRoot, [tag '_binned.mat']), 'spk', 'predFull', 'predVS', 'predPos', 'occ', 'speedTB', 'trialCond', 'cellNum');
    [tf, row] = ismember([repmat(C.sessions{s,3}, numel(B.cellNum), 1) B.cellNum], [T.sessionid T.cellNum], 'rows');
    assert(all(tf), '%s: cells missing from the table', tag);
    base = find(B.trialCond == 1); bA = base(1:2:end); bB = base(2:2:end);
    mp = @(Y, tr) squeeze(sum(Y(:, tr, :), 2, 'omitnan')) ./ sum(B.occ(tr, :), 1);
    P.mapObs(row,:) = mp(B.spk, base); P.mapObsA(row,:) = mp(B.spk, bA); P.mapObsB(row,:) = mp(B.spk, bB);
    P.mapFull(row,:) = mp(B.predFull, base); P.mapVS(row,:) = mp(B.predVS, base); P.mapPos(row,:) = mp(B.predPos, base);
    occB = B.occ(base, :); spB = B.speedTB(base, :);
    P.speedMap(row,:) = repmat(single(sum(spB .* occB, 1, 'omitnan') ./ sum(occB .* ~isnan(spB), 1)), numel(row), 1);
    % speed: trial x bin samples with enough occupancy
    ok = B.occ(:) >= 3 & ~isnan(B.speedTB(:)); sp = B.speedTB(ok); w = B.occ(ok);
    ed = [-inf prctile(sp, 12.5:12.5:87.5) inf]; [~, ~, sb] = histcounts(sp, ed);
    sx = accumarray(sb, sp .* w) ./ accumarray(sb, w);
    R = reshape(B.spk, numel(B.cellNum), []) ./ reshape(B.occ, 1, []); R = R(:, ok);
    for c = 1:numel(row)
        r = double(R(c,:))'; P.speedCorr(row(c)) = corr(r, sp, 'rows', 'complete');
        P.speedTune(row(c),:) = accumarray(sb, r .* w) ./ accumarray(sb, w);
        P.speedTuneX(row(c),:) = sx;
        a = P.mapObsA(row(c),:); b = P.mapObsB(row(c),:); P.relMap(row(c)) = corr(a(:), b(:), 'rows', 'complete');
    end
    fprintf('%s: %d cells\n', tag, numel(row));
end
save(f, '-struct', 'P', '-v7.3');
end
