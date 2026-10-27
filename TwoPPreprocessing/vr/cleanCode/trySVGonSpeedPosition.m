%% SVD_separability_session_withShuffle.m
% Population-level SVD separability analysis for one session, with a
% shuffle test for significance, and a summary histogram matching the
% visual style of the other population figures in this analysis.
%
% METHODOLOGY NOTE: the speed x position matrix used here is rebuilt
% directly from rawROIData/speedData (using the SAME speedEdges as your
% saved continuous.matrix), rather than using the officially saved
% continuous.matrix directly. This is necessary for the shuffle test to
% be valid — the real statistic and the null statistic must be computed
% by the identical procedure, and getPositionSpeedMatrix's internals
% aren't available to call repeatedly inside a shuffle loop. Real
% separability values here may differ slightly from a version computed
% directly on continuous.matrix as a result.
%
% Test: is a cell's separability significantly HIGHER than a null built
% from shuffling lap-speed labels? A cell passing this test shows a
% genuine, coherent, single-shape (gain-like) speed effect beyond what
% chance/noise alone would produce. Failing does NOT mean "no speed
% effect" — it means a single coherent shape can't be confirmed.
%
% Requires (same variables used elsewhere in the full-corridor script,
% for ONE session): rawROIData [nROIs x nLaps x nPosBins], speedData
% [nLaps x nPosBins], FilteredROIS_idx, and
% response.speedPositionActivity.continuous.speedEdges.

nShuffles = 1000; % same convention as the existing modulation-index shuffle test
pThresh   = 0.05; % uncorrected, matching convention used elsewhere in this thesis

speedEdges = response.speedPositionActivity.continuous.speedEdges;
nSpeedBins = length(speedEdges) - 1;
nPosBins   = size(speedData, 2);
nLaps      = size(speedData, 1);
nROIs_test = length(FilteredROIS_idx);

realSep = nan(nROIs_test, 1);
pVal    = nan(nROIs_test, 1);

fprintf('Running SVD separability + shuffle test on %d ROIs...\n', nROIs_test);

for i = 1:nROIs_test
    roi = FilteredROIS_idx(i);
    roiActivity = squeeze(rawROIData(roi, :, :)); % [nLaps x nPosBins]

    %% Real statistic
    realMatrix = buildContinuousMatrix(roiActivity, speedData, speedEdges);
    realMatrix(isnan(realMatrix)) = 0;
    centeredReal = realMatrix - mean(realMatrix, 1);
    [~, S, ~] = svd(centeredReal);
    sv = diag(S);
    realSep(i) = sv(1)^2 / sum(sv.^2);

    %% Null distribution — shuffle lap-speed labels, same trick as the
    % existing modulation-index shuffle test (reassign which lap's speed
    % profile pairs with which lap's activity, whole profile at once)
    nullSep = nan(nShuffles, 1);
    for sh = 1:nShuffles
        shuffIdx = randperm(nLaps);
        shuffledSpeedData = speedData(shuffIdx, :);

        shuffMatrix = buildContinuousMatrix(roiActivity, shuffledSpeedData, speedEdges);
        shuffMatrix(isnan(shuffMatrix)) = 0;
        centeredShuff = shuffMatrix - mean(shuffMatrix, 1);
        [~, Ss, ~] = svd(centeredShuff);
        svs = diag(Ss);
        nullSep(sh) = svs(1)^2 / sum(svs.^2);
    end

    pVal(i) = mean(nullSep >= realSep(i));

    if mod(i, 20) == 0
        fprintf('  %d / %d ROIs done\n', i, nROIs_test);
    end
end

sigIdx = pVal < pThresh;
fprintf('\nSVD separability results:\n');
fprintf('  n tested       = %d\n', nROIs_test);
fprintf('  n significant  = %d (%.1f%%) [uncorrected for multiple comparisons]\n', ...
        sum(sigIdx), sum(sigIdx)/nROIs_test*100);
fprintf('  mean sep (all)        = %.3f\n', mean(realSep, 'omitnan'));
fprintf('  mean sep (significant) = %.3f\n', mean(realSep(sigIdx), 'omitnan'));
fprintf('  mean sep (non-sig)     = %.3f\n', mean(realSep(~sigIdx), 'omitnan'));

%% Summary figure — same visual convention as the other population plots
figure('Position', [100 100 500 350], 'Color', 'w');
hold on
histogram(realSep(~sigIdx), 30, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
histogram(realSep(sigIdx),  30, 'FaceColor', [0.0 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9)
xline(mean(realSep, 'omitnan'), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
xlabel('Separability index (fraction of speed-dependent variation, rank-1)')
ylabel('Number of ROIs')
legend('Non-sig (vs. shuffle)', 'Sig — coherent gain-like effect', 'Box', 'off')
title(sprintf('SVD separability, %s | n=%d | sig: %d (%.1f%%)', ...
    strrep(sprintf('%s_%s', mousenumber, sessionName), '_', '\_'), ...
    nROIs_test, sum(sigIdx), sum(sigIdx)/nROIs_test*100))
set(gca, 'Box', 'off', 'TickDir', 'out')

%% Example ROIs — chosen FROM the population results, not decided upfront
% One clearly significant (most separable), one clearly non-significant
% (least separable, closest to noise-like), and one borderline case near
% the significance threshold. This mirrors the population histogram
% rather than illustrating cases picked before seeing the data.
sigROIpositions    = find(sigIdx);
nonSigROIpositions = find(~sigIdx);

exampleList = struct('pos', {}, 'label', {});

if ~isempty(sigROIpositions)
    [~, ord] = sort(realSep(sigROIpositions), 'descend');
    exampleList(end+1) = struct('pos', sigROIpositions(ord(1)), 'label', 'Most separable (significant)');
end
if ~isempty(nonSigROIpositions)
    [~, ord] = sort(realSep(nonSigROIpositions), 'ascend');
    exampleList(end+1) = struct('pos', nonSigROIpositions(ord(1)), 'label', 'Least separable (non-sig, noise-like)');
end
[~, ord] = sort(abs(pVal - pThresh));
exampleList(end+1) = struct('pos', ord(1), 'label', 'Borderline (closest to threshold)');

nExamples = length(exampleList);
figExamples = figure('Color', 'w', 'Position', [100 100 1000 nExamples*180]);

for i = 1:nExamples
    pos = exampleList(i).pos;
    roi = FilteredROIS_idx(pos);
    roiActivity = squeeze(rawROIData(roi, :, :));

    matrix = buildContinuousMatrix(roiActivity, speedData, speedEdges);
    matrix(isnan(matrix)) = 0;
    centered = matrix - mean(matrix, 1);
    [U, S, V] = svd(centered);
    sv = diag(S);
    rank1Recon = U(:,1) * S(1,1) * V(:,1)';
    residual   = centered - rank1Recon;

    sharedMax = max(abs([centered(:); rank1Recon(:)]));
    if sharedMax == 0, sharedMax = 1; end

    subplot(nExamples, 3, (i-1)*3 + 1)
    imagesc(centered, [-sharedMax, sharedMax])
    title(sprintf('ROI %d — Centered\n%s', roi, exampleList(i).label), 'FontSize', 9)
    if i == nExamples, xlabel('Position bin'); end
    ylabel('Speed bin')
    colormap(gca, redbluecmap_fallback()); colorbar
    set(gca, 'Box', 'off', 'TickDir', 'out')

    subplot(nExamples, 3, (i-1)*3 + 2)
    imagesc(rank1Recon, [-sharedMax, sharedMax])
    title(sprintf('Rank-1 (sep=%.2f, p=%.3f)', realSep(pos), pVal(pos)))
    if i == nExamples, xlabel('Position bin'); end
    colormap(gca, redbluecmap_fallback()); colorbar
    set(gca, 'Box', 'off', 'TickDir', 'out')

    subplot(nExamples, 3, (i-1)*3 + 3)
    residMax = max(abs(residual(:)));
    if residMax == 0, residMax = 1; end
    imagesc(residual, [-residMax, residMax])
    title('Residual')
    if i == nExamples, xlabel('Position bin'); end
    colormap(gca, redbluecmap_fallback()); colorbar
    set(gca, 'Box', 'off', 'TickDir', 'out')
end

sgtitle('Example ROIs, chosen from population results (significant / non-sig / borderline)')

%% Helper function — rebuilds the speed x position matrix from raw lap
% data, using the SAME speedEdges as the officially saved
% continuous.matrix, so real and shuffled versions are computed
% identically.
function matrix = buildContinuousMatrix(roiActivity, speedData, speedEdges)
    nSpeedBins = length(speedEdges) - 1;
    nPosBins   = size(speedData, 2);
    matrix = nan(nSpeedBins, nPosBins);

    for b = 1:nPosBins
        binSpeeds = speedData(:, b);
        validLaps = binSpeeds > 1 & ~isnan(binSpeeds); % exclude stationary, same convention as elsewhere

        for k = 1:nSpeedBins
            if k < nSpeedBins
                inThisBin = validLaps & binSpeeds >= speedEdges(k) & binSpeeds < speedEdges(k+1);
            else
                inThisBin = validLaps & binSpeeds >= speedEdges(k) & binSpeeds <= speedEdges(k+1);
            end
            if sum(inThisBin) >= 2
                matrix(k, b) = mean(roiActivity(inThisBin, b), 'omitnan');
            end
        end
    end
end

%% Local fallback colormap if you don't have redbluecmap (Bioinformatics Toolbox)
function cmap = redbluecmap_fallback()
    n = 128;
    cmap = [linspace(0.1,1,n)', linspace(0.1,1,n)', linspace(0.8,1,n)'; ...
            linspace(1,0.8,n)', linspace(1,0.1,n)', linspace(1,0.1,n)'];
end