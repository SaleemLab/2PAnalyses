%% ================== COMPARE DATASETS: RESPONSIVE + TRUSTED-FIT BAR CHART ==================
% Run this AFTER save_RF_summary.m has been run once per dataset (V1_somas, V1_boutons,
% RSP_boutons). Loads each summary_*.mat, builds a stacked bar (%responsive, with the
% %trusted-fit portion stacked underneath/shaded differently), and overlays per-session
% dots for any dataset that has multiple sessions (e.g. RSP boutons).

summaryDir  = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Supp_Section1_Fig4_1_VISp\summaries';
files       = {'V1_somas_summary.mat', 'V1_boutons_summary.mat', 'RSP_boutons_summary.mat'};
groupNames  = {'V1 somas', 'V1 boutons', 'RSP boutons'};
nGroups     = numel(files);

allSummaries  = cell(1, nGroups);
pctResponsive = nan(1, nGroups);
pctTrusted    = nan(1, nGroups);
nTotalArr     = nan(1, nGroups);

for g = 1:nGroups
    fpath = fullfile(summaryDir, files{g});
    if ~isfile(fpath)
        warning('Missing summary file for %s (%s) -- skipping this group.', groupNames{g}, fpath);
        continue;
    end
    S = load(fpath, 'summary');
    allSummaries{g} = S.summary;
    nTotalArr(g)    = S.summary.nTotal;

    % Handle either format: full logical vectors, or hand-entered scalar counts
    if isfield(S.summary, 'isResponsive') && ~isempty(S.summary.isResponsive)
        pctResponsive(g) = 100 * sum(S.summary.isResponsive) / S.summary.nTotal;
    elseif isfield(S.summary, 'nResponsive')
        pctResponsive(g) = 100 * S.summary.nResponsive / S.summary.nTotal;
    end

    if isfield(S.summary, 'isTrustedFull') && ~isempty(S.summary.isTrustedFull)
        pctTrusted(g) = 100 * sum(S.summary.isTrustedFull) / S.summary.nTotal;
    elseif isfield(S.summary, 'nTrusted')
        pctTrusted(g) = 100 * S.summary.nTrusted / S.summary.nTotal;
    end
end

fprintf('\n--- Summary ---\n');
for g = 1:nGroups
    if isnan(nTotalArr(g)), continue; end
    fprintf('%-14s n=%-5d  responsive=%.1f%%  trusted fit=%.1f%%\n', ...
        groupNames{g}, nTotalArr(g), pctResponsive(g), pctTrusted(g));
end

%% ---- Stacked bar: trusted fit (bottom) + responsive-but-not-trusted (top) ----
%% ---- Stacked bar: trusted fit (bottom) + responsive-but-not-trusted (top) ----
figBar = figure('Position', [100 100 420 440], 'Color', 'w');
axB = axes; hold(axB, 'on');

remainder = pctResponsive - pctTrusted;               % responsive but fit not trusted
barMatrix = [pctTrusted(:), remainder(:)];

xPos = (1:nGroups) * 0.35;   % <-- closer group spacing (was 1:nGroups)

hBar = bar(axB, xPos, barMatrix, 'stacked', 'BarWidth', 0.5);   % <-- narrower bars (was 0.25... wider % of slot but slot is now smaller)
hBar(1).FaceColor = [0.20 0.40 0.70];   % trusted good fit
hBar(2).FaceColor = [0.75 0.85 0.95];   % responsive, fit not trusted
set(hBar, 'EdgeColor', 'k', 'LineWidth', 0.5);

set(axB, 'XTick', xPos, 'XTickLabel', groupNames, ...
    'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 10);
xlim(axB, [xPos(1)-0.25, xPos(end)+0.25]);
ylabel(axB, '% of ROIs', 'FontName', 'Arial', 'FontSize', 10);
ylim(axB, [0 100]);
legend(axB, {'Trusted Gaussian fit', 'Responsive (fit not trusted)'}, ...
    'Location', 'northoutside', 'Box', 'off', 'FontSize', 9);
title(axB, 'RF responsiveness and fit quality by population', 'FontName', 'Arial', 'FontSize', 11);

% n= labels above each bar
for g = 1:nGroups
    if isnan(nTotalArr(g)), continue; end
    text(axB, xPos(g), pctResponsive(g) + 3, sprintf('n=%d', nTotalArr(g)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontName', 'Arial');
end
%% ---- Overlay per-session scatter for any group with session labels ----
% Points = % responsive within that single session (dots let you see spread vs. the
% pooled bar height, e.g. for RSP boutons where you have many sessions).
jitterWidth = 0.28;
for g = 1:nGroups
    if isnan(nTotalArr(g)) || isempty(allSummaries{g}), continue; end
    S = allSummaries{g};
    if ~isfield(S, 'sessionLabels') || isempty(S.sessionLabels)
        continue;   % no session breakdown available for this group (e.g. hand-entered V1 numbers)
    end

    uniqueSess = unique(S.sessionLabels);
    nSess = numel(uniqueSess);
    sessPctResp = nan(nSess, 1);

    for s = 1:nSess
        m = strcmp(S.sessionLabels, uniqueSess{s});
        sessPctResp(s) = 100 * sum(S.isResponsive(m)) / sum(m);
    end

    jx = xPos(g) + (rand(nSess, 1) - 0.5) * jitterWidth;
    scatter(axB, jx, sessPctResp, 40, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
end
defaultAxesProperties(gca, true);
drawnow;

%% ---- Save ----
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Supp_Section1_Fig4_1_VISp\comparison';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figBar, fullfile(outputDir, 'responsive_trusted_comparison_barplot'));
