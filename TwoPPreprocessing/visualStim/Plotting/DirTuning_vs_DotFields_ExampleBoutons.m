% DirTuning_vs_DotFields_ExampleBoutons.m
%
% Picks a handful of example boutons from a chosen tuning category
% (both/dir-only/dot-only) and plots DirTuning (polar tuning + mini
% panels) side by side with DotFields (tuning curve by speed index),
% one page per bouton, appended to a single PDF.
%
% Requires, all still in the workspace:
%   comboUnits    -- from DirTuning_vs_DotFields_CombinedAnalysis.m
%                    (sessionLabel, roiIdx, isTuned_dir, isTuned_dot)
%   allDirTuning  -- from DirTuning_IncTuningCurveAnalysis_CompareStates_2PData.m
%                    (sessionLabel, roiIdx, fullTrace_stat, fullTrace_run,
%                    timeVec, stimOnDuration, stimValues)
%   allDotUnits   -- from DotFields_IncTuningCurveAnalysis_compareStatesV2_2PData.m
%                    (sessionLabel, roiIdx, alltraces [nSpeeds x 2])

%% choose which category to plot
category = 'both'; % 'both', 'dirOnly', 'dotOnly'
nToPlot  = 6;

switch category
    case 'both',    catIdx = find([comboUnits.isTuned_dir] & [comboUnits.isTuned_dot]);
    case 'dirOnly', catIdx = find([comboUnits.isTuned_dir] & ~[comboUnits.isTuned_dot]);
    case 'dotOnly', catIdx = find(~[comboUnits.isTuned_dir] & [comboUnits.isTuned_dot]);
    otherwise, error('Unknown category: %s', category);
end

fprintf('%d boutons in category "%s".\n', numel(catIdx), category);
if isempty(catIdx)
    error('No boutons in this category -- nothing to plot.');
end

nToPlot = min(nToPlot, numel(catIdx));
plotIdx = catIdx(randperm(numel(catIdx), nToPlot));

%% build lookup indices: (sessionLabel, roiIdx) -> index into allDirTuning / allDotUnits
dirTuningSessionLabels = {allDirTuning.sessionLabel};
dirTuningRoiIdx        = [allDirTuning.roiIdx];
dotUnitsSessionLabels  = {allDotUnits.sessionLabel};
dotUnitsRoiIdx         = [allDotUnits.roiIdx];

%% output
outputFolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section2_Fig4_2\isTunedBoutonsOnly';
pdfPath = fullfile(outputFolder, sprintf('DirTuning_vs_DotFields_Examples_%s.pdf', category));
if exist(pdfPath, 'file')
    delete(pdfPath);
end

hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 550]);
nSuccess = 0;

for i = 1:numel(plotIdx)
    c = plotIdx(i);
    thisSessionLabel = comboUnits(c).sessionLabel;
    thisRoiIdx        = comboUnits(c).roiIdx;

    dirMatch = find(strcmp(dirTuningSessionLabels, thisSessionLabel) & dirTuningRoiIdx == thisRoiIdx, 1);
    dotMatch = find(strcmp(dotUnitsSessionLabels, thisSessionLabel) & dotUnitsRoiIdx == thisRoiIdx, 1);

    if isempty(dirMatch) || isempty(dotMatch)
        warning('Could not find matching entry for %s roiIdx %d -- skipping.', thisSessionLabel, thisRoiIdx);
        continue;
    end

    if ~isvalid(hFig)
        hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 550]);
    end
    clf(hFig);

    try
        plotComparisonPage(allDirTuning(dirMatch), allDotUnits(dotMatch), ...
            comboUnits(c), thisSessionLabel, thisRoiIdx, hFig);
        exportgraphics(hFig, pdfPath, 'Append', true);
        nSuccess = nSuccess + 1;
    catch ME
        warning('Skipping %s roiIdx %d due to error: %s', thisSessionLabel, thisRoiIdx, ME.message);
    end
end

if isvalid(hFig)
    close(hFig);
end
fprintf('\nDone. %d / %d bouton pages written to:\n%s\n', nSuccess, numel(plotIdx), pdfPath);

%% ===================== local function =====================
function plotComparisonPage(sDir, sDot, comboRow, sessionLabel, roiIdx, hFig)

figure(hFig);

%% ---- LEFT: DirTuning (polar + mini panels), combined trace ----
thetaDeg = sDir.stimValues(:)';
[thetaSorted, sortIdx] = sort(mod(round(thetaDeg), 360));
nDir = numel(thetaSorted);

fullTrace_comb = cell(nDir, 1);
for d = 1:nDir
    fullTrace_comb{d} = [sDir.fullTrace_stat{d}; sDir.fullTrace_run{d}];
end
meanDirResponse_comb = cellfun(@(x) mean(x(:), 'omitnan'), fullTrace_comb); % rough combined mean, display only
R_plot = max(meanDirResponse_comb(sortIdx)', 0);

arrowMap = containers.Map( ...
    {0, 45, 90, 135, 180, 225, 270, 315}, ...
    {char(8594), char(8599), char(8593), char(8598), char(8592), char(8601), char(8595), char(8600)});
arrowLabels = cell(1, nDir);
for dIdx = 1:nDir
    if isKey(arrowMap, thetaSorted(dIdx)), arrowLabels{dIdx} = arrowMap(thetaSorted(dIdx)); else, arrowLabels{dIdx} = '?'; end
end

axPolar = polaraxes('Position', [0.02 0.30 0.16 0.62]);
thetaRadClosed = deg2rad([thetaSorted, thetaSorted(1)]);
Rclosed = [R_plot, R_plot(1)];
polarplot(axPolar, thetaRadClosed, Rclosed, 'k-o', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'MarkerSize', 3);
axPolar.ThetaZeroLocation = 'right'; axPolar.ThetaDir = 'counterclockwise';
axPolar.ThetaTick = thetaSorted;
axPolar.ThetaTickLabel = arrayfun(@(x) sprintf('%d%s', x, char(176)), thetaSorted, 'UniformOutput', false);
title(axPolar, 'DirTuning', 'FontSize', 10);

timeVec = sDir.timeVec(:)';
stimOnDuration = sDir.stimOnDuration;
displayWindow = [-0.5, 4];
displayMask = timeVec >= displayWindow(1) & timeVec <= displayWindow(2);
timeVecDisp = timeVec(displayMask);

allBounds = [];
for dIdx = 1:nDir
    tm = fullTrace_comb{sortIdx(dIdx)}(:, displayMask);
    if isempty(tm), continue; end
    allBounds = [allBounds, tm(:)']; %#ok<AGROW>
end
yLims = [min(allBounds, [], 'omitnan'), max(allBounds, [], 'omitnan')];
if any(isnan(yLims)) || diff(yLims) == 0, yLims = [-0.1 0.1]; end

panelLeft = 0.20; panelWidth = 0.55 / nDir; panelBottom = 0.12; panelHeight = 0.28;
stimBarHeight = yLims(1) + 0.08 * diff(yLims);

for d = 1:nDir
    axP = axes('Position', [panelLeft + (d-1)*panelWidth, panelBottom, panelWidth*0.9, panelHeight]);
    hold(axP, 'on');
    if ~isnan(stimOnDuration)
        patch(axP, [0 stimOnDuration stimOnDuration 0], [yLims(1) yLims(1) stimBarHeight stimBarHeight], [0.9 0.9 0.9], 'EdgeColor', 'none');
    end
    traceMat = fullTrace_comb{sortIdx(d)}(:, displayMask);
    if ~isempty(traceMat)
        hTrials = plot(axP, timeVecDisp, traceMat', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.4);
        for k = 1:numel(hTrials), hTrials(k).Color(4) = 0.25; end
        plot(axP, timeVecDisp, mean(traceMat, 1, 'omitnan'), 'k-', 'LineWidth', 1.2);
    end
    xline(axP, 0, 'k:', 'LineWidth', 0.6);
    ylim(axP, yLims); xlim(axP, displayWindow); axis(axP, 'off');
    text(axP, 0.5, -0.22, arrowLabels{d}, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
end

%% ---- RIGHT: DotFields tuning curve (combined stat+run), by speed index ----
axDot = axes('Position', [0.80 0.15 0.18 0.7]);
hold(axDot, 'on');

alltraces = sDot.alltraces; % [nSpeeds x 2]
nSpeeds = size(alltraces, 1);
meanResp = nan(1, nSpeeds);
semResp  = nan(1, nSpeeds);
for sp = 1:nSpeeds
    valsStat = alltraces{sp, 1}; valsStat = valsStat(~isnan(valsStat));
    valsRun  = alltraces{sp, 2}; valsRun  = valsRun(~isnan(valsRun));
    valsComb = [valsStat, valsRun];
    meanResp(sp) = mean(valsComb, 'omitnan');
    semResp(sp)  = std(valsComb, 'omitnan') / sqrt(max(numel(valsComb), 1));
end

errorbar(axDot, 1:nSpeeds, meanResp, semResp, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1.2, 'MarkerSize', 5);
xlim(axDot, [0.5, nSpeeds+0.5]); xticks(axDot, 1:nSpeeds);
xlabel(axDot, 'Visual Speed Index'); ylabel(axDot, '\DeltaF/F');
title(axDot, 'DotFields', 'FontSize', 10);

%% title
sgtitle(sprintf('%s | ROI %d | DirTuning cvR^2=%.3f | DotFields cvR^2=%.3f', ...
    sessionLabel, roiIdx, comboRow.cvR2_dir, comboRow.cvR2_dot), 'Interpreter', 'none');

end
