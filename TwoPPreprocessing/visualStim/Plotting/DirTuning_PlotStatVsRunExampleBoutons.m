%DirTuning_PlotStatVsRunExampleBoutons.m 
% Produces TWO separate PDFs from the compare-states allDirTuning struct:
%   - one page per bouton passing isTunedCVR2_stat  -> ..._Stationary.pdf
%   - one page per bouton passing isTunedCVR2_run   -> ..._Running.pdf
%
% Each page: polar tuning plot + 8 mini-panel raw trial traces (+ mean),
% same style as the other example-bouton plots. No OSI/DSI needed here
% (not computed on this struct) -- annotated with cvR2/p-value instead.

outputFolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section2_Fig4_2\isTunedBoutonsOnly';

%% ===================== stationary =====================
selectedIdx_stat = find([allDirTuning.isTunedCVR2_stat]);
fprintf('Stationary: %d / %d boutons pass isTunedCVR2_stat.\n', numel(selectedIdx_stat), numel(allDirTuning));

pdfPath_stat = fullfile(outputFolder, 'DirTuning_SelectedBoutons_Stationary.pdf');
if exist(pdfPath_stat, 'file'), delete(pdfPath_stat); end

if isempty(selectedIdx_stat)
    warning('No boutons pass isTunedCVR2_stat -- skipping stationary PDF.');
else
    hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
    nSuccess = 0;
    for i = 1:numel(selectedIdx_stat)
        b = selectedIdx_stat(i);
        if ~isvalid(hFig)
            hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
        end
        clf(hFig);
        try
            s = allDirTuning(b);
            plotSingleStateBouton(s, s.fullTrace_stat, s.meanDirResponse_stat, ...
                s.cvR2_stat, s.cvPval_stat, 'Stationary', b, hFig);
            exportgraphics(hFig, pdfPath_stat, 'Append', true);
            nSuccess = nSuccess + 1;
        catch ME
            warning('Skipping bouton %d (stationary) due to error: %s', b, ME.message);
        end
        if mod(i, 25) == 0
            fprintf('Stationary: processed %d / %d (%d successfully plotted)...\n', i, numel(selectedIdx_stat), nSuccess);
        end
    end
    if isvalid(hFig), close(hFig); end
    fprintf('Stationary: done. %d / %d bouton pages written to:\n%s\n', nSuccess, numel(selectedIdx_stat), pdfPath_stat);
end

%% ===================== running =====================
selectedIdx_run = find([allDirTuning.isTunedCVR2_run]);
fprintf('\nRunning: %d / %d boutons pass isTunedCVR2_run.\n', numel(selectedIdx_run), numel(allDirTuning));

pdfPath_run = fullfile(outputFolder, 'DirTuning_SelectedBoutons_Running.pdf');
if exist(pdfPath_run, 'file'), delete(pdfPath_run); end

if isempty(selectedIdx_run)
    warning('No boutons pass isTunedCVR2_run -- skipping running PDF.');
else
    hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
    nSuccess = 0;
    for i = 1:numel(selectedIdx_run)
        b = selectedIdx_run(i);
        if ~isvalid(hFig)
            hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
        end
        clf(hFig);
        try
            s = allDirTuning(b);
            plotSingleStateBouton(s, s.fullTrace_run, s.meanDirResponse_run, ...
                s.cvR2_run, s.cvPval_run, 'Running', b, hFig);
            exportgraphics(hFig, pdfPath_run, 'Append', true);
            nSuccess = nSuccess + 1;
        catch ME
            warning('Skipping bouton %d (running) due to error: %s', b, ME.message);
        end
        if mod(i, 25) == 0
            fprintf('Running: processed %d / %d (%d successfully plotted)...\n', i, numel(selectedIdx_run), nSuccess);
        end
    end
    if isvalid(hFig), close(hFig); end
    fprintf('Running: done. %d / %d bouton pages written to:\n%s\n', nSuccess, numel(selectedIdx_run), pdfPath_run);
end

%% ===================== local function =====================
function plotSingleStateBouton(s, fullTrace, meanDirResponse, cvR2, cvPval, stateLabel, boutonIdx, hFig)
% Draws a full-page figure (polar plot + 8 mini panels) for one bouton,
% one state, using RAW traces (no baseline subtraction).

figure(hFig);

thetaDeg = s.stimValues(:)';
[thetaSorted, sortIdx] = sort(mod(round(thetaDeg), 360));
nDir = numel(thetaSorted);
R_plot = max(meanDirResponse(sortIdx)', 0);

arrowMap = containers.Map( ...
    {0, 45, 90, 135, 180, 225, 270, 315}, ...
    {char(8594), char(8599), char(8593), char(8598), char(8592), char(8601), char(8595), char(8600)});
arrowLabels = cell(1, nDir);
for dIdx = 1:nDir
    if isKey(arrowMap, thetaSorted(dIdx))
        arrowLabels{dIdx} = arrowMap(thetaSorted(dIdx));
    else
        arrowLabels{dIdx} = '?';
    end
end

%% polar plot
axPolar = polaraxes('Position', [0.04 0.30 0.24 0.62]);
thetaRadClosed = deg2rad([thetaSorted, thetaSorted(1)]);
Rclosed = [R_plot, R_plot(1)];
polarplot(axPolar, thetaRadClosed, Rclosed, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 4);
axPolar.ThetaZeroLocation = 'right';
axPolar.ThetaDir = 'counterclockwise';
axPolar.ThetaTick = thetaSorted;
axPolar.ThetaTickLabel = arrayfun(@(x) sprintf('%d%s', x, char(176)), thetaSorted, 'UniformOutput', false);
title(axPolar, sprintf('Bouton %d', boutonIdx));

%% text annotation
annotStr = sprintf('%s\nCross-val R^2 = %.3f\np = %.3g', stateLabel, cvR2, cvPval);
annotation(hFig, 'textbox', [0.30 0.55 0.28 0.35], 'String', annotStr, ...
    'EdgeColor', 'none', 'FontSize', 10, 'VerticalAlignment', 'top');

%% mini PSTH panels
timeVec = s.timeVec(:)';
stimOnDuration = s.stimOnDuration;
displayWindow = [-0.5, 4];
displayMask   = timeVec >= displayWindow(1) & timeVec <= displayWindow(2);
timeVecDisp   = timeVec(displayMask);

allBounds = [];
for dIdx = 1:nDir
    tm = fullTrace{sortIdx(dIdx)};
    tm = tm(:, displayMask);
    if isempty(tm), continue; end
    allBounds = [allBounds, tm(:)']; %#ok<AGROW>
end
yLims = [min(allBounds, [], 'omitnan'), max(allBounds, [], 'omitnan')];
if any(isnan(yLims)) || diff(yLims) == 0
    yLims = [-0.1 0.1];
end

panelLeft   = 0.04;
panelWidth  = 0.90 / nDir;
panelBottom = 0.12;
panelHeight = 0.28;
stimBarHeight = yLims(1) + 0.08 * diff(yLims);

for d = 1:nDir
    axP = axes('Position', [panelLeft + (d-1)*panelWidth, panelBottom, panelWidth*0.9, panelHeight]);
    hold(axP, 'on');

    if ~isnan(stimOnDuration)
        patch(axP, [0 stimOnDuration stimOnDuration 0], ...
            [yLims(1) yLims(1) stimBarHeight stimBarHeight], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
    end

    traceMat = fullTrace{sortIdx(d)};
    traceMat = traceMat(:, displayMask);

    if ~isempty(traceMat)
        hTrials = plot(axP, timeVecDisp, traceMat', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
        for k = 1:numel(hTrials)
            hTrials(k).Color(4) = 0.3;
        end
        plot(axP, timeVecDisp, mean(traceMat, 1, 'omitnan'), 'k-', 'LineWidth', 1.5);
    end

    xline(axP, 0, 'k:', 'LineWidth', 0.75);
    ylim(axP, yLims);
    xlim(axP, displayWindow);
    axis(axP, 'off');

    text(axP, 0.5, -0.18, arrowLabels{d}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 14);
end

% scale bar
targetFrac = 0.4;
rawScaleBarHeight = targetFrac * diff(yLims);
niceSteps = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10];
[~, niceIdx] = min(abs(niceSteps - rawScaleBarHeight));
scaleBarHeight = niceSteps(niceIdx);

axScale = axes('Position', [panelLeft + nDir*panelWidth, panelBottom, 0.05, panelHeight], 'Visible', 'off');
xlim(axScale, [0 1]); ylim(axScale, yLims);
line(axScale, [0.2 0.2], [yLims(1) yLims(1)+scaleBarHeight], 'Color', 'k', 'LineWidth', 1.5);
text(axScale, 0.3, yLims(1) + scaleBarHeight/2, sprintf('%.3g \\DeltaF/F', scaleBarHeight), ...
    'FontSize', 8, 'Units', 'data');

sgtitle(sprintf('Bouton %d | %s | cvR^2=%.3f', boutonIdx, stateLabel, cvR2));

end