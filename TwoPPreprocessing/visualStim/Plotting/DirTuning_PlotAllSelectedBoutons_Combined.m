% DirTuning_PlotAllSelectedBoutons_Combined.m
%
% Loops through every bouton passing isTunedCVR2_comb (cross-validated
% R^2, stationary+running trials pooled), plots each one (polar tuning
% plot + 8 mini-panel raw trial traces), and appends one page per bouton
% to a single PDF.
%
% Requires allDirTuning from
% DirTuning_IncTuningCurveAnalysis_CompareStates_2PData.m (has
% fullTrace_stat, fullTrace_run, timeVec, stimOnDuration, stimValues,
% cvR2_comb, cvPval_comb, isTunedCVR2_comb).
%
% The combined trace shown here is built on the fly by concatenating
% fullTrace_stat{d} and fullTrace_run{d} (trials) for each direction --
% matching exactly what fed into the combined cvR2 computation.

%%
selectedIdx = find([allDirTuning.isTunedCVR2_comb]);
fprintf('%d / %d boutons pass isTunedCVR2_comb.\n', numel(selectedIdx), numel(allDirTuning));

if isempty(selectedIdx)
    error('No boutons pass isTunedCVR2_comb -- nothing to plot.');
end

%%
outputFolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section2_Fig4_2\isTunedBoutonsOnly';
pdfPath = fullfile(outputFolder, 'DirTuning_SelectedBoutons_Combined.pdf');
if exist(pdfPath, 'file')
    delete(pdfPath);
end

hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
nSuccess = 0;

for i = 1:numel(selectedIdx)
    b = selectedIdx(i);
    if ~isvalid(hFig)
        hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
    end
    clf(hFig);

    try
        s = allDirTuning(b);
        nDir = numel(s.fullTrace_stat);

        % build the combined trace on the fly: concatenate stat+run
        % trials per direction (same trials that fed cvR2_comb)
        fullTrace_comb = cell(nDir, 1);
        meanDirResponse_comb = nan(nDir, 1);
        respIdxLocal = true(1, numel(s.timeVec)); % placeholder, replaced below using stimOnDuration-based window
        for d = 1:nDir
            fullTrace_comb{d} = [s.fullTrace_stat{d}; s.fullTrace_run{d}];
        end

        % recompute meanDirResponse_comb using the same resp window
        % logic (respWin) is not stored on s directly, so approximate
        % it using the same [0.1, stimOnDuration+... ] convention is not
        % available here -- instead just average meanDirResponse_stat
        % and meanDirResponse_run weighted by trial count, which is
        % algebraically equivalent to averaging over the pooled trials.
        nStat = cellfun(@(x) size(x,1), s.fullTrace_stat);
        nRun  = cellfun(@(x) size(x,1), s.fullTrace_run);
        meanDirResponse_comb = (s.meanDirResponse_stat .* nStat(:) + s.meanDirResponse_run .* nRun(:)) ./ (nStat(:) + nRun(:));

        plotSingleStateBouton(s, fullTrace_comb, meanDirResponse_comb, ...
            s.cvR2_comb, s.cvPval_comb, 'Combined (stat+run)', b, hFig);

        exportgraphics(hFig, pdfPath, 'Append', true);
        nSuccess = nSuccess + 1;
    catch ME
        warning('Skipping bouton %d due to error: %s', b, ME.message);
    end

    if mod(i, 25) == 0
        fprintf('Processed %d / %d boutons (%d successfully plotted so far)...\n', i, numel(selectedIdx), nSuccess);
    end
end

if isvalid(hFig)
    close(hFig);
end
fprintf('\nDone. %d / %d bouton pages successfully written to:\n%s\n', nSuccess, numel(selectedIdx), pdfPath);

%% ===================== local function =====================
function plotSingleStateBouton(s, fullTrace, meanDirResponse, cvR2, cvPval, stateLabel, boutonIdx, hFig)
% Draws a full-page figure (polar plot + 8 mini panels) for one bouton,
% using RAW traces (no baseline subtraction).

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

    if d == 1
        axP.XColor = 'none';
        axP.YColor = 'k';
        axP.Box    = 'off';
        axP.TickDir = 'out';
        axP.FontSize = 7;
        ylabel(axP, '\DeltaF/F', 'FontSize', 8);
    else
        axis(axP, 'off');
    end

    text(axP, 0.5, -0.18, arrowLabels{d}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 14);
end

sgtitle(sprintf('Bouton %d | %s | cvR^2=%.3f', boutonIdx, stateLabel, cvR2));

end
