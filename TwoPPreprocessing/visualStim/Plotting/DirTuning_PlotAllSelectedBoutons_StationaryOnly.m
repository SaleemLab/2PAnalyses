% DirTuning_PlotAllSelectedBoutons_StationaryOnly.m
%
% Loops through every bouton passing isTunedCVR2_stat (cross-validated
% R^2, stationary trials only), plots each one using
% plotExampleDirTuningBouton's style (polar tuning plot + mini
% per-direction PSTH traces), and appends one page per bouton to a
% single PDF.
%
% Requires: allDirTuning from DirTuning_IncTuningCurveAnalysis_StationaryOnly_2PData.m
% (already has isTunedCVR2_stat, fullTraceSub, timeVec, stimOnDuration).
% Requires plotExampleDirTuningBouton.m, computeDirTuningOSI.m,
% computeDirTuningDSI.m on your MATLAB path.

%%
getFieldSafe = @(s, f) isfield(s, f) && s.(f) == 1;

selectedIdx = find([allDirTuning.isTunedCVR2_stat]);
fprintf('%d / %d boutons pass the cross validated r2 (stationary).\n', ...
    numel(selectedIdx), numel(allDirTuning));

if isempty(selectedIdx)
    error('No boutons pass isTunedCVR2_stat -- nothing to plot.');
end

%%
% Gates OSI/DSI computation on isTunedCVR2_stat -- same criterion used to
% build selectedIdx above, so every bouton in the PDF has a valid,
% non-NaN OSI_simple/DSI_simple to display.
allDirTuning = computeDirTuningOSI(allDirTuning, 'isTunedCVR2_stat', false);
allDirTuning = computeDirTuningDSI(allDirTuning, 'isTunedCVR2_stat', false);

%%
outputFolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section2_Fig4_2\isTunedBoutonsOnly';
pdfPath = fullfile(outputFolder, 'DirTuning_SelectedBoutons_StationaryOnly.pdf'); % separate filename -- won't overwrite the running-only or all-trials PDFs
if exist(pdfPath, 'file')
    delete(pdfPath);
end

hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
nSuccess = 0;

for i = 1:numel(selectedIdx)
    b = selectedIdx(i);

    if ~isvalid(hFig)
        warning('hFig became invalid before bouton %d -- recreating figure.', b);
        hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
    end

    try
        plotExampleDirTuningBouton(allDirTuning, b, hFig);

        % annotate which criteria this bouton passed, for reference
        % (STATIONARY-ONLY field names)
        critStr = {};
        if getFieldSafe(allDirTuning(b), 'isResponsive_stat'),         critStr{end+1} = 'SD-heuristic (stat)'; end
        if getFieldSafe(allDirTuning(b), 'isResponsive_ttest_stat'),   critStr{end+1} = 't-test (stat)';        end
        if getFieldSafe(allDirTuning(b), 'isResponsive_ranksum_stat'), critStr{end+1} = 'Wilcoxon (stat)';       end
        if getFieldSafe(allDirTuning(b), 'isTunedCVR2_stat'),          critStr{end+1} = 'CV_r2 (stat)';          end
        annotation(hFig, 'textbox', [0.05 0.95 0.9 0.04], 'String', ...
            sprintf('Passed: %s', strjoin(critStr, ', ')), ...
            'EdgeColor', 'none', 'FontSize', 8, 'Interpreter', 'none');

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
