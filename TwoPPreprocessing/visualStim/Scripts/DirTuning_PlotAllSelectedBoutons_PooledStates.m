% DirTuning_PlotAllSelectedBoutons_PooledStates.m
%
% Loops through every bouton passing isTunedCVR2_comb (cross-validated
% R^2, stationary + running trials pooled), plots each one using
% plotExampleDirTuningBouton's style (polar tuning plot + mini
% per-direction PSTH traces), and appends one page per bouton to a
% single PDF.
%
% Requires: allDirTuning from DirTuning_IncTuningCurveAnalysis_PoolStates_2PData.m
% (already has isTunedCVR2_comb, meanDirResponse, timeVec, stimOnDuration).
% Requires plotExampleDirTuningBouton.m, computeDirTuningOSI.m,
% computeDirTuningDSI.m on your MATLAB path.

%%
getFieldSafe = @(s, f) isfield(s, f) && s.(f) == 1;
selectedIdx = find([allDirTuning.isTunedCVR2_comb]);
fprintf('%d / %d boutons pass the cross validated r2 (pooled states).\n', ...
    numel(selectedIdx), numel(allDirTuning));
if isempty(selectedIdx)
    error('No boutons pass isTunedCVR2_comb -- nothing to plot.');
end

%%
% Gates OSI/DSI computation on isTunedCVR2_comb -- same criterion used to
% build selectedIdx above, so every bouton in the PDF has a valid,
% non-NaN OSI_simple/DSI_simple to display.
allDirTuning = computeDirTuningOSI(allDirTuning, 'isTunedCVR2_comb', false);
allDirTuning = computeDirTuningDSI(allDirTuning, 'isTunedCVR2_comb', false);

%%
outputFolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section2_Fig4_2\isTunedBoutonsOnly';
pdfPath = fullfile(outputFolder, 'DirTuning_SelectedBoutons_PooledStates.pdf'); % separate filename -- won't overwrite the stationary-only, running-only, or all-trials PDFs
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
        plotExampleDirTuningBouton(allDirTuning, b, hFig, 15);

        % annotate which criteria this bouton passed, for reference
        % (POOLED-STATES field names)
        critStr = {};
        if getFieldSafe(allDirTuning(b), 'isResponsive_comb'),         critStr{end+1} = 'SD-heuristic (pooled)'; end
        if getFieldSafe(allDirTuning(b), 'isResponsive_ttest_comb'),   critStr{end+1} = 't-test (pooled)';        end
        if getFieldSafe(allDirTuning(b), 'isResponsive_ranksum_comb'), critStr{end+1} = 'Wilcoxon (pooled)';       end
        if getFieldSafe(allDirTuning(b), 'isTunedCVR2_comb'),          critStr{end+1} = 'CV_r2 (pooled)';          end
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