% DirTuning_PlotAllSelectedBoutons.m
%
% Loops through every bouton that passed ANY of the built
% responsiveness/tuning criteria (SD-heuristic, ANOVA-protected t-test,
% ANOVA-protected Wilcoxon, cross-validated R^2), plots each one using
% plotExampleDirTuningBouton's style (polar tuning plot + mini
% per-direction PSTH traces), and appends one page per bouton to a
% single PDF.
%
% Requires: allDirTuning already has isResponsive, isResponsive_ttest,
% isResponsive_ranksum, isTunedCV computed (from the main pooling +
% responsiveness script), and fullTraceSub/timeVec/stimOnDuration
% (updated pooling script). Requires plotExampleDirTuningBouton.m,
% computeDirTuningOSI.m, computeDirTuningDSI.m on your MATLAB path.

%%
getFieldSafe = @(s, f) isfield(s, f) && s.(f) == 1;

% isSelectedAny = arrayfun(@(s) getFieldSafe(s, 'isResponsive') || ...
%                               getFieldSafe(s, 'isResponsive_ttest') || ...
%                               getFieldSafe(s, 'isResponsive_ranksum'), ...
%                           allDirTuning);
% NOTE: cross-val R^2 (isTunedCV) deliberately excluded -- too many
% boutons pass it with only a marginal, non-substantive effect size
% (no minimum R^2 magnitude requirement, only a p-value), producing
% low-quality example pages.


% for b = 1:numel(allDirTuning)
%     allDirTuning(b).isTunedCVR2 = isSelectedAny(b);
% end

selectedIdx = find([allDirTuning.isTunedCVR2_run]);
fprintf('%d / %d boutons pass the cross validated r2).\n', ...
    numel(selectedIdx), numel(allDirTuning));

if isempty(selectedIdx)
    error('No boutons pass any criterion -- nothing to plot.');
end

%%
% Gates OSI/DSI computation on isTunedCVR2_run (cross-validated R^2,
% running trials only) -- same criterion used to build selectedIdx
% above, so every bouton in the PDF has a valid, non-NaN OSI_simple/
% DSI_simple to display.
allDirTuning = computeDirTuningOSI(allDirTuning, 'isTunedCVR2_run', false);
allDirTuning = computeDirTuningDSI(allDirTuning, 'isTunedCVR2_run', false);

%% 
outputFolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section2_Fig4_2\isTunedBoutonsOnly'; 
pdfPath = fullfile(outputFolder, 'DirTuning_SelectedBoutons_RunningOnly.pdf'); % RUNNING-ONLY -- separate filename so it doesn't overwrite the all-trials PDF
if exist(pdfPath, 'file')
    delete(pdfPath);
end

hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
nSuccess = 0;

for i = 1:numel(selectedIdx)
    b = selectedIdx(i);

    %  if hFig got deleted/corrupted by a previous error,
    % recreate it rather than letting every subsequent iteration fail
    % against the same broken handle
    if ~isvalid(hFig)
        warning('hFig became invalid before bouton %d -- recreating figure.', b);
        hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
    end

    try
        plotExampleDirTuningBouton(allDirTuning, b, hFig, 5);

        % annotate which criteria this bouton passed, for reference
        % (RUNNING-ONLY field names -- the non-running-suffixed fields
        % don't exist on a running-only allDirTuning, so checking those
        % here would silently always be false)
        critStr = {};
        if getFieldSafe(allDirTuning(b), 'isResponsive_run'),         critStr{end+1} = 'SD-heuristic (run)'; end
        if getFieldSafe(allDirTuning(b), 'isResponsive_ttest_run'),   critStr{end+1} = 't-test (run)';        end
        if getFieldSafe(allDirTuning(b), 'isResponsive_ranksum_run'), critStr{end+1} = 'Wilcoxon (run)';       end
        if getFieldSafe(allDirTuning(b), 'isTunedCVR2_run'),          critStr{end+1} = 'CV_r2 (run)';          end
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
