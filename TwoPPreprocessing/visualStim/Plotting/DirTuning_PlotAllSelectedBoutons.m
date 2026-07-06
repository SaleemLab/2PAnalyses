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

isSelectedAny = arrayfun(@(s) getFieldSafe(s, 'isResponsive') || ...
                              getFieldSafe(s, 'isResponsive_ttest') || ...
                              getFieldSafe(s, 'isResponsive_ranksum'), ...
                          allDirTuning);
% NOTE: cross-val R^2 (isTunedCV) deliberately excluded -- too many
% boutons pass it with only a marginal, non-substantive effect size
% (no minimum R^2 magnitude requirement, only a p-value), producing
% low-quality example pages.

for b = 1:numel(allDirTuning)
    allDirTuning(b).isSelectedAny = isSelectedAny(b);
end

selectedIdx = find(isSelectedAny);
fprintf('%d / %d boutons pass at least one criterion (SD-heuristic OR t-test OR Wilcoxon).\n', ...
    numel(selectedIdx), numel(allDirTuning));

if isempty(selectedIdx)
    error('No boutons pass any criterion -- nothing to plot.');
end

%%
% Uses isSelectedAny (not just one specific criterion) so every bouton
% in the PDF has a valid, non-NaN OSI/DSI to display, regardless of
% which specific criterion it passed.
allDirTuning = computeDirTuningOSI(allDirTuning, 'isSelectedAny', false);
allDirTuning = computeDirTuningDSI(allDirTuning, 'isSelectedAny', false);

%% 
outputFolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section2_Fig4_2\selectedBoutons'; 
pdfPath = fullfile(outputFolder, 'DirTuning_SelectedBoutons.pdf');
if exist(pdfPath, 'file')
    delete(pdfPath);
end

hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
nSuccess = 0;

for i = 1:numel(selectedIdx)
    b = selectedIdx(i);

    % self-healing: if hFig got deleted/corrupted by a previous error,
    % recreate it rather than letting every subsequent iteration fail
    % against the same broken handle
    if ~isvalid(hFig)
        warning('hFig became invalid before bouton %d -- recreating figure.', b);
        hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 550]);
    end

    try
        plotExampleDirTuningBouton(allDirTuning, b, hFig);

        % annotate which criteria this bouton passed, for reference
        critStr = {};
        if getFieldSafe(allDirTuning(b), 'isResponsive'),         critStr{end+1} = 'SD-heuristic'; end 
        if getFieldSafe(allDirTuning(b), 'isResponsive_ttest'),   critStr{end+1} = 't-test';        end 
        if getFieldSafe(allDirTuning(b), 'isResponsive_ranksum'), critStr{end+1} = 'Wilcoxon';       end 
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
