function pdfPath = plotAllROIs_ConditionStack_PDF(sessionFileInfo, response, signalToUse, roiList)
% plotAllROIs_ConditionStack_PDF: Loops over ROIs, generates the horizontal
% condition-stack figure for each one (via plotSingleROI_ConditionStack),
% and appends every figure as a single page into one multi-page PDF

if nargin < 3 || isempty(signalToUse), signalToUse = 'dFFNeuropilCorrected'; end

data = response.lapPositionActivity.(signalToUse);
numROIsTotal = size(data, 1);

if nargin < 4 || isempty(roiList)
    roiList = 1:numROIsTotal;
end

saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ROI_Summaries');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
pdfPath = fullfile(saveDir, sprintf('%s_AllROIs_ConditionStack_%s.pdf', ...
    sessionFileInfo.session_name, signalToUse));

% Start fresh: exportgraphics 'Append' just keeps adding pages, so delete
% any pre-existing file with the same name first.
if exist(pdfPath, 'file')
    delete(pdfPath);
end

for idx = 1:length(roiList)
    neuronIdx = roiList(idx);
    fprintf('Plotting ROI %d (%d/%d)...\n', neuronIdx, idx, length(roiList));

    % saveIndividualPNG = false -> skip per-ROI PNG, we only want the PDF
    fig = plotSingleROI_ConditionStack(sessionFileInfo, response, neuronIdx, signalToUse);

    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType', 'vector');

    close(fig);
end

fprintf('Done. Saved %d-page PDF to:\n%s\n', length(roiList), pdfPath);

end