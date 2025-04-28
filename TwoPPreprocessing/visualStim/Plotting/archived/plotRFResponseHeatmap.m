function plotRFResponseHeatmap(sessionFileInfo, RFMappingStimName, azLim, elLim)
% Generates a heatmap (ROI x time) of RF responses averaged over az/el window
% Inputs:
%   - sessionFileInfo : struct
%   - RFMappingStimName : string (e.g. 'M25041_RFMapping_20250417_00001')
%   - azLim : [min max] azimuth window (e.g. [-80 -25])
%   - elLim : [min max] elevation window (e.g. [-30 30])

if nargin < 3
    azLim = [-80 -25];
end
if nargin < 4
    elLim = [-30 30];
end

% --- Load RF response
iStim = find(strcmp(RFMappingStimName, {sessionFileInfo.stimFiles.name}), 1);
RFResponseStruct = load(sessionFileInfo.stimFiles(iStim).Response, 'response');
RFResponse = RFResponseStruct.response;
psthData = RFResponse.psthData;

tVec = psthData(1).timeVector;
nTime = length(tVec);
nROIs = size(psthData(1).alignedResponses, 1);

% --- Collect RF responses within desired az/el
RF_matrix = NaN(nROIs, nTime);

for roiIdx = 1:nROIs
    roiTraces = [];
    for i = 1:length(psthData)
        stimPos = psthData(i).stimValue;
        if stimPos(1) >= azLim(1) && stimPos(1) <= azLim(2) && ...
           stimPos(2) >= elLim(1) && stimPos(2) <= elLim(2)

            aligned = squeeze(psthData(i).alignedResponses(roiIdx, :, :));
            if size(aligned, 1) > size(aligned, 2), aligned = aligned'; end
            roiTraces = cat(1, roiTraces, aligned); % trials x time
        end
    end
    if ~isempty(roiTraces)
        RF_matrix(roiIdx, :) = mean(roiTraces, 1, 'omitnan');
    end
end

% --- Normalize (optional)
RF_matrix_z = normalize(RF_matrix, 2);  % z-score or range normalisation

% --- Sort by peak response time
[~, peakIdx] = max(RF_matrix_z, [], 2, 'omitnan');
[~, sortOrder] = sort(peakIdx);
RF_sorted = RF_matrix_z(sortOrder, :);

% --- Plot heatmap
figure('Position', [100 100 1000 600]);
imagesc(tVec, 1:nROIs, RF_sorted);
colormap(redWhiteBlue);
colorbar;
xline(0, 'k--', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('ROIs (sorted by peak time)');
title(sprintf('RF Responses | Az: [%d %d], El: [%d %d]', azLim(1), azLim(2), elLim(1), elLim(2)));
set(gca, 'FontSize', 12);

end
