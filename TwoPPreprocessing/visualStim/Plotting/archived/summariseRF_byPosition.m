function summariseRF_byPosition(sessionFileInfo, RFMappingStimName, prestimWindow, poststimWindow, threshold)
% summarizeRF_byPosition: summary of RF responses across population per stimulus position
%   Computes mean and count of responsive ROIs at each azimuth/elevation.
% Inputs:
%   sessionFileInfo - struct with .stimFiles and .Directories.save_folder
%   RFMappingStimName - string for RF mapping stimulus name
%   prestimWindow - [t0 t1] baseline period (default [-0.5 0])
%   poststimWindow - [t0 t1] response window (default [0 0.5])
%   threshold - ΔF/F threshold to count responsive ROIs (default 0.1)

if nargin<3, prestimWindow = [-0.5 0]; end
if nargin<4, poststimWindow = [0 0.5]; end

% Interpret scalar windows: allow scalar specifying duration
if isscalar(prestimWindow)
    prestimWindow = [-abs(prestimWindow), 0];
end
if isscalar(poststimWindow)
    poststimWindow = [0, abs(poststimWindow)];
end
if nargin<5, threshold = 0.1; end

% Load psthData
iStim = find(strcmp(RFMappingStimName,{sessionFileInfo.stimFiles.name}),1);
assert(~isempty(iStim),'Stimulus not found');
resp = load(sessionFileInfo.stimFiles(iStim).Response,'response');
psthData = resp.response.psthData;

% Extract unique azimuth and elevation
stimVals = vertcat(psthData.stimValue);
az = unique(stimVals(:,1));
el = unique(stimVals(:,2));

% Initialize summary matrices
meanRespMap = nan(numel(el),numel(az));
countRespMap = zeros(numel(el),numel(az));

% Time vector and masks
tVec = psthData(1).timeVector;
postMask = tVec >= poststimWindow(1) & tVec <= poststimWindow(2);
baselineMask = tVec >= prestimWindow(1) & tVec <= prestimWindow(2);

nROI = size(psthData(1).alignedResponses,1);

% Loop over positions
for k = 1:numel(psthData)
    pos = psthData(k).stimValue;
    aIdx = find(az==pos(1));
    eIdx = find(el==pos(2));
    if isempty(aIdx)||isempty(eIdx), continue; end
    
    % Extract data: ROI x time x trials
    data = squeeze(psthData(k).alignedResponses); % [ROI x time x trials]
    if size(data,1)==numel(tVec) && size(data,2)==nROI
        data = permute(data,[2 1 3]); % ensure [ROI x time x trials]
    end
    % Compute ΔF/F relative to baseline for each ROI/trial
    baseline = mean(data(:,baselineMask,:),2); % ROI x 1 x trials
    df_f = (data - baseline) ./ baseline;           
    % Average across time in post window, flatten trials
    roiResp = squeeze(mean(df_f(:,postMask,:),2)); % ROI x trials
    % Compute mean response across all ROI & all trials
    meanRespMap(eIdx,aIdx) = mean(roiResp(:),'omitnan');
    % Count responsive ROI-trials > threshold
    countRespMap(eIdx,aIdx) = sum(roiResp(:) > threshold);
end

% Figure 1: Mean response heatmap
figure('Position',[100 100 800 600]);
imagesc(az,el,meanRespMap);
axis xy;
xlabel('Azimuth (°)'); ylabel('Elevation (°)');
title('Mean RF response across population');
colormap(redWhiteBlue); colorbar; grid off;

% Figure 2: Count heatmap
figure('Position',[100 100 800 600]);
imagesc(az,el,countRespMap);
axis xy;
xlabel('Azimuth (°)'); ylabel('Elevation (°)');
title(sprintf('Count of ROI-trials > %.2f',threshold));
colormap(parula); colorbar; grid off;

end
