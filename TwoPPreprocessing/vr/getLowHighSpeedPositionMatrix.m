function response = getLowHighSpeedPositionMatrix(sessionFileInfo, VRStimName, options)
% Computes per-position low and high speed activity matrices for each ROI.
% Low/high split is defined by the per-position occupancy-based median —
% guaranteeing equal frame counts per group at every position bin.
%
% Output saved under: speedPositionActivity.lowHigh
%   .matrixLow        [numROIs x numPosBins] — mean activity in low speed frames
%   .matrixHigh       [numROIs x numPosBins] — mean activity in high speed frames
%   .matrixDiff       [numROIs x numPosBins] — matrixLow - matrixHigh
%   .medianThreshLine [1 x numPosBins]       — per-position median speed threshold
%   .sessionMedian    scalar                 — median speed across all valid frames
%   .sessionIQR       scalar                 — IQR across all valid frames
%
% Options:
%   .signalToUse : 'dFFNeuropilCorrected' (default)
%   .minSpeed    : 1 (default, exclude stationary frames >1 cm/s)
%
% Example usage:
%   options.signalToUse = 'dFFNeuropilCorrected';
%   response = getLowHighSpeedPositionMatrix(sessionFileInfo, VRStimName, options);

%% Defaults
if nargin < 3, options = struct(); end
if ~isfield(options, 'signalToUse'), options.signalToUse = 'dFFNeuropilCorrected'; end
if ~isfield(options, 'minSpeed'),    options.minSpeed    = 1; end

%% Find stimulus
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found.'); end
filePath = sessionFileInfo.stimFiles(stimIdx).Response;

%% Load data
data       = load(filePath, 'lapPositionActivity', 'lapPositionRunningSpeed', 'trialIndicesByCondition');
rawROIData = data.lapPositionActivity.(options.signalToUse); % [numROIs x numTrials x numPosBins]
speedData  = data.lapPositionRunningSpeed;                   % [numTrials x numPosBins]

%% Baseline trials only
baselineTrials = data.trialIndicesByCondition.Baseline;
rawROIData     = rawROIData(:, baselineTrials, :);
speedData      = speedData(baselineTrials, :);

[numROIs, ~, numPosBins] = size(rawROIData);

%% Session-level speed summary
allSpeeds     = speedData(:);
validSpeeds   = allSpeeds(allSpeeds > options.minSpeed & ~isnan(allSpeeds));
sessionMedian = median(validSpeeds);
sessionIQR    = iqr(validSpeeds);
fprintf('  Session median speed: %.1f cm/s | IQR: %.1f cm/s\n', sessionMedian, sessionIQR);

%% Per-position occupancy-based median threshold
medianThreshLine = nan(1, numPosBins);

for b = 1:numPosBins
    binSpeeds = speedData(:, b);
    validBin  = binSpeeds(binSpeeds > options.minSpeed & ~isnan(binSpeeds));
    if ~isempty(validBin)
        medianThreshLine(b) = median(validBin);
    end
end

%% Build low / high matrices
% 
matrixLow  = nan(numROIs, numPosBins);
matrixHigh = nan(numROIs, numPosBins);

for b = 1:numPosBins
    if isnan(medianThreshLine(b)), continue; end

    binSpeeds = speedData(:, b);
    lowIdx    = binSpeeds > options.minSpeed & binSpeeds <= medianThreshLine(b) & ~isnan(binSpeeds);
    highIdx   = binSpeeds > medianThreshLine(b) & ~isnan(binSpeeds);

    if sum(lowIdx) >= 2
        matrixLow(:, b)  = mean(squeeze(rawROIData(:, lowIdx, b)), 2, 'omitnan');
    end
    if sum(highIdx) >= 2
        matrixHigh(:, b) = mean(squeeze(rawROIData(:, highIdx, b)), 2, 'omitnan');
    end
end

%% Difference matrix
matrixDiff = matrixLow - matrixHigh;

%% Package output
lh.matrixLow        = matrixLow;
lh.matrixHigh       = matrixHigh;
lh.matrixDiff       = matrixDiff;
lh.medianThreshLine = medianThreshLine;
lh.sessionMedian    = sessionMedian;
lh.sessionIQR       = sessionIQR;

%% Save to file — appends to existing speedPositionActivity struct
if isfile(filePath)
    existing = load(filePath, 'speedPositionActivity');
    if isfield(existing, 'speedPositionActivity')
        spa = existing.speedPositionActivity;
    else
        spa = struct();
    end
else
    spa = struct();
end

spa.lowHigh           = lh;
speedPositionActivity = spa;

fprintf('  Saving lowHigh matrices to %s\n', filePath);
save(filePath, 'speedPositionActivity', '-append');

response.speedPositionActivity = spa;
end
