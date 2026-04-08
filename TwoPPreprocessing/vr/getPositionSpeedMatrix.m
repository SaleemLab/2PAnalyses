function [matrix, speedBinCenters, speedEdges, response] = getPositionSpeedMatrix(sessionFileInfo, VRStimName, options)
% Shapes neural data into a Position-Speed matrix.
% Output: response.speedPositionActivity.matrix [SpeedBins x PositionBins x ROIs]
%       : response.speedEdges
%       : response.speedBinCenters
%       : response.response [this will only contain the above three
%       variables] 
% Options:
% .signalToUse      : 'dFFNeuropilCorrected' (default)
% .OnlyRunningIdx   : true (default, >1cm/s)
% .useQuantileBins  : true (Equal data per bin) or false (Fixed speed width)
% .numBins          : Number of speed rows to create (default 10)
% .applySmoothing   : false (default)
% .smoothSigma      : [1, 1.2] (default)


% Example usage: [response] = getPositionSpeedMatrix(sessionFileInfo, VRStimName, options);
% Example usage with options: 
% options.signalToUse      = 'dFFNeuropilCorrected'; 
% options.useQuantileBins  = false;  % Matches the paper's linear speed ruler
% options.numBins          = 10;     % Create 10 rows of speed
% options.applySmoothing   = true;   % Smooth the "Pillar"
% options.smoothSigma      = [1, 1.2]; % [Speed, Position]



%% Defaults
if nargin < 3, options = struct(); end

if ~isfield(options, 'signalToUse'),     options.signalToUse = 'dFFNeuropilCorrected'; end
if ~isfield(options, 'OnlyRunningIdx'),  options.OnlyRunningIdx = true; end
if ~isfield(options, 'useQuantileBins'), options.useQuantileBins = true; end
if ~isfield(options, 'numBins'),         options.numBins = 10; end
if ~isfield(options, 'applySmoothing'),  options.applySmoothing = false; end
if ~isfield(options, 'smoothSigma'),     options.smoothSigma = [1, 1.2]; end

%% load data 
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found.'); end

filePath = sessionFileInfo.stimFiles(stimIdx).Response;
data = load(filePath, 'lapPositionActivity', 'lapPositionRunningSpeed');

rawROIData = data.lapPositionActivity.(options.signalToUse);
speedData  = data.lapPositionRunningSpeed;

[numROIs, ~, numPosBins] = size(rawROIData);

%% Speed binning 
allSpeeds = speedData(:);
if options.OnlyRunningIdx
    % running index only 
    runningIdx = allSpeeds > 1 & ~isnan(allSpeeds);
else
    % remove nans if present 
    runningIdx = ~isnan(allSpeeds);
end
runningSpeeds = allSpeeds(runningIdx);

if options.useQuantileBins
    % Logic: Divides speed into bins with equal numbers of data points.
    % This "fills in" speed gaps by grouping rare speeds together.
    targetSamplesPerBox = 5;
    nSpeedBins = floor(length(runningSpeeds) / (numPosBins * targetSamplesPerBox));
    nSpeedBins = max(3, min(nSpeedBins, 10));

    if isfield(options, 'numBins') && ~isempty(options.numBins)
        nSpeedBins = options.numBins;
    end
    speedEdges = quantile(runningSpeeds, linspace(0, 1, nSpeedBins + 1));
    fprintf(' quantileBins (%d bins with equal data counts)\n', nSpeedBins);
else
    % Logic: Standard linear spacing 
    % Every bin covers the same number of cm/s (e.g., 5cm/s, 10cm/s).
    nSpeedBins = options.numBins;
    minS = min(runningSpeeds);
    maxS = max(runningSpeeds);
    speedEdges = linspace(minS, maxS, nSpeedBins + 1);
    fprintf('Mode: fixedBins (%d bins with equal speed widths)\n', nSpeedBins);
end

% midpoints for plotting
speedBinCenters = (speedEdges(1:end-1) + speedEdges(2:end)) / 2;

%% create speed-posiiton matrix
matrix = nan(nSpeedBins, numPosBins, numROIs);
for b = 1:numPosBins
    binSpeeds = speedData(:, b);
    for s = 1:nSpeedBins
        idx = binSpeeds >= speedEdges(s) & binSpeeds < speedEdges(s+1);
        if sum(idx) >= 2
            matrix(s, b, :) = mean(rawROIData(:, idx, b), 2, 'omitnan');
        end
    end
end

%% optional smoothning
if options.applySmoothing
    for r = 1:numROIs
        currROI = matrix(:, :, r);
        mask = ~isnan(currROI);
        dataZeroed = currROI;
        dataZeroed(~mask) = 0;

        % This is to deal with nans [gemini]
        blurredData = imgaussfilt(dataZeroed, options.smoothSigma, 'Padding', 'replicate');
        blurredMask = imgaussfilt(double(mask), options.smoothSigma, 'Padding', 'replicate');
        matrix(:, :, r) = blurredData ./ blurredMask;
    end
end

%% append matrix to response 
response.speedPositionActivity.matrix = matrix;
response.speedPositionActivity.speedBinCenters = speedBinCenters;
response.speedPositionActivity.speedEdges = speedEdges;

disp(['Saving speed-posiiton-activity matrix to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
end