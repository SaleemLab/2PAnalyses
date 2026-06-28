function response = getLowMedHigh_SpeedPosActivityMatrix(sessionFileInfo, VRStimName, options)
% Shapes neural data into a [3 x numPosBins x numROIs] Position-Speed matrix
% by pooling frames into Low / Med / High speed categories.
% Output saved under: speedPositionActivity.stratified
%   .speedRanges.low/med/high  [numPosBins x numROIs]
%   .speedEdges                [1 x 4]
%   .method                    'fixed' or 'percentile'
%   .qualityMetrics.low/med/high
%
% Options:
%   .signalToUse        : 'dFFNeuropilCorrected' (default)
%   .useFixedThresholds : true (default) -> edges [1, 15, 27, Inf]
%   .fixedEdges         : [1, 15, 27, Inf] (default)
%   .minFramesPerBin    : 5 (default, used for quality metrics coverage only)

%% Defaults
if nargin < 3, options = struct(); end
if ~isfield(options, 'signalToUse'),        options.signalToUse        = 'dFFNeuropilCorrected'; end
if ~isfield(options, 'useFixedThresholds'), options.useFixedThresholds = true; end
if ~isfield(options, 'fixedEdges'),         options.fixedEdges         = [1, 15, 27, Inf]; end
if ~isfield(options, 'minFramesPerBin'),    options.minFramesPerBin    = 5; end

%% Load data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found.'); end
filePath = sessionFileInfo.stimFiles(stimIdx).Response;

data       = load(filePath, 'lapPositionActivity', 'lapPositionRunningSpeed', 'trialIndicesByCondition');
rawROIData = data.lapPositionActivity.(options.signalToUse); % [numROIs x numTrials x numPosBins]
speedData  = data.lapPositionRunningSpeed;                   % [numTrials x numPosBins]

%% Baseline trials only
baselineTrials = data.trialIndicesByCondition.Baseline;
rawROIData     = rawROIData(:, baselineTrials, :);
speedData      = speedData(baselineTrials, :);

[numROIs, ~, numPosBins] = size(rawROIData);

%% Speed edges
allSpeeds     = speedData(:);
runningMask   = allSpeeds > 1 & ~isnan(allSpeeds);
runningSpeeds = allSpeeds(runningMask);

if options.useFixedThresholds
    speedEdges = options.fixedEdges;
    method     = 'fixed';
else
    speedEdges = prctile(runningSpeeds, [0, 33.33, 66.67, 100]);
    method     = 'percentile';
end

fprintf('Method: %s | Edges: %.1f  %.1f  %.1f  %.1f\n', method, speedEdges(1), speedEdges(2), speedEdges(3), speedEdges(4));

%% Per-frame speed label [numTrials x numPosBins]
speedLabel = zeros(size(speedData));
for s = 1:3
    inBin = speedData >= speedEdges(s) & speedData < speedEdges(s+1);
    speedLabel(inBin) = s;
end
% inclusive upper edge for bin 3
speedLabel(speedData >= speedEdges(3) & speedData <= speedEdges(4)) = 3;

%% Build matrix + quality metrics
binNames       = {'low', 'med', 'high'};
matrix         = nan(3, numPosBins, numROIs);
qualityMetrics = struct();

for s = 1:3
    binMask2D    = speedLabel == s;
    binSpeeds    = speedData(binMask2D);
    nFramesTotal = sum(binMask2D(:));
    posBinCounts = sum(binMask2D, 1);
    coverage     = mean(posBinCounts >= options.minFramesPerBin);

    qm.nFrames      = nFramesTotal;
    qm.coverage     = coverage;
    qm.speed.mean   = mean(binSpeeds,   'omitnan');
    qm.speed.median = median(binSpeeds, 'omitnan');
    qm.speed.std    = std(binSpeeds,    'omitnan');
    qm.speed.min    = min(binSpeeds);
    qm.speed.max    = max(binSpeeds);

    qualityMetrics.(binNames{s}) = qm;

    for b = 1:numPosBins
        trialMask = speedLabel(:, b) == s;
        if sum(trialMask) >= 1
            frameData       = rawROIData(:, trialMask, b);
            matrix(s, b, :) = mean(frameData, 2, 'omitnan');
        end
    end
end

%% Print quality metrics
fprintf('\n--- Quality Metrics ---\n');
for s = 1:3
    qm = qualityMetrics.(binNames{s});
    fprintf('%s | nFrames: %d | coverage: %.0f%% | speed mean: %.1f median: %.1f std: %.1f (min: %.1f max: %.1f)\n', ...
        binNames{s}, qm.nFrames, qm.coverage*100, ...
        qm.speed.mean, qm.speed.median, qm.speed.std, qm.speed.min, qm.speed.max);
end
fprintf('-----------------------\n');

%% Save
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

spa.stratified.speedRanges.low  = squeeze(matrix(1, :, :)); % [numPosBins x numROIs]
spa.stratified.speedRanges.med  = squeeze(matrix(2, :, :));
spa.stratified.speedRanges.high = squeeze(matrix(3, :, :));
spa.stratified.speedEdges       = speedEdges;
spa.stratified.method           = method;
spa.stratified.qualityMetrics   = qualityMetrics;

speedPositionActivity = spa;
disp(['Saving stratified speed-position matrix to ', filePath]);
save(filePath, 'speedPositionActivity', '-append');

response.speedPositionActivity = spa;
end