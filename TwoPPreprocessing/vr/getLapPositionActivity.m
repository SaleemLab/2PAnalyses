function [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, VRStimName, overwrite, onlyIncludeROIs, useZScoredProcessedSignals)

% Calculates binned lap activity (Neurons x Laps x Bins) for a single stimulus.
% Returns the unshifted signal and its time indices for stitching/shuffling later.
% Using modified version Dec 2025. 

%% Handle optional inputs
if nargin < 3, overwrite = true; end
if nargin < 4, onlyIncludeROIs = false; end
if nargin < 5, useZScoredProcessedSignals = true; end

%% Find VR stimulus and load data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
disp(['Loading data for stimulus: ', VRStimName]);
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');
load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');

%% Overwrite check (only remove activity fields, keep shuffle matrix if it exists)
if overwrite && isfield(response, 'lapPositionActivity')
    fieldsToRemove = {'lapPositionActivity', 'lapPositionActivity_meanShift', 'cellROIs'};
    response = rmfield(response, intersect(fieldsToRemove, fieldnames(response)));
end

%% Get cell ROIs
if onlyIncludeROIs
    ROIs = find(processedTwoPData.iscell(:, 1));
else
    ROIs = 1:size(processedTwoPData.iscell, 1);
end
response.cellROIs = ROIs;

%% Select appropriate signal and extract matrix
if useZScoredProcessedSignals
    signals = processedTwoPData.zScoredProcessedSignals;
    response.signalsZScored = true;
else
    signals = processedTwoPData.processedSignals;
    response.signalsZScored = false;
end
signalNames = fieldnames(signals);
numSignals = length(signalNames);
numBins = 140;

% Initialize output for signal matrix (Cell X Time) and time indices
signalMatrix = struct();
timeIndices = processedTwoPData.TwoPFrameTime;
% totalFrames = length(timeIndices);
nLaps = size(response.lapPosition2PFrameIdx, 1);

lapPositionActivity = struct(); % To store the real (unshuffled) activity (Neurons x Laps x Bins)

%% Calculate Real Activity (Neurons x Laps x Bins)
disp('Calculating Real Activity (Neurons x Laps x Bins)...');
for iSignal = 1:numSignals
    currentSignalName = signalNames{iSignal};
    currentSignalMatrix = signals.(currentSignalName)(ROIs, :); % Select only the ROIs
    numROIs = size(currentSignalMatrix, 1);

    % Store the full signal matrix (Cell x Time) for later stitching
    signalMatrix.(currentSignalName) = currentSignalMatrix;

    % Initialize storage
    lapPositionActivity.(currentSignalName) = nan(numROIs, nLaps, numBins);

    % Calculate activity
    for thisLap = 1:nLaps
        for thisBin = 1:numBins
            frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
            if ~isempty(frameIdx)
                meanActivityReal = mean(currentSignalMatrix(:, frameIdx), 2, 'omitnan');
                lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = meanActivityReal;
            end
        end
    end
end

response.lapPositionActivity = lapPositionActivity;

disp(['Saving updated Response struct (Lap Activity) to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-append');
end