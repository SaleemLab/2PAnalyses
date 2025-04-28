function psthData = getTrialResponsePSTHsV3(sessionFileInfo, stimName, method, interpRate, frameRate)
%
% Extracts trial-aligned neural responses for different stimulus conditions.
%
% Inputs:
%   - sessionFileInfo (struct): Contains session/sessionFile info.
%   - stimName (string): The name of the stimulus condition.
%   - method (int): PSTH calculation method:
%       1: Interpolate & average
%       2: Temporal smoothing
%       3: Direct mean response (no interpolation)
%   - interpRate (float): Interpolation rate in Hz (default: 60 Hz)
%
% Output:
%   - psthData (struct): One entry per unique stimulus value, each with:
%       - stimValue, alignedResponses, meanResponse, stdResponse,
%         semResponse, timeVector, responseType, neuronPlaneIndex

if nargin < 4
    interpRate = 60;  % Default interpolation rate
end

if nargin < 5
    frameRate = 7.28;  % Default interpolation rate
end

if nargin <6 
    computeDeltaFOverF = true;
end

if nargin <7 
    neuropilCorrection = true;
end
% Locate the requested stimulus
isStim = strcmp(stimName, {sessionFileInfo.stimFiles.name});
iStim = find(isStim, 1);

if isempty(iStim)
    error('Stimulus name "%s" not found in sessionFileInfo.', stimName);
end

% Load relevant data
if exist(sessionFileInfo.stimFiles(iStim).Response, 'file') && ...
   exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file')
    load(sessionFileInfo.stimFiles(iStim).Response, 'response');
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
else
    error('Missing Response or BonsaiData for stimulus "%s".', stimName);
end

if exist(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'file')
    load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData');
else
    error('Missing twoPData for stimulus "%s".', stimName);
end

% Build full neural matrix
neuronPlaneIndex = [];
allNeurons = [];

for thisPlane = 1:length(twoPData)
    isRoiArray = twoPData(thisPlane).iscell;
    isRoiIdx = find(isRoiArray(:, 1) == 1);
    rawFRois = twoPData(thisPlane).F(isRoiIdx, :);
    numNeurons = size(rawFRois, 1);
    allNeurons = [allNeurons; rawFRois];
    neuronPlaneIndex = [neuronPlaneIndex; repmat(thisPlane, numNeurons, 1)];
end

totalNeurons = size(allNeurons, 1);
numStimuli = length(bonsaiData.trialGroups);

% Preallocate output struct with consistent fields
psthData = repmat(struct('stimValue', [], ...
                         'alignedResponses', [], ...
                         'meanResponse', [], ...
                         'stdResponse', [], ...
                         'semResponse', [], ...
                         'timeVector', [], ...
                         'responseType', [], ...
                         'neuronPlaneIndex', neuronPlaneIndex), numStimuli, 1);

for i = 1:numStimuli
    stimValue = bonsaiData.trialGroups(i).value;
    trialIndices = bonsaiData.trialGroups(i).trials;
    trialIndices = trialIndices(~isnan(trialIndices));
    numTrials = length(trialIndices);

    % Determine global trial time range
    allMaxTimes = [];
    allMinTimes = [];
    for thisPlane = 1:length(response)
        maxTime = max(cellfun(@max, response(thisPlane).responseFrameRelTimes(trialIndices)));
        minTime = min(cellfun(@min, response(thisPlane).responseFrameRelTimes(trialIndices)));
        allMaxTimes = [allMaxTimes, maxTime];
        allMinTimes = [allMinTimes, minTime];
    end
    maxTime = max(allMaxTimes);
    minTime = min(allMinTimes);
    timeVector = linspace(minTime, maxTime, interpRate);

    alignedResponses = NaN(totalNeurons, length(timeVector), numTrials);
    neuronOffset = 0;

    for thisPlane = 1:length(twoPData)
        isRoiArray = twoPData(thisPlane).iscell;
        isRoiIdx = find(isRoiArray(:, 1) == 1);
        
        rawFRois = twoPData(thisPlane).F(isRoiIdx, :);
        f0 = get_F0(rawFRois', frameRate, 10);
        deltaF_F = get_delta_F_over_F(rawFRois', f0)';

        numNeurons = size(deltaF_F, 1);

        for t = 1:numTrials
            trialIdx = trialIndices(t);
            selectedFrames = response(thisPlane).responseFrameIdx{trialIdx};
            trialData = deltaF_F(:, selectedFrames);
            trialTimes = response(thisPlane).responseFrameRelTimes{trialIdx};

            for neuronIdx = 1:numNeurons
                globalNeuronIdx = neuronOffset + neuronIdx;
                alignedResponses(globalNeuronIdx, :, t) = interp1(trialTimes, ...
                    trialData(neuronIdx, :), timeVector, 'linear', NaN);
            end
        end
        neuronOffset = neuronOffset + numNeurons;
    end

    % Compute PSTH
    switch method
        case 1  % Interpolate & average
            psthMean = nanmean(alignedResponses, 3);
            psthStd = nanstd(alignedResponses, 0, 3);

        case 2  % Temporal smoothing
            w = gausswin(100); w = w / sum(w);
            psthMean = nanmean(alignedResponses, 3);
            nanMeanValue = nanmean(psthMean(:));
            psthMean(isnan(psthMean)) = nanMeanValue;
            psthMean = filtfilt(w, 1, psthMean')';
            psthStd = nanstd(alignedResponses, 0, 3);
            psthStd = filtfilt(w, 1, psthStd')';

        case 3  % Direct mean (no interpolation)
            maxFrames = max(cellfun(@(x) sum(x), response(1).responseFrameIdx(trialIndices)));
            directMeanResponse = NaN(totalNeurons, maxFrames, numTrials);

            for t = 1:numTrials
                trialIdx = trialIndices(t);
                selectedFrames = response(1).responseFrameIdx{trialIdx};
                trialData = allNeurons(:, selectedFrames);
                directMeanResponse(:, 1:size(trialData, 2), t) = trialData;
            end

            psthMean = nanmean(directMeanResponse, 3);
            psthStd = nanstd(directMeanResponse, 0, 3);
            timeVector = 1:size(psthMean, 2);  % crude frame-based time

        otherwise
            error('Invalid method. Use 1 (Interpolate), 2 (Smooth), or 3 (Direct)');
    end

    psthSEM = psthStd ./ sqrt(sum(~isnan(alignedResponses), 3));

    % Store all in preallocated struct
    psthData(i).stimValue = stimValue;
    psthData(i).alignedResponses = alignedResponses;
    psthData(i).meanResponse = psthMean;
    psthData(i).stdResponse = psthStd;
    psthData(i).semResponse = psthSEM;
    psthData(i).timeVector = timeVector;
    psthData(i).responseType = method;
end
end
