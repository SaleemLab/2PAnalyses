function [response, sessionFileInfo] = computeShuffleMatrixForSession(sessionFileInfo, response, vrStimNamesToStitch, useZScoredProcessedSignals)

%
% Computes and adds the shuffle matrix to the input 'response' structure.
%
% vrStimNamesToStitch: Cell array of VR stimulus names whose RAW signals (processedTwoPData) 
%                      will be stitched together to form the continuous time series for shuffling.

%% Handle optional inputs and basic checks

if nargin < 4, useZScoredProcessedSignals = true; end

if isempty(vrStimNamesToStitch)
    error('vrStimNamesToStitch list cannot be empty for shuffle computation.');
end

% Check required fields
if  ~isfield(response, 'lapPosition2PFrameIdx')
    error('Input response structure is missing required fields (lapPosition2PFrameIdx).');
end

%% 
[numROIs, ~, numBins] = size(response.lapPositionActivity.dFF); 
ROIs = 1:numROIs;
signalNames = {}; % Will be set during the first load

% if overwrite && isfield(response, 'lapPositionActivity_ShuffleMatrix')
%     disp('Overwrite is true. Removing old shuffle analysis field...');
%     response = rmfield(response, 'lapPositionActivity_ShuffleMatrix');
% end

signalVarName = iif(useZScoredProcessedSignals, 'zScoredProcessedSignals', 'processedSignals');

%% Stitch RAW Signals Only

disp('Stitching RAW signals from processedTwoPData across all specified VR runs...');
stitchedRawSignalMatrices = struct();
totalCombinedFrames = 0;

for iStim = 1:length(vrStimNamesToStitch)
    stimName = vrStimNamesToStitch{iStim};
    stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
    
    if isempty(stimIdx), error(['Raw signal file for ', stimName, ' not found.']); end
    
    % Pull only the struct indicated by signalVarName
    % This ignores massive F and Fneu matrices to save RAM.
    filePath = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;
    dataChunk = load(filePath, signalVarName);
    signals = dataChunk.(signalVarName);
    
    if isempty(signalNames)
        signalNames = fieldnames(signals);
    end
    
    % Stitch the signal matrices (Cell x Time)
    for iSignal = 1:length(signalNames)
        currentSignalName = signalNames{iSignal};
        
        % Select ROIs from the current signal chunk
        currentChunk = signals.(currentSignalName)(ROIs, :);
        
        if iStim == 1
            stitchedRawSignalMatrices.(currentSignalName) = currentChunk;
        else
            stitchedRawSignalMatrices.(currentSignalName) = [stitchedRawSignalMatrices.(currentSignalName), currentChunk];
        end
    end
    totalCombinedFrames = totalCombinedFrames + size(currentChunk, 2);
end
clear dataChunk signals currentChunk; 

%% Create Global Frame-to-Bin Map (Uses input response.lapPosition2PFrameIdx; this has already been combined across runs where appropriate)

% Map indices directly, as they are assumed to be global time stamps
lapPosition2PFrameIdx = response.lapPosition2PFrameIdx;
if ndims(lapPosition2PFrameIdx) == 3, lapPosition2PFrameIdx = squeeze(lapPosition2PFrameIdx); end % little bug fix if present..

nLaps = size(lapPosition2PFrameIdx, 1);

frameToBinMap_Global = nan(1, totalCombinedFrames);
% frameToLapMap has been removed as I am not computing the lap-to-lap
% activity 
disp('Creating global frame-to-bin map...');
for thisLap = 1:nLaps
    for thisBin = 1:numBins
        frameIdx = lapPosition2PFrameIdx{thisLap, thisBin};
        if ~isempty(frameIdx)
            % Use indices directly; they should already point into the stitched signal length
            frameToBinMap_Global(frameIdx) = thisBin;
        end
    end
end
% Frame map is now complete and reflects the stitching order 

%% Calculate Shuffle Matrix 
totalTimes = totalCombinedFrames; % changed from the length of F to total frames combined across VRRuns 
maxShift = round(totalTimes./2); % total time divided by 2 
numShifts = 1000;
rng(1) % will only work once 
randShifts = randi(maxShift,[1 numShifts]);

lapPositionActivity_ShuffleMatrix = struct();



disp(['Calculating Shuffle Matrix across ', num2str(totalTimes), ' frames...']);
tl = tic;

for iSignal = 1:length(signalNames)
    currentSignalName = signalNames{iSignal};
    currentSignalMatrix_Combined = stitchedRawSignalMatrices.(currentSignalName);

    lapPositionActivity_ShuffleMatrix.(currentSignalName) = nan(numROIs, numBins, numShifts);

    disp(['Processing Shuffle Matrix for: ', currentSignalName]);
    for thisShift = 1:numShifts
        if mod(thisShift, 100) == 0 || thisShift == 1
            fprintf('  -> Processing Shuffle %d of %d... (Shift amount: %d)\n', thisShift, numShifts, randShifts(thisShift));
        end
        thisShiftAmount = randShifts(thisShift);

        % Circularly shift the whole combined signal
        shiftedSignalMatrix = circshift(currentSignalMatrix_Combined, [0 thisShiftAmount]);

        medianLapActivityForShift = nan(numROIs, numBins);

        % Loop over position bins
        for thisBin = 1:numBins
            % Find all frames belonging to THIS position bin (across ALL combined laps)
            allFramesInBin = find(frameToBinMap_Global == thisBin);

            if ~isempty(allFramesInBin)
                activityInBin = shiftedSignalMatrix(:, allFramesInBin);
                medianActivityInBin = median(activityInBin, 2, 'omitnan');
                medianLapActivityForShift(:, thisBin) = medianActivityInBin;
            end
        end

        lapPositionActivity_ShuffleMatrix.(currentSignalName)(:, :, thisShift) = medianLapActivityForShift;
    end
end
toc(tl)

%% 5. Return Results
response.lapPositionActivity_ShuffleMatrix = lapPositionActivity_ShuffleMatrix; 
stimIdx = find(strcmp(response.stimName, {sessionFileInfo.stimFiles.name}));

disp(['Saving updated Response struct to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end

% Helper function (for conditional assignment) @gemini and stupid.. 
function out = iif(condition, true_value, false_value)
    if condition
        out = true_value;
    else
        out = false_value;
    end
end