function [response, sessionFileInfo] = computeShuffleMatrixForSession(sessionFileInfo, response, vrStimNamesToStitch, useZScoredProcessedSignals)
%
% Computes and adds the shuffle matrix to the input 'response' structure.
% Handles standalone 'spks' variables with temporal smoothing.
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

%% Setup dimensions and filters
[numROIs, ~, numBins] = size(response.lapPositionActivity.dFF); 
ROIs = 1:numROIs;
signalNames = {}; 
signalVarName = iif(useZScoredProcessedSignals, 'zScoredProcessedSignals', 'processedSignals');

% Setup temporal smoothing window matching getLapPositionActivity
w = gausswin(15); 
w = w / sum(w);

%% Stitch RAW Signals Only
disp('Stitching RAW signals from processedTwoPData across all specified VR runs...');
stitchedRawSignalMatrices = struct();
totalCombinedFrames = 0;

for iStim = 1:length(vrStimNamesToStitch)
    stimName = vrStimNamesToStitch{iStim};
    stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
    
    if isempty(stimIdx), error(['Raw signal file for ', stimName, ' not found.']); end
    
    filePath = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;
    varsInFile = who('-file', filePath);
    
    % Pull the main signals struct
    dataChunk = load(filePath, signalVarName);
    signals = dataChunk.(signalVarName);
    
    if isempty(signalNames)
        signalNames = fieldnames(signals);
    end
    
    % Check for standalone spikes exactly like the previous function
    hasSpks = false;
    if ismember('spks', varsInFile)
        spksData = load(filePath, 'spks');
        if ~isempty(spksData.spks) && any(spksData.spks(:))
            hasSpks = true;
            if ~ismember('spks', signalNames)
                signalNames{end+1} = 'spks'; % Append to the processing track
            end
        end
    end
    
    % Stitch the signal matrices (Cell x Time)
    for iSignal = 1:length(signalNames)
        currentSignalName = signalNames{iSignal};
        
        if strcmp(currentSignalName, 'spks')
            if hasSpks
                % Filter/Smooth the spikes matrix prior to stitching
                currentChunk = filtfilt(w, 1, spksData.spks')';
                currentChunk = currentChunk(ROIs, :);
            else
                % If this specific run chunk is missing spikes, fill with NaNs to keep matrix aligned
                arbitraryFields = fieldnames(signals);
                firstField = arbitraryFields{1};
                currentChunk = nan(length(ROIs), size(signals.(firstField), 2));
            end
        else
            % Pull and smooth standard signals (dFF, etc.)
            rawChunk = signals.(currentSignalName)(ROIs, :);
            currentChunk = filtfilt(w, 1, rawChunk')';
        end
        
        if iStim == 1
            stitchedRawSignalMatrices.(currentSignalName) = currentChunk;
        else
            stitchedRawSignalMatrices.(currentSignalName) = [stitchedRawSignalMatrices.(currentSignalName), currentChunk];
        end
    end
    
    % Fixed: Extract the character array from the cell array first before querying size
    arbitraryFields = fieldnames(signals);
    firstField = arbitraryFields{1};
    totalCombinedFrames = totalCombinedFrames + size(signals.(firstField), 2);
end

clear dataChunk signals currentChunk rawChunk spksData; 

%% Create Global Frame-to-Bin Map
lapPosition2PFrameIdx = response.lapPosition2PFrameIdx;
if ndims(lapPosition2PFrameIdx) == 3, lapPosition2PFrameIdx = squeeze(lapPosition2PFrameIdx); end 
nLaps = size(lapPosition2PFrameIdx, 1);
frameToBinMap_Global = nan(1, totalCombinedFrames);

disp('Creating global frame-to-bin map...');
for thisLap = 1:nLaps
    for thisBin = 1:numBins
        frameIdx = lapPosition2PFrameIdx{thisLap, thisBin};
        cleanFrameIdx = frameIdx(~isnan(frameIdx));
        
        if ~isempty(cleanFrameIdx)
            cleanFrameIdx = round(cleanFrameIdx); 
            frameToBinMap_Global(cleanFrameIdx) = thisBin;
        end
    end
end

%% Calculate Shuffle Matrix 
totalTimes = totalCombinedFrames; 
maxShift = round(totalTimes./2); 
numShifts = 1000;
rng(1) 
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

%% Return Results
response.lapPositionActivity_ShuffleMatrix = lapPositionActivity_ShuffleMatrix; 
stimIdx = find(strcmp(response.stimName, {sessionFileInfo.stimFiles.name}));
disp(['Saving updated Response struct to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end

% Helper function (for conditional assignment)
function out = iif(condition, true_value, false_value)
    if condition
        out = true_value;
    else
        out = false_value;
    end
end