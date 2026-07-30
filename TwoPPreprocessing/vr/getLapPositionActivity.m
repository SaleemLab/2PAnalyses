function [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, VRStimName, useZScoredProcessedSignals, onlyIncludeROIs)
% Calculates binned lap activity (Neurons x Laps x Bins) for a single stimulus.
% Returns the unshifted signal and its time indices for stitching/shuffling later.
% Handles cases where 'spks' may be missing, empty, or all zeros.
%
% Applies temporal smoothing (gausswin 15) to ALL signals prior to spatial
% binning. Default is 250ms smoothning 

if nargin < 3, useZScoredProcessedSignals = true; end
if nargin < 4, onlyIncludeROIs = false; end

stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
disp(['Loading data for stimulus: ', VRStimName]);

hasSpks = false; 
filePath = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;
varsInFile = who('-file', filePath);

if useZScoredProcessedSignals 
    load(filePath, 'zScoredProcessedSignals', 'iscell');
    signals = zScoredProcessedSignals;
else
    load(filePath, 'processedSignals', 'iscell');
    signals = processedSignals;
end

if ismember('spks', varsInFile)
    load(filePath, 'spks');
    if ~isempty(spks) && any(spks(:))
        hasSpks = true;
    else
        warning('Variable "spks" found but is empty or all zeros for %s. Skipping.', VRStimName);
    end
else
    warning('Variable "spks" not found in file: %s. Skipping.', filePath);
end

% Ensure flaggedLaps field exists to prevent crashes down stream
response = load(sessionFileInfo.stimFiles(stimIdx).Response);
if ~isfield(response, 'flaggedLaps')
    response.flaggedLaps = [];
end

w = gausswin(15); 
w = w / sum(w);
sigFields = fieldnames(signals);
for f = 1:length(sigFields)
    fieldName = sigFields{f};
    if isnumeric(signals.(fieldName)) && ~isempty(signals.(fieldName))
        signals.(fieldName) = filtfilt(w, 1, signals.(fieldName)')';
    end
end

spks_smoothed = []; 
if hasSpks
    spks_smoothed = filtfilt(w, 1, spks')'; 
end

if onlyIncludeROIs
    ROIs = find(iscell(:, 1));
else
    ROIs = (1:size(iscell, 1))';
end

response.lapPositionActivityZScored = useZScoredProcessedSignals;

signalNames = {'dFF', 'dFFNeuropilCorrected', 'spks'}; 
numSignals = length(signalNames);

if contains(VRStimName, 'Baseline') || contains(VRStimName, 'LandManipCorridor')
    numBins = 200; % was 200 when the previous version was shifted forward.. 
elseif contains(VRStimName, 'VRCorr')
    numBins = 140;
end

%% --- RESTRICT TO COMPLETED LAPS ONLY ---
if isfield(response, 'completedLaps_AbsoluteIdx') && ~isempty(response.completedLaps_AbsoluteIdx)
    completedLapIndices = response.completedLaps_AbsoluteIdx;
else
    warning('completedLaps_AbsoluteIdx field not found or empty. Defaulting to all rows in lapPosition2PFrameIdx.');
    completedLapIndices = (1:size(response.lapPosition2PFrameIdx, 1))';
end

% Convert logical masks to index list if necessary
if islogical(completedLapIndices)
    completedLapIndices = find(completedLapIndices);
end

numCompletedLaps = length(completedLapIndices);
lapPositionActivity = struct(); 

for iSignal = 1:numSignals
    currentSignalName = signalNames{iSignal};
    
    if strcmp(currentSignalName, 'spks') && ~hasSpks
        continue; 
    end
    
    if strcmp(currentSignalName, 'spks')
        currentSignalMatrix = spks_smoothed(ROIs, :);
    else 
        if isfield(signals, currentSignalName)
            currentSignalMatrix = signals.(currentSignalName)(ROIs, :); 
        else
            warning('Field %s not found in signals struct. Skipping.', currentSignalName);
            continue;
        end
    end 
    
    numROIs = size(currentSignalMatrix, 1);
    
    % Preallocate matrix using completed lap count
    lapPositionActivity.(currentSignalName) = nan(numROIs, numCompletedLaps, numBins);
    
    for iLap = 1:numCompletedLaps
        % Map the sequential index loop to the target lap row
        actualLapRow = completedLapIndices(iLap);
        
        for thisBin = 1:numBins
            frameIdx = response.lapPosition2PFrameIdx{actualLapRow, thisBin};
            
            validFrameMask = ~isnan(frameIdx);
            cleanFrameIdx = frameIdx(validFrameMask);
            
            if ~isempty(cleanFrameIdx)
                meanActivity = mean(currentSignalMatrix(:, cleanFrameIdx), 2, 'omitnan');
                lapPositionActivity.(currentSignalName)(:, iLap, thisBin) = meanActivity;
            end
        end
    end
end

response.lapPositionActivity = lapPositionActivity;

disp(['Saving updated Response struct (Lap Activity) to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Done.');
end