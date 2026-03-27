% function [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, VRStimName, useZScoredProcessedSignals, onlyIncludeROIs)
% 
% % Calculates binned lap activity (Neurons x Laps x Bins) for a single stimulus.
% % Returns the unshifted signal and its time indices for stitching/shuffling later.
% % Using modified version Dec 2025. 
% 
% %% Handle optional inputs
% if nargin < 3, useZScoredProcessedSignals = true; end
% if nargin < 4, onlyIncludeROIs = false; end
% 
% 
% %% Find VR stimulus and load data
% stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
% if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
% 
% filePath = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;
% varsInFile = who('-file', filePath); 
% 
% if useZScoredProcessedSignals
%     load(filePath, 'zScoredProcessedSignals', 'iscell');
%     signals = zScoredProcessedSignals;
% else
%     % Load processedSignals and iscell 
%     load(filePath, 'processedSignals', 'iscell');
%     signals = processedSignals;
% 
%     % Only load spks if it actually exists in the file
%     if ismember('spks', varsInFile)
%         load(filePath, 'spks');
%         hasSpks = true;
%     else
%         hasSpks = false;
%         warning('Spikes (spks) not found for this session. Skipping spike processing.');
%     end
% end
% 
% response = load(sessionFileInfo.stimFiles(stimIdx).Response, 'lapPosition2PFrameIdx');
% 
% %% Apply temporal smoothning to spikes (similar to mika elife paper)
% 
% spks_smoothed = []; 
% if ~useZScoredProcessedSignals && hasSpks
%     disp('Applying temporal smoothing to spikes (gausswin 15)...');
%     w = gausswin(15); 
%     w = w / sum(w);
%     spks_smoothed = filtfilt(w, 1, spks')'; 
% end
% 
% %% Get cell ROIs
% if onlyIncludeROIs
%     ROIs = find(iscell(:, 1));
% else
%     ROIs = (1:size(iscell, 1))';
% end
% % response.numROIs = ROIs;
% response.lapPositionActivityZScored = useZScoredProcessedSignals;
% 
% %% Calculate Real Activity (Neurons x Laps x Bins)
% % signalNames = fieldnames(signals);
% signalNames = {'dFF', 'dFFNeuropilCorrected', 'spks'}; 
% numSignals = length(signalNames);
% if contains(VRStimName, 'Baseline') || contains(VRStimName, 'LandManipCorridor')
%     numBins = 200;
% elseif contains(VRStimName, 'VRCorr')
%     numBins = 140;
% end
% nLaps = size(response.lapPosition2PFrameIdx, 1);
% lapPositionActivity = struct(); 
% 
% disp('Calculating Real Activity (Neurons x Laps x Bins)...');
% for iSignal = 1:numSignals
%     currentSignalName = signalNames{iSignal};
% 
%     if strcmp(currentSignalName, 'spks') && (~hasSpks || useZScoredProcessedSignals)
%         continue; 
%     end
% 
%     if strcmp(currentSignalName, 'spks')
%         currentSignalMatrix = spks_smoothed(ROIs,:);
%     else 
%         currentSignalMatrix = signals.(currentSignalName)(ROIs, :); 
%     end
% 
%     numROIs = size(currentSignalMatrix, 1);
%     % Initialize storage
%     lapPositionActivity.(currentSignalName) = nan(numROIs, nLaps, numBins);
%     % frameIdx = [];
%     % for thisLap = 1:nLaps
%     %     for thisBin = 1:numBins
%     %         frameIdx(thisLap,thisBin) = length(response.lapPosition2PFrameIdx{thisLap, thisBin});
%     %     end
%     % end
%     for thisLap = 1:nLaps
%         for thisBin = 1:numBins
%             frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
%             if ~isempty(frameIdx)
%                 % Compute mean activity for this bin
%                 lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = ...
%                     mean(currentSignalMatrix(:, frameIdx), 2, 'omitnan');
%             end
%         end
%     end
% end
% 
% %% Now spks 
% 
% response.lapPositionActivity = lapPositionActivity;
% 
% disp(['Saving updated Response struct (Lap Activity) to ', sessionFileInfo.stimFiles(stimIdx).Response]);
% save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
% save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
% % save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-append');
% end

function [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, VRStimName, useZScoredProcessedSignals, onlyIncludeROIs)
% Calculates binned lap activity (Neurons x Laps x Bins) for a single stimulus.
% Returns the unshifted signal and its time indices for stitching/shuffling later.
% Handles cases where 'spks' may be missing, empty, or all zeros.

%% Handle optional inputs
if nargin < 3, useZScoredProcessedSignals = true; end
if nargin < 4, onlyIncludeROIs = false; end

%% Find VR stimulus and load data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end

disp(['Loading data for stimulus: ', VRStimName]);


hasSpks = false; 
filePath = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;

% Check what variables exist in the .mat file before loading
varsInFile = who('-file', filePath);

% 
if useZScoredProcessedSignals 
    load(filePath, 'zScoredProcessedSignals', 'iscell');
    signals = zScoredProcessedSignals;
else
    load(filePath, 'processedSignals', 'iscell');
    signals = processedSignals;
end


% Check if 'spks' exists, is not empty, and contains actual data (not all
% zeros); this needs to be fixed in the 2p data 
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

% Load response indexing
response = load(sessionFileInfo.stimFiles(stimIdx).Response, 'lapPosition2PFrameIdx');


%% Apply temporal smoothing to spikes
spks_smoothed = []; 
if hasSpks
    disp('Applying temporal smoothing to spikes (gausswin 15)...');
    w = gausswin(15); 
    w = w / sum(w);
    % Vectorised smoothing along time dimension
    spks_smoothed = filtfilt(w, 1, spks')'; 
end

%% Get cell ROIs
if onlyIncludeROIs
    ROIs = find(iscell(:, 1));
else
    ROIs = (1:size(iscell, 1))';
end

response.lapPositionActivityZScored = useZScoredProcessedSignals;

%% Calculate Real Activity (Neurons x Laps x Bins)
signalNames = {'dFF', 'dFFNeuropilCorrected', 'spks'}; 
numSignals = length(signalNames);

% Define binning based on Stimulus Name
if contains(VRStimName, 'Baseline') || contains(VRStimName, 'LandManipCorridor')
    numBins = 200;
elseif contains(VRStimName, 'VRCorr')
    numBins = 140;
else
    numBins = 100; % Default fallback
end

nLaps = size(response.lapPosition2PFrameIdx, 1);
lapPositionActivity = struct(); 

disp('Calculating Real Activity (Neurons x Laps x Bins)...');

for iSignal = 1:numSignals
    currentSignalName = signalNames{iSignal};
    
    % SKIP logic: If this is spks but we don't have valid data, skip to next signal
    if strcmp(currentSignalName, 'spks') && ~hasSpks
        continue; 
    end

    % Select the correct matrix to process
    if strcmp(currentSignalName, 'spks')
        currentSignalMatrix = spks_smoothed(ROIs, :);
    else 
        % Check if the field exists in the signals struct (e.g. dFF)
        if isfield(signals, currentSignalName)
            currentSignalMatrix = signals.(currentSignalName)(ROIs, :); 
        else
            warning('Field %s not found in signals struct. Skipping.', currentSignalName);
            continue;
        end
    end 
    
    numROIs = size(currentSignalMatrix, 1);
    

    lapPositionActivity.(currentSignalName) = nan(numROIs, nLaps, numBins);
    
    for thisLap = 1:nLaps
        for thisBin = 1:numBins
            frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
            if ~isempty(frameIdx)
                % Compute mean activity for this bin
                lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = ...
                    mean(currentSignalMatrix(:, frameIdx), 2, 'omitnan');
            end
        end
    end
end

%% save
response.lapPositionActivity = lapPositionActivity;

disp(['Saving updated Response struct (Lap Activity) to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

end