function [response, sessionFileInfo] = getLapPositionActivity_Flatten(sessionFileInfo, VRStimName, overwrite, onlyIncludeROIs)
% Calculates and saves binned lap activity for four signal types.
% Modified to only load necessary variables from disk to save memory/time.
%
% Aman and Sonali February 2025
% Modified Oct 2025

%% Handle optional inputs
if nargin < 3, overwrite = false; end 
if nargin < 4, onlyIncludeROIs = false; end

%% Find VR stimulus
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end

%% Load data (Selective Method)
disp('Loading Response and necessary processedTwoPData fields...');

% --- 1. Load Response ---
load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');

% --- 2. Ensure lapPosition2PFrameIdx is 2D ---
if ndims(response.lapPosition2PFrameIdx) == 3
     % If it was saved as 1 x Laps x Bins, squeeze it to Laps x Bins
     response.lapPosition2PFrameIdx = squeeze(response.lapPosition2PFrameIdx(1,:,:));
end

% --- 3. Load ONLY necessary parts of processedTwoPData ---
p2pFile = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;

% Try to load ALL potential necessary root variables at once.
% MATLAB will only load the ones that actually exist in the file.
loadedSubset = load(p2pFile, 'iscell', 'processedSignals', 'processedTwoPData');

if isfield(loadedSubset, 'processedTwoPData')
    % OLD nested format detected
    processedTwoPData = loadedSubset.processedTwoPData;
elseif isfield(loadedSubset, 'iscell') && isfield(loadedSubset, 'processedSignals')
    % NEW flattened format detected
    % Reconstruct a minimal struct so downstream code doesn't need changes
    processedTwoPData.iscell = loadedSubset.iscell;
    processedTwoPData.processedSignals = loadedSubset.processedSignals;
else
    error('File does not contain required variables (iscell/processedSignals OR processedTwoPData).');
end
clear loadedSubset;

%% Overwrite check
if overwrite && isfield(response, 'lapPositionActivity')
    disp('Overwrite is true. Removing old analysis fields...');
    fieldsToRemove = {'lapPositionActivity', 'cellROIs'}; 
    response = rmfield(response, intersect(fieldsToRemove, fieldnames(response)));
end

%% Get cell ROIs if needed
if onlyIncludeROIs
    ROIs = find(processedTwoPData.iscell(:, 1));
else
    ROIs = 1:size(processedTwoPData.iscell, 1);
end
numCells = length(ROIs);

%% Binning parameters and initialisation
numBins = 140; % 1cm bins for a 140cm track
nLaps = size(response.lapPosition2PFrameIdx, 1); 

lapPositionActivity = struct();

% Initialise storage for all signal types found
signalNames = fieldnames(processedTwoPData.processedSignals);
for thisSignal = 1:length(signalNames)
    lapPositionActivity.(signalNames{thisSignal}) = nan(numCells, nLaps, numBins);
end

%% --- EFFICIENT LOOP ---
disp('Binning lap-position-activity for all signals...');
for thisLap = 1:nLaps
    for thisBin = 1:numBins
        frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
        
        if ~isempty(frameIdx)
            for iSignal = 1:length(signalNames)
                currentSignalName = signalNames{iSignal};
                currentSignalMatrix = processedTwoPData.processedSignals.(currentSignalName);
                
                % Vectorized mean across time dimension
                lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = ...
                    mean(currentSignalMatrix(ROIs, frameIdx), 2, 'omitnan');
            end
        end
    end
end

%% Save results to the response struct
response.lapPositionActivity = lapPositionActivity;
if onlyIncludeROIs
    response.cellROIs = ROIs;
end

disp('Saving response with updated lapPositionActivities...');
save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');
end