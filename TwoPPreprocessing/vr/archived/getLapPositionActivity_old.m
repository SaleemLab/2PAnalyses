function [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, VRStimName, overwrite, onlyIncludeROIs, useZScoredProcessedSignals)
% Calculates and saves binned lap activity for four signal types.
% This version uses an efficient, vectorized approach assuming that
% lapPosition2PFrameIdx is a 2D cell array (nLaps x nBins).
%
% Aman and Sonali February 2025
% Modified Oct 2025

%% Handle optional inputs
if nargin < 3, overwrite = false; end % Default overwrite to false
if nargin < 4, onlyIncludeROIs = false; end
if nargin < 5, useZScoredProcessedSignals = true; end


%% Find VR stimulus and load data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end

disp('Loading processedTwoPData and Response structs...');
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');
load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');

%% Overwrite check
if overwrite && isfield(response, 'lapPositionActivity')
    disp('Overwrite is true. Removing old analysis fields...');
    fieldsToRemove = {'lapPositionActivity', 'cellROIs'}; % Add any other related fields
    response = rmfield(response, intersect(fieldsToRemove, fieldnames(response)));
end

%% Get cell ROIs if needed
if onlyIncludeROIs
    ROIs = find(processedTwoPData.iscell(:, 1));
else
    ROIs = 1:size(processedTwoPData.iscell, 1);
end
numCells = length(ROIs);

%% Select appropriate signal to use

if useZScoredProcessedSignals
    disp('Using zScored dFF and dFFNeuropilCorrected for spatial tuning curves..')
    signals = processedTwoPData.zScoredProcessedSignals; 
    response.signalsZScored = true; 
else  
    signals = processedTwoPData.processedSignals; 
    response.signalsZScored = false; 
    disp('Using dFF and dFFNeuropilCorrected (without zScoring) for spatial tuning curves..')
end 
%% Binning parameters and initialisation
if ndims(response.lapPosition2PFrameIdx) == 3
    response.lapPosition2PFrameIdx = squeeze(response.lapPosition2PFrameIdx);
end
numBins = 140; % 1cm bins for a 140cm track
nLaps = size(response.lapPosition2PFrameIdx, 1); % Get nLaps from the 2D cell array
lapPositionActivity = struct();


% Initialise storage for all signal types
signalNames = fieldnames(signals);
for thisSignal = 1:length(signalNames)
    lapPositionActivity.(signalNames{thisSignal}) = nan(numCells, nLaps, numBins);
end

%% 
disp('Binning lap-position-activity for all signals...');

% Loop over laps and bins first
tl = tic;
for thisLap = 1:nLaps
    for thisBin = 1:numBins
        % 1. Get frame indices ONCE for this lap and bin.
        % lapPosition2PFrameIdx is now {nLaps, nBins}
        frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
        
        % If frames exist in this bin, process all signals and cells
        if ~isempty(frameIdx)
            for iSignal = 1:length(signalNames)
                currentSignalName = signalNames{iSignal};
                currentSignalMatrix = signals.(currentSignalName);
                
                % Vectorised calculation for ALL cells at once.
                % This is much faster than an inner for-loop.
                % It takes the mean across the time dimension (dim 2).
                meanActivity = mean(currentSignalMatrix(ROIs, frameIdx), 2, 'omitnan');
                
                % Store the resulting vector of activities for all cells.
                lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = meanActivity;
            end
        end
    end
end
toc(tl)

maxShift = 2000;
numShifts = 20;
randShifts = randi(maxShift,[1 numShifts]);
for iShift = 1:length(randShifts)
    newSignal = circshift(currentSignalMatrix,[0 randShifts(iShift)]);
    meanActivity = mean(newSignal(ROIs, frameIdx), 2, 'omitnan');
    lapPositionActivity_randShift.(currentSignalName)(:, thisLap, thisBin, iShift) = meanActivity;
end
%% Save results to the response struct
response.lapPositionActivity = lapPositionActivity;
if onlyIncludeROIs
    response.cellROIs = ROIs;
end

disp('Saving response with updated lapPositionActivities...');
save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');
end