function [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, signalField, applySmoothing, VRStimName, onlyIncludeROIs)
%   Extracts mean binned fluorescence/activity values per lap from 2P data.
%   Only includes ROIs labeled as "cells" and only for completed laps.
%   Optionally applies Gaussian temporal smoothing across position bins.
%
% Inputs:
%   sessionFileInfo : struct
%
%   signalField : string (optional, default = 'F')
%       The signal to extract (e.g., 'F', 'Fneu', 'spks', 'deltaF_over_F')
%
%   applySmoothing : logical (optional, default = false)
%       Whether to smooth activity across position bins using a Gaussian filter
%
% Output:
%   response : struct (updated)
%       Adds the following fields:
%         - lapPositionActivity : [nCells x nLaps x nBins] mean activity per bin
%         - cellROIs             : indices of cell ROIs
%         - signalUsed           : which signal type was used
%         - smoothingApplied     : whether smoothing was applied
%
% Example usage:
%   response = getLapPositionActivity(sessionFileInfo, 'F', false);

% Handle optional inputs
if nargin < 2 || isempty(signalField)
    signalField = 'F';
end
if nargin < 3
    applySmoothing = false;
end

if nargin < 4
    onlyIncludeROIs = true;
end 

%% Find the VRCorr stimulus index
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName, sessionFileInfo.stimFiles(iStim).name);
end
iStim = find(bonsaiData.isVRstim == 1);

if isempty(iStim)
    error('No VRCorr stimulus found in sessionFileInfo.');
end

%% Load required data
if exist(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'file') && ...
   exist(sessionFileInfo.stimFiles(iStim).Response, 'file')
    disp('Loading processedTwoPData..')
    load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'processedTwoPData')
    disp('Loaded response..')
    load(sessionFileInfo.stimFiles(iStim).Response, 'response')
else
    error('Missing processed files: Response or processedTwoPData.');
end

%% Validate and pull signal
if ~isfield(processedTwoPData, signalField)
    error(['Signal field "' signalField '" not found in twoPData.']);
end
signalMatrix = processedTwoPData.(signalField);  % Size: ROI x time

if onlyIncludeROIs
    % Get ROIs that are cells
    isCell = logical(processedTwoPData.iscell(:, 1));
    cellROIs = find(isCell);
    numCells = length(cellROIs);
else
    numCells = size(processedTwoPData.F, 1);
end 
% Lap + bin setup
binCentres = 0.5:1:139.5;
numBins = length(binCentres);
nLaps = length(response.completedStartTimes);

% Init output matrix
lapPositionActivity = nan(numCells, nLaps, numBins);

% Compute mean signal per bin per lap
for thisCell = 1:numCells
    if onlyIncludeROIs
        roiIdx = cellROIs(thisCell);
    else
        roiIdx = thisCell;
    end
    for lapIdx = 1:nLaps
        for binIdx = 1:numBins
            frameIdx = response.lapPosition2PFrameIdx{roiIdx, lapIdx, binIdx};
            if ~isempty(frameIdx)
                lapPositionActivity(thisCell, lapIdx, binIdx) = ...
                    mean(signalMatrix(roiIdx, frameIdx));
            end
        end
    end
end

%% Optional Gaussian smoothing across bins
if applySmoothing
    w = gausswin(9); w = w / sum(w);
    for thisCell = 1:numCells
        for lapIdx = 1:nLaps
            signal = squeeze(lapPositionActivity(thisCell, lapIdx, :));
            if all(isnan(signal))
                continue 
            end
            nanMask = isnan(signal);F
            signal(nanMask) = 0; % Turn nans to 0 temporarily 
            smoothed = filtfilt(w, 1, signal); % Smoothen
            smoothed(nanMask) = NaN; % Turn back to nans 
            lapPositionActivity(thisCell, lapIdx, :) = smoothed;
        end
    end
end

%% Store in response struct and save
response.lapPositionActivity = lapPositionActivity;
if onlyIncludeROIs
    response.cellROIs = cellROIs;
end
response.signalUsed = signalField;
response.smoothingApplied = applySmoothing;

% Save updated response
disp('Saving updated response with lapPositionActivity...');
save(sessionFileInfo.stimFiles(iStim).Response, 'response', '-v7.3');

end
