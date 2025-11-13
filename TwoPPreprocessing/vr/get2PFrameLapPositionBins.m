function [response, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, VRStimName)
%   For each lap and spatial bin (1 cm), this function finds the
%   corresponding two-photon (2P) frame indices. This index is the same
%   for all ROIs and is stored once per lap/bin.
%
%   Only includes frames where wheel speed > 1 cm/s.
%
% Inputs:
%   sessionFileInfo : struct
%       Contains paths to all processed data files for the session.
%
% Output:
%   response : struct (updated)
%       Adds:
%         - lapPosition2PFrameIdx{lap, bin}   : 2P frame indices per lap per bin.
%         - lapPositionRelativeTime{lap, bin} : Time relative to lap start per lap per bin.
%
% Aman and Sonali - April 2025
% Optimized to save lapPosition2PFrameIdx as a 2D array - October 2025

%% Load processed data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx)
    error('Specified VRStimName ''%s'' not found in sessionFileInfo.', VRStimName);
end

disp('Loading data...');
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');
load(sessionFileInfo.stimFiles(stimIdx).BonsaiData, 'bonsaiData');
load(sessionFileInfo.stimFiles(stimIdx).processedPeripheralData, 'peripheralData');
load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');

%% Define bins and setup
posBinEdges = 0:140;
numPosBins = length(posBinEdges) - 1;
nLaps = length(response.completedStartTimes);
timeVec = processedTwoPData.(processedTwoPData.resample2PTimeUsed);
speedFilter = response.wheelSpeed > 1 & response.wheelSpeed < 100;

% Initialise as a 2D cell array 
% No need for the ROI dimension.
lapPosition2PFrameIdx = cell(nLaps, numPosBins);
lapPositionRelativeTime = cell(nLaps, numPosBins);

%% Loop over laps to find frame indices for each spatial bin
for thisLap = 1:nLaps
    disp(['Processing Lap: ' num2str(thisLap)]);
    
    lapStart = response.completedStartTimes(thisLap);
    lapEnd = response.completedEndTimes(thisLap);
    
    % Find all 2P frames that occurred within this lap's time window
    lapFrameIdx = find(timeVec >= lapStart & timeVec <= lapEnd);
    
    if isempty(lapFrameIdx)
        disp(['No 2P frames found for lap ' num2str(thisLap) '. Continuing...']);
        continue;
    end
    
    lapPosition = response.mouseVirtualPosition(lapFrameIdx);
    
    % Assign each frame in the lap to a position bin
    positionIdx = discretize(lapPosition, posBinEdges);
    
    for thisBin = 1:numPosBins
        % Create a mask for frames that are in the current bin AND pass the speed filter
        binMask = (positionIdx == thisBin) & speedFilter(lapFrameIdx);
        frameIdxInBin = lapFrameIdx(binMask);
        
        % --- CHANGE 2: Store the indices directly into the 2D array ---
        % The unnecessary loop over ROIs is completely removed.
        lapPosition2PFrameIdx{thisLap, thisBin} = frameIdxInBin;

        if ~isempty(frameIdxInBin)
            frameTimeInBin = timeVec(frameIdxInBin);
            lapPositionRelativeTime{thisLap, thisBin} = frameTimeInBin - lapStart;
        else
            lapPositionRelativeTime{thisLap, thisBin} = [];
        end
    end
end

%% Save the updated 2D cell arrays to the response struct
response.lapPosition2PFrameIdx = lapPosition2PFrameIdx;
response.lapPositionRelativeTime = lapPositionRelativeTime;

disp('Saving updated response struct...');
save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Done.');
end