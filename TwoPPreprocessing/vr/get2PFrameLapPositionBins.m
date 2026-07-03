%% 
function [response, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, VRStimName, turnBadFramesToNans)
%   For each lap and spatial bin (1 cm), this function finds the
%   corresponding two-photon (2P) frame indices. This index is the same
%   for all ROIs and is stored once per lap/bin.
%
%   Only includes frames where wheel speed > 1 cm/s and less than 100 cm/s
%   Bad frames are replaced with NaN to preserve timing structure.
%
% Inputs:
%   sessionFileInfo : struct
%       Contains paths to all processed data files for the session.
%   turnBadFramesToNans : logical (optional, default = false)
%       If true, filters out frames marked in processedTwoPData.badFrames.
%
% Output:
%   response : struct (updated)
%       Adds:
%         - lapPosition2PFrameIdx{lap, bin}   : 2P frame indices per lap per bin (NaN if invalid).
%         - lapPositionRelativeTime{lap, bin} : Time relative to lap start per lap per bin (NaN if invalid).
%
% Aman and Sonali - April 2025
% Optimized to save lapPosition2PFrameIdx as a 2D array - October 2025
% Modified May 2026 to exclude bad frames 
if nargin < 3, turnBadFramesToNans = true; end
%% Load processed data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx)
    error('Specified VRStimName ''%s'' not found in sessionFileInfo.', VRStimName);
end
disp('Loading data...');
% TwoPdata timeVec
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'badFrames', 'resample2PTimeUsed');
timeData = load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, resample2PTimeUsed);
timeVec = timeData.(resample2PTimeUsed);
% Ensure timeVec is a row vector
if size(timeVec, 1) > size(timeVec, 2), timeVec = timeVec'; end
% Load required response fields
response = load(sessionFileInfo.stimFiles(stimIdx).Response, 'completedStartTimes', 'completedEndTimes', ...
                         'wheelSpeed', 'mouseVirtualPosition', 'movementVisualGain', 'startTimeAll', 'endTimeAll');
load(sessionFileInfo.stimFiles(stimIdx).BonsaiData, 'bonsaiData');
% load(sessionFileInfo.stimFiles(stimIdx).processedPeripheralData, 'peripheralData');

% create a copy of the interpolated mouse position value 
bonsaiData.MousePos.Value_nanCorrected = bonsaiData.MousePos.Value;
%% Define filters and setup
minSpeedBin = 1;
maxSpeedBin = 100;
% Initial speed filter 
speedFilter = logical(response.wheelSpeed > minSpeedBin & response.wheelSpeed < maxSpeedBin);
if size(speedFilter, 1) > size(speedFilter, 2), speedFilter = speedFilter'; end

% fallback if bad frames mask are all false good by default
isBadFrame = false(1, length(timeVec));

% Bad frame filter logic
if turnBadFramesToNans && exist('badFrames', 'var') && ~isempty(badFrames)
    for p = 1:numel(badFrames)
        % Pull plane mask and force to logical
        thisPlaneMask = logical(badFrames{p});
        if size(thisPlaneMask, 1) > size(thisPlaneMask, 2)
            thisPlaneMask = thisPlaneMask';
        end
        % Combine using element-wise OR
        isBadFrame = isBadFrame | thisPlaneMask;
    end
    combinedFilter = speedFilter & ~isBadFrame;
    disp(['Excluding ' num2str(sum(isBadFrame)) ' bad frames from analysis via NaN masking.']);
    clear badFrames; % Free up memory
else
    combinedFilter = speedFilter;
end
if contains(VRStimName, 'Baseline') || contains(VRStimName, 'LandManipCorridor')
    % Track [1 to 200], ITI is 201.
    posBinEdges = 1:201; 
elseif contains(VRStimName, 'VRCorr')
    % Track [0 to 139], ITI is 140.
    posBinEdges = 0:140;
end 
numPosBins = length(posBinEdges) - 1;
nLaps = length(response.startTimeAll);
lapPosition2PFrameIdx = cell(nLaps, numPosBins);
lapPositionRelativeTime = cell(nLaps, numPosBins);

lapPosition_speedDerived = cell(nLaps,1);
%% Loop over laps to find frame indices for each spatial bin
flaggedLaps = []; % Initialize empty array to track bad laps
for thisLap = 1:nLaps
    disp(['Processing Lap: ' num2str(thisLap)]);
    lapStart = response.startTimeAll(thisLap);
    lapEnd = response.endTimeAll(thisLap);
    
    % Find indices for frames in this lap
    lapFrameIdx = find(timeVec >= lapStart & timeVec <= lapEnd);
    if isempty(lapFrameIdx)
        disp(['No 2P frames found for lap ' num2str(thisLap) '. Continuing...']);
        continue;
    end
    
    % Enforce lapFrameIdx to be a row vector to prevent dimension mismatch errors
    if size(lapFrameIdx, 1) > size(lapFrameIdx, 2), lapFrameIdx = lapFrameIdx'; end
    
    % Data loss check (Evaluates strictly based on bad frames)
    totalFramesInLap = length(lapFrameIdx);
    badFramesInLap = sum(isBadFrame(lapFrameIdx));
    if (badFramesInLap / totalFramesInLap) > 0.5
        flaggedLaps = [flaggedLaps, thisLap]; % Append the bad lap index to the array
        warning('Lap %d: More than 50%% of frames excluded due to BadFrames \n.', thisLap);
    end

    % Get position and force to row vector alignment
    lapPosition = response.mouseVirtualPosition(lapFrameIdx);
    if size(lapPosition, 1) > size(lapPosition, 2), lapPosition = lapPosition'; end

    % check position of the first frame; this can happen due to the
    % interpolation. For example if the first position is 159cm (due to
    % linear interpolation) and the second position is 1.x then change the
    % first position to match the second one.
    %     if length(lapPosition) > 1 && lapPosition(1) > lapPosition(2)
    %
    %         fprintf('  [PATCH] Lap %-2d: Shifting F1 Position from %6.2f cm -> %6.2f cm\n', ...
    %             thisLap, lapPosition(1), lapPosition(2));
    %         lapPosition(1) = lapPosition(2);
    %     end

   if length(lapPosition) > 1
       if thisLap <nLaps
           nextStart = response.startTimeAll(thisLap+1);
       else
           nextStart = response.endTimeAll(thisLap);
       end
       % Find indices for frames from start of this lap to the start of
       % next lap
       FrameIdx = find(timeVec >= lapStart & timeVec <= nextStart);
       lapPosition_speedDerived{thisLap} = [FrameIdx' cumsum(response.movementVisualGain/0.0612*response.wheelSpeed(FrameIdx)/60)];% Wheel speed also scaled to match VR speed
      
       rawFirstFive = nan(1,5);
       rawFirstFive(1:min(5, length(lapPosition))) = lapPosition(1:min(5, length(lapPosition)));
        lapPosition_original = lapPosition;
        lapPosition(find(abs(diff(rawFirstFive))>3)) =nan;

         bonsaiData.MousePos.Value_nanCorrected(lapFrameIdx) = lapPosition;

        fprintf(' Lap %-2d | Raw First 5: [%.2f, %.2f, %.2f, %.2f, %.2f] -> [%.2f, %.2f, %.2f, %.2f, %.2f]', ...
            thisLap, rawFirstFive(1), rawFirstFive(2), rawFirstFive(3), rawFirstFive(4), rawFirstFive(5),...
            lapPosition(1), lapPosition(2), lapPosition(3), lapPosition(4), lapPosition(5));
    end

    % Assign each frame in the lap to a position bin
    positionIdx = discretize(lapPosition, posBinEdges);
    if size(positionIdx, 1) > size(positionIdx, 2), positionIdx = positionIdx'; end

    for thisBin = 1:numPosBins
        % binMask captures all frames that  occurred in this position bin
        binMask = (positionIdx == thisBin);

        % Extract raw global frame indices as integers first
        rawFrameIdxInBin = lapFrameIdx(binMask);

        if ~isempty(rawFrameIdxInBin)
            % Identify frames that fail our filter criteria
            badFrameMask = ~combinedFilter(rawFrameIdxInBin);
            
            % Pull the times using the valid, raw integers 
            frameTimeInBin = timeVec(rawFrameIdxInBin); 
            relativeTime = frameTimeInBin - lapStart;
            
            % Convert indices to double so they can accept NaN
            frameIdxInBin = double(rawFrameIdxInBin);
            
            % Mask bad frames out with NaN for both outputs
            frameIdxInBin(badFrameMask) = NaN;
            relativeTime(badFrameMask) = NaN; 
            
            % Store the masked arrays
            lapPosition2PFrameIdx{thisLap, thisBin} = frameIdxInBin;
            lapPositionRelativeTime{thisLap, thisBin} = relativeTime;
        else
            % If the mouse skipped or never occupied this spatial bin during the lap
            lapPosition2PFrameIdx{thisLap, thisBin} = [];
            lapPositionRelativeTime{thisLap, thisBin} = [];
        end
    end
end
fprintf('%d out of %d laps flagged due to over 50%% loss of frames [bad-frames only] \n', length(flaggedLaps) ,nLaps)
%% Save the updated 2D cell arrays to the response struct
response.lapPosition2PFrameIdx = lapPosition2PFrameIdx;
response.lapPositionRelativeTime = lapPositionRelativeTime;
response.flaggedLaps= flaggedLaps; 
response.lapPosition_speedDerived = lapPosition_speedDerived; 

% bonsaiData.MousePos.Value_nanCorrected = Value_nanCorrected; 
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append'); 
save(sessionFileInfo.stimFiles(stimIdx).BonsaiData, 'bonsaiData', '-append'); % added a new corrected mousepos.value variable that turns uncertains positions to nans (due to interp)
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Done.');
end