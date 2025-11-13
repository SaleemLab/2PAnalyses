function [responseCombined, sessionFileInfo] = combineResponseForVRRuns(VRruns, sessionFileInfo, overwrite)
% Combines VR experiment 'response' structures.
% Includes an 'overwrite' flag and robustly handles array reshaping.

%% Handle optional 'overwrite' input
if nargin < 3, overwrite = false; end

%% Define combined stimulus name and file path
combinedStimName = [sessionFileInfo.animal_name '_VRCorr_' sessionFileInfo.session_name '_CombinedRuns'];
combinedResponseFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_Response_' combinedStimName '.mat'];
combinedResponseFilePath = fullfile(sessionFileInfo.Directories.save_folder, combinedResponseFileName);
allStimNames = {sessionFileInfo.stimFiles.name};
existingStimIndex = find(strcmp(allStimNames, combinedStimName), 1);

%% Overwrite Check
if ~overwrite && ~isempty(existingStimIndex) && exist(combinedResponseFilePath, 'file')
    disp('Existing combined run found and overwrite is false. Loading existing file.');
    load(combinedResponseFilePath, 'response');
    responseCombined = response;
    return;
end
disp('Proceeding with combination (either no existing file found or overwrite is true).');

%% Handle sessionFileInfo Entry (Create or Update)
if ~isempty(existingStimIndex)
    sessionFileInfo.stimFiles(existingStimIndex).Response = combinedResponseFilePath;
    combinedStimIndex = existingStimIndex;
else
    disp('No existing combined stimulus found in sessionFileInfo. Creating a new entry.');
    newStimEntry = struct('name', combinedStimName, 'Response', combinedResponseFilePath, ...
        'bonsai_filepaths', [], 'tif_filepaths', [], 'BonsaiData', [], ...
        'TwoPMetaData', [], 'mergedBonsai2PSuite2pData', [], ...
        'processedMergedBonsaiSuite2pData', [], 'processedPeripheralData', []);
    sessionFileInfo.stimFiles(end+1) = newStimEntry;
    combinedStimIndex = length(sessionFileInfo.stimFiles);
end

%% Handle case of single run
numRuns = length(VRruns);
if numRuns < 2
    responseCombined = VRruns{1};
    disp('Only one VR run provided. No combination performed.');
    return;
end

%% Initialise Combination
responseCombined = VRruns{1};
numLapsInFirstRun = size(responseCombined.lapCount, 1);
responseCombined.VRLapRun = ones(numLapsInFirstRun, 1);

%% Main Combination Loop
for thisRun = 2:numRuns
    responseCurrent = VRruns{thisRun};
    lastEndTime = responseCombined.endTimeAll(end);
    totalLapsSoFar = size(responseCombined.lapCount, 1);

    % Adjust time, lap counters, and indices
    responseCurrent.startTimeAll = responseCurrent.startTimeAll + lastEndTime;
    responseCurrent.endTimeAll = responseCurrent.endTimeAll + lastEndTime;
    responseCurrent.completedStartTimes = responseCurrent.completedStartTimes + lastEndTime;
    responseCurrent.completedEndTimes = responseCurrent.completedEndTimes + lastEndTime;
    responseCurrent.lapCount = responseCurrent.lapCount + totalLapsSoFar;
    if ~isempty(responseCurrent.completedLaps), responseCurrent.completedLaps = responseCurrent.completedLaps + totalLapsSoFar; end
    if ~isempty(responseCurrent.abortedLaps), responseCurrent.abortedLaps = responseCurrent.abortedLaps + totalLapsSoFar; end

    % Combine 3D activity data
    responseCombined.lapPositionActivity.dFF = cat(2, responseCombined.lapPositionActivity.dFF, responseCurrent.lapPositionActivity.dFF);
    responseCombined.lapPositionActivity.dFFNeuropilCorrected = cat(2, responseCombined.lapPositionActivity.dFFNeuropilCorrected, responseCurrent.lapPositionActivity.dFFNeuropilCorrected);

    %% 
    % Enforce a consistent 2D shape (laps x bins) before concatenating.
    nBins = 140; % Assuming the number of bins is constant

    % For the already combined data:
    % Calculate number of laps from the total elements divided by number of bins.
    nLapsCombined = numel(responseCombined.lapPosition2PFrameIdx) / nBins;
    % Reshape it into the correct 2D format. This works even if it's already 2D.
    reshapedCombined = reshape(responseCombined.lapPosition2PFrameIdx, nLapsCombined, nBins);

    % For the new data to be added:
    nLapsCurrent = numel(responseCurrent.lapPosition2PFrameIdx) / nBins;
    reshapedCurrent = reshape(responseCurrent.lapPosition2PFrameIdx, nLapsCurrent, nBins);

    % Now that both arrays are guaranteed to be 2D with the same number of columns,
    % vertical concatenation is safe.
    responseCombined.lapPosition2PFrameIdx = [reshapedCombined; reshapedCurrent];
    % --- END OF FIX ---
    
    % Combine all other fields
    numLapsInCurrentRun = size(responseCurrent.lapCount, 1);
    currentRunLapTracker = ones(numLapsInCurrentRun, 1) * thisRun;
    responseCombined.wheelSpeed = [responseCombined.wheelSpeed; responseCurrent.wheelSpeed];
    responseCombined.mouseVirtualPosition = [responseCombined.mouseVirtualPosition; responseCurrent.mouseVirtualPosition];
    responseCombined.trackIDFromMousePosition = [responseCombined.trackIDFromMousePosition; responseCurrent.trackIDFromMousePosition];
    responseCombined.mouseRecordedPosition = [responseCombined.mouseRecordedPosition; responseCurrent.mouseRecordedPosition];
    responseCombined.trackIDs = [responseCombined.trackIDs; responseCurrent.trackIDs];
    responseCombined.lapCount = [responseCombined.lapCount; responseCurrent.lapCount];
    responseCombined.blockIDs = [responseCombined.blockIDs; responseCurrent.blockIDs];
    responseCombined.trialType = [responseCombined.trialType; responseCurrent.trialType];
    responseCombined.completedLaps = [responseCombined.completedLaps; responseCurrent.completedLaps];
    responseCombined.abortedLaps = [responseCombined.abortedLaps; responseCurrent.abortedLaps];
    responseCombined.endTimeAll = [responseCombined.endTimeAll; responseCurrent.endTimeAll];
    responseCombined.startTimeAll = [responseCombined.startTimeAll; responseCurrent.startTimeAll];
    responseCombined.completedStartTimes = [responseCombined.completedStartTimes; responseCurrent.completedStartTimes];
    responseCombined.completedEndTimes = [responseCombined.completedEndTimes; responseCurrent.completedEndTimes];
    responseCombined.stimName = [responseCombined.stimName; responseCurrent.stimName];
    responseCombined.lapPositionRelativeTime = [responseCombined.lapPositionRelativeTime; responseCurrent.lapPositionRelativeTime];
    responseCombined.VRLapRun = [responseCombined.VRLapRun; currentRunLapTracker];
end
disp(['Successfully combined ' num2str(numRuns) ' experiments.']);

%% Save final files
response = responseCombined;
save(combinedResponseFilePath, 'response', '-v7.3');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Saved combined response and updated sessionFileInfo.');

end