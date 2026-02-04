function [responseCombined, sessionFileInfo] = combineResponseForVRRuns(VRruns, sessionFileInfo, overwrite)
% Combines VR experiment 'response' structures into a flattened .mat file.

%% Handle optional 'overwrite' input
if nargin < 3, overwrite = true; end

%% Define combined stimulus name and file path
combinedStimName = [sessionFileInfo.animal_name '_VRCorr_' sessionFileInfo.session_name '_CombinedRuns'];
combinedResponseFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_Response_' combinedStimName '.mat'];
combinedResponseFilePath = fullfile(sessionFileInfo.Directories.save_folder, combinedResponseFileName);

% stimFiles is a struct array and get existing names
if isempty(sessionFileInfo.stimFiles)
    allStimNames = {};
else
    allStimNames = {sessionFileInfo.stimFiles.name};
end
existingStimIndex = find(strcmp(allStimNames, combinedStimName), 1);

%% Overwrite 
if ~overwrite && ~isempty(existingStimIndex) && exist(combinedResponseFilePath, 'file')
    disp('Existing combined run found and overwrite is false. Loading flattened file...');
    responseCombined = load(combinedResponseFilePath);
    return;
end

%% Handle sessionFileInfo Entry (Create or Update)
if ~isempty(existingStimIndex)
    % Update existing entry
    sessionFileInfo.stimFiles(existingStimIndex).Response = combinedResponseFilePath;
else
    % Create new entry using the existing struct as a template to avoid field mismatches
    if ~isempty(sessionFileInfo.stimFiles)
        % Copy the first element to get all field names exactly right
        newStimEntry = sessionFileInfo.stimFiles(1);
        % Reset all fields to empty/default values
        fnames = fieldnames(newStimEntry);
        for i = 1:length(fnames)
            newStimEntry.(fnames{i}) = [];
        end
    else
        % Fallback: manually define if stimFiles is totally empty
        newStimEntry = struct('name', [], 'Response', [], 'bonsai_filepaths', [], ...
            'tif_filepaths', [], 'BonsaiData', [], 'TwoPMetaData', [], ...
            'mergedBonsai2PSuite2pData', [], 'processedMergedBonsaiSuite2pData', [], ...
            'processedPeripheralData', []);
    end
    
    % Assign specific values
    newStimEntry.name = combinedStimName;
    newStimEntry.Response = combinedResponseFilePath;
    
    % Append to array
    sessionFileInfo.stimFiles(end+1) = newStimEntry;
end

%% Handle case of single run
numRuns = length(VRruns);
if numRuns < 2
    responseCombined = VRruns{1};
    disp('Only one VR run provided. No combination performed.');
else
    %% Initialize Combination
    responseCombined = VRruns{1};
    numLapsInFirstRun = size(responseCombined.lapCount, 1);
    responseCombined.VRLapRun = ones(numLapsInFirstRun, 1);
    responseCombined.stimName = combinedStimName; 
    responseCombined.AllstimName = VRruns{1}.stimName;
    
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
        sigNames = fieldnames(responseCombined.lapPositionActivity);
        for s = 1:length(sigNames)
            sn = sigNames{s};
            responseCombined.lapPositionActivity.(sn) = cat(2, responseCombined.lapPositionActivity.(sn), responseCurrent.lapPositionActivity.(sn));
        end
        
        % Concatenate lapPosition2PFrameIdx (flattened 2D handling)
        nBins = 140; 
        nLapsCombined = numel(responseCombined.lapPosition2PFrameIdx) / nBins;
        reshapedCombined = reshape(responseCombined.lapPosition2PFrameIdx, nLapsCombined, nBins);
        nLapsCurrent = numel(responseCurrent.lapPosition2PFrameIdx) / nBins;
        reshapedCurrent = reshape(responseCurrent.lapPosition2PFrameIdx, nLapsCurrent, nBins);
        responseCombined.lapPosition2PFrameIdx = [reshapedCombined; reshapedCurrent];
        
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
        responseCombined.AllstimName = [responseCombined.AllstimName, ' | ', responseCurrent.stimName];
        responseCombined.lapPositionRelativeTime = [responseCombined.lapPositionRelativeTime; responseCurrent.lapPositionRelativeTime];
        responseCombined.VRLapRun = [responseCombined.VRLapRun; currentRunLapTracker];
    end
end

%% Save final flattened file
disp(['Saving combined response (Flattened) to ', combinedResponseFilePath]);
save(combinedResponseFilePath, '-struct', 'responseCombined', '-v7.3');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Done.');
end