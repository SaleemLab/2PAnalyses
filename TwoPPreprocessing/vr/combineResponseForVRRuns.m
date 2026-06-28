function [responseCombined, sessionFileInfo] = combineResponseForVRRuns(VRruns, sessionFileInfo, overwrite)
% combineResponseForVRRuns: Concatenates multiple VR sessions into one.
% Updated to ensure trial indices align with matrix rows for unified plotting.
if nargin < 3, overwrite = true; end

% Generate naming for the combined file
parts = split(VRruns{2}.stimName, '_');
corridorType = 'VRRun';
if length(parts) >= 2, corridorType = parts{2}; end
combinedStimName = [sessionFileInfo.animal_name '_' corridorType '_' sessionFileInfo.session_name '_CombinedRuns'];
combinedResponseFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_Response_' combinedStimName '.mat'];
combinedResponseFilePath = fullfile(sessionFileInfo.Directories.save_folder, combinedResponseFileName);

if isempty(sessionFileInfo.stimFiles)
    allStimNames = {};
else
    allStimNames = {sessionFileInfo.stimFiles.name};
end
existingStimIndex = find(strcmp(allStimNames, combinedStimName), 1);

if ~overwrite && ~isempty(existingStimIndex) && exist(combinedResponseFilePath, 'file')
    disp('Existing combined run found. Loading flattened file...');
    responseCombined = load(combinedResponseFilePath);
    return;
end

% Add or update fields in sessionFileInfo 
if ~isempty(existingStimIndex)
    sessionFileInfo.stimFiles(existingStimIndex).Response = combinedResponseFilePath;
else
    if ~isempty(sessionFileInfo.stimFiles)
        newStimEntry = sessionFileInfo.stimFiles(1);
        fnames = fieldnames(newStimEntry);
        for i = 1:length(fnames)
            newStimEntry.(fnames{i}) = [];
        end
    else
        newStimEntry = struct('name', [], 'Response', [], 'bonsai_filepaths', [], ...
            'tif_filepaths', [], 'BonsaiData', [], 'TwoPMetaData', [], ...
            'mergedBonsai2PSuite2pData', [], 'processedMergedBonsaiSuite2pData', [], ...
            'processedPeripheralData', []);
    end
    newStimEntry.name = combinedStimName;
    newStimEntry.Response = combinedResponseFilePath;
    sessionFileInfo.stimFiles(end+1) = newStimEntry;
end

numRuns = length(VRruns);
if numRuns < 2
    responseCombined = VRruns{1};
    disp('Only one VR run provided.');
else
    % Initialise with the first run
    responseCombined = VRruns{1};
    responseCombined.stimName = combinedStimName; 
    responseCombined.AllstimName = VRruns{1}.stimName;
    
    % Ensure flaggedLaps field exists in the initial struct (even if empty)
    if ~isfield(responseCombined, 'flaggedLaps')
        responseCombined.flaggedLaps = [];
    end
    
    % Force flaggedLaps from Run 1 to be a row vector for consistent appending
    if size(responseCombined.flaggedLaps, 1) > size(responseCombined.flaggedLaps, 2)
        responseCombined.flaggedLaps = responseCombined.flaggedLaps';
    end
    
    % Track how many completed laps are stored in the matrix columns
    sigNames = fieldnames(responseCombined.lapPositionActivity);
    numCompletedSoFar = size(responseCombined.lapPositionActivity.(sigNames{1}), 2);
    
    for thisRun = 2:numRuns
        responseCurrent = VRruns{thisRun};
        lastEndTime = responseCombined.endTimeAll(end);
        
        % Adjust timestamps
        timeFields = {'startTimeAll', 'endTimeAll', 'completedStartTimes', 'completedEndTimes'};
        for tf = 1:length(timeFields)
            fieldName = timeFields{tf};
            responseCurrent.(fieldName) = responseCurrent.(fieldName) + lastEndTime;
        end
        
        % Shift trial indices
        if isfield(responseCurrent, 'trialIndicesByCondition')
            condNames = fieldnames(responseCurrent.trialIndicesByCondition);
            for c = 1:length(condNames)
                cn = condNames{c};
                currIdx = responseCurrent.trialIndicesByCondition.(cn);
                
                if ~isempty(currIdx)
                    shiftedIdx = currIdx + numCompletedSoFar;
                    if isfield(responseCombined.trialIndicesByCondition, cn)
                        responseCombined.trialIndicesByCondition.(cn) = ...
                            [responseCombined.trialIndicesByCondition.(cn); shiftedIdx];
                    else
                        responseCombined.trialIndicesByCondition.(cn) = shiftedIdx;
                    end
                end
            end
        end
        
        %% Process and Shift newly added flaggedLaps variable
        if isfield(responseCurrent, 'flaggedLaps') && ~isempty(responseCurrent.flaggedLaps)
            currFlaggedLaps = responseCurrent.flaggedLaps;
            if size(currFlaggedLaps, 1) > size(currFlaggedLaps, 2), currFlaggedLaps = currFlaggedLaps'; end
            shiftedFlaggedLaps = currFlaggedLaps + numCompletedSoFar;
            responseCombined.flaggedLaps = [responseCombined.flaggedLaps, shiftedFlaggedLaps];
        end
        
        % Concatenate lapPositionActivity (Matrix Columns = Laps)
        for s = 1:length(sigNames)
            sn = sigNames{s};
            responseCombined.lapPositionActivity.(sn) = cat(2, ...
                responseCombined.lapPositionActivity.(sn), ...
                responseCurrent.lapPositionActivity.(sn));
        end
        
        % Concatenate LapRunningSpeed related data 
        responseCombined.lapPositionRunningSpeed = cat(1, ...
            responseCombined.lapPositionRunningSpeed, ...
            responseCurrent.lapPositionRunningSpeed); 
        
        %% --- INTEGRATION OF lapPosition_speedDerived (CELL ARRAY) ---
        if isfield(responseCurrent, 'lapPosition_speedDerived')
            responseCombined.lapPosition_speedDerived = [responseCombined.lapPosition_speedDerived; ...
                                                         responseCurrent.lapPosition_speedDerived];
        end
        %% -------------------------------------------------------------
        
        % Concatenate lap position data 
        responseCombined.lapPosition2PFrameIdx = [responseCombined.lapPosition2PFrameIdx; ...
                                                  responseCurrent.lapPosition2PFrameIdx];
        responseCombined.lapPositionRelativeTime = [responseCombined.lapPositionRelativeTime; ...
                                                    responseCurrent.lapPositionRelativeTime];
        
        landmarkFields = {'completedLandmarkNames', 'completedLandmarkPositions', ...
                          'completedLandmarkSizes', 'completedNumLandmarks', 'completedLandmarkCenterOffsets', 'completedLandmarkReward'};
        for thisLandField = 1:length(landmarkFields)
            fieldName = landmarkFields{thisLandField};
            if isfield(responseCurrent, fieldName)
                responseCombined.(fieldName) = [responseCombined.(fieldName); responseCurrent.(fieldName)];
            end
        end
        
        % concatenation of other variables 
        responseCombined.wheelSpeed = [responseCombined.wheelSpeed; responseCurrent.wheelSpeed];
        responseCombined.mouseVirtualPosition = [responseCombined.mouseVirtualPosition; responseCurrent.mouseVirtualPosition];
        responseCombined.lapCountAll = [responseCombined.lapCountAll; responseCurrent.lapCountAll];
        responseCombined.startTimeAll = [responseCombined.startTimeAll; responseCurrent.startTimeAll];
        responseCombined.endTimeAll = [responseCombined.endTimeAll; responseCurrent.endTimeAll];
        responseCombined.completedStartTimes = [responseCombined.completedStartTimes; responseCurrent.completedStartTimes];
        responseCombined.completedEndTimes = [responseCombined.completedEndTimes; responseCurrent.completedEndTimes];
        responseCombined.trialTypeAll = [responseCombined.trialTypeAll; responseCurrent.trialTypeAll];
        
        if ~isfield(responseCombined, 'AllstimName')
            responseCombined.AllstimName = {};
        end
        
        if isfield(responseCombined, 'movementVisualGain')
             responseCombined.movementVisualGain = [responseCombined.movementVisualGain; responseCurrent.movementVisualGain]; 
        end 
        responseCombined.AllstimName = [responseCombined.AllstimName; {responseCurrent.stimName}];
       
        % Mean running speed 
        responseCombined.lapRunningSpeed = [responseCombined.lapRunningSpeed; responseCurrent.lapRunningSpeed]; 
        
        % Update the offset for the next run in the loop
        numCompletedSoFar = size(responseCombined.lapPositionActivity.(sigNames{1}), 2);
    end
end

if numRuns >= 2
    responseCombined.completedLaps_AbsoluteIdx = [responseCombined.completedLaps_AbsoluteIdx; responseCurrent.completedLaps_AbsoluteIdx];
    responseCombined.abortedLaps_AbsoluteIdx = [responseCombined.abortedLaps_AbsoluteIdx; responseCurrent.abortedLaps_AbsoluteIdx];
end

% Save results
disp(['Saving combined response to ', combinedResponseFilePath]);
save(combinedResponseFilePath, '-struct', 'responseCombined', '-v7.3');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Done.');
end