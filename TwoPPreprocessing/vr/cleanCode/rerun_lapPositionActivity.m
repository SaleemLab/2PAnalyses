% Error logging setup
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'errorLog_rerun0722.csv');
logHeaders = {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'};

% If the log file doesn't exist, create it with headers
if ~exist(logFilePath, 'file')
    logTable = table('Size', [0, numel(logHeaders)], ...
        'VariableTypes', repmat("string", 1, numel(logHeaders)), ...
        'VariableNames', logHeaders);
    writetable(logTable, logFilePath);
end
fprintf('Error logging enabled. Log will be saved to: %s\n', logFilePath);

vrKeywords = {'VRCorr', 'BaselineCorridor', 'LandManipCorridor'};
rfKeywords = {'RFMapping'};
doNotCombine = {'M25040_VRCorr_20250507_00001', 'M25040_VRCorr_20250507_00002', 'M25057_VRCorr_20250526_00001', 'M25057_VRCorr_20250526_00002', 'M25126_VRCorr_20260123_00001', 'M25126_VRCorrBaseline_20260123_00002', 'M25126_VRCorrWithManipulations_20260123_00003', 'M25132_BaselineCorridor_20260219_00001', 'M25132_BaselineCorridor_20260219_00002', ...
    'M26003_BaselineCorridor_20260322_00001',  'M26003_BaselineCorridor_20260322_00002','M26003_BaselineCorridor_20260324_00001', 'M25131_BaselineCorridor_20260421_00001', 'M26005_BaselineCorridor_20260421_00001'};

%  % ran these sessions: on May 28th 
% % pairs = struct; 
% pairs.M26005 = ['20260305', '20260306', '20260311']; % ,
% pairs.M26004 = ['20260305', '20260307', '20260312', '20260313', '20260314']; % '
% % % pairs.M25131 = ['20260312', '20260313', '20260314'] done 
% % %pairs.M25126 = ['20260311', '20260312', '20260313'];
% pairs.M25132 = ['20260221', '20260223', '20260226', '20260228', '20260313']; % done '20260219','20260220',
% pairs.M25133 = ['20260219','20260220','20260221','20260223','20260224'];
% pairs.M26003 = ['20260316','20260317', '20260320','20260321','20260322', '20260324', '20260325'];

pairs = struct; 
pairs.M26003 = ['20260324'];

useZScoredProcessedSignals=true;
% filteredTable = filterMasterTable_usingNameSessionPairs('MouseID', {'M26003'}, 'DayOfExperience',300,'Exclude', 0);
filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs,'Exclude', 0);
mouseInfo = sessionsToProcess(filteredTable);

totalSessionsToProcess = 0;
sessionsProcessedCount = 0;
for i = 1:size(mouseInfo, 1)
    totalSessionsToProcess = totalSessionsToProcess + length(mouseInfo{i, 2});
end

for thisMouse = 1:size(mouseInfo, 1)
    mousenumber = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};
    
    fprintf('Processing Mouse: %s\n', mousenumber);
    
    % Loop through each session for the current mouse
    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        sessionsProcessedCount = sessionsProcessedCount + 1;
        fprintf('  Processing Session: %s\n', sessionName);
        
        try
            % Add your session processing code here
            infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
            if ~isfile(infoPath)
                 error('Info Missing'); % Changed from return to error so it logs and continues
            end
            
            loadedInfo = load(infoPath, 'sessionFileInfo');
            sessionFileInfo = loadedInfo.sessionFileInfo;
            stimNames = {sessionFileInfo.stimFiles.name};
            
            % select combined runs if present 
            targetIdx = find(contains(stimNames, "Corridor") & ~contains(stimNames, "CombinedRuns"));
            VRStimNames = stimNames(targetIdx);
            
            if isempty(targetIdx)
                error('No valid LandManip or Baseline corridor found.');
            end
        
            for thisTargetIdx = 1:length(targetIdx)
                idx = targetIdx(thisTargetIdx);
                VRStimName = sessionFileInfo.stimFiles(targetIdx).name;
                % including 2p frame lap position activity
                [~, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, VRStimName);
                [~, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, VRStimName, useZScoredProcessedSignals,false); %this still contains the flaggedLaps
                [~, sessionFileInfo] = getRunningSpeedAcrossLaps(sessionFileInfo, VRStimName);
            end
            
            shouldProcessSeparately = any(ismember(VRStimNames, doNotCombine));
            runAsIndividualLoop = (isscalar(VRStimNames)) || (length(VRStimNames) > 1 && shouldProcessSeparately);
            
            if runAsIndividualLoop
                fprintf('Processing via INDIVIDUAL LOOP path.\n');
                for thisVRStim = 1:length(VRStimNames)
                    vrStimName = VRStimNames{thisVRStim};
                    % stimIdx = find(strcmp(vrStimName, {sessionFileInfo.stimFiles.name}), 1);
                            % load the individual response structure into the environment
                    % response = load(sessionFileInfo.stimFiles(stimIdx).Response, 'lapPositionActivity', 'trialIndicesByCondition' , 'flaggedLaps');
                    % excluding all other fields temporary 
                    crossValExpVar = getCrossValidatedExplainedVariance(sessionFileInfo, vrStimName);
                    [~, ~] = computeSpatialModulationIdex(sessionFileInfo, vrStimName,true,true); % exclude laps; apply smoothning
                    % [~,~,~] = checkOddEvenCorrelation(sessionFileInfo, response,signalName, applyTemporalSmoothing,true, false); % Run on final response
                    % [~,~,~] = checkHalvesCorrelation(sessionFileInfo, response, signalName, applyTemporalSmoothing, true, false);
                end
            else
                fprintf('Processing via COMBINED PATH\n');
                % Combine Response sturctures
                fprintf('Checking for and combining Response across multiple VR runs\n')
                responseVRRuns = loadDataStructuresByKeyword(sessionFileInfo, 'Corridor', 'Response');
                [response, sessionFileInfo] = combineResponseForVRRuns(responseVRRuns, sessionFileInfo);
                crossValExpVar = getCrossValidatedExplainedVariance(sessionFileInfo, response.stimName);
                [~, ~] = computeSpatialModulationIdex(sessionFileInfo,response.stimName,true,true); 
            end
            
        catch ME
            % --- ERROR LOGGING ---
            fprintf('    --> ERROR encountered. Logging to CSV...\n');
            
            % Get current time
            timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
            
            % Get error details
            errMsg = string(ME.message);
            if isempty(ME.stack)
                errFunc = "Unknown";
                errLine = "Unknown";
            else
                errFunc = string(ME.stack(1).name);
                errLine = string(ME.stack(1).line);
            end
            
            % Create table row for the error
            errorData = table(timestamp, string(mousenumber), string(sessionName), ...
                              errMsg, errFunc, errLine, ...
                              'VariableNames', logHeaders);
                          
            % Append to the CSV
            writetable(errorData, logFilePath, 'WriteMode', 'append');
            
            % Print to console and continue to next session
            fprintf('    --> Logged error: %s\n', errMsg);
            continue;
        end
    end 
end