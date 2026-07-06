% rerun the rf mapping and direction tuning sessions

%% INITIALISE ERROR LOG
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'TuningRerun_20260705_RSP_dotfields.csv');
logHeaders = {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'};
if ~exist(logFilePath, 'file')
    logTable = table('Size', [0, numel(logHeaders)], 'VariableTypes', repmat("string", 1, numel(logHeaders)), 'VariableNames', logHeaders);
    writetable(logTable, logFilePath);
end
fprintf('Error logging enabled. Log will be saved to: %s\n', logFilePath);

%% rerun the rf mapping and direction tuning sessions
filteredTable = filterMasterTable_usingNameSessionPairs('MouseID', {'M26003', 'M25133', 'M25132'},'HasStimulus',{'DotMotion_SpeedTuning'},'Exclude', 0);

% filteredTable = filterMasterTable_usingNameSessionPairs('MouseID', {'M26004', 'M26005', 'M25131', 'M25126'},'HasStimulus',{'DotMotion_SpeedTuning','DirTuning', 'RFMapping'},'Exclude', 0);
mouseInfo = sessionsToProcess(filteredTable);
signalName = 'dFFNeuropilCorrected';
totalSessionsToProcess = 0;
sessionsProcessedCount = 0;

for i = 1:size(mouseInfo, 1)
    totalSessionsToProcess = totalSessionsToProcess + length(mouseInfo{i, 2});
end

for thisMouse = 1:size(mouseInfo, 1)
    mousenumber = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};
    fprintf('Processing Mouse: %s\n', mousenumber);

    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        sessionsProcessedCount = sessionsProcessedCount + 1;
        fprintf('  Processing Session: %s\n', sessionName);

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if ~isfile(infoPath)
            warning('  sfi missing for %s — skipping.', sessionName);
            logError(logFilePath, mousenumber, sessionName, 'sessionFileInfo file missing', 'findSessionFileInfoFilePath', NaN);
            continue;
        end
        loadedInfo = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames = {sessionFileInfo.stimFiles.name};
        
        tuningStimIdx = find(contains(stimNames, {'DotMotion_SpeedTuning'}, 'IgnoreCase', true));
        % tuningStimIdx = find(contains(stimNames, {'RFMapping', 'DotMotion_SpeedTuning', 'DirTuning'}, 'IgnoreCase', true));
        if isempty(tuningStimIdx)
            disp('  RFMapping, DotMotion_SpeedTuning, or DirTuning not found for this session.');
            continue;
        end

        for i = 1:length(tuningStimIdx)
            thisName = stimNames{tuningStimIdx(i)};
            fprintf('Processing stimulus: %s\n', thisName);

            [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, thisName, 1,false, false, 5,5,20); %zscore; do not smooth or do red channel correction 
            [~, sessionFileInfo] = getStimTimes(sessionFileInfo, thisName, pdThresholdForStimEvents);

            try
                if contains(thisName, 'RFMapping')
                    fprintf('  Running RFMapping Analysis...\n');
                    [~, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, thisName, true);
                    [~, sessionFileInfo] = getRFMappingTrialGroups(sessionFileInfo, thisName);
                    [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, thisName, signalName);
                    [sessionFileInfo, RFMapping, RFMappingMetadata, allCenters] = analyseRFMapping(sessionFileInfo, thisName);
                    [~] = getRunningAndPupilDataByTrials(sessionFileInfo, thisName);

                elseif contains(thisName, 'DotMotion_SpeedTuning')
                    fprintf('  Running DotMotion Analysis...\n');
                    [~, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, thisName, true);
                    [bonsaiData, sessionFileInfo] = getTrialGroups(sessionFileInfo, thisName);
                    [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, thisName, signalName);
                    [~] = getRunningAndPupilDataByTrials(sessionFileInfo, thisName);

                elseif contains(thisName, 'DirTuning')
                    fprintf('  Running Direction Tuning Analysis...\n');
                    [~, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, thisName, true);
                    [~, sessionFileInfo] = getTrialGroups(sessionFileInfo, thisName);
                    [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, thisName, signalName);
                    [~] = getRunningAndPupilDataByTrials(sessionFileInfo, thisName);
                    disp('  Direction tuning processed.');
                end

                fprintf('  Completed: %s\n', thisName);

            catch ME
                warning('  Error processing %s / %s / %s: %s', mousenumber, sessionName, thisName, ME.message);
                if ~isempty(ME.stack)
                    logError(logFilePath, mousenumber, sessionName, ME.message, ME.stack(1).name, ME.stack(1).line);
                else
                    logError(logFilePath, mousenumber, sessionName, ME.message, thisName, NaN);
                end
                continue; % move on to the next stimulus/session rather than stopping the whole run
            end
        end
    end
end

fprintf('\nDone. %d/%d sessions processed. See %s for any errors.\n', ...
    sessionsProcessedCount, totalSessionsToProcess, logFilePath);


%% local error logging helper
function logError(logFilePath, mouseID, sessionName, errMsg, funcName, lineNum)
    newRow = table(string(datetime('now')), string(mouseID), string(sessionName), ...
        string(errMsg), string(funcName), string(lineNum), ...
        'VariableNames', {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'});
    writetable(newRow, logFilePath, 'WriteMode', 'Append');
end