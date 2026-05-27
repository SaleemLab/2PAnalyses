pairs.M26005 = ['20260305', '20260306', '20260311', '20260318', '20260321', '20260322'];
pairs.M26004 = ['20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322'];
pairs.M25131 = ['20260312', '20260313', '20260314', '20260318', '20260321', '20260322'];
pairs.M25126 = ['20260311', '20260312', '20260313'];
pairs.M25132 = ['20260228', '20260313'];
pairs.M26003 = ['20260324', '20260325'];
pairs.M25133 = '20260224';

mouseInfo = sessionsToProcess(filteredTable);


signalToUse='dFFNeuropilCorrected';
nFolds =5; 

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

        % Process the current session
        sessionsProcessedCount = sessionsProcessedCount + 1;
        fprintf('  Processing Session: %s\n', sessionName);
        % Add your session processing code here

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if ~isfile(infoPath)
             msg = 'Info Missing'; return; 
        end
        loadedInfo = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;

        stimNames = string({sessionFileInfo.stimFiles.name});
        % select combined runs if present 
        targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "CombinedRuns"), 1);
        
        if isempty(targetIdx)
            targetIdx = find(contains(stimNames, "Corridor"), 1);
        end
        
        if length(targetIdx)>1
            targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "00002"), 1);
        end 

        if isempty(targetIdx)
            error('No valid LandManip or Baseline corridor found.');
        end
        fprintf('Loading response from %s \n', sessionFileInfo.stimFiles(targetIdx).name)
        response = load(sessionFileInfo.stimFiles(targetIdx).Response);
        
        if ~isfield(response.lapPositionActivity, signalToUse)
            error('Signal %s missing from data.', signalToUse);
        end
        

        % [~] = getCrossValidatedExplainedVariance(sessionFileInfo, response, signalToUse);
        % rerun the correaltion checks using baseline conditions only 
        [~,~,~] = checkOddEvenCorrelation(sessionFileInfo, response, signalToUse, true, false, false);

        [~,~,~] = checkHalvesCorrelation(sessionFileInfo, response, signalToUse, true, false, false);



    end 
end