
% day of experience 300 should be the open-closed loop experiments 

% RSP_OL_CL = filterMasterTable_usingNameSessionPairs('MouseID', {'M25132', 'M26003'}, 'DayOfExperience',300,'Exclude', 0);
%stimName = if mouse 'M25132' stim name ''M25132_BaselineCorridor_20260423_00001' is closed loop and 'M25132_BaselineCorridor_OpenLoop_20260423_00001' is open loop

%if mouseID is 'M26003' then use
%'M26003_BaselineCorridor_20260421_CombinedRuns' for open-loop and
%'M26003_BaselineCorridor_20260421_00003' for closed loop 



RSPData = appendFilteredROIs(RSPData,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);


mouseInfo = sessionsToProcess(RSP_OL_CL);

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

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if ~isfile(infoPath)
             msg = 'Info Missing'; return; 
        end
        loadedInfo = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;

        sessionROIData = sessionFileInfo.otherSessFilePaths




        stimNames = {sessionFileInfo.stimFiles.name};

        % create a filepath for closed loop 
        closedLoopName = sprintf(['%s_BaselineCorridor_%s_00001', ])




    end

end 