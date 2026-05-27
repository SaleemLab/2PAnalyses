% pairs=struct; % RSP 
% pairs.M25132 = ['20260226', '20260228', '20260313'];
% pairs.M25133 = '20260224';
% pairs.M26003 = ['20260322', '20260324', '20260325'];

pairs=struct; 
 
pairs.M26005 = ['20260305', '20260306', '20260311', '20260318', '20260321', '20260322'];
pairs.M26004 = ['20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322'];
pairs.M25131 = ['20260312', '20260313', '20260314', '20260318', '20260321', '20260322'];
pairs.M25126 = ['20260311', '20260312', '20260313'];
pairs.M25132 = ['20260228', '20260313'];
pairs.M26003 = ['20260324', '20260325'];
pairs.M25133 = '20260224';

filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

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

        VRStimName = sessionFileInfo.stimFiles(targetIdx).name;

        

        [~, ~] = computeSpatialModulationIdex(sessionFileInfo, VRStimName,true,false);
        crossValExpVar = getCrossValidatedExplainedVariance(sessionFileInfo, VRStimName);

    end 
end


% newStim.name                             = 'M26003_20260324_Response_M26003_LandManipCorridor_20260324_CombinedRuns';
% newStim.bonsai_filepaths                 = [];
% newStim.eyetracking_filepaths            = [];
% newStim.tif_filepaths                    = [];
% newStim.TwoPMetaData                     = [];
% newStim.processedPeripheralData          = [];
% newStim.mergedBonsai2PSuite2pData        = [];
% newStim.processedMergedBonsaiSuite2pData = [];
% newStim.BonsaiData                       = [];
% newStim.Response                         = 'Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_Response_M26003_LandManipCorridor_20260324_CombinedRuns.mat'; 
% 
% % 2. Append the empty entry directly to the end of your struct array
% sessionFileInfo.stimFiles(end + 1) = newStim;
% 
% % Display the updated size
% disp(sessionFileInfo.stimFiles);