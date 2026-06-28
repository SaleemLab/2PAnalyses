clear; clc;
% load sessions with gray screen and darkness 
% rsp-post learning and visp mice 
pairs = struct; 

% % Use cell arrays {} to prevent MATLAB from flattening distinct dates into a single long string
% pairs.M26005 = {'20260305', '20260306', '20260311', '20260318', '20260321', '20260322'};
% pairs.M26004 = {'20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322'};
% pairs.M25131 = {'20260312', '20260313', '20260314', '20260318', '20260321', '20260322'};
% pairs.M25126 = {'20260311', '20260312', '20260313'};
% pairs.M25132 = {'20260226', '20260228', '20260313'};
% pairs.M26003 = {'20260322', '20260324', '20260325'};
% pairs.M25133 = {'20260224'};

% missing unique rsp fovs
% pairs.M25132 = {'20260219', '20260220', '20260221', '20260223', '20260228', '20260303', '20260306', '20260312'};
% pairs.M26003 = {'20260316', '20260317', '20260320', '20260321'};
% pairs.M25133 = {'20260219', '20260220', '20260221', '20260223'};
pairs.M25132 = {'20260306', '20260219'};

paramsBoutons.signalName = 'dFFNeuropilCorrected';
paramsBoutons.useZScoredProcessedSignals = false; 

paramsSomas.signalName = 'spks'; 
paramsSomas.useZScoredProcessedSignals = false; 

filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'Darkness', 'GrayScreen'});
mouseInfo = sessionsToProcess(filteredTable);

logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'DarkGary_20260617_RSPmissingsessions.csv');
logHeaders = {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'};

% If the log file doesn't exist, create it with headers
if ~exist(logFilePath, 'file')
    logTable = table('Size', [0, numel(logHeaders)], 'VariableTypes', repmat("string", 1, numel(logHeaders)), 'VariableNames', logHeaders);
    writetable(logTable, logFilePath);
end
fprintf('Error logging enabled. Log will be saved to: %s\n', logFilePath);

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
        
        fprintf('  Processing Session [%d/%d]: %s  %s \n', sessionsProcessedCount, totalSessionsToProcess, mousenumber, sessionName);
        
        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if isempty(infoPath) || ~isfile(infoPath)
             fprintf('    WARNING: File path missing or invalid for %s - %s. Skipping.\n', mousenumber, sessionName);
             continue; % Using continue instead of return ensures your loop keeps iterating over remaining sessions
        end
        
  
        matchRow = find(strcmp(filteredTable.MouseID, mousenumber) & strcmp(filteredTable.Session, sessionName), 1);
        if isempty(matchRow)
            fprintf('    WARNING: Could not resolve matching row inside metadata table for %s - %s. Skipping.\n', mousenumber, sessionName);
            continue;
        end
        
        typeImaged = filteredTable.TypeImaged{matchRow};
         fprintf('   Selecting perams for imaged type %s .\n', typeImaged);
        
        % Condition parameter routing assignment block
        if strcmpi(typeImaged, 'Boutons')
            signalName        = paramsBoutons.signalName;
            useZScoredSignals = paramsBoutons.useZScoredProcessedSignals;
        else
            signalName        = paramsSomas.signalName;
            useZScoredSignals = paramsSomas.useZScoredProcessedSignals;
        end
        
        try
            loadedInfo = load(infoPath, 'sessionFileInfo');
            sessionFileInfo = loadedInfo.sessionFileInfo;
            stimNames = {sessionFileInfo.stimFiles.name}; 
            
            % Select combined runs if present 
            targetIdx = find(contains(stimNames, {'Darkness', 'GrayScreen'}));
            
            if any(targetIdx)
                fprintf('    Found %d target stimulus track file(s).\n', length(targetIdx));
                for thisStim = 1:length(targetIdx)
                    thisStimName = stimNames{targetIdx(thisStim)}; % Single braces access contents cleanly
                    
                    fprintf('    Executing Pipeline Processing on: %s\n', thisStimName);
                    
                    %  Compute speed traces and variance shuffles (now updates sessionFileInfo natively)
                    [response, sessionFileInfo] = getRunningSpeedTuningCurves(sessionFileInfo, thisStimName, useZScoredSignals);  %shuffle = true
                    
                    % Classify templates across both structures and save internally
                    response = classifySpeedTuningFromCurves(sessionFileInfo, response, signalName); 

                
                    % Define thresholds 
                    p_thresh  = 0.01; 
                    r2_thresh = 0.1;
                    
                    % Find indices that pass  tuning criteria
                    isTuned = (response.tuningCurve.(signalName).pValFull <= p_thresh) & ...
                              (response.tuningCurve.(signalName).classification.R2 > r2_thresh);
                    
                    % ROIs that fail this criteria are your biologically 'untuned' cohort!
                    untunedIdx = find(~isTuned);
                    tunedIdx   = find(isTuned);
                    
                    fprintf('True Tuned Boutons: %d\n', length(tunedIdx));
                    fprintf('Untuned/Noisy Boutons: %d\n', length(untunedIdx));
                    % example plotting function 
                    plotAllSpeedTuningCategories(response,signalName)
                    
                    % example plotting code 

                    % targetStruct = 'tuningCurve';              % or 'tuningCurveFixedBins'
                    % useField = signalName;         % change if needed
                    % plotOpts = struct();
                    % plotOpts.maxPlots = 50;                    % how many examples to show
                    % plotOpts.nCols = 10;
                    % plotOpts.showSEM = true;
                    % plotOpts.onlySignificantMoving = true;    % set true if you want shuffle-significant only
                    % plotOpts.sortBy = 'R2';                    % 'R2', 'preferredSpeed'
                    % plotSpeedTuningCategoryExamples(response, targetStruct, useField, 'lowpass', plotOpts);
                    % plotSpeedTuningCategoryExamples(response, targetStruct, useField, 'highpass', plotOpts);
                    % plotSpeedTuningCategoryExamples(response, targetStruct, useField, 'bandpass', plotOpts);
                    % plotSpeedTuningCategoryExamples(response, targetStruct, useField, 'trough_inverted', plotOpts);


                    clear response 
                end 
            else
                fprintf('    No Darkness or GrayScreen stim files detected in this session layout.\n');
            end
            
        catch ME
            % Error Catch Block: Formats and logs any unexpected pipeline failures to your target CSV path
            timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
            errLine = "NaN"; if ~isempty(ME.stack), errLine = string(ME.stack(1).line); end
            errFunc = "Unknown"; if ~isempty(ME.stack), errFunc = string(ME.stack(1).name); end
            
            newErrorRow = table(timestamp, string(mousenumber), string(sessionName), ...
                                string(ME.message), errFunc, errLine, 'VariableNames', logHeaders);
            writetable(newErrorRow, logFilePath, 'WriteMode', 'append');
            fprintf('    ERROR LOGGED: %s (Line: %s)\n', ME.message, errLine);
        end
                       
    end 
end
fprintf('\n=== Batch Run Finished! Successfully completed processing loop. ===\n');