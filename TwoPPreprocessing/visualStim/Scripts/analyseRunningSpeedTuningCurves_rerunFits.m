% used this script to rerun the gaussian template fits 

clear; clc;

pairs = struct; 
% pairs.M26005 = {'20260305', '20260306', '20260311', '20260318', '20260321', '20260322'};
% pairs.M26004 = {'20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322'};
% pairs.M25131 = {'20260312', '20260313', '20260314', '20260318', '20260321', '20260322'};
% pairs.M25126 = {'20260311', '20260312', '20260313'};
% pairs.M25132 = {'20260226', '20260228', '20260313'};
% pairs.M26003 = {'20260322', '20260324', '20260325'};
% pairs.M25133 = {'20260224'};

% pairs.M25132 = {'20260219', '20260220', '20260221', '20260223', '20260228', '20260303', '20260306', '20260312', '20260226', '20260228', '20260313'};
% pairs.M26003 = {'20260316', '20260317', '20260320', '20260321', '20260322', '20260324', '20260325'};
% pairs.M25133 = {'20260219', '20260220', '20260221', '20260223', '20260224'};

pairs.M25132 = {'20260306', '20260219'};

paramsBoutons.signalName = 'dFFNeuropilCorrected';
paramsBoutons.useZScoredProcessedSignals = false; 
paramsSomas.signalName   = 'spks'; 
paramsSomas.useZScoredProcessedSignals = false; 

filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'Darkness', 'GrayScreen'});
mouseInfo = sessionsToProcess(filteredTable);

logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'DarkGary_20260620_VISpAndRSpDay5_200.csv');
logHeaders = {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'};

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

targetStruct = 'tuningCurve';

for thisMouse = 1:size(mouseInfo, 1)
    mousenumber = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};
    
    fprintf('Processing Mouse: %s\n', mousenumber);
    
    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        sessionsProcessedCount = sessionsProcessedCount + 1;
        
        fprintf('  Processing Session [%d/%d]: %s  %s \n', sessionsProcessedCount, totalSessionsToProcess, mousenumber, sessionName);
        
        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if isempty(infoPath) || ~isfile(infoPath)
             fprintf('    WARNING: File path missing or invalid for %s - %s. Skipping.\n', mousenumber, sessionName);
             continue; 
        end
        
        matchRow = find(strcmp(filteredTable.MouseID, mousenumber) & strcmp(filteredTable.Session, sessionName), 1);
        if isempty(matchRow)
            fprintf('    WARNING: Could not resolve matching row inside metadata table for %s - %s. Skipping.\n', mousenumber, sessionName);
            continue;
        end
        
        typeImaged = filteredTable.TypeImaged{matchRow};
        fprintf('   Selecting params for imaged type %s .\n', typeImaged);
        
        if strcmpi(typeImaged, 'Boutons')
            signalName = paramsBoutons.signalName;
        else
            signalName = paramsSomas.signalName;
        end
        
        % Check variable field name saved inside the file structure to handle standard/split versions
        try
            loadedInfo = load(infoPath, 'sessionFileInfo');
            sessionFileInfo = loadedInfo.sessionFileInfo;
            stimNames = {sessionFileInfo.stimFiles.name}; 
            
            targetIdx = find(contains(stimNames, {'Darkness', 'GrayScreen'}));
            
            if any(targetIdx)
                fprintf('    Found %d target stimulus track file(s).\n', length(targetIdx));
                for thisStim = 1:length(targetIdx)
                    thisStimName = stimNames{targetIdx(thisStim)}; 
                    
                    stimFileName = sprintf('%s_%s_Response_%s.mat', mousenumber, sessionName, thisStimName);
                    fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
                    
                    if ~isfile(fileFullPath)
                        fprintf('    WARNING: Response file missing on disk: %s. Skipping.\n', stimFileName);
                        continue;
                    end
                    
                    fprintf('    Loading Response File: %s\n', stimFileName);
                    loadedData = load(fileFullPath, 'response');
                    response = loadedData.response;
                    
                    fprintf('    Executing Pipeline Processing...\n');
                    response = classifySpeedTuningFromCurves(sessionFileInfo, response, signalName); 
                    if strcmpi(typeImaged, 'Somas')
                        response = classifySpeedTuningFromCurves(sessionFileInfo, response, 'dFFNeuropilCorrected'); 
                        plotAllSpeedTuningCategories(response, 'dFFNeuropilCorrected')
                    end 
                    
%                     plotAllSpeedTuningCategories(response, signalName)
                    % 
                    % plotOpts = struct();
                    % plotOpts.maxPlots = 24;                    
                    % plotOpts.nCols = 4;
                    % plotOpts.showSEM = true;
                    % plotOpts.onlySignificantMoving = true;    
                    % plotOpts.sortBy = 'R2';                    
                    
                    % % fig1 = plotSpeedTuningCategoryExamples(response, targetStruct, signalName, 'lowpass', plotOpts);
                    % % fig2 = plotSpeedTuningCategoryExamples(response, targetStruct, signalName, 'highpass', plotOpts);
                    % % fig3 = plotSpeedTuningCategoryExamples(response, targetStruct, signalName, 'bandpass', plotOpts);
                    % % fig4 = plotSpeedTuningCategoryExamples(response, targetStruct, signalName, 'trough_inverted', plotOpts);
                    
                    % saveFolder = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
                    % if ~exist(saveFolder, 'dir'), mkdir(saveFolder); end
                    % 
                    % pdfPath = fullfile(saveFolder, sprintf('%s_%s_Response_%s_SpeedTuningExamples.pdf', mousenumber, sessionName, thisStimName));
                    % if exist(pdfPath, 'file'), delete(pdfPath); end
                    % 
                    % figHandles = {fig1, fig2, fig3, fig4};
                    % for fH = 1:length(figHandles)
                    %     if ~isempty(figHandles{fH}) && ishandle(figHandles{fH})
                    %         exportgraphics(figHandles{fH}, pdfPath, 'Append', true);
                    %         close(figHandles{fH});
                    %     end
                    % end
                    % 
                    % figPop = plotAllSpeedTuningCategories(response, signalName);
                    % if ~isempty(figPop) && ishandle(figPop)
                    %     popFigPath = fullfile(saveFolder, sprintf('%s_%s_Response_%s_BestExamples.png', mousenumber, sessionName, thisStimName));
                    %     exportgraphics(figPop, popFigPath, 'Resolution', 300);
                    %     close(figPop);
                    % end
                    
                    clear response 
                end 
            else
                fprintf('    No Darkness or GrayScreen stim files detected in this session layout.\n');
            end
            
        catch ME
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