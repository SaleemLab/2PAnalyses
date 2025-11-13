% Script to find and combine multiple VR sessions 
clear; clc;

%% Define mouse information
filteredTable = filterMasterTable('Exclude', 0, 'Suite2PPreprocessing', 1);

%mouseInfo = sessionsToProcess(filteredTable);
mouseInfo = {
    'M25058', {'20250529'};...
    % Add new mice here, e.g.: 'M25042', {'YYYYMMDD', 'YYYYMMDD'};
};
%% Define error log 
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'CombineVRResponseRuns_20251016.csv');
logHeaders = {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'};

% If the log file doesn't exist, create it with headers
if ~exist(logFilePath, 'file')
    logTable = table('Size', [0, numel(logHeaders)], 'VariableTypes', repmat("string", 1, numel(logHeaders)), 'VariableNames', logHeaders);
    writetable(logTable, logFilePath);
end
fprintf('Error logging enabled. Log will be saved to: %s\n', logFilePath);

%% START PROCESSING LOOP
% Loop through each mouse defined in mouseInfo
for thisMouse = 1:size(mouseInfo, 1)
    
    mousenumber = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};
    
    fprintf('Processing Mouse: %s\n', mousenumber);
    % Loop through each session for the current mouse
    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        
        fprintf('\n-- Processing Session: %s --\n', sessionName);
        
        try
            % Dynamically load sessionFileInfo file to run proceeding
            % analyses 
            sessionFileInfoFilePath = findSessionFileInfoFilePath(mousenumber,sessionName);
            disp('Loading sessionFileInfo...');
            load(sessionFileInfoFilePath)
          
            % Combine multiple VR runs within sessions 
      
            fprintf('Checking for and combining Response across multiple VR runs if present')
            responseVRRuns = loadDataStructuresByKeyword(sessionFileInfo, 'VRCorr', 'Response');
            processedTwoPDataVRRuns = loadDataStructuresByKeyword(sessionFileInfo, 'VRCorr', 'processedMergedBonsaiSuite2pData');
            response = combineResponseForVRRuns(responseVRRuns, sessionFileInfo);
            plotSortedPopulationResponse_OddEven(sessionFileInfo, response, 'dFFNeuropilCorrected', true);
            fprintf('  Session %s processed successfully!\n', sessionName);

        catch ME
          
            % Still display the warning for immediate feedback
            warning('    Error processing %s, session %s: %s', mousenumber, sessionName, ME.message);
            fprintf(2, '    Error in function %s, line %d.\n', ME.stack(1).name, ME.stack(1).line);
            
            % Prepare the data for the log file
            timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
            mouseStr = string(mousenumber);
            sessionStr = string(sessionName);
            % Clean up the error message to remove newlines for better CSV formatting
            errorMessage = string(strrep(ME.message, newline, ' '));
            errorFunction = string(ME.stack(1).name);
            errorLine = ME.stack(1).line;
            
            % Create a table with the new error information
            errorData = table(timestamp, mouseStr, sessionStr, errorMessage, errorFunction, errorLine, ...
                'VariableNames', logHeaders);
            
            % Append the error data to the CSV file
            writetable(errorData, logFilePath, 'WriteMode', 'append', 'WriteVariableNames', false);
        end
    end
end
disp('All mice and sessions VR-ResonseRuns processed!');