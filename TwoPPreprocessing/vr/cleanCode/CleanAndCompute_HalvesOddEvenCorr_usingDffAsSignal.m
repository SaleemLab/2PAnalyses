% process_all_mice_sessions_combined_final.m
% Combined script with dynamic parameter loading based on TypeImaged from filteredTable.
clear; clc;
%% DEFINE MICE, SESSIONS, AND KEYWORDS
% Define all mice and their sessions to be processed
% The mouseInfo is generated from the filteredTable
vrKeywords = {'VRCorr'};
rfKeywords = {'RFMapping'};
signalName = 'dFF'; %changed to dff 2026 jan
% filteredTable now holds the key metadata (TypeImaged) needed for parameter setting.

filteredTable = filterMasterTable('Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'TargetArea', 'ENTm'); % 'TypeImaged', 'Boutons'

mouseInfo = sessionsToProcess(filteredTable);

%% FOR SANITY 

totalSessionsToProcess = 0;
sessionsProcessedCount = 0;

for i = 1:size(mouseInfo, 1)
    totalSessionsToProcess = totalSessionsToProcess + length(mouseInfo{i, 2});
end
%% DEFINE PARAMETER SETS BASED ON IMAGING TYPE
% Parameters for Somas Imaging;


%% INITIALISE ERROR LOG
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'clearCode_rerunWithDffAsSignalVarianceTuning_20260113.csv');
logHeaders = {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'};
% If the log file doesn't exist, create it with headers
if ~exist(logFilePath, 'file')
    logTable = table('Size', [0, numel(logHeaders)], 'VariableTypes', repmat("string", 1, numel(logHeaders)), 'VariableNames', logHeaders);
    writetable(logTable, logFilePath);
end
fprintf('Error logging enabled. Log will be saved to: %s\n', logFilePath);

%% START CLEANING LOOP
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
            % Find the path to the session's info file
        sessionInfoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if ~isfile(sessionInfoPath)
            errorMessage = 'InfoMissing';
            return;
        end
        
        
        % Load the session file info
        D = load(sessionInfoPath, 'sessionFileInfo');
        sfi = D.sessionFileInfo; % <-- sfi is assigned here
        
        % Find the correct VR stimulus file
        stimFileNames = string({sfi.stimFiles.name});
        vrFileIndex = find(contains(stimFileNames, "VRCorr") & contains(stimFileNames, "CombinedRuns"), 1);
        
        if isempty(vrFileIndex) % Fallback: Find last run if no 'CombinedRuns'
            vrFileIndex = find(contains(stimFileNames, "VRCorr") & ~contains(stimFileNames, "CombinedRuns"), 1, 'last'); 
        end
        
        if isempty(vrFileIndex)
            errorMessage = 'NoVRCorr';
            return; 
        end
        
        % Get the response file path
        responseFilePath = sfi.stimFiles(vrFileIndex).Response;
        if isempty(responseFilePath) || ~isfile(responseFilePath)
            errorMessage = 'RespMissing'; 
            return; 
        end
        
        % Load the response data
        response = load(responseFilePath, 'response');
        
        % Check if the signal exists in the lap data structure
        if ~isfield(response.response, 'lapPositionActivity') || ~isfield(response.response.lapPositionActivity, signalName)
            errorMessage = 'SignalMissing'; 
            return;
        end

        % rerun these
        % [~,~,~] = checkOddEvenCorrelation(sfi, response.response, signalName); % Run on final response
        % [~,~,~] = checkHalvesCorrelation(sfi, response.response, signalName);
        [~,~,~] = computeVarianceAcrossPositionBins(sfi, response.response, signalName);

        catch ME
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
disp('All mice and sessions processed! Check error log file...');