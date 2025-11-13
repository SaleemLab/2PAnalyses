% process_all_mice_sessions.m
% Script to automate processing for multiple mice and sessions with error logging.
clear; clc;

%% DEFINE MICE, SESSIONS, AND KEYWORDS
% Define all mice and their sessions to be processed
mouseInfo = {

    'M25058', {'20250529'};...
    % 'M25040', {'20250511A', '20250507'};...
    % Add new mice here, e.g.: 'M25042', {'YYYYMMDD', 'YYYYMMDD'};
};
% Keywords to identify stimulus types from filenames
vrKeywords = {'VRCorr'};
% rfKeywords = {'RFMapping'};

% filteredTable = filterMasterTable('Exclude', 0, 'Suite2PPreprocessing', 1, 'TargetArea', {'V1'});
% mouseInfo = sessionsToProcess(filteredTable);
%% SET PROCESSING PARAMETERS:


interpRate = 60;        % Hz, for resampling
frameRate = 7.50;       % Frame rate 
pdthreshold = 10;       % Photodiode threshold
planeNums = 8;          % Number of planes in recording
channelsSaved = 1;      % Number of channels saved
method = 2;  
isZcorrected = true;
% Method for PSTH extraction (e.g., 2 for mean)
% RF Mapping PSTH parameters
% rfpreStimTime = 0.5;    % seconds
% rfpostStimTime = 3;     % seconds
% applyNeuropilCorrection = true; 
% calculateDFF = true;            

%% INITIALISE ERROR LOG 
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'runVRPipeline_20251110.csv');
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
            % Dynamically get the list of stimuli for this specific session
            stimList = getStimList(mousenumber, sessionName);
            fprintf('  Found stimuli: %s\n', strjoin(stimList, ', '));
            
            % Identify which stimuli are VR and which are RF based on keywords
            vrStimNames = {};
            % rfStimNames = {};
            for k = 1:length(stimList)
                if any(contains(stimList{k}, vrKeywords, 'IgnoreCase', true))
                    vrStimNames{end+1} = stimList{k};
                % elseif any(contains(stimList{k}, rfKeywords, 'IgnoreCase', true))
                %     rfStimNames{end+1} = stimList{k};
                end
            end
            
            % General file processing for the entire session
            sessionFileInfo = get2PsessionFilePaths(mousenumber, sessionName, stimList);
            sessionFileInfo = get2PMetadata(sessionFileInfo);
            [sessionFileInfo] = get2PFrameTimes_SpeedyVersion(sessionFileInfo, planeNums, isZcorrected);
            sessionFileInfo = processPeripheralFiles(sessionFileInfo);
            sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo);
            
            % B. Process all VR Stimuli
            if ~isempty(vrStimNames)
                fprintf('Found %d VR stimulus file(s).\n', length(vrStimNames));
                for thisVRStim = 1:length(vrStimNames)
                    vrStimName = vrStimNames{thisVRStim};
                    fprintf('Processing VR Stim: %s\n', vrStimName);

                    [~, sessionFileInfo]        = getVRBonsaiFiles(sessionFileInfo, vrStimName);
                    [~, sessionFileInfo]        = findBonsaiPeripheralLag(sessionFileInfo, 1, 60, vrStimName);
                    [~, sessionFileInfo]        = alignVRBonsaiToPeripheralData(sessionFileInfo,vrStimName);
                    [~, ~, ~, sessionFileInfo]  = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo,60,'TwoPFrameTime', vrStimName, true, true);
                    [~, sessionFileInfo]        = extractVRAndPeripheralData(sessionFileInfo, vrStimName, true);
                    [~, sessionFileInfo]        = get2PFrameLapPositionBins(sessionFileInfo, vrStimName);
                    [~, sessionFileInfo]        = computeNeuropilCorrectionAndDFF(sessionFileInfo, vrStimName); % default is overwrite
                    [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, vrStimName, true);
                    % Get roi stability matrix

                    plotSortedPopulationResponse_OddEven(sessionFileInfo, response, 'dFFNeuropilCorrected', true);
                end
            end

            % 
            % % B.2. Combine multiple VR runs within sessions 
            % if ~isempty(vrStimName)
            %     fprintf('Checking for and combining Response across multiple VR runs if present')
            % 
            % 
            %     responseVRRuns = loadDataStructuresByKeyword(sessionFileInfo, 'VRCorr', 'Response');
            %     response = combineResponseForVRRuns(responseVRRuns, sessionFileInfo);
            %     plotSortedPopulationResponse_OddEven(sessionFileInfo, response, 'dFFNeuropilCorrected', true);
            % 
            % end
            % fprintf('  Session %s processed successfully!\n', sessionName);

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
            % --- END OF MODIFICATION ---
        end
    end
end
disp('All mice and sessions processed!');