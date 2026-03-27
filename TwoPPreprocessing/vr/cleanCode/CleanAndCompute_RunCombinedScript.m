clear; clc;
%% DEFINE MICE, SESSIONS, AND KEYWORDS
% Define all mice and their sessions to be processed
% The mouseInfo is generated from the filteredTable
vrKeywords = {'VRCorr', 'BaselineCorridor', 'LandManipCorridor'};
rfKeywords = {'RFMapping'};
doNotCombine = {'M25040_VRCorr_20250507_00001', 'M25040_VRCorr_20250507_00002', 'M25057_VRCorr_20250526_00001', 'M25057_VRCorr_20250526_00002', 'M25126_VRCorr_20260123_00001', 'M25126_VRCorrBaseline_20260123_00002', 'M25126_VRCorrWithManipulations_20260123_00003', 'M25132_BaselineCorridor_20260219_00001', 'M25132_BaselineCorridor_20260219_00002'};
signalName = 'dFF'; %changed to dff 2026 jan 
% filteredTable now holds the key metadata (TypeImaged) needed for parameter setting.

filteredTable = filterMasterTable('Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'MouseID', {'M25133', 'M25132'}, ...
     'Session', {'20260219','20260220','20260221','20260223','20260224', '20260226'}); 
mouseInfo = sessionsToProcess(filteredTable);

%% FOR SANITY 

totalSessionsToProcess = 0;
sessionsProcessedCount = 0;

for i = 1:size(mouseInfo, 1)
    totalSessionsToProcess = totalSessionsToProcess + length(mouseInfo{i, 2});
end
%% DEFINE PARAMETER SETS BASED ON IMAGING TYPE
% Parameters for Somas Imaging;
paramsSomas.interpRate = 60;
%paramsSomas.frameRate = 30; %changed
paramsSomas.pdthreshold = 10;
paramsSomas.isZcorrected = false;
paramsSomas.zScoreProcessedSignals = true;
paramsSomas.applyTemporalSmoothing = true;
paramsSomas.prctlF = 8; % The percentile from which to take F0 (baseline F).
paramsSomas.windowSize = 60; % The rolling window over which to calculate F0.



% Parameters for Bouton Imaging
paramsBoutons.interpRate = 60;
%paramsBoutons.frameRate = 7.28;
paramsBoutons.pdthreshold = 10;
paramsBoutons.isZcorrected = true;
paramsBoutons.zScoreProcessedSignals = true;
paramsBoutons.applyTemporalSmoothing = true;
paramsBoutons.prctlF = 8; % The percentile from which to take F0 (baseline F).
paramsBoutons.windowSize = 60; % The rolling window over which to calculate F0.

% Common processing parameters
rfpreStimTime = 2;    % seconds
rfpostStimTime = 3;     % seconds
method = 2;             % Method for PSTH extraction (e.g., 2 for mean)

%% INITIALISE ERROR LOG
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'runNewVRCorr_rerunCombinedRuns_20260302.csv');
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
            
            sessionsProcessedCount = sessionsProcessedCount + 1;
            fprintf(' --------------- Session: %d of %d ------------ \n', ...
                 sessionsProcessedCount, totalSessionsToProcess);
            
            % Find the row in filteredTable corresponding to the current mouse and session
            sessionIdx = find(strcmp(filteredTable.MouseID, mousenumber) & strcmp(filteredTable.Session, sessionName), 1);

            if isempty(sessionIdx)
                error('Session info not found in filteredTable for %s, %s.', mousenumber, sessionName);
            end

            typeImaged = filteredTable.TypeImaged{sessionIdx}; % Get the imaging type



            fprintf('  Type Imaged: %s. Loading specific processing parameters.\n', typeImaged);

            % Assign parameters based on the imaging type
            if strcmpi(typeImaged, 'Somas')
                currentParams = paramsSomas;
            elseif strcmpi(typeImaged, 'Boutons')
                currentParams = paramsBoutons;
            else
                error('Unknown TypeImaged: %s. Cannot set processing parameters.', typeImaged);
            end


            interpRate = currentParams.interpRate;
            zScoreProcessedSignals = currentParams.zScoreProcessedSignals;
            isZcorrected = currentParams.isZcorrected;
            applyTemporalSmoothing = currentParams.applyTemporalSmoothing;
            prctlF = currentParams.prctlF;
            windowSize = currentParams.windowSize;
            pdThresholdForStimEvents = currentParams.pdthreshold;

            stimList = getStimList(mousenumber, sessionName); % If the names of tif files have been changed after running through suite2p this will not work! 

            fprintf('  Found stimuli: %s\n', strjoin(stimList, ', '));

            vrStimNames = {};
            otherVisualStim = {};
            rfStimNames = {};
            % Will currently pull out all Corr stim names; baseline and
            % manipulations grouped together 
            for thisStim = 1:length(stimList)
                if any(contains(stimList{thisStim}, vrKeywords, 'IgnoreCase', true))
                    vrStimNames{end+1} = stimList{thisStim};
                else
                    otherVisualStim{end+1} = stimList{thisStim};
                    if any(contains(stimList{thisStim}, rfKeywords, 'IgnoreCase', true))
                        rfStimNames{end+1} = stimList{thisStim};
                    end
                end
            end

            sessionFileInfoFilePath = getSessionFileInfoFilePath(mousenumber, sessionName);
            load(sessionFileInfoFilePath)
           
            % C. Process all VR Stimuli
            if ~isempty(vrStimNames)
                fprintf('Found %d VR stimulus file(s).\n', length(vrStimNames));

                shouldProcessSeparately = any(ismember(vrStimNames, doNotCombine));


                % Path condition: TRUE if it's a single run OR if it's a multi-run exception.
                runAsIndividualLoop = (isscalar(vrStimNames)) || (length(vrStimNames) > 1 && shouldProcessSeparately);
                % Execute Shuffle and Analysis

                

                if ~shouldProcessSeparately

                    % Combination Required (N>1 and NOT an exception) ---
                    fprintf('Processing via COMBINED PATH...This is to keep Sonali sane\n');
                    % Combine Response sturctures
                    fprintf('Checking for and combining Response across multiple VR runs...\n')
                    
                    responseVRRuns = loadDataStructuresByKeyword(sessionFileInfo, 'BaselineCorridor', 'Response');
    
                    [response, sessionFileInfo] = combineResponseForVRRuns(responseVRRuns, sessionFileInfo);

                    % Compute Shuffle Matrix (Stitch ALL raw signals together)
                    fprintf('Computing shuffle matrix for COMBINED signals...\n')

                    [response, sessionFileInfo] = computeShuffleMatrixForSession(sessionFileInfo, response, vrStimNames);
                    plotSortedPopulationResponse_OddEven(sessionFileInfo, response, signalName, true)

                    % Run checks
                    [~,~] = getRangeSignificance_fromShuffle(sessionFileInfo, response); % Run on final response
                    [~, ~] = getPeakSignificance_fromShuffle(sessionFileInfo, response); % Run on final response
                    [~,~,~]= computeVarianceAcrossPositionBins(sessionFileInfo, response); % Run on final response
                    [~,~,~] = checkOddEvenCorrelation(sessionFileInfo, response); % Run on final response
                    [~,~,~] = checkHalvesCorrelation(sessionFileInfo, response);

                    % For boutons specifically
                    if strcmpi(typeImaged, 'Boutons')
                        disp('Finding highly correlated boutons for this session')
                        [~,~,~,~] = findHightlyCorrelatedROIs(sessionFileInfo); %TODO: flatten!!! 
                    end

                end
            end

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
disp('All mice and sessions processed!');