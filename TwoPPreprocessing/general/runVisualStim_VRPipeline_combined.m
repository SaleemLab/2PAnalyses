% process_all_mice_sessions_combined_final.m
% Combined script with dynamic parameter loading based on TypeImaged from filteredTable.
clear; clc;
%% DEFINE MICE, SESSIONS, AND KEYWORDS
% Define all mice and their sessions to be processed
% The mouseInfo is generated from the filteredTable
vrKeywords = {'VRCorr'};
rfKeywords = {'RFMapping'};
doNotCombine = {'M25040_VRCorr_20250507_00001', 'M25040_VRCorr_20250507_00002', 'M25057_VRCorr_20250526_00001', 'M25057_VRCorr_20250526_00002'};
signalName = 'dFFNeuropilCorrected';
% filteredTable now holds the key metadata (TypeImaged) needed for parameter setting.

filteredTable = filterMasterTable('Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'MouseID', {'M25057', 'M25058'}); %, 'TypeImaged', 'Boutons'

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
paramsSomas.frameRate = 7.5;
paramsSomas.pdthreshold = 10;
paramsSomas.isZcorrected = false;
paramsSomas.zScoreProcessedSignals = true;
paramsSomas.applyTemporalSmoothing = true;
paramsSomas.prctlF = 8; % The percentile from which to take F0 (baseline F).
paramsSomas.windowSize = 60; % The rolling window over which to calculate F0.



% Parameters for Bouton Imaging
paramsBoutons.interpRate = 60;
paramsBoutons.frameRate = 7.28;
paramsBoutons.pdthreshold = 10;
paramsBoutons.isZcorrected = true;
paramsBoutons.zScoreProcessedSignals = true;
paramsBoutons.applyTemporalSmoothing = true;
paramsBoutons.prctlF = 8; % The percentile from which to take F0 (baseline F).
paramsBoutons.windowSize = 60; % The rolling window over which to calculate F0.

% Common processing parameters
rfpreStimTime = 0.5;    % seconds
rfpostStimTime = 3;     % seconds
method = 2;             % Method for PSTH extraction (e.g., 2 for mean)

%% INITIALISE ERROR LOG
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'runVR_VisualStimPipeline_20251203.csv');
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

            % Get the list of stimuli for this specific session from
            % Suite2p; These stim list will be extracted
            % the order of stimuli that stimuli were concatenated before
            % feeding into Suite2p.
            stimList = getStimList(mousenumber, sessionName);
            fprintf('  Found stimuli: %s\n', strjoin(stimList, ', '));

            % Identify which stimuli are VR and which are others (including RF)
            vrStimNames = {};
            otherVisualStim = {};
            rfStimNames = {};

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

            % A. General file processing for the entire session.
            sessionFileInfo = get2PsessionFilePaths(mousenumber, sessionName, stimList, 1); % Overwrite is a must
            sessionFileInfo = get2PMetadata(sessionFileInfo); % This function will not run as the tif files have been moved to a different repo.
            [sessionFileInfo] = get2PFrameTimes_SpeedyVersion(sessionFileInfo, isZcorrected); % Uses dynamic planeNums, isZcorrected
            sessionFileInfo = processPeripheralFiles(sessionFileInfo);
            sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo);
            [sessionFileInfo] = createSessionROIData(sessionFileInfo);

            % B. Process all other visual stimuli; Process Visual stimuli
            % first as it is needed to check for highly correlated boutons
            if ~isempty(otherVisualStim)
                fprintf('Found %d other visual stimulus file(s).\n', length(otherVisualStim));
                for thisOtherVisualStim = 1:length(otherVisualStim)
                    otherVisualStimName = otherVisualStim{thisOtherVisualStim};
                    fprintf('Processing Stimulus file: %s\n', otherVisualStimName);

                    try
                        % Loads Bonsai data files where appropriate
                        [~, sessionFileInfo]       = getTuningStimEventsBonsaiFile(sessionFileInfo, otherVisualStimName, true);
                        [~, ~, ~, sessionFileInfo] = ...
                            resamplAndAlignVisualStim_BonsaiPeripheralSuite2P(sessionFileInfo,interpRate,'TwoPFrameTime', otherVisualStimName);  % Plotflag false, trimNaNs true and overwrite true;
                        [~, sessionFileInfo]       = computeNeuropilCorrectionAndDFF(sessionFileInfo, otherVisualStimName, zScoreProcessedSignals, applyTemporalSmoothing, prctlF, windowSize); % Overwrite is true
                    catch
                        fprintf('Missing bonsai data structure for stimulus file: %s\n', otherVisualStimName)
                    end
                end
            end

            % C. Process all VR Stimuli
            if ~isempty(vrStimNames)
                fprintf('Found %d VR stimulus file(s).\n', length(vrStimNames));

                % C.1 Preprocessing Loop (Generates individual response.mat files for ALL runs)
                for thisVRStim = 1:length(vrStimNames)
                    vrStimName = vrStimNames{thisVRStim};
                    fprintf('Processing VR Stim: %s\n', vrStimName);

                    % --- Preprocessing Steps (Assumed correct) ---
                    [~, sessionFileInfo] = getVRBonsaiFiles(sessionFileInfo, vrStimName);
                    [~, sessionFileInfo] = findBonsaiPeripheralLag(sessionFileInfo, vrStimName, 1, interpRate);
                    [~, sessionFileInfo] = alignVRBonsaiToPeripheralData(sessionFileInfo,vrStimName);
                    [~, ~, ~, sessionFileInfo] = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo,interpRate,'TwoPFrameTime', vrStimName);
                    [~, sessionFileInfo] = extractVRAndPeripheralData(sessionFileInfo, vrStimName);
                    [~, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, vrStimName);
                    [~, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, vrStimName, zScoreProcessedSignals, applyTemporalSmoothing, prctlF, windowSize);

                    % Calculates and saves LapPositionActivity (un-shuffled)
                    [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, vrStimName);
                    clear response;
                end

                % Determine Path
                shouldProcessSeparately = any(ismember(vrStimNames, doNotCombine));


                % Path condition: TRUE if it's a single run OR if it's a multi-run exception.
                runAsIndividualLoop = (isscalar(vrStimNames)) || (length(vrStimNames) > 1 && shouldProcessSeparately);
                % Execute Shuffle and Analysis

                if runAsIndividualLoop
                    fprintf('Processing via INDIVIDUAL LOOP path.\n');

                    for thisVRStim = 1:length(vrStimNames)

                        vrStimName = vrStimNames{thisVRStim};
                        stimIdx = find(strcmp(vrStimName, {sessionFileInfo.stimFiles.name}), 1);
                        % load the individual response structure into the environment
                        load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');
                        % Compute Shuffle Matrix (Stitch to itself)
                        fprintf(' -> Computing shuffle matrix for run: %s\n', vrStimName);
                        [response, sessionFileInfo] = computeShuffleMatrixForSession(sessionFileInfo,response,{vrStimName});
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
                            [~,~,~,~] = findHightlyCorrelatedROIs(sessionFileInfo);
                        end
                    end

                else

                    % Combination Required (N>1 and NOT an exception) ---
                    fprintf('Processing via COMBINED PATH...This is to keep Sonali sane\n');
                    % Combine Response sturctures
                    fprintf('Checking for and combining Response across multiple VR runs...\n')
                    responseVRRuns = loadDataStructuresByKeyword(sessionFileInfo, 'VRCorr', 'Response');
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
                        [~,~,~,~] = findHightlyCorrelatedROIs(sessionFileInfo);
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