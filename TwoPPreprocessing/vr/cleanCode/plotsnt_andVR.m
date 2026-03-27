% process_all_mice_sessions_combined_final.m
clear; clc;

%% 1. DEFINE MICE, SESSIONS, AND KEYWORDS
vrKeywords = {'VRCorr', 'BaselineCorridor', 'LandManipCorridor'};
rfKeywords = {'RFMapping'};
doNotCombine = {'M25040_VRCorr_20250507_00001', 'M25040_VRCorr_20250507_00002', 'M25057_VRCorr_20250526_00001', 'M25057_VRCorr_20250526_00002', 'M25126_VRCorr_20260123_00001', 'M25126_VRCorrBaseline_20260123_00002', 'M25126_VRCorrWithManipulations_20260123_00003', 'M25132_BaselineCorridor_20260219_00001', 'M25132_BaselineCorridor_20260219_00002'};

filteredTable = filterMasterTable('Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'MouseID', {'M25131', 'M25126', 'M26004'}, ...
    'Session', {'20260312'}); 

mouseInfo = sessionsToProcess(filteredTable);

%% 2. INITIALIZE COUNTERS
totalSessionsToProcess = 0;
for i = 1:size(mouseInfo, 1)
    totalSessionsToProcess = totalSessionsToProcess + length(mouseInfo{i, 2});
end
sessionsProcessedCount = 0;

%% 3. DEFINE PARAMETER SETS (Original formatting)
% Parameters for Somas Imaging
paramsSomas.interpRate = 60;
paramsSomas.pdthreshold = 10;
paramsSomas.isZcorrected = false;
paramsSomas.zScoreProcessedSignals = true;
paramsSomas.applyTemporalSmoothing = true;
paramsSomas.prctlF = 8; 
paramsSomas.windowSize = 60; 

% Parameters for Bouton Imaging
paramsBoutons.interpRate = 60;
paramsBoutons.pdthreshold = 10;
paramsBoutons.isZcorrected = true;
paramsBoutons.zScoreProcessedSignals = true;
paramsBoutons.applyTemporalSmoothing = true;
paramsBoutons.prctlF = 8; 
paramsBoutons.windowSize = 60; 

% Common processing parameters
rfpreStimTime = 2;    
rfpostStimTime = 3;     
method = 2;             

%% 4. INITIALIZE ERROR LOG
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'runNewVRCorr_plotRfs_20260302.csv');
logHeaders = {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'};

if ~exist(logFilePath, 'file')
    logTable = table('Size', [0, numel(logHeaders)], 'VariableTypes', repmat("string", 1, numel(logHeaders)), 'VariableNames', logHeaders);
    writetable(logTable, logFilePath);
end

%% 5. START PROCESSING LOOP
for thisMouse = 1:size(mouseInfo, 1)
    mousenumber = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};
    
    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        sessionsProcessedCount = sessionsProcessedCount + 1;
        
        fprintf('\n-- Processing Session: %s (%d of %d) --\n', sessionName, sessionsProcessedCount, totalSessionsToProcess);
        
        try
            % Get session metadata
            sessionIdx = find(strcmp(filteredTable.MouseID, mousenumber) & strcmp(filteredTable.Session, sessionName), 1);
            if isempty(sessionIdx), error('Session info not found in filteredTable.'); end
            
            typeImaged = filteredTable.TypeImaged{sessionIdx};
            
            % Assign parameters based on type
            if strcmpi(typeImaged, 'Somas')
                currentParams = paramsSomas;
            elseif strcmpi(typeImaged, 'Boutons')
                currentParams = paramsBoutons;
            else
                error('Unknown TypeImaged: %s', typeImaged);
            end
            
            % Extract parameters for use
            applyTemporalSmoothing = currentParams.applyTemporalSmoothing;
            
            % Load file info
            sessionFileInfoFilePath = getSessionFileInfoFilePath(mousenumber, sessionName);
            load(sessionFileInfoFilePath); % Loads 'sessionFileInfo'

            % --- PROCESS VR STIM ---
            stimNames = {sessionFileInfo.stimFiles.name};
            
            % Find all VR-related stimulus names
            vrStimNames = stimNames(contains(stimNames, vrKeywords, 'IgnoreCase', true));
            
            if ~isempty(vrStimNames)
                % Check for CombinedRuns
                combinedIdx = find(contains(vrStimNames, 'CombinedRuns', 'IgnoreCase', true), 1);
                
                if ~isempty(combinedIdx)
                    % If "Combined" exists, plot ONLY that
                    finalVrList = vrStimNames(combinedIdx);
                    fprintf('  Found CombinedRuns. Plotting merged file only.\n');
                else
                    % Otherwise, plot all individual VR files found
                    finalVrList = vrStimNames;
                    fprintf('  No CombinedRuns found. Processing %d VR file(s).\n', length(finalVrList));
                end
                
                % Process the decided list
                for i = 1:length(finalVrList)
                    vrStimName = finalVrList{i};
                    
                    try
                        stimIdx = find(strcmp(vrStimName, {sessionFileInfo.stimFiles.name}), 1);
                        
                        if isempty(stimIdx)
                            fprintf('  WARNING: Stimulus %s not in sessionFileInfo. Skipping.\n', vrStimName);
                            continue;
                        end

                        % Load data and plot
                        response = load(sessionFileInfo.stimFiles(stimIdx).Response);
                        
                        plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing, 'dFF')
                        plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing, 'dFFNeuropilCorrected')
                        
                        fprintf('  Successfully processed: %s\n', vrStimName);
                        clear response;
                        
                    catch ME_inner
                        fprintf('  WARNING: Could not process %s. Error: %s\n', vrStimName, ME_inner.message);
                    end
                end
            end
           
        catch ME
            % Main Loop Error Logging
            warning('Error processing %s, session %s: %s', mousenumber, sessionName, ME.message);
            timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
            errData = table(timestamp, string(mousenumber), string(sessionName), ...
                            string(strrep(ME.message, newline, ' ')), string(ME.stack(1).name), ME.stack(1).line, ...
                            'VariableNames', logHeaders);
            writetable(errData, logFilePath, 'WriteMode', 'append', 'WriteVariableNames', false);
        end
    end
end