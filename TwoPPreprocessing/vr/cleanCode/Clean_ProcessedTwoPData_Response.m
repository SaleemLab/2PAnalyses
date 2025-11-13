% =========================================================================
% Loads processedTwoPData data for 

% Sonali Oct 2025, 
% =========================================================================
clear; clc; close all;

%% --- Define filepaths and criteria ---
masterTableFilePath = "Z:\ibn-vision\USERS\Sonali\datatable\MasterExpDatatable.csv";
rootDir = 'Z:\ibn-vision\DATA\SUBJECTS';
% DEFINE where to save the log file for any unprocessed stimuli
logFilePath = "Z:\ibn-vision\USERS\Sonali\errorLogs\unprocessed_log20251013.csv";


%% --- Initialize an error log to track unprocessed files ---
unprocessedStimuli = struct('MouseID', {}, 'Session', {}, 'StimulusFile', {}, 'Reason', {});

%% --- Load, parse and filter master table ---

filteredTable = filterMasterTable('TypeImaged', 'Boutons', ...
    'HasStimulus', 'VRCorr');

if isempty(filteredTable)
    disp('No sessions matched the criteria. Stopping script.');
    return;
end
disp([num2str(height(filteredTable)) ' sessions found. Starting analysis loop...']);

%% --- Iterate through each session and rerun scripts ---
for thisRow = 1:height(filteredTable)
    currentRow = filteredTable(thisRow, :);
    
    % Extract the MouseID and Session for the current row
    mouseID = currentRow.MouseID;
    session = currentRow.Session;
    
    fprintf('\nProcessing Mouse: %s, Session: %s (%d of %d)\n', mouseID, session, thisRow, height(filteredTable));
    
    % Construct the expected filename and the full path
    fileName = sprintf('%s_%s_sessionFileInfo.mat', mouseID, session);
    sessionFilePath = fullfile(rootDir, mouseID, 'Analysis', session, fileName);
    
    % Check if the sessionFileInfo file actually exists before trying to load
    if isfile(sessionFilePath)
        fprintf('Loading session file: %s\n', sessionFilePath);
        load(sessionFilePath, 'sessionFileInfo');

        [sessionFileInfo, allFilesOK] = findAndAppendMissingFilePaths(sessionFileInfo);

        % If any files were updated, you should save the changes
        if allFilesOK
            disp('All files are present. Saving updated session info.');
            save(sessionFilePath, 'sessionFileInfo');
        else
            disp('Some files were missing. Please check the warnings.');
            save(sessionFilePath, 'sessionFileInfo');
        end
        
        % Find all stimulus files containing 'VRCorr'
        allStimNames = {sessionFileInfo.stimFiles.name};
        isVRCorr = contains(allStimNames, 'VRCorr');
        vrCorrStimFiles = allStimNames(isVRCorr);
  %%      
        % --- Iterate through each found VRCorr stimulus file ---
        if ~isempty(vrCorrStimFiles)
            fprintf(' Found %d VRCorr files. Iterating through each...\n', length(vrCorrStimFiles));
            
            for thisVRIdx = 1:length(vrCorrStimFiles)
                currentVRFile = vrCorrStimFiles{thisVRIdx};
                fprintf('    -- Analysing VR file: %s\n', currentVRFile);
                
                % Use a try-catch block to  handle errors for any stimulus
                try
                    % Rerun neuropil correction and DFF; Append to processedTwoPData
                    [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, vrStimName);
                    % Rerun lapPositionActivity
                    [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, currentVRFile, true); % overwrite = true; 
                    % Save one copy of 2P frame index 
                    response.lapPosition2PFrameIdx = response.lapPosition2PFrameIdx(1,:,:);
                    % Save response 
                    stimIdx = find(strcmp(currentVRFile, {sessionFileInfo.stimFiles.name}));
                    save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');

                    
                catch ME % ME is an object containing information about the error
                    % If an error occurs, log it and continue to the next file
                    fprintf('ERROR processing %s. Logging and skipping.\n', currentVRFile);
                    
                    newEntry.MouseID = mouseID;
                    newEntry.Session = session;
                    newEntry.StimulusFile = currentVRFile;
                    newEntry.Reason = ME.message; % Get the specific error message
                    unprocessedStimuli(end+1) = newEntry;
                end
            end
        else
            disp('No VRCorr stimulus files found in this session.');
        end
        
    else
        warning('Session file not found, skipping session: %s', sessionFilePath);
        % Log this session-level failure
        newEntry.MouseID = mouseID;
        newEntry.Session = session;
        newEntry.StimulusFile = 'N/A (Session File Missing)';
        newEntry.Reason = 'sessionFileInfo.mat not found at expected path.';
        unprocessedStimuli(end+1) = newEntry;
    end
end

%% --- FINAL REPORT: Display and save any stimuli that were not processed ---

disp('All sessions processed.');

if isempty(unprocessedStimuli)
    disp('All stimuli were processed successfully.');
else
    disp('The following stimuli were NOT processed due to errors:');
    % Convert the struct to a table for a clean display and saving
    errorTable = struct2table(unprocessedStimuli);
    disp(errorTable);
    
    % Save the error log to a permanent file
    try
        writetable(errorTable, logFilePath);
        fprintf('Error log has been saved to: %s\n', logFilePath);
    catch saveME
        warning('Could not save the error log file. Please check permissions or the path.');
        disp(saveME.message);
    end
end

