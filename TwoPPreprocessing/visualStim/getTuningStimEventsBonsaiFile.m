function [bonsaiData, sessionFileInfo] = getTuningStimEventsBonsaiFile(sessionFileInfo, stimName, stimTypeTableName)
    % Get the Bonsai stim-events/Trial Parameters table (if saved) and add them to a .mat file.
    % CAUTION! New stimuli will need to be defined as a new case. 
    % Add all stim identities to bonsaiData.stimType 
    % 
    % Inputs:
    %   - sessionFileInfo (.mat) : File info .mat file
    %   - stimName (str) : Name of the stimulus to process, e.g., 'DirTuning' 
    %   - stimTypeTableName (str) : Name of CSV containing the stim specific information (e.g., 'StimEvents')  
    %  
    % Output:
    %   - bonsaiData (struct) : Collated stim-events table; stimulus identity always saved as StimType 
    %   
    % Aman and Sonali - Feb 2025

    if nargin < 3
        stimTypeTableName = 'StimEvents';
    end

    % Find matching stimulus files
    isTuningStim = strcmp({sessionFileInfo.stimFiles.name}, stimName);
    iStim = find(isTuningStim, 1);  % Get first match

    if isempty(iStim)
        error('Stimulus "%s" not found in sessionFileInfo.stimFiles.', stimName);
    end

    stimFileName = sprintf('%s_%s_BonsaiData_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
    
    sessionFileInfo.stimFiles(iStim).BonsaiData = ...
        fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

    tuningFilePath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, stimTypeTableName);
    stimEventsTable = readtable(tuningFilePath);
    
        % Normalize stimName into a known label for switch
    if contains(stimName, 'DirTuning')
        stimTypeKey = 'DirTuning';
    elseif contains(stimName, 'DotMotion_SpeedTuning')
        stimTypeKey = 'DotMotion_SpeedTuning';
    elseif contains(stimName, 'RFMapping')
        stimTypeKey = 'RFMapping';
    else
        stimTypeKey = 'Unknown';
    end

    % Extract relevant columns based on stimulus type; If adding new cases 
    % save stimulus identity in bonsaiData.stimType 
    switch stimTypeKey
        case 'DirTuning'
            bonsaiData.bonsaiStimOnset = stimEventsTable.Var2;
            bonsaiData.stimOnsetRenderFrameIdx = stimEventsTable.Var1;
            bonsaiData.stimType = round(rad2deg(stimEventsTable.Var5));
        
        case 'DotMotion_SpeedTuning'
            bonsaiData.bonsaiStimOnset = stimEventsTable.BonsaiTime;
            bonsaiData.ArduinoTime = stimEventsTable.ArduinoTime;
            bonsaiData.stimType = stimEventsTable.VelX1;
            bonsaiData.stimID = stimEventsTable.Id;
        case 'RFMapping'
            bonsaiData.stimID = stimEventsTable.Var2;
            bonsaiData.delay = stimEventsTable.Var4; 
            bonsaiData.duration = stimEventsTable.Var6; 
            bonsaiData.diameter = stimEventsTable.Var8;
            bonsaiData.positionX = stimEventsTable.Var10;
            bonsaiData.positionY = stimEventsTable.Var12;
            bonsaiData.contrast = stimEventsTable.Var14;
            bonsaiData.spatialFrequency = stimEventsTable.Var16;
            bonsaiData.temporalFrequency = stimEventsTable.Var18;
            bonsaiData.orientation = stimEventsTable.Var20;
            bonsaiData.bonsaiStimOnset = stimEventsTable.Var21;
            bonsaiData.RenderFrameCount = stimEventsTable.Var22;
            bonsaiData.LastSyncPulseTime = stimEventsTable.Var23;
            bonsaiData.ArduinoTime = stimEventsTable.Var24;

        otherwise
            error('Unknown stimulus type: %s', stimName);
    end

    bonsaiData.stimEventsTable = stimEventsTable; 

    % Save the extracted Bonsai data
    save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

end
