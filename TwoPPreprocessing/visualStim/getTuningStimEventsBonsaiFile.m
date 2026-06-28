function [bonsaiData, sessionFileInfo] = getTuningStimEventsBonsaiFile(sessionFileInfo, stimName, overwrite)
    % Get the Bonsai stim-events/Trial Parameters table (if saved) and add them to a .mat file.
    % CAUTION! New stimuli will need to be defined as a new case. 
    % Add all stim identities to bonsaiData.stimType 
    % 
    % Inputs:
    %   - sessionFileInfo (.mat) : File info .mat file
    %   - stimName (str) : Name of the stimulus to process, e.g., 'DirTuning' 
    %   - overwrite (bool) : (Optional) If true, overwrites existing file. Default false.
    %  
    % Output:
    %   - bonsaiData (struct) : Collated stim-events table; stimulus identity always saved as StimType 
    %   
    % Aman and Sonali - Feb 2025

    if nargin < 3
        overwrite = false;
    end
    
    % Normalise stimName into a known label for switch and define
    % stimTypeTbaleName to load bonsai file. 
    bonsaiData = [];
    stimTypeKey = 'Unknown';
    stimTypeTableName = 'Unknown';
    if contains(stimName, 'DirTuning')
        stimTypeKey = 'DirTuning';
        stimTypeTableName = 'StimEvents';
    elseif contains(stimName, 'DotMotion_SpeedTuning')
        stimTypeKey = 'DotMotion_SpeedTuning';
        stimTypeTableName = 'TrialParams';
    elseif contains(stimName, 'DotMotion_RFMapping')
        stimTypeKey = 'DotMotion_RFMapping';
        stimTypeTableName = 'TrialParams';
    elseif contains(stimName, 'RFMapping')
        stimTypeKey = 'RFMapping';
        stimTypeTableName = 'StimulusParams'; 
    elseif contains(stimName, 'SparseNoiseTexture') || contains(stimName, 'SparseNoise')
        stimTypeKey = 'SparseNoiseTexture';
        stimTypeTableName = 'Log';
    elseif contains(stimName, 'Contrast')
        stimTypeKey = 'Contrast';
        stimTypeTableName = 'stimEvents';
    elseif contains(stimName, 'Position')
        stimTypeKey = 'Position';
        stimTypeTableName = 'stimEvents';
    elseif contains(stimName, {'CheckerBoard', 'DriftingBar', 'FullFieldFlash', 'GrayScreen', 'Darkness', 'NoDisplay'}, "IgnoreCase",true)
        fprintf('%s stimulus does not contain BonsaiFile to load. Skipping this file..\n', stimName);
        stimTypeKey = 'N/A';
        stimTypeTableName = 'N/A';
        return
    end 
    
    % Find matching stimulus files
    isTuningStim = strcmp({sessionFileInfo.stimFiles.name}, stimName);
    iStim = find(isTuningStim, 1);  % Get first match
    if isempty(iStim)
        error('Stimulus "%s" not found in sessionFileInfo.stimFiles.', stimName);
    end
    
    stimFileName = sprintf('%s_%s_BonsaiData_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
    
    saveFilePath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
    sessionFileInfo.stimFiles(iStim).BonsaiData = saveFilePath;

    % Check if file exists and handle overwrite logic
    if exist(saveFilePath, 'file') && ~overwrite
        fprintf('Loading existing Bonsai data for %s...\n', stimName);
        loadedData = load(saveFilePath, 'bonsaiData');
        bonsaiData = loadedData.bonsaiData;
        return;
    end

    tuningFilePath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, stimTypeTableName);
    
    
    % Extract relevant columns based on stimulus type; If adding new cases 
    % save stimulus identity in bonsaiData.stimType 
    switch stimTypeKey
        case 'DirTuning'
            stimEventsTable = readtable(tuningFilePath);
            bonsaiData.bonsaiStimOnsetRaw = stimEventsTable.Timestamp; % Not sure if this is bonsaiStimOnset! 
            bonsaiData.stimOnsetRenderFrameIdx = stimEventsTable.Frame;
            bonsaiData.stimType = round(rad2deg(stimEventsTable.Value));
        
        case 'DotMotion_SpeedTuning'
            stimEventsTable = readtable(tuningFilePath);
            bonsaiData.bonsaiStimOnsetRaw = stimEventsTable.BonsaiTime;
            bonsaiData.ArduinoTimeRaw = stimEventsTable.ArduinoTime./1000;
            bonsaiData.stimType = stimEventsTable.VelX1;
            bonsaiData.stimID = stimEventsTable.Id;
            

        case 'DotMotion_RFMapping'
            stimEventsTable = readtable(tuningFilePath);
            bonsaiData.bonsaiStimOnsetRaw = stimEventsTable.BonsaiTime;
            bonsaiData.ArduinoTimeRaw = stimEventsTable.ArduinoTime./1000;
            bonsaiData.stimType = stimEventsTable.VelX1;
            bonsaiData.stimID = stimEventsTable.Id;
        case 'RFMapping'
            stimEventsTable = readtable(tuningFilePath);
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
            bonsaiData.bonsaiStimOnsetRaw = stimEventsTable.Var21;
            bonsaiData.RenderFrameCount = stimEventsTable.Var22;
            bonsaiData.LastSyncPulseTime = stimEventsTable.Var23;
            bonsaiData.ArduinoTimeRaw = stimEventsTable.Var24./1000;
        
        case 'SparseNoiseTexture' %change to a function instead..TODO @sonali


            bonsaiData.gridSize = [8 12];
            fprintf('Using defualt grid size %d for SparseNoiseTexture..\n', bonsaiData.gridSize)
            fileID=fopen(tuningFilePath);
            thisBinFile=fread(fileID);
            fclose(fileID);

            % Translate stimulus into -1:1 scale
            stimMatrix = zeros(1,length(thisBinFile));
            stimMatrix(thisBinFile==0)=-1;
            stimMatrix(thisBinFile==255)=1;
            stimMatrix(thisBinFile==128)=0;

            % Make a NxM grid from the stimulus log
            stimMatrix = reshape(stimMatrix, [bonsaiData.gridSize(1), bonsaiData.gridSize(2), length(thisBinFile)/bonsaiData.gridSize(1)/bonsaiData.gridSize(2)]);
            stimMatrix = stimMatrix(:,:,1:end-1); % The last 'stimulus'

            for thisTrial = 1:size(stimMatrix,3)
                bonsaiData.stimMatrix{thisTrial,1} = squeeze(stimMatrix(:,:,thisTrial));
            end


        case 'Contrast'
            stimEventsTable = readtable(tuningFilePath);
            bonsaiData.bonsaiStimOnsetRaw = stimEventsTable.Timestamp; % Not sure if this is bonsaiStimOnset!
            bonsaiData.stimOnsetRenderFrameIdx = stimEventsTable.Frame;
            bonsaiData.stimType = stimEventsTable.Var5;

        case 'Position'
            stimEventsTable = readtable(tuningFilePath, 'VariableNamingRule', 'preserve');
            
           
           
            bonsaiData.bonsaiStimOnsetRaw = stimEventsTable.Var2; 
            bonsaiData.stimOnsetRenderFrameIdx = stimEventsTable.Var1;
            
            idxLoc = strcmp(stimEventsTable.Var4, 'LocationX');
            bonsaiData.locationX = stimEventsTable.Var5(idxLoc);
            
            idxExt = strcmp(stimEventsTable.Var4, 'ExtentX');
            bonsaiData.extentX = stimEventsTable.Var5(idxExt);
            
            idxTF = strcmp(stimEventsTable.Var4, 'TF');
            bonsaiData.TF = stimEventsTable.Var5(idxTF);
            % bonsaiData.stimType = string(bonsaiData.locationX) + "deg_" + string(bonsaiData.extentX) + "width";

            
        otherwise
            error('Unknown stimulus type: %s', stimName);
    end
    
    try % SparseNoiseTexture will not have stimEventsTable 
        bonsaiData.stimEventsTable = stimEventsTable;
    catch
        fprintf('Could not save stimEvents raw table in Bonsai.mat for %s...\n', stimName)
    end    
    
    % Save the extracted Bonsai data (overwrites file if it exists)
    save(saveFilePath, 'bonsaiData');
    save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end