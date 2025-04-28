function [bonsaiData] = getTrialGroups(sessionFileInfo, stimName)
    % Groups trials based on stimulus type from Bonsai data.
    %
    % Extracts stimulus-related trials by identifying unique stimulus
    % types (`stimType`) recorded in Bonsai data. If `stimType` is unavailable, it assigns 
    % all trials to a general category using `stimName` as the group identifier.
    %
    % Inputs:
    %   sessionFileInfo (struct) - Structure containing file paths and session details. 
    %                              Used to locate the corresponding Bonsai data file.
    %   stimName (char) - Name of the stimulus for which trial groups need to be extracted.
    %
    % Outputs:
    %   stimGroups (struct array) - Structure with fields:
    %       .value        - Numeric identifier for each stimulus type (e.g., 0, 45, 90, etc.).
    %       .stimTypeName - Name of the stimulus category (uses `stimName` if not in Bonsai data).
    %       .trials       - Indices of trials corresponding to each stimulus type.
    %
    % Notes/Possible future fixes or improvements.. 
    %   - The function loads Bonsai-generated data if available.
    %   - If `stimType` does not exist in Bonsai data, it assigns all trials to a generic group
    %     and labels it using `stimName`.
    %   - If Bonsai data files are missing, the function throws an error.
    %
    % Aman and Sonali - Feb 2025
    
    % Locate the requested stimulus in sessionFileInfo
    isStim = strcmp(stimName, {sessionFileInfo.stimFiles.name}); 
    iStim = find(isStim, 1);  % Find first matching index
    
    % Check if the requested stimulus exists
    if isempty(iStim)
        error('Stimulus name "%s" not found in sessionFileInfo.', stimName);
    end
    
    % Check and load Bonsai data and response file
    if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
       exist(sessionFileInfo.stimFiles(iStim).Response, 'file')
        load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData'); % Load stimulus data
        load(sessionFileInfo.stimFiles(iStim).Response, 'response');     % Load response data
    else
        error('Missing BonsaiData and/or response files for stimulus "%s".', stimName);
    end


    % Check if Bonsai data contains the `stimType` field
    
    trialGroups = struct();
    if isfield(bonsaiData, 'stimType')
        uniqueStimTypes = unique(bonsaiData.stimType); % Get unique stimulus types
        
        % Group trials by each stimulus type
        for thisStimType = 1:length(uniqueStimTypes)
            stimValue = uniqueStimTypes(thisStimType);
            trialGroups(thisStimType).value = stimValue;  % e.g., 0, 45, 90, etc.
            trialGroups(thisStimType).stimTypeName = stimName; % Use stimName if stimTypeName is missing
            trialGroups(thisStimType).trials = find(bonsaiData.stimType == stimValue); % Get trial indices
        end
    else
        % If `stimType` is not available, assign all trials to a single group using `stimName`
        trialGroups(1).value = nan;  % Undefined group value
        trialGroups(1).stimTypeName = stimName; % Use stimName as identifier
        trialGroups(1).trials = 1:length(response(1).responseFrameIdx); % Include all trials
    end
    bonsaiData.trialGroups = trialGroups;
    save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData')
end
