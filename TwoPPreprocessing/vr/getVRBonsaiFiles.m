function [bonsaiData, sessionFileInfo] = getVRBonsaiFiles(sessionFileInfo, VRStimName)
% Get the bonsai VR mouse position and trial info tables and saves in a .mat file
% Aman and Sonali - Jan 2025
% Modified Jan 2026 to integrate UCL-open VRCorridor and landmark reformatting

%% Load VR stimulus 
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName, sessionFileInfo.stimFiles(iStim).name);
end

iStim = find(bonsaiData.isVRstim == 1);
if isempty(iStim)
    error('VRStimName %s not found in sessionFileInfo.', VRStimName);
end

% Construct file path for saving Bonsai data
stimFileName = sprintf('%s_%s_BonsaiData_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
sessionFileInfo.stimFiles(iStim).BonsaiData = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

%% Look for all files 
mouseposFilePath       = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'MousePos');
trialInfoFilePath      = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'TrialInfoLog');
VrPositionFilePath     = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'VrPosition'); 
LandmarkChoiceFilePath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'LandmarkChoice');

%% VR mouse position 
if exist(mouseposFilePath, 'file')
    
    mousePosTable = readtable(mouseposFilePath);
    if ismember('ArduinoTime', mousePosTable.Properties.VariableNames)
        rowsToRemove = (mousePosTable.ArduinoTime == 0) | (mousePosTable.MousePosition == 0);
        mousePosTable(rowsToRemove, :) = [];   
        
        [~, keep_idx, ~] = unique(mousePosTable.ArduinoTime);
        
        mousepos.rawArduinoTime       = mousePosTable.ArduinoTime(keep_idx) ./ 1000;
        mousepos.rawBonsaiTime        = mousePosTable.BonsaiTime(keep_idx);
        mousepos.rawValue             = mousePosTable.MousePosition(keep_idx);
        mousepos.rawLastSyncPulseTime = mousePosTable.LastSyncPulseTime(keep_idx) ./ 1000;
        mousepos.rawRenderFrameCount  = mousePosTable.RenderFrameCount(keep_idx); 
        
        bonsaiData.MousePos = mousepos;
    end
elseif exist(VrPositionFilePath, 'file')
    % UCL-Open format
    vrPositionTable = readtable(VrPositionFilePath, 'VariableNamingRule', 'preserve');
    if ismember('Seconds', vrPositionTable.Properties.VariableNames)
        rows_to_remove = (vrPositionTable.Seconds == 0); % Only using seconds here @Aman 
        vrPositionTable(rows_to_remove, :) = [];
        
        [~, keep_idx, ~] = unique(vrPositionTable.Seconds);
        mousepos.rawArduinoTime   = vrPositionTable.Seconds(keep_idx);
        mousepos.rawValue         = vrPositionTable.("Value.Z")(keep_idx); % Primary VR depth
        mousepos.rawValueX        = vrPositionTable.("Value.X")(keep_idx);
        mousepos.rawValueY        = vrPositionTable.("Value.Y")(keep_idx);
        mousepos.rawLength        = vrPositionTable.("Value.Length")(keep_idx);
        mousepos.rawLengthFast    = vrPositionTable.("Value.LengthFast")(keep_idx);
        
        % Removed and redundant here 
        % if ismember('Value.LastSyncPulseTime', vrPositionTable.Properties.VariableNames)
        %     mousepos.rawLastSyncPulseTime = vrPositionTable.("Value.LastSyncPulseTime")(keep_idx) ./ 1000;
        % end
        % 
        bonsaiData.MousePos = mousepos;
    end
else
    bonsaiData.MousePos = [];
    warning('CRITICAL: Mouse/VR position missing for session: %s', VRStimName);
end    

%% Trial/Landmark info
if exist(trialInfoFilePath, 'file')
    
    trialInfoTable = readtable(trialInfoFilePath);
    trialInfo.rawArduinoTime       = trialInfoTable.ArduinoTime ./ 1000;
    trialInfo.rawBonsaiTime        = trialInfoTable.BonsaiTime;
    trialInfo.rawLapCount          = trialInfoTable.LapCount;
    trialInfo.rawLastSyncPulseTime = trialInfoTable.LastSyncPulseTime ./ 1000;
    trialInfo.rawRenderFrameCount  = trialInfoTable.RenderFrameCount; 
    
    bonsaiData.TrialInfo = trialInfo;
    
elseif exist(LandmarkChoiceFilePath, 'file')
    % UCL-Open Format with reformatting @Aman 
    LandmarkChoiceTable = readtable(LandmarkChoiceFilePath, 'VariableNamingRule', 'preserve');
    
    % Reformat landmarks (per row) into trials (one row per unique timestamp)
    reformattedTrialTable = reformatLandmarkTableToTrialTable(LandmarkChoiceTable);
    
    trialInfo.rawArduinoTime        = reformattedTrialTable.Seconds;
    trialInfo.rawLapCount           = reformattedTrialTable.Lap; % New unique counter
    trialInfo.LandmarkNames         = reformattedTrialTable.LandmarkNames;
    trialInfo.LandmarkPositions     = reformattedTrialTable.Positions;
    trialInfo.LandmarkSizes         = reformattedTrialTable.Sizes;
    trialInfo.LandmarkCenterOffsets = reformattedTrialTable.CenterOffsets;
    trialInfo.LandmarkRewardValence = reformattedTrialTable.RewardValence;
    trialInfo.NumLandmarks          = reformattedTrialTable.NumLandmarks;
    
    bonsaiData.TrialInfo = trialInfo;
else 
    bonsaiData.TrialInfo = [];
    warning('Trial/Landmark info missing for stimulus: %s', VRStimName);
end 

%% Save 
save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
fprintf('Bonsai VR data saved for: %s\n', VRStimName);

end

%% Sub-function for landmark reformatting @Could move this further up the pipeline TODO CHECK 
function trialTable = reformatLandmarkTableToTrialTable(landmarkTable)
    % Groups landmark entries by unique timestamp to define discrete trials/laps
    [uniqueTimestamps, ~, trialIdx] = unique(landmarkTable.Seconds);
    numTrials = length(uniqueTimestamps);
    
    % Pre-allocate output
    trialTable = table('Size', [numTrials, 8], ...
        'VariableTypes', {'double', 'double', 'cell', 'cell', 'cell', 'cell', 'cell', 'double'}, ...
        'VariableNames', {'Seconds', 'Lap', 'LandmarkNames', 'Positions', 'Sizes', 'CenterOffsets', 'RewardValence', 'NumLandmarks'});

    for thisNum = 1:numTrials
        thisTrial = landmarkTable(trialIdx == thisNum, :);
        
        trialTable.Seconds(thisNum)      = uniqueTimestamps(thisNum);
        trialTable.Lap(thisNum)          = thisNum; % Lap number based on unique occurrence
        trialTable.LandmarkNames{thisNum}    = thisTrial.("Value.Value");
        trialTable.Positions{thisNum}        = thisTrial.("Value.Value.Position");
        trialTable.Sizes{thisNum}            = thisTrial.("Value.Value.Size");
        trialTable.CenterOffsets{thisNum}    = thisTrial.("Value.Value.CenterOffset");
        trialTable.RewardValence{thisNum}    = thisTrial.("Value.Value.RewardValence");
        trialTable.NumLandmarks(thisNum)     = height(thisTrial);
    end
end