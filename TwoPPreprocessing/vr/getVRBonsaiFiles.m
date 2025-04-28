function [bonsaiData, sessionFileInfo] = getVRBonsaiFiles(sessionFileInfo, VRStimName)
% Get the bonsai VR mouse position and trial info tables and saves in a 
% Usage: [bonsaiData] = getVRBonsaiFiles(sessionFileInfo)
% Output: bonsaiData.mousePos; bonsaiData.trialInfo; bonsaiData.isVRStim
% Aman and Sonali - Jan 2025
% Change to include multi-VR stim in one session 

%% Load VR stimulus 
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName,sessionFileInfo.stimFiles(iStim).name);
end

% numVRstim = sum(VRInfo.isVRstim);
% for iVR = 1:numVRstim
iStim = find(bonsaiData.isVRstim==1);

% Construct file path for saving Bonsai data
stimFileName = sprintf('%s_%s_BonsaiData_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
% Save filepath to sessionFileInfo     
sessionFileInfo.stimFiles(iStim).BonsaiData = ...
    fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

%% VR mouse position 
mouseposFilePath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'MousePos');
mousePosTable = readtable(mouseposFilePath);
if ismember('ArduinoTime', mousePosTable.Properties.VariableNames) && ismember('MousePosition', mousePosTable.Properties.VariableNames)
    rowsToRemove = (mousePosTable.ArduinoTime == 0) | (mousePosTable.MousePosition == 0);
    mousePosTable(rowsToRemove, :) = [];   
end

% Remove repeated time measures
[~,keep_idx,~] = unique(mousePosTable.ArduinoTime);

mousepos.rawArduinoTime  = mousePosTable.ArduinoTime(keep_idx)./1000;
mousepos.rawBonsaiTime   = mousePosTable.BonsaiTime(keep_idx);
mousepos.rawValue  = mousePosTable.MousePosition(keep_idx);
mousepos.rawLastSyncPulseTime = mousePosTable.LastSyncPulseTime(keep_idx);
mousepos.rawRenderFrameCount = mousePosTable.RenderFrameCount(keep_idx); 

bonsaiData.MousePos = mousepos;
%% Load and save all the trial info variables from the bonsai files into a mat file. 
trialInfoFilePath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths,'TrialInfoLog');
trialInfoTable = readtable(trialInfoFilePath);

% Discrete data associated with trial start times 
trialInfo.rawArduinoTime  = trialInfoTable.ArduinoTime./1000;
trialInfo.rawBonsaiTime   = trialInfoTable.BonsaiTime;
trialInfo.rawLapCount  = trialInfoTable.LapCount;
trialInfo.rawLastSyncPulseTime = trialInfoTable.LastSyncPulseTime;
trialInfo.rawRenderFrameCount = trialInfoTable.RenderFrameCount; 
bonsaiData.TrialInfo = trialInfo;

%% Save filepath in SessionFileInfo and the new Bonsai Data file
% trial_path = findFile(sessionFileInfo.stimbonsai_filepaths.VRCorr, 'TrialInfoLog');
save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData')
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');