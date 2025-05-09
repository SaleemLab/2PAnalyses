% process_all_mice_sessions.m
% Script to process all mice and sessions

clear; clc;

% Define mouse and session information 
mouseInfo = {
    'M25040', {'20250506'}, {'M25040_RFMapping_20250506_00001'};
%     'M25038', {'20250423A','20250423B', '20250423C', '20250423D', '20250426A', '20250426B' }, {'M25038_RFMapping_20250423_00001', 'M25038_RFMapping_20250423_00003','M25038_RFMapping_20250423_00004' ,'M25038_RFMapping_20250423_00005', 'M25038_RFMapping_20250426_00001','M25038_RFMapping_20250426_00002'};
};

% Processing parameters 
interpRate = 60;        % Hz
preStimTime = 0.5;     % seconds
planeNums = 8;      % Total channels in recording
channelsSaved = 2;       % Channel to align frame times to
method = 2;            % Method for PSTH extraction
frameRate = 7.28;       % Degrees/bin or other relevant stimulus unit
applyNeuropilCorrection = true;
calculateDFF = true;
pdthreshold = 8;
postStimTime = 3; 

% === Loop through each mouse ===
for thisMouse = 1:size(mouseInfo,1)
    
    mousenumber = mouseInfo{thisMouse,1};
    sessionNames = mouseInfo{thisMouse,2};
    stimlists = mouseInfo{thisMouse,3};
    fprintf('Processing mouse: %s\n', mousenumber);
    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        stimName = stimlists{thisSession};
        fprintf('  Session: %s | Stimulus: %s\n', sessionName, stimName);
        try
            % Get session file paths
            sessionFileInfo = get2PsessionFilePaths(mousenumber, sessionName, {stimName});
            % Load metadata
            sessionFileInfo = get2PMetadata(sessionFileInfo);
            % Get frame times for two channels
            sessionFileInfo = get2PFrameTimes_TwoChannels(sessionFileInfo, planeNums, channelsSaved);
            % Process peripheral files (wheel, position etc)
            sessionFileInfo = processPeripheralFiles(sessionFileInfo);
            % Merge Bonsai and Suite2p files
            sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo);
            % Get stimulus events from Bonsai file
            [bonsaiData, sessionFileInfo] = getTuningStimEventsBonsaiFile(sessionFileInfo, stimName, 'StimulusParams');
            % Get stimulus times
            [bonsaiData, sessionFileInfo] = getStimTimes(sessionFileInfo, stimName, pdthreshold);
            % Get 2P frames by trial
            [response, sessionFileInfo] = get2PFramesByTrialV3(sessionFileInfo, stimName, postStimTime, preStimTime);
            % Group trials
            [bonsaiData, sessionFileInfo] = getTrialGroupsV2(sessionFileInfo, stimName);
            % PSTHs
            [response, sessionFileInfo] = getTrialResponsePSTHsV4(sessionFileInfo, stimName, method, interpRate, frameRate, applyNeuropilCorrection, calculateDFF);
            %
            plotRFGrid_byPosition_ROIs(sessionFileInfo, stimName);
            %
            plotAllROIPSTHsByPosition(sessionFileInfo, response);
            % 
            
            %
            fprintf('    Done!\n'); 
        catch ME
            warning('    Error processing %s %s: %s', mousenumber, sessionName, ME.message);
        end
        
    end
end

disp('All mice and sessions processed and figures saved in session-wise analysis folders');
