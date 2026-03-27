%% temp test
function [response, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, stimName, excludeBadFrames, preStimTime, postStimTime, useoffARDTimes)
% Extract frames by directly filtering timestamps.
% If excludeBadFrames is true, relative times for artifact frames are set to NaN
% to preserve the temporal grid of the trial.
%
% Aman and Sonali - Feb 2025

%% Handle Input Defaults
if nargin < 3, excludeBadFrames = false; end
if nargin < 5, useoffARDTimes = 'false'; end
if nargin < 4, preStimTime = 0; end

% Stimulus specific overrides
if contains(stimName, 'SparseNoiseTexture', 'IgnoreCase',true)
    postStimTime = 3; useoffARDTimes = 'true'; 
elseif contains(stimName, 'RFMapping', 'IgnoreCase',true)
    postStimTime = 4; preStimTime = 2; useoffARDTimes = 'false';
elseif contains(stimName, 'DotMotion_SpeedTuning','IgnoreCase',true)
    postStimTime = 4; preStimTime = 2; useoffARDTimes = 'false'; 
elseif contains(stimName, 'DirTuning','IgnoreCase',true)
    postStimTime = 4; preStimTime = 2; useoffARDTimes = 'false'; 
end

% Locate the current stimulus
isStim = strcmp(stimName, {sessionFileInfo.stimFiles.name});
iStim = find(isStim, 1);
if isempty(iStim)
    error('Stimulus name not found in sessionFileInfo');
end

%% Load flattened data
if exist(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData,'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file')
    processedData = load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData);
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
else
    error('Missing files.');
end

stimFileName = sprintf('%s_%s_Response_%s.mat', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
sessionFileInfo.stimFiles(iStim).Response = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

if strcmpi(useoffARDTimes, 'true') || (islogical(useoffARDTimes) && useoffARDTimes)
    combinedStimARDTimes = sort([bonsaiData.onARDTimes; bonsaiData.offARDTimes]);
else
    combinedStimARDTimes = bonsaiData.onARDTimes;
end

frameTimes = processedData.TwoPFrameTime;
if size(frameTimes, 1) > size(frameTimes, 2), frameTimes = frameTimes'; end

%% Handle bad frames 
isBadFrameGlobal = false(size(frameTimes));
if excludeBadFrames && isfield(processedData, 'badFrames') && ~isempty(processedData.badFrames)
    tempMask = processedData.badFrames{1};
    if size(tempMask, 1) > size(tempMask, 2)
        tempMask = tempMask'; 
    end
    len = min(length(isBadFrameGlobal), length(tempMask));
    isBadFrameGlobal(1:len) = logical(tempMask(1:len));
end

%% Process Trials
nTrials = length(combinedStimARDTimes);
response.responseFrameIdx = cell(nTrials, 1);
response.responseFrameRelTimes = cell(nTrials, 1);
response.badTrialMask = false(nTrials, 1); 

for iTrial = 1:nTrials
    t_start = combinedStimARDTimes(iTrial) - preStimTime;
    t_end   = combinedStimARDTimes(iTrial) + postStimTime;
    
    % Identify all frames in the window 
    timeWindowMask = (frameTimes >= t_start) & (frameTimes <= t_end);
    
    % Check excluded frames to display
    totalFramesInWindow = sum(timeWindowMask);
    badFramesInWindow = sum(timeWindowMask & isBadFrameGlobal);
    
    if totalFramesInWindow > 0 && ( (totalFramesInWindow - badFramesInWindow) / totalFramesInWindow < 0.5 )
        response.badTrialMask(iTrial) = true;
    end
    
   
    % We store the FULL window mask so the indexing into F/spks remains consistent
    response.responseFrameIdx{iTrial} = timeWindowMask;
    
    % Relative Times with NaNs for bad frames
    if any(timeWindowMask)
        relTimes = frameTimes(timeWindowMask) - combinedStimARDTimes(iTrial);
        
        if excludeBadFrames
            % 
            localBadMask = isBadFrameGlobal(timeWindowMask);
            % Set bad frame relative times to NaN
            relTimes(localBadMask) = NaN;
        end
        
        response.responseFrameRelTimes{iTrial} = relTimes;
    else
        response.responseFrameRelTimes{iTrial} = [];
    end
end

response.preStimTime = preStimTime;
response.postStimTime = postStimTime;
response.excludeBadFramesUsed = excludeBadFrames;

%% Save
save(sessionFileInfo.stimFiles(iStim).Response, 'response');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp(['Done. ' num2str(sum(response.badTrialMask)) ' trials flagged for data loss.']);
end