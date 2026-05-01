%% temp test
function [response, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, stimName, excludeBadFrames, preStimTime, postStimTime)
%% Handle Input Defaults and Stimulus Overrides
if nargin < 3, excludeBadFrames = true; end % Default to true to clean data
if nargin < 4, preStimTime = 2; end
if nargin < 5, postStimTime = 4; end

% (Your existing stimulus overrides logic here...)
if contains(stimName, 'RFMapping', 'IgnoreCase',true)
    postStimTime = 4; preStimTime = 2; 
end

if contains(stimName, 'DotMotion_SpeedTuning_Contrast')
    postStimTime = 6; preStimTime = 2;
end 

if contains(stimName, 'DotMotion_RFMapping')
    postStimTime = 6; preStimTime = 2;
end 
%% Load Data
isStim = strcmp(stimName, {sessionFileInfo.stimFiles.name});
iStim = find(isStim, 1);
processedData = load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData);
load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');

% Set output path
stimFileName = sprintf('%s_%s_Response_%s.mat', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, stimName);
sessionFileInfo.stimFiles(iStim).Response = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

%% Setup Timebase
frameTimes = processedData.TwoPFrameTime(:)'; % Ensure row vector
combinedStimARDTimes = bonsaiData.onARDTimes;
nTrials = length(combinedStimARDTimes);

% We create a fixed vector from -pre to +post
Fs = 60; 
perfectTime = -preStimTime : (1/Fs) : postStimTime;
nFramesInPerfect = length(perfectTime);

%% Handle Bad Frames
isBadFrameGlobal = false(size(frameTimes));
if excludeBadFrames && isfield(processedData, 'badFrames')
    % Suite2p badframes are usually saved as a cell array per plane
    % We assume they were already interpolated to the 60Hz grid in Fun 2
    tempMask = processedData.badFrames{1}; 
    isBadFrameGlobal(1:length(tempMask)) = logical(tempMask);
end

%% Process Trials
response.alignedTimes = perfectTime; % The "Ground Truth" x-axis
response.badTrialMask = false(nTrials, 1);
response.responseFrameIdx = cell(nTrials, 1);
response.responseFrameRelTimes = cell(nTrials, 1);

for iTrial = 1:nTrials
    % Define the window for this specific trial
    t_start = combinedStimARDTimes(iTrial) - preStimTime;
    t_end   = combinedStimARDTimes(iTrial) + postStimTime;

    % Find frames within the window
    timeWindowMask = (frameTimes >= t_start - 0.1) & (frameTimes <= t_end + 0.1); 

    % --- THE QUALITY CHECK ---
    totalFramesInWindow = sum(timeWindowMask);
    badFramesInWindow = sum(timeWindowMask & isBadFrameGlobal);

    % If > 50% of frames are bad, mark the whole trial
    if totalFramesInWindow > 0 && ((totalFramesInWindow - badFramesInWindow) / totalFramesInWindow < 0.5)
        response.badTrialMask(iTrial) = true;
    end

    % Store the mask for later indexing
    response.responseFrameIdx{iTrial} = timeWindowMask;

    % --- THE ALIGNMENT FIX ---
    % Calculate exact relative times for the frames we found
    relTimes = frameTimes(timeWindowMask) - combinedStimARDTimes(iTrial);

    if excludeBadFrames
        % Set bad frames to NaN so they are ignored during later interpolation
        localBadMask = isBadFrameGlobal(timeWindowMask);
        relTimes(localBadMask) = NaN;
    end

    response.responseFrameRelTimes{iTrial} = relTimes;
end

%% Save
response.preStimTime = preStimTime;
response.postStimTime = postStimTime;
save(sessionFileInfo.stimFiles(iStim).Response, 'response');
fprintf('Done. %d trials flagged. Alignment grid set to 60Hz.\n', sum(response.badTrialMask));
end


