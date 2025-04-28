function [response] = get2PFramesByTrialV4(sessionFileInfo, stimName,  postStimTime, preStimTime, useoffARDTimes)
% Extract frames by directly filtering timestamps based on a time window around stimulus events.
% If off transitions are set as the 'start' of the next stimulus (i.e.,
% no gray screen) -- include offARDTimes.
%
% @Aman :
% (0) Do we want to save this as a new .mat file called response?
% (1) Check if its necessary to loop through planes; if all grabs per
% (2) Trials will have different length
%
% Inputs:
%   - sessionFileInfo: Structure containing stimFiles and associated data file paths.
%   - stimName: Name of the stimulus to use for frame extraction.
%   - useoffARDTimes: Boolean flag; if true, use offARDTimes along with onARDTimes.
%   - preStimTime: Time before stimulus onset to include.
%   - postStimTime: Time after stimulus onset to include.
%
% Outputs:
%   - response: A struct array (one per plane) with fields:
%         responseFrameIdx         - Cell array of logical arrays indicating frames within the stimulus window.
%         responseFrameRelTimes - Cell array of relative frame times.

%
% Example:
% [response] = get2PFramesByTrialV3(sessionFileInfo, 'SparseNoise', true, 0, 0.7)
%
% Aman and Sonali - Feb 2025

if nargin<5
    useoffARDTimes = 'false';
end

if nargin<4
    preStimTime = 0;
end

% Locate the current stimulus in sessionFileInfo
isStim = false(1, length(sessionFileInfo.stimFiles));
for iStim = 1:length(sessionFileInfo.stimFiles)
    isStim(iStim) = strcmp(stimName, sessionFileInfo.stimFiles(iStim).name);
end
iStim = find(isStim, 1);  % take the first match
if isempty(iStim)
    error('Stimulus name not found in sessionFileInfo');
end

%% Check for existence of required files and load them
if exist(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file')
    load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData');
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
else
    error('Missing: TwoPData and/or bonsaiData files.');
end

% Create response.mat file and file path to sessionFileInfo
stimFileName = sprintf('%s_%s_Response_%s.mat', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
sessionFileInfo.stimFiles(iStim).Response = ...
    fullfile(sessionFileInfo.Directories.save_folder, stimFileName);


% Determine the combined event times based on the useoffARDTimes flag
if useoffARDTimes
    combinedStimARDTimes = sort([bonsaiData.onARDTimes; bonsaiData.offARDTimes]);
else
    combinedStimARDTimes = bonsaiData.onARDTimes;
end

% Process each plane in twoPData / Possibly not necessary..
for thisPlane = 1:length(twoPData)
    frameTimes = twoPData(thisPlane).TwoPFrameTime;
    response(thisPlane).planeName = twoPData(thisPlane).planeName;
    
    % Preallocate cell arrays for each trial in this plane
    response(thisPlane).responseFrameIdx = nan(length(combinedStimARDTimes), length(frameTimes));
    response(thisPlane).responseFrameRelTimes = nan(length(combinedStimARDTimes), length(frameTimes));

    % Loop over each stimulus event
    for iTrial = 1:length(combinedStimARDTimes)
        % Define the time window around the stimulus event
        %             startTimes = combinedTimes(iTrial) - preStimTime;
        %             endTimes   = combinedTimes(iTrial) + postStimTime;

        % Logical: true/1 for frames within the window
        frameIdxToAnalyse = (frameTimes >= (combinedStimARDTimes(iTrial) - preStimTime)) & (frameTimes <= (combinedStimARDTimes(iTrial) + postStimTime));
        response(thisPlane).responseFrameIdx(iTrial,:) = frameIdxToAnalyse;

        % Compute relative frame times with respect to the event time
        response(thisPlane).responseFrameRelTimes(iTrial,:) = frameTimes - combinedStimARDTimes(iTrial);
        %             response(thisPlane).responseFrameRealTimes{iTrial} = frameTimes(frameIdxToAnalyse);
    end
end
save(sessionFileInfo.stimFiles(iStim).Response, 'response');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

end
