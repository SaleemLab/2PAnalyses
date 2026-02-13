function [response, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, stimName,  postStimTime, preStimTime, useoffARDTimes)
% Extract frames by directly filtering timestamps based on a time window around stimulus events.
% If off transitions are set as the 'start' of the next stimulus (i.e.,
% no gray screen) -- include offARDTimes.
%
% @Aman :
% (0) Do we want to save this as a new .mat file called response?
% (1) Check if its necessary to loop through planes; 
% (2) Trials will have different length..fill them up with nans 
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
%         .responseFrameIdx         - Cell array of logical arrays indicating frames within the stimulus window.
%         .responseFrameRelTimes - Cell array of relative frame times.
%         .preStimTime - PreStimulus time used
%         .postStimTime - PostStimulus time used


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

if contains(stimName, 'SparseNoise', 'IgnoreCase',true)
    postStimTime = 0.7; 
end 

if contains(stimName, 'RFMapping', 'IgnoreCase',true)
    postStimTime = 3; 
    preStimTime = 0.5;
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
if exist(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData,'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file')
    processedTwoPData = load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData);
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
else
    error('Missing: TwoPData and/or bonsaiData files.');
end

% Create response.mat file and file path to sessionFileInfo
stimFileName = sprintf('%s_%s_Response_%s.mat', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
sessionFileInfo.stimFiles(iStim).Response = fullfile(sessionFileInfo.Directories.save_folder,stimFileName);


% Determine the combined event times based on the useoffARDTimes flag
if useoffARDTimes
    combinedStimARDTimes = sort([bonsaiData.onARDTimes; bonsaiData.offARDTimes]);
else
    combinedStimARDTimes = bonsaiData.onARDTimes;
end

% Process each plane in twoPData / Possibly not necessary..CHANGE; not
% iterating through planes anymore.. 

frameTimes = processedTwoPData.TwoPFrameTime;
% response = twoPData(thisPlane).planeName;

% Preallocate cell arrays for each trial in this plane
response.responseFrameIdx = cell(length(combinedStimARDTimes), 1);
response.responseFrameRelTimes = cell(length(combinedStimARDTimes), 1);

% Loop over each stimulus event
for iTrial = 1:length(combinedStimARDTimes)
    % Define the time window around the stimulus event
    %             startTimes = combinedTimes(iTrial) - preStimTime;
    %             endTimes   = combinedTimes(iTrial) + postStimTime;

    % Logical: true/1 for frames within the window
    frameIdxToAnalyse = (frameTimes >= (combinedStimARDTimes(iTrial) - preStimTime)) & (frameTimes <= (combinedStimARDTimes(iTrial) + postStimTime));
    response.responseFrameIdx{iTrial} = frameIdxToAnalyse;

    % Compute relative frame times with respect to the event time
    response.responseFrameRelTimes{iTrial} = frameTimes(frameIdxToAnalyse) - combinedStimARDTimes(iTrial);
    response.preStimTime = preStimTime;
    response.postStimTime = postStimTime;
end

save(sessionFileInfo.stimFiles(iStim).Response, 'response');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end





