function [responseFrameIdx, responseFrameRelTimesIdx] = get2PFramesByTrial(onARDTimes, twoPData, preStimTime, postStimTime, offARDTimes)
    % Extract frames by directly filtering timestamps based on a time window around stimulus events.
    % If off transitions are set as the 'start' of the next stimulus (i.e.,
    % no gray screen) -- include offARDTimes.
    %
    % Inputs:
    %   - onARDTimes: Arduino times for detected stimulus ON events
    %   - twoPData: Suite2p and Bonsai two-photon data structure for all planes (trimmed to have the same length)  
    %   - preStimTime: Time before stimulus onset to include
    %   - postStimTime: Time after stimulus onset to include
    %   - offARDTimes: (Optional) Arduino times for detected stimulus OFF events
    %
    % Outputs:
    %   - responseFrameIdx: Logical array indicating frames within stimulus periods
    %   - responseFrameRelTimesIdx: Relative timestamps of the selected frames
    %
    % Aman and Sonali - Feb 2025

    if nargin < 5 || isempty(offARDTimes)
        combinedTimes = onARDTimes; % Only use on-times if off-times are not provided
    else
        combinedTimes = sort([onARDTimes; offARDTimes]); % Merge and sort event times
    end

    frameTimes = twoPData(1).TwoPFrameTime; 

    responseFrameIdx = cell(length(combinedTimes), 1);
    responseFrameRelTimesIdx = cell(length(combinedTimes), 1);

    for iTrial = 1:length(combinedTimes)
        % Define the time window
        startTime = combinedTimes(iTrial) - preStimTime;
        endTime = combinedTimes(iTrial) + postStimTime;

        % Logical array: true for frames within the window
        responseFrameIdx{iTrial} = (frameTimes >= startTime) & (frameTimes <= endTime);

        % Extract relative frame times (only for selected frames)
        responseFrameRelTimesIdx{iTrial} = frameTimes(responseFrameIdx{iTrial}) - combinedTimes(iTrial);
    end    
end
