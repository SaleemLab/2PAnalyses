function [responseFrameIdx, responseFrameRelTimes] = get2PFramesByTrialV2(onARDTimes, twoPData, preStimTime, postStimTime, offARDTimes, nExpected)
    % Extract frames by filtering timestamps around stimulus events and force
    % every trial to have the same number of frames (indices) by cropping or padding.
    %
    % Inputs:
    %   - onARDTimes (1D array): Arduino times for detected stimulus ON events
    %   - twoPData (struct): Suite2p and Bonsai two-photon data structure for all planes
    %   - preStimTime (float): Time before stimulus onset to include
    %   - postStimTime (float): Time after stimulus onset to include
    %   - offARDTimes (1D array): Optional, Arduino times for detected
    %     stimulus OFF events. Use with Stimuli that do not have a ITI
    %     (e.g., SpareNoise) 
    %   - nExpected (1D array): Optional, Desired number of frames per trial in the output
    %
    % Outputs:
    %   - responseFrameIdx: Matrix (nTrials x nFrames) of frame indices from twoPData
    %   - responseFrameRelTimesIdx: Matrix (nTrials x nFrames) of relative frame times (frameTime - stimulusTime)
    %
    % Aman and Sonali - Feb 2025

    % Use only on-times if off-times are not provided
    if nargin < 5 || isempty(offARDTimes)
        combinedTimes = onARDTimes;
    else
        combinedTimes = sort([onARDTimes; offARDTimes]);
    end

    % Get two-photon frame timestamps
    frameTimes = twoPData(1).TwoPFrameTime; 
    nTrials = length(combinedTimes);

    % If nExpected is not provided, determine it as the median number of frames
    % across all trials (using the actual frames found in the window).
    % Alternatively, could use fs but fs will not be uniformely sampled for
    % z-motion corrected datasets, unless interpolated beforehand? (@Aman?)
    if nargin < 6 || isempty(nExpected)
        framesPerTrial = zeros(nTrials,1);
        for iTrial = 1:nTrials
            startTime = combinedTimes(iTrial) - preStimTime;
            endTime   = combinedTimes(iTrial) + postStimTime;
            framesPerTrial(iTrial) = sum(frameTimes >= startTime & frameTimes <= endTime);
        end
        nExpected = round(median(framesPerTrial));
        if nExpected < 1
            error('No frames found in at least one trial window.');
        end
    end

    % Preallocate outputs with NaNs (in case some trials have fewer frames)
    responseFrameIdx = NaN(nTrials, nExpected);
    responseFrameRelTimes = NaN(nTrials, nExpected);

    for iTrial = 1:nTrials
        % Define the time window for this trial
        startTime = combinedTimes(iTrial) - preStimTime;
        endTime   = combinedTimes(iTrial) + postStimTime;
        
        % Get indices of frames within the window
        trialFrameIndices = find(frameTimes >= startTime & frameTimes <= endTime);
        % Relative times for these frames (relative to stimulus onset)
        trialRelTimes = frameTimes(trialFrameIndices) - combinedTimes(iTrial);

        % If there are more frames than expected, crop (e.g., take the first nExpected)
        if length(trialFrameIndices) >= nExpected
            responseFrameIdx(iTrial, :) = trialFrameIndices(1:nExpected);
            responseFrameRelTimes(iTrial, :) = trialRelTimes(1:nExpected);
        else
            % Otherwise, store what you have (the remaining columns will stay as NaN)
            nFound = length(trialFrameIndices);
            responseFrameIdx(iTrial, 1:nFound) = trialFrameIndices;
            responseFrameRelTimes(iTrial, 1:nFound) = trialRelTimes;
        end
    end
    
end
