function [response, sessionFileInfo] = extractVRAndPeripheralData(sessionFileInfo,  VRStimName, plotFlag)
%   Extracts wheel speed, virtual position, and lap-related info
%   from aligned Bonsai and peripheral data during 2P-VR Aman's Classical Corridor.
%   Handles lap classification (completed/aborted) and optionally plots lap timing.
%
% Inputs:
%   sessionFileInfo : struct
%       Contains file paths and metadata for the session. 
%
%   plotFlag : logical (optional)
%       If true, plots a visual summary of lap start/end times. Default: false.
%
% Outputs:
%   response : struct
%       Contains key behavioral and timing variables for downstream lap-by-lap analysis:
%       - wheelSpeed             : real-time wheel speed (cm/s)
%       - mouseVirtualPosition   : virtual track position (1–140 cm)
%       - trackIDFromMousePosition : track IDs inferred from mouse position (usually 1)
%       - mouseRecordedPosition  : raw position signal (-1141 to -1000)
%       - trackIDs               : track ID per lap (all 1s in classical VR corridor)
%       - lapCount               : unified lap index
%       - blockIDs               : cumulative index of track switches (e.g., block transitions; NA)
%       - trialType              : task/trial type ID per lap 
%                                  (0-NoTask; 1-Passive; 2-Hybrid)
%       - completedLaps          : indices of completed laps
%       - abortedLaps            : indices of aborted laps
%       - startTimeAll           : Bonsai start times for all laps
%       - endTimeAll             : parsed lap end time (based on position trace)
%       - completedStartTimes    : start times for completed laps
%       - completedEndTimes      : end times for completed laps
%
% Example usage:
%   response = extractVRAndPeripheralData(sessionFileInfo, true);
%
% Authors: Sonali & Aman (based on Diao/Masa code), March 2025

    if nargin < 3, plotFlag = false; end

    %% Load data files 
    for iStim = 1:length(sessionFileInfo.stimFiles) % Find VR Stimulus
        bonsaiData.isVRstim(iStim) = strcmp(VRStimName, sessionFileInfo.stimFiles(iStim).name);
    end
    iStim = find(bonsaiData.isVRstim==1); 
    
    if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
            exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file') && ...
            exist(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'file')
        load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
        load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');     
        load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'processedTwoPData');
    else
        error('PeriphralData and/or BonsaiData missing for VR recording');
    end
    
    %% Create response data file 
    stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_Response_' sessionFileInfo.stimFiles(iStim).name '.mat'];
    sessionFileInfo.stimFiles(iStim).Response = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

    %% 'Real' wheel position and real speed
    tickToCmConversion = 3.1415 * 20 / 1024;  % Wheel radius 20 cm, 1024 ticks per revolution
    displacement = [0; diff(peripheralData.Wheel.Value * tickToCmConversion)];
    
    % Handle unrealistic large changes (e.g., due to teleportation or resets)
    displacement(displacement < -100) = 0;  % Negative large jumps
    displacement(displacement > 100) = 0;   % Positive large jumps
    
    % Calculate speed (in cm/s)
    response.wheelSpeed = displacement ./ [0; diff(peripheralData.Wheel.sampleTimes)]; % Change to peripheralData.Wheel.ArduinoTime
    
    %% Virtual position and virtual speed
    mouseVirtualPosition = nan(1,length(bonsaiData.MousePos.Value));
    trackIDFromPosition = nan(1,length(bonsaiData.MousePos.Value));
    
    % Diao's track 1 excluding the contextual 
    % Convert raw mouse positions between -1141 and -1000 into virtual positions.
    % The conversion involves adding 1140 to the raw mouse position and taking the absolute value.
    mouseVirtualPosition(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -1000)) ...
        = abs(bonsaiData.MousePos.Value(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -1000))+1140);
    trackIDFromPosition(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -1000)) = 1; % Only one track (in Sonali's exps)
    mouseVirtualPosition(mouseVirtualPosition>140) = 140;
    
    response.mouseVirtualPosition = mouseVirtualPosition';
    response.trackIDFromMousePosition = trackIDFromPosition';
    response.mouseRecordedPosition = bonsaiData.MousePos.Value;

    
    %% Lap track Info
    % Save track ID as 1 for all the laps. 
    response.trackIDs = ones(1, length(bonsaiData.TrialInfo.StartTimeAll))';
    % LapCounts
    response.lapCount = (1:length(bonsaiData.TrialInfo.StartTimeAll))';  % Unified lap numbering
    
    % Block ID of each lap; same for all tracks
    blockTransition = [1; diff(response.trackIDs)];
    blockTransition(blockTransition~=0) = 1;
    response.blockIDs = cumsum(blockTransition);
    
    % Trial type for each lap
    if isfield(bonsaiData.TrialInfo, 'Trial_type')
        response.trialType = bonsaiData.TrialInfo.trialType; 
    else 
        % calling it 0 i.e., no task component; Masa - 1 is active only and
        % 2 is hybrid(?)
        response.trialType = zeros(1, length(bonsaiData.TrialInfo.StartTimeAll))';
    end
    

   %% Find completed and aborted laps 
    completedLaps = [];
    abortedLaps = [];
    
    lapStartTimeAll = bonsaiData.TrialInfo.StartTimeAll;
    trackIDs = response.trackIDs;
    EndTimeAll = NaN(length(trackIDs), 1); % Preallocate with NaNs for safety
    
    if ~isempty(lapStartTimeAll)
        x = response.mouseVirtualPosition;  % Virtual Position trace
        t = processedTwoPData.(processedTwoPData.resample2PTimeUsed); % Find the closest twoP time? Mouse position time? @Aman 
        startIdx = zeros(length(trackIDs), 1); % Index into time vector for each lap start
    
        for nlap = 1:length(trackIDs)
            % Find the time index closest to each lap start time
            [~, startIdx(nlap)] = min(abs(t - lapStartTimeAll(nlap)));
        end
    
        for nlap = 1:length(startIdx)
            % Extract position and time for current lap
            if nlap < length(startIdx)
                currentLapX = x(startIdx(nlap):startIdx(nlap+1));
                currentLapT = t(startIdx(nlap):startIdx(nlap+1));
            else
                currentLapX = x(startIdx(nlap):end);
                currentLapT = t(startIdx(nlap):end);
            end
    
            % Only proceed if there’s more than 1 non-NaN datapoint (avoid  missing data)
            if length(currentLapX) > 1 && sum(~isnan(currentLapX)) > 1
                onTrackX = currentLapX(~isnan(currentLapX));
                onTrackT = currentLapT(~isnan(currentLapX));
    
                % If large jumps in the initial part of the lap (likely due to lag), remove them
                if length(onTrackX) > 30
                    endFrame = 30;
    
                    % Positive jump
                    if ~isempty(find(diff(onTrackX(1:endFrame)) > 5, 1))
                        jumpIndex = find(diff(onTrackX(1:endFrame)) > 5, 1, 'last');
                        onTrackX(1:jumpIndex) = [];
                        onTrackT(1:jumpIndex) = [];
                    end
    
                    % Negative jump
                    if length(onTrackX) > 1 && ~isempty(find(diff(onTrackX(1:endFrame)) < -5, 1))
                        jumpIndex = find(diff(onTrackX(1:endFrame)) < -5, 1, 'last');
                        onTrackX(1:jumpIndex) = [];
                        onTrackT(1:jumpIndex) = [];
                    end
                end
    
                % If lap starts somewhere mid-track, align to the first time position = 0
                if sum(onTrackX == 0) > 0
                    startPosition = find(onTrackX == 0, 1);
                    if startPosition < length(onTrackX) - 10
                        onTrackX = onTrackX(startPosition:end);
                        onTrackT = onTrackT(startPosition:end);
                    end
                end
    
                % Sometimes final point jumps incorrectly — fix by matching to second-last point
                if length(onTrackX) >= 2 && onTrackX(end) ~= onTrackX(end-1)
                    onTrackX(end) = onTrackX(end-1);
                    onTrackT(end) = onTrackT(end-1);
                end
    
                % Get final position and time index
                [lastPosition, lastPositionIndex] = max(onTrackX);
    
                % If track end (140 cm) is reached almost instantly remove early part
                if lastPositionIndex * mean(diff(onTrackT)) < 0.1
                    onTrackX(1:lastPositionIndex) = [];
                    onTrackT(1:lastPositionIndex) = [];
    
                    if isempty(onTrackX)
                        fprintf('Lap %d aborted: only fast 140cm jump found.\n', nlap);
                        abortedLaps = [abortedLaps; nlap];
                        EndTimeAll(nlap) = NaN;
                        continue
                    end
    
                    [lastPosition, lastPositionIndex] = max(onTrackX);
                end
    
                % If the end of track was reached properly
                if lastPosition >= 139 % sometimes last lap ends before 140 cm
                    EndTimeAll(nlap) = onTrackT(lastPositionIndex); % End time when track completed
                    completedLaps = [completedLaps; nlap]; % Save lap number as completed
                else
                    % If end of track was not reached, use the last recorded time
                    EndTimeAll(nlap) = onTrackT(end); 
                    abortedLaps = [abortedLaps; nlap]; % Save lap number as aborted
                end
            else
                % Not enough data points to evaluate this lap
                fprintf('Lap %d aborted: insufficient valid position data.\n', nlap);
                abortedLaps = [abortedLaps; nlap];
                EndTimeAll(nlap) = NaN;
            end
        end
    end
    
    response.completedLaps = completedLaps;
    response.abortedLaps = abortedLaps;
    response.endTimeAll = EndTimeAll;
    response.startTimeAll = lapStartTimeAll;
    % Extract lap-wise variables for only completed laps
    response.completedStartTimes = response.startTimeAll(response.completedLaps);
    response.completedEndTimes   = response.endTimeAll(response.completedLaps);
    
    % Final check to ensure all laps were accounted for
    assert(length(completedLaps) + length(abortedLaps) == length(trackIDs), ...
        'This is to keep Sonali sane: Some laps were not classified into completed or aborted.');

    %% Sanity check plot: Lap start and end times across session
    if nargin < 2 || plotFlag
        figure('Name', 'Lap Start and End Times'); clf;
    
        startTimes = response.startTimeAll;
        endTimes = response.endTimeAll;
        nLaps = min(length(startTimes), length(endTimes));
        lapIDs = 1:nLaps;
    
        completedLaps = response.completedLaps;
        abortedLaps = response.abortedLaps;
    
        hold on;
    
        % --- Plot lap connectors ---
        for i = 1:nLaps
            if ismember(i, completedLaps)
                plot([startTimes(i), endTimes(i)], [lapIDs(i), lapIDs(i)], 'k-', 'LineWidth', 1); % black line
            elseif ismember(i, abortedLaps)
                plot([startTimes(i), endTimes(i)], [lapIDs(i), lapIDs(i)], 'r-', 'LineWidth', 1); % red line
            end
        end
    
        % --- Start and End markers ---
        plot(startTimes(completedLaps), completedLaps, 'go', 'MarkerFaceColor', 'k'); % completed starts
        plot(endTimes(completedLaps), completedLaps, 'ko', 'MarkerFaceColor', 'k');   % completed ends
    
        plot(startTimes(abortedLaps), abortedLaps, 'ro', 'MarkerFaceColor', 'r');     % aborted starts
        plot(endTimes(abortedLaps), abortedLaps, 'ro', 'MarkerFaceColor', 'r');       % aborted ends
    
        % --- Labels and axis formatting ---
        xlabel('Time (s)');
        ylabel('Lap #');
        title('Lap Start and End Timeline');
        legend({'Aborted Lap', 'Completed Lap'}, 'Location', 'southeast');

       
        summaryText = sprintf('Total laps: %d\nCompleted: %d\nAborted: %d', ...
            length(trackIDs), length(completedLaps), length(abortedLaps));
        
        % Add it to the upper-right corner of the axes
        xPos = max(endTimes) + 1;
        yPos = nLaps;
        
        text(xPos, yPos, summaryText, ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 10, ...
            'FontWeight', 'bold');

    end

    
%% Saving 
save(sessionFileInfo.stimFiles(iStim).Response, 'response');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end
