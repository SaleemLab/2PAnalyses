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



%% Defualt parameters 
if nargin < 3, plotFlag = false; end

%% Load data files
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName, sessionFileInfo.stimFiles(iStim).name);
end
iStim = find(bonsaiData.isVRstim==1);

if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
   exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file') && ...
   exist(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'file')

    % Load peripheral and bonsai structs
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');

    % "Map" string to see which timebase was used
    load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'resample2PTimeUsed');

    % This avoids loading massive F or spks matrices
    timeData = load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, resample2PTimeUsed);
    t = timeData.(resample2PTimeUsed); 

else
    error('Peripheral, Bonsai, or Processed 2P Data missing for VR recording');
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
if isfield(peripheralData, 'Pupil')
    response.pupilArea = peripheralData.Pupil.Value.Area;
end 

%% Virtual position and virtual speed
mouseVirtualPosition = nan(1,length(bonsaiData.MousePos.Value));
trackIDFromPosition = nan(1,length(bonsaiData.MousePos.Value));

% Diao's track 1 excluding the contextual
% Convert raw mouse positions between -1141 and -1000 into virtual positions.
% The conversion involves adding 1140 to the raw mouse position and taking the absolute value.
% adds +1140 to these values
% -1140 becomes 0
% -1141 becomes -1
% -1000 becomes 140 

if contains(VRStimName, 'VRCorr', 'IgnoreCase', true)
    % condition one 
    mouseVirtualPosition(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -990)) ...
        = floor(abs(bonsaiData.MousePos.Value(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -990))+1140)) + 1;

    trackIDFromPosition(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -990)) = 1; % Only one track (in Sonali's exps)


    mouseVirtualPosition(mouseVirtualPosition > 140) = 140;

    response.mouseVirtualPosition = mouseVirtualPosition';


elseif contains(VRStimName, 'Baseline') || contains(VRStimName, 'LandManipCorridor') % this is temporary; might move to a new function @aman 
    % condition two 
    % % bonsaiData.MousePos.Value(bonsaiData.MousePos.Value < 0) = 0; % Positions less than 0 are all assigned as 0 
    % % mouseVirtualPosition = bonsaiData.MousePos.Value;
    % % 
    % % 
    % % % mouseVirtualPosition(mouseVirtualPosition > 140) = 140; % Change in future.. [currently the location goes beyond 140 because end of track is gray screen and the animal cannot see past this..] 
    % % mouseVirtualPosition(mouseVirtualPosition > 200) = 200; % The last position in vr is 200 but it seems to wobbly between 200.1xx and can also drop down to 199? 
    % % % mouseVirtualPosition(mouseVirtualPosition < 1)   = 1;
    % % 
    % % trackIDFromPosition(:) = 2; %@AMAN 
    % % response.mouseVirtualPosition = mouseVirtualPosition;

    rawPos = bonsaiData.MousePos.Value;
    % ensure no values are below 0 or above 200 before we apply the offset.
    % cleaning up 'wobble' before shifting 
    rawPos(rawPos < 0) = 0;
    rawPos(rawPos > 200) = 200;

    % original 0-199 (Track) is now 1-200
    % original 200 (ITI) is now 201
    mouseVirtualPosition = rawPos; %+ 1;
    if contains(VRStimName, 'Baseline')
        trackIDFromPosition(:) = 2;
    elseif contains(VRStimName, 'LandManipCorridor')
        trackIDFromPosition(:) = 3; 
    end 

    response.mouseVirtualPosition = mouseVirtualPosition;
end


% removed because I dont think this makes sense anymore; all coordintes
% occupy the same virtual positions 1 to 200 unlike masa's previous
% version. 
%response.trackIDFromMousePosition = trackIDFromPosition';
%response.mouseRecordedPosition = bonsaiData.MousePos.Value;

%% Lap track Info TODO: change to include block structure 
% Save track ID as 1 for all the laps.
%response.trackIDs = ones(1, length(bonsaiData.TrialInfo.StartTimeAll))';
% LapCounts
response.lapCountAll = (1:length(bonsaiData.TrialInfo.StartTimeAll))';  % Unified lap numbering

% Block ID of each lap; same for all tracks
% blockTransition = [1; diff(response.trackIDs)];
% blockTransition(blockTransition~=0) = 1;
% response.blockIDs = cumsum(blockTransition);

% Trial type for each lap
if isfield(bonsaiData.TrialInfo, 'Trial_type')
    response.trialTypeAll = bonsaiData.TrialInfo.trialType;
else
    % calling it 0 i.e., no task component; Masa - 1 is active only and
    % 2 is hybrid(?)
    response.trialTypeAll = zeros(1, length(bonsaiData.TrialInfo.StartTimeAll))';
end


%% Find completed and aborted laps
completedLaps_AbsoluteIdx = [];
abortedLaps_AbsoluteIdx = [];

lapStartTimeAll = bonsaiData.TrialInfo.StartTimeAll;
% trackIDs = response.trackIDAlls;
EndTimeAll = NaN(length(lapStartTimeAll), 1); % Preallocate with NaNs for safety

if ~isempty(lapStartTimeAll)
    x = response.mouseVirtualPosition;  % Virtual Position trace
    % t = processedTwoPData.(processedTwoPData.resample2PTimeUsed); % t was
    % loaded in the load data sections above.. 
    startIdx = zeros(length(response.lapCountAll), 1); % Index into time vector for each lap start

    for nlap = 1:length(response.lapCountAll)
        % Find the time index closest to each lap start time
        [~, startIdx(nlap)] = min(abs(t - lapStartTimeAll(nlap)));
    end

    for nlap = 1:length(startIdx)
        % Extract position and time for current lap; changed this 
        if nlap < length(startIdx)
            currentLapX = x(startIdx(nlap):startIdx(nlap+1)); % x(startIdx(nlap):startIdx(nlap+1));
            currentLapT = t(startIdx(nlap):startIdx(nlap+1)); % t(startIdx(nlap):startIdx(nlap+1));
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
            % For when the animal moves very fast and the position starts
            % at 5cm for example. 
            if sum(onTrackX == 0) > 0 %  sum(onTrackX == 0) > 0
                startPosition = find(onTrackX == 0, 1); %  find(onTrackX == 0, 1);
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
                    abortedLaps_AbsoluteIdx = [abortedLaps_AbsoluteIdx; nlap];
                    EndTimeAll(nlap) = NaN;
                    continue
                end

                [lastPosition, lastPositionIndex] = max(onTrackX);
            end

            % If the end of track was reached properly
            % if lastPosition >= 139
            if contains(VRStimName, 'Baseline') || contains(VRStimName, 'LandManipCorridor')
                %track ends at 199 and ITI is 200, 199.5 is a safe "finish line"
                condition = lastPosition >= 199.5;
            elseif contains(VRStimName, 'VRCorr')
                condition = lastPosition >= 139.5;
            end
            % lastPosition >= 199
            if condition % sometimes last lap ends before 140 cm; eg 139.99 
                EndTimeAll(nlap) = onTrackT(lastPositionIndex); % End time when track completed
                completedLaps_AbsoluteIdx = [completedLaps_AbsoluteIdx; nlap]; % Save lap number as completed
            else
                % If end of track was not reached, use the last recorded time
                EndTimeAll(nlap) = onTrackT(end);
                abortedLaps_AbsoluteIdx = [abortedLaps_AbsoluteIdx; nlap]; % Save lap number as aborted
            end
        else
            % Not enough data points to evaluate this lap
            fprintf('Lap %d aborted: insufficient valid position data.\n', nlap);
            abortedLaps_AbsoluteIdx = [abortedLaps_AbsoluteIdx; nlap];
            EndTimeAll(nlap) = NaN;
        end
    end
end

response.completedLaps_AbsoluteIdx = completedLaps_AbsoluteIdx;
response.abortedLaps_AbsoluteIdx = abortedLaps_AbsoluteIdx;
response.endTimeAll = EndTimeAll;
response.startTimeAll = lapStartTimeAll;
% Extract lap-wise variables for only completed laps
response.completedStartTimes = response.startTimeAll(completedLaps_AbsoluteIdx);
response.completedEndTimes   = response.endTimeAll(completedLaps_AbsoluteIdx);
response.stimName = VRStimName;
% Final check to ensure all laps were accounted for
assert(length(completedLaps_AbsoluteIdx) + length(abortedLaps_AbsoluteIdx) == length(response.lapCountAll), ...
    'This is to keep Sonali sane: Some laps were not classified into completed or aborted. If you see this something is work and requires debugging..');

%% Save and Filter Landmark Information (Baseline/LandManip only)
if contains(VRStimName, 'BaselineCorridor') || contains(VRStimName, 'LandManipCorridor')
    % Renaming with 'completed' prefix to distinguish from raw Bonsai data
    response.completedLandmarkNames         = bonsaiData.TrialInfo.LandmarkNames(completedLaps_AbsoluteIdx);
    response.completedLandmarkPositions     = bonsaiData.TrialInfo.LandmarkPositions(completedLaps_AbsoluteIdx);
    response.completedLandmarkSizes         = bonsaiData.TrialInfo.LandmarkSizes(completedLaps_AbsoluteIdx);
    response.completedLandmarkCenterOffsets = bonsaiData.TrialInfo.LandmarkCenterOffsets(completedLaps_AbsoluteIdx);
    response.completedLandmarkRewardValence = bonsaiData.TrialInfo.LandmarkRewardValence(completedLaps_AbsoluteIdx);
    response.completedNumLandmarks = bonsaiData.TrialInfo.NumLandmarks(completedLaps_AbsoluteIdx);
    % included movementvisual gain here 
    response.movementVisualGain = bonsaiData.movementVisualGain;
end

%%
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
save(sessionFileInfo.stimFiles(iStim).Response,'-struct', 'response', '-v7.3'); %using -v7.3 as this is first save
disp('Saved Response');
%% Sanity check plot: Lap start and end times across session
if nargin < 2 || plotFlag
    figure('Name', 'Lap Start and End Times'); clf;

    startTimes = response.startTimeAll;
    endTimes = response.endTimeAll;
    nLaps = min(length(startTimes), length(endTimes));
    lapIDs = 1:nLaps;

    completedLaps_AbsoluteIdx = response.completedLaps_AbsoluteIdx;
    abortedLaps_AbsoluteIdx = response.abortedLaps_AbsoluteIdx;

    hold on;

    % --- Plot lap connectors ---
    for i = 1:nLaps
        if ismember(i, completedLaps_AbsoluteIdx)
            plot([startTimes(i), endTimes(i)], [lapIDs(i), lapIDs(i)], 'k-', 'LineWidth', 1); % black line
        elseif ismember(i, abortedLaps_AbsoluteIdx)
            plot([startTimes(i), endTimes(i)], [lapIDs(i), lapIDs(i)], 'r-', 'LineWidth', 1); % red line
        end
    end

    % --- Start and End markers ---
    plot(startTimes(completedLaps_AbsoluteIdx), completedLaps_AbsoluteIdx, 'go', 'MarkerFaceColor', 'k'); % completed starts
    plot(endTimes(completedLaps_AbsoluteIdx), completedLaps_AbsoluteIdx, 'ko', 'MarkerFaceColor', 'k');   % completed ends

    plot(startTimes(abortedLaps_AbsoluteIdx), abortedLaps_AbsoluteIdx, 'ro', 'MarkerFaceColor', 'r');     % aborted starts
    plot(endTimes(abortedLaps_AbsoluteIdx), abortedLaps_AbsoluteIdx, 'ro', 'MarkerFaceColor', 'r');       % aborted ends

    % --- Labels and axis formatting ---
    xlabel('Time (s)');
    ylabel('Lap #');
    title('Lap Start and End Timeline');
    legend({'Aborted Lap', 'Completed Lap'}, 'Location', 'southeast');


    summaryText = sprintf('Total laps: %d\nCompleted: %d\nAborted: %d', ...
        length(response.lapCountAll), length(completedLaps_AbsoluteIdx), length(abortedLaps_AbsoluteIdx));

    % Add it to the upper-right corner of the axes
    xPos = max(endTimes) + 1;
    yPos = nLaps;

    text(xPos, yPos, summaryText, ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

end

end


% function [response, sessionFileInfo] = extractVRAndPeripheralData(sessionFileInfo,  VRStimName, plotFlag)
% %   Extracts wheel speed, virtual position, and lap-related info
% %   from aligned Bonsai and peripheral data during 2P-VR Aman's Classical Corridor.
% %   Handles lap classification (completed/aborted) and optionally plots lap timing.
% %
% % Inputs:
% %   sessionFileInfo : struct
% %       Contains file paths and metadata for the session.
% %
% %   plotFlag : logical (optional)
% %       If true, plots a visual summary of lap start/end times. Default: false.
% %
% % Outputs:
% %   response : struct
% %       Contains key behavioral and timing variables for downstream lap-by-lap analysis:
% %       - wheelSpeed             : real-time wheel speed (cm/s)
% %       - mouseVirtualPosition   : virtual track position (1–140 cm)
% %       - trackIDFromMousePosition : track IDs inferred from mouse position (usually 1)
% %       - mouseRecordedPosition  : raw position signal (-1141 to -1000)
% %       - trackIDs               : track ID per lap (all 1s in classical VR corridor)
% %       - lapCount               : unified lap index
% %       - blockIDs               : cumulative index of track switches (e.g., block transitions; NA)
% %       - trialType              : task/trial type ID per lap
% %                                  (0-NoTask; 1-Passive; 2-Hybrid)
% %       - completedLaps          : indices of completed laps
% %       - abortedLaps            : indices of aborted laps
% %       - startTimeAll           : Bonsai start times for all laps
% %       - endTimeAll             : parsed lap end time (based on position trace)
% %       - completedStartTimes    : start times for completed laps
% %       - completedEndTimes      : end times for completed laps
% %
% % Example usage:
% %   response = extractVRAndPeripheralData(sessionFileInfo, true);
% %
% % Authors: Sonali & Aman (based on Diao/Masa code), Updated June 2026
% 
% %% Default parameters 
% if nargin < 3, plotFlag = false; end
% 
% %% Load data files
% for iStim = 1:length(sessionFileInfo.stimFiles)
%     bonsaiData.isVRstim(iStim) = strcmp(VRStimName, sessionFileInfo.stimFiles(iStim).name);
% end
% iStim = find(bonsaiData.isVRstim==1);
% if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
%    exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file') && ...
%    exist(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'file')
%     % Load peripheral and bonsai structs
%     load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
%     load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
%     
%     % Direct variable name extraction targeting TwoPFrameTimes explicitly
%     %     % "Map" string to see which timebase was used
%     load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'resample2PTimeUsed');
% 
%     % This avoids loading massive F or spks matrices
%     timeData = load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, resample2PTimeUsed);
%     t = timeData.(resample2PTimeUsed); 
% else
%     error('Peripheral, Bonsai, or Processed 2P Data missing for VR recording');
% end
% 
% %% Create response data file
% stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_Response_' sessionFileInfo.stimFiles(iStim).name '.mat'];
% sessionFileInfo.stimFiles(iStim).Response = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
% 
% %% 'Real' wheel position and real speed
% tickToCmConversion = 3.1415 * 20 / 1024;  % Wheel radius 20 cm, 1024 ticks per revolution
% displacement = [0; diff(peripheralData.Wheel.Value * tickToCmConversion)];
% % Handle unrealistic large changes (e.g., due to teleportation or resets)
% displacement(displacement < -100) = 0;  % Negative large jumps
% displacement(displacement > 100) = 0;   % Positive large jumps
% % Calculate speed (in cm/s)
% response.wheelSpeed = displacement ./ [0; diff(peripheralData.Wheel.sampleTimes)]; 
% if isfield(peripheralData, 'Pupil')
%     response.pupilArea = peripheralData.Pupil.Value.Area;
% end 
% 
% %% Virtual position and virtual speed
% mouseVirtualPosition = nan(1,length(bonsaiData.MousePos.Value));
% trackIDFromPosition = nan(1,length(bonsaiData.MousePos.Value));
% 
% if contains(VRStimName, 'VRCorr', 'IgnoreCase', true)
%     mouseVirtualPosition(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -990)) ...
%         = floor(abs(bonsaiData.MousePos.Value(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -990))+1140)) + 1;
%     trackIDFromPosition(find(bonsaiData.MousePos.Value >= -1141 & bonsaiData.MousePos.Value < -990)) = 1; 
%     mouseVirtualPosition(mouseVirtualPosition > 140) = 140;
%     response.mouseVirtualPosition = mouseVirtualPosition';
% elseif contains(VRStimName, 'Baseline', 'IgnoreCase', true) || contains(VRStimName, 'LandManipCorridor', 'IgnoreCase', true) 
%     rawPos = bonsaiData.MousePos.Value;
%     % Clean up hardware wobble before shifting
%     rawPos(rawPos < 0) = 0;
%     rawPos(rawPos > 200) = 200;
%     % Shifting strategy: Track (0-199) maps to (1-200), ITI (200) maps to 201
%     mouseVirtualPosition = rawPos + 1;
%     if contains(VRStimName, 'Baseline', 'IgnoreCase', true)
%         trackIDFromPosition(:) = 2;
%     elseif contains(VRStimName, 'LandManipCorridor', 'IgnoreCase', true)
%         trackIDFromPosition(:) = 3; 
%     end 
%     response.mouseVirtualPosition = mouseVirtualPosition;
% end
% 
% % LapCounts
% response.lapCountAll = (1:length(bonsaiData.TrialInfo.StartTimeAll))';  
% 
% % Trial type for each lap
% if isfield(bonsaiData.TrialInfo, 'Trial_type')
%     response.trialTypeAll = bonsaiData.TrialInfo.trialType;
% else
%     response.trialTypeAll = zeros(1, length(bonsaiData.TrialInfo.StartTimeAll))';
% end
% 
% %% Find completed and aborted laps
% completedLaps_AbsoluteIdx = [];
% abortedLaps_AbsoluteIdx = [];
% lapStartTimeAll = bonsaiData.TrialInfo.StartTimeAll;
% EndTimeAll = NaN(length(lapStartTimeAll), 1); 
% 
% if ~isempty(lapStartTimeAll)
%     x = response.mouseVirtualPosition;  
%     
%     % --- SECTION 2: MATCH INDICES WITH AUTO-FORWARD COMPENSATING JUMP ---
%     startIdx = zeros(length(response.lapCountAll), 1); 
%     for nlap = 1:length(response.lapCountAll)
%         matchIdx = find(t == lapStartTimeAll(nlap), 1);
%         if isempty(matchIdx)
%             [~, matchIdx] = min(abs(t - lapStartTimeAll(nlap)));
%         end
%         
%         % If matched frame catches the old track or an interpolation mid-point slide,
%         % but the next frame is firmly at baseline, advance start index 1 frame.
%         if matchIdx < length(x) && x(matchIdx) > 5 && x(matchIdx+1) < 5
%             startIdx(nlap) = matchIdx + 1;
%         else
%             startIdx(nlap) = matchIdx;
%         end
%     end
%     
%     % --- DIAGNOSTIC PRINT: REAL DATA INITIAL SLICE TRANSITION INSPECTION ---
%     fprintf('\n===================================================================================================\n');
%     fprintf('  DIAGNOSTIC REPORT: F0 LOOKBACK WITH INITIAL 5 SLICE FRAMES (ALIGNMENT CORRECTED)\n');
%     fprintf('===================================================================================================\n');
%     for dlap = 1:length(startIdx)
%         stRow = startIdx(dlap);
%         F0_val = NaN; if stRow > 1, F0_val = x(stRow - 1); end
%         
%         if dlap < length(startIdx)
%             tempSlice = x(stRow:startIdx(dlap+1)-1);
%         else
%             tempSlice = x(stRow:end);
%         end
%         paddedPrint = nan(1,5);
%         paddedPrint(1:min(5, length(tempSlice))) = tempSlice(1:min(5, length(tempSlice)));
%         
%         fprintf('Lap %-3d | Start Row: %-6d | F0 (Prev): %-7.2f | F1: %-7.2f | F2: %-7.2f | F3: %-7.2f | F4: %-7.2f | F5: %-7.2f\n', ...
%             dlap, stRow, F0_val, paddedPrint(1), paddedPrint(2), paddedPrint(3), paddedPrint(4), paddedPrint(5));
%     end
%     fprintf('===================================================================================================\n\n');
% 
%     % --- SECTION 3: LAP SLICING LOOP ---
%     for nlap = 1:length(startIdx)
%         if nlap < length(startIdx)
%             currentLapX = x(startIdx(nlap):startIdx(nlap+1)-1);
%             currentLapT = t(startIdx(nlap):startIdx(nlap+1)-1);
%         else
%             currentLapX = x(startIdx(nlap):end);
%             currentLapT = t(startIdx(nlap):end);
%         end
%         
%         % --- THE COSMETIC CLIFF CLEANER ---
%         % Overwrites lingering 60Hz interpolation slopes at F1 and F2 with 
%         % the true baseline position established at F3.
%         if length(currentLapX) > 3 && currentLapX(3) < 5
%             if currentLapX(1) > 5, currentLapX(1) = currentLapX(3); end
%             if currentLapX(2) > 5, currentLapX(2) = currentLapX(3); end
%         end
%         
%         if length(currentLapX) > 1 && sum(~isnan(currentLapX)) > 1
%             onTrackX = currentLapX(~isnan(currentLapX));
%             onTrackT = currentLapT(~isnan(currentLapX));
%             
%             % Initial lag-induced artifact filtration
%             if length(onTrackX) > 30
%                 endFrame = 30;
%                 if ~isempty(find(diff(onTrackX(1:endFrame)) > 5, 1))
%                     jumpIndex = find(diff(onTrackX(1:endFrame)) > 5, 1, 'last');
%                     onTrackX(1:jumpIndex) = []; onTrackT(1:jumpIndex) = [];
%                 end
%                 if length(onTrackX) > 1 && ~isempty(find(diff(onTrackX(1:endFrame)) < -5, 1))
%                     jumpIndex = find(diff(onTrackX(1:endFrame)) < -5, 1, 'last');
%                     onTrackX(1:jumpIndex) = []; onTrackT(1:jumpIndex) = [];
%                 end
%             end
%             
%             % Mid-track entry alignment rule (maps to 1 cm due to offsets)
%             if sum(onTrackX == 1) > 0
%                 startPosition = find(onTrackX == 1, 1);
%                 if startPosition < length(onTrackX) - 10
%                     onTrackX = onTrackX(startPosition:end);
%                     onTrackT = onTrackT(startPosition:end);
%                 end
%             end
%             
%             if length(onTrackX) >= 2 && onTrackX(end) ~= onTrackX(end-1)
%                 onTrackX(end) = onTrackX(end-1);
%                 onTrackT(end) = onTrackT(end-1);
%             end
%             
%             [lastPosition, lastPositionIndex] = max(onTrackX);
%            
%             if lastPositionIndex * mean(diff(onTrackT)) < 0.1
%                 onTrackX(1:lastPositionIndex) = [];
%                 onTrackT(1:lastPositionIndex) = [];
%                 if isempty(onTrackX)
%                     fprintf('Lap %d aborted: only fast track jump found.\n', nlap);
%                     abortedLaps_AbsoluteIdx = [abortedLaps_AbsoluteIdx; nlap];
%                     EndTimeAll(nlap) = NaN;
%                     continue
%                 end
%                 [lastPosition, lastPositionIndex] = max(onTrackX);
%             end
%             
%             if contains(VRStimName, 'Baseline', 'IgnoreCase', true) || contains(VRStimName, 'LandManipCorridor', 'IgnoreCase', true)
%                 condition = lastPosition >= 199.5;
%             elseif contains(VRStimName, 'VRCorr', 'IgnoreCase', true)
%                 condition = lastPosition >= 139.5;
%             end
%             
%             if condition 
%                 EndTimeAll(nlap) = onTrackT(lastPositionIndex); 
%                 completedLaps_AbsoluteIdx = [completedLaps_AbsoluteIdx; nlap]; 
%             else
%                 EndTimeAll(nlap) = onTrackT(end);
%                 abortedLaps_AbsoluteIdx = [abortedLaps_AbsoluteIdx; nlap]; 
%             end
%         else
%             fprintf('Lap %d aborted: insufficient valid position data.\n', nlap);
%             abortedLaps_AbsoluteIdx = [abortedLaps_AbsoluteIdx; nlap];
%             EndTimeAll(nlap) = NaN;
%         end
%     end
% end
% 
% response.completedLaps_AbsoluteIdx = completedLaps_AbsoluteIdx;
% response.abortedLaps_AbsoluteIdx = abortedLaps_AbsoluteIdx;
% response.endTimeAll = EndTimeAll;
% response.startTimeAll = lapStartTimeAll;
% response.completedStartTimes = response.startTimeAll(completedLaps_AbsoluteIdx);
% response.completedEndTimes   = response.endTimeAll(completedLaps_AbsoluteIdx);
% response.stimName = VRStimName;
% 
% assert(length(completedLaps_AbsoluteIdx) + length(abortedLaps_AbsoluteIdx) == length(response.lapCountAll), ...
%     'Error: Some laps were not classified into completed or aborted.');
% 
% %% Save and Filter Landmark Information (Baseline/LandManip only)
% if contains(VRStimName, 'Baseline', 'IgnoreCase', true) || contains(VRStimName, 'LandManipCorridor', 'IgnoreCase', true)
%     response.completedLandmarkNames         = bonsaiData.TrialInfo.LandmarkNames(completedLaps_AbsoluteIdx);
%     response.completedLandmarkPositions     = bonsaiData.TrialInfo.LandmarkPositions(completedLaps_AbsoluteIdx);
%     response.completedLandmarkSizes         = bonsaiData.TrialInfo.LandmarkSizes(completedLaps_AbsoluteIdx);
%     response.completedLandmarkCenterOffsets = bonsaiData.TrialInfo.LandmarkCenterOffsets(completedLaps_AbsoluteIdx);
%     response.completedLandmarkRewardValence = bonsaiData.TrialInfo.LandmarkRewardValence(completedLaps_AbsoluteIdx);
%     response.completedNumLandmarks           = bonsaiData.TrialInfo.NumLandmarks(completedLaps_AbsoluteIdx);
%     response.movementVisualGain             = bonsaiData.movementVisualGain;
% end
% 
% %% Save Data Matrices
% save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
% save(sessionFileInfo.stimFiles(iStim).Response,'-struct', 'response', '-v7.3'); 
% disp('Saved Response File Successfully.');
% 
% %% Sanity check plot
% if plotFlag
%     figure('Name', 'Lap Start and End Times'); clf;
%     startTimes = response.startTimeAll;
%     endTimes = response.endTimeAll;
%     nLaps = min(length(startTimes), length(endTimes));
%     lapIDs = 1:nLaps;
%     hold on;
%     
%     for i = 1:nLaps
%         if ismember(i, completedLaps_AbsoluteIdx)
%             plot([startTimes(i), endTimes(i)], [lapIDs(i), lapIDs(i)], 'k-', 'LineWidth', 1); 
%         elseif ismember(i, abortedLaps_AbsoluteIdx)
%             plot([startTimes(i), endTimes(i)], [lapIDs(i), lapIDs(i)], 'r-', 'LineWidth', 1); 
%         end
%     end
%     
%     plot(startTimes(completedLaps_AbsoluteIdx), completedLaps_AbsoluteIdx, 'go', 'MarkerFaceColor', 'k'); 
%     plot(endTimes(completedLaps_AbsoluteIdx), completedLaps_AbsoluteIdx, 'ko', 'MarkerFaceColor', 'k');   
%     plot(startTimes(abortedLaps_AbsoluteIdx), abortedLaps_AbsoluteIdx, 'ro', 'MarkerFaceColor', 'r');     
%     plot(endTimes(abortedLaps_AbsoluteIdx), abortedLaps_AbsoluteIdx, 'ro', 'MarkerFaceColor', 'r');       
%     
%     xlabel('Time (s)'); ylabel('Lap #'); title('Lap Start and End Timeline');
%     legend({'Aborted Lap', 'Completed Lap'}, 'Location', 'southeast');
%     summaryText = sprintf('Total laps: %d\nCompleted: %d\nAborted: %d', ...
%         length(response.lapCountAll), length(completedLaps_AbsoluteIdx), length(abortedLaps_AbsoluteIdx));
%     text(max(endTimes) + 1, nLaps, summaryText, 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');
% end
% end