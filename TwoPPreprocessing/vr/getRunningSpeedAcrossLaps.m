function [response, sessionFileInfo] = getRunningSpeedAcrossLaps(sessionFileInfo, VRStimName, processedTwoPData, doPlot)
% getRunningSpeedAcrossLaps Loads running speed data, splits it into laps and 
% position bins, and prepares the data for distribution plots.
%
%   sessionFileInfo: Struct containing information about the session.
%   VRStimName: String specifying the name of the stimulus file to load.
%   processedTwoPData: Struct containing the two-photon and time data.
%
%   response.lapRunningSpeed: Cell array of speed vectors for each lap (time-based, unfiltered).
%   response.lapPositionRunningSpeed: Matrix of mean speed per lap and position bin (calculated on the fly).

if nargin < 4, doPlot = false; end 
% Load data 
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}), 1);
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
disp('Loading Response and processedTwoPData struct...');


% Load the 'response' struct
response = load(sessionFileInfo.stimFiles(stimIdx).Response);

% Load 'processedTwoPData' if not passed as an argument
if nargin < 3
    % "Map" string to see which timebase was used
    load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'resample2PTimeUsed');
    % This avoids loading massive F or spks matrices
    timeData = load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, resample2PTimeUsed);
    timeVec = timeData.(resample2PTimeUsed); 
else 
    timeVec = processedTwoPData.(processedTwoPData.resample2PTimeUsed); 
end 

% Get the main variables
sessionWheelSpeed = response.wheelSpeed;
if isempty(sessionWheelSpeed)
    error('response.wheelSpeed is empty');
end
if ~isfield(response, 'mouseVirtualPosition')
    error('response.mouseVirtualPosition is required for position binning but is missing.');
end

% Define Variables for Lap/Binning
nLaps = length(response.completedStartTimes);
if contains(VRStimName, 'Baseline') || contains(VRStimName, 'LandManipCorridor')
    posBinEdges = 1:201; 
elseif contains(VRStimName, 'VRCorr')
    posBinEdges = 1:141; 
end

numBins = length(posBinEdges) - 1;

% Initialize cell array to temporarily store indices for both lap speed and binning
% Each cell will contain the [speed vector, index vector] for the lap
lapData = cell(1, nLaps);
%%
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

pngFilePath = fullfile(figSaveDir, ...
    [sessionFileInfo.stimFiles(stimIdx).name '_RunningSpeedAcrossLaps.png']);
%% Calculate Full Lap Speed (Time-Based, Unfiltered)
if nLaps > 0
    for thisLap = 1:nLaps
        disp(['Processing Lap: ' num2str(thisLap)]);
        
        lapStart = response.completedStartTimes(thisLap);
        lapEnd = response.completedEndTimes(thisLap);
        
        % Find all 2P frame indices that occurred within this lap's time window
        lapWheelSpeedIdx = find(timeVec >= lapStart & timeVec <= lapEnd);
        
        if isempty(lapWheelSpeedIdx)
            disp(['No wheel speed indices found for lap ' num2str(thisLap) '. Assigning empty array.']);
            lapData{thisLap} = {[], []}; % {speed, indices}
            continue;
        end
        
        % Store speed and indices for later use
        lapSpeed = sessionWheelSpeed(lapWheelSpeedIdx);
        lapData{thisLap} = {lapSpeed, lapWheelSpeedIdx}; 
    end
else
    warning('Lap timing information not found. Skipping lap-wise analysis.');
end

%% Calculate Position-Binned Speed 
if nLaps > 0
    binnedSpeed = nan(nLaps, numBins); 
    
    % Loop over laps using the indices stored in lapData
    for thisLap = 1:nLaps
        [~, lapWheelSpeedIdx] = lapData{thisLap}{:};
        
        if isempty(lapWheelSpeedIdx)
            continue;
        end

        % Get position and speed data for the current lap
        lapPosition = response.mouseVirtualPosition(lapWheelSpeedIdx);
        lapSpeedData = sessionWheelSpeed(lapWheelSpeedIdx);
        
        % Assign each frame's position to a position bin
        % discretize returns the bin number (1 to numBins) for each position point
        binID = discretize(lapPosition, posBinEdges); 
        
        % Here we enforce valid bins (1 to numBins).
        validIndices = ~isnan(binID);
        validBinID = binID(validIndices);
        validSpeed = lapSpeedData(validIndices);

        % Loop over bins to calculate the mean speed
        for thisBin = 1:numBins
            % Find all speed values that fall into the current bin
            binSpeedValues = validSpeed(validBinID == thisBin);
            
            if ~isempty(binSpeedValues)
                % Calculate the mean speed for this bin
                meanSpeed = mean(binSpeedValues, 'omitnan');
                
                % Store the resulting mean speed
                binnedSpeed(thisLap, thisBin) = meanSpeed;
            end
        end
    end
else
    binnedSpeed = [];
end

% Extract final lapSpeed cell array from lapData for saving
lapSpeedFinal = cellfun(@(x) x{1}, lapData, 'UniformOutput', false);
response.lapRunningSpeed = lapSpeedFinal';
response.lapPositionRunningSpeed = binnedSpeed;

disp('Running speed data has been successfully processed into session, lap, and binned structures.');
disp('Saving response with updated lapRunningSpeed and lapPositionRunningSpeed...');
disp(['Saving updated Response struct (Lap Activity) to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');

%% 
if doPlot

    if nLaps > 0 && ~isempty(binnedSpeed)
        figure('Name', 'Running Speed Analysis');
        subplot(121);

        %mean speed profile across all laps
        meanSpeedProfile = mean(binnedSpeed, 1, 'omitnan');

        if any(~isnan(meanSpeedProfile))
            plot(meanSpeedProfile, 'LineWidth', 2);

            % Set X-axis to match position (1-140)
            xlim([1 numBins]);

            title('Mean Speed');
            ylabel('Mean Speed (cm/s)');
            xlabel('Position on Track (cm)');
            xticks([1 40 80 120 160 200]);
            xticklabels({'1', '40', '80', '120', '160', '200'});
            xline(40, 'k--', 'LineWidth', 2.5);
            xline(80, 'k--', 'LineWidth', 2.5);
            xline(120, 'k--', 'LineWidth', 2.5);
            xline(160, 'k--', 'LineWidth', 2.5);
            grid on;
        else
            text(0.5, 0.5, 'Mean speed profile is empty.', 'HorizontalAlignment', 'center');
        end

        %  Heatmap of Binned Speed (Lap vs. Position)
        subplot(122);
        imagesc(binnedSpeed);
        xticks([1 40 80 120 160 200]);
        xticklabels({'1', '40', '80', '120', '160', '200'});

        yTickPositions = 1:5:nLaps; % Show every 5th lap
        set(gca, 'YTick', yTickPositions);
        colorbar;
        xline(40, 'k--', 'LineWidth', 2.5);
        xline(80, 'k--', 'LineWidth', 2.5);
        xline(120, 'k--', 'LineWidth', 2.5);
        xline(160, 'k--', 'LineWidth', 2.5);
        title('Binned Running Speed Across Laps');
        xlabel('Position on Track (cm)');
        ylabel('Lap Number');
        axis tight;

        set(gcf, 'PaperUnits', 'inches', ...
            'PaperPosition', [0 0 11 8.5], ...
            'PaperOrientation', 'landscape');
        print(gcf, pngFilePath, '-dpng', '-r300');

    end

end
end