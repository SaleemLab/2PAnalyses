function [response, sessionFileInfo] = getRunningSpeedAcrossLaps(sessionFileInfo, VRStimName, processedTwoPData)
% getRunningSpeedAcrossLaps Loads running speed data, splits it into laps and 
% position bins, and prepares the data for distribution plots.
%
%   sessionFileInfo: Struct containing information about the session.
%   VRStimName: String specifying the name of the stimulus file to load.
%   processedTwoPData: Struct containing the two-photon and time data.
%
%   response.lapRunningSpeed: Cell array of speed vectors for each lap (time-based, unfiltered).
%   response.lapPositionRunningSpeed: Matrix of mean speed per lap and position bin (calculated on the fly).

% Load data 
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}), 1);
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
disp('Loading Response and processedTwoPData struct...');

% Load the 'response' struct
load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');

% Load 'processedTwoPData' if not passed as an argument
if nargin < 3
    load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');
    disp('Loading processedTwoPData Struct from SessionFileInfo ')
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
posBinEdges = 0:140;
numBins = length(posBinEdges) - 1;
timeVec = processedTwoPData.(processedTwoPData.resample2PTimeUsed);

% Initialize cell array to temporarily store indices for both lap speed and binning
% Each cell will contain the [speed vector, index vector] for the lap
lapData = cell(1, nLaps);


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


%% Calculate Position-Binned Speed (Calculated On-the-Fly)

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
        [~, binID] = histc(lapPosition, posBinEdges); 
        
        % The last bin returned by histc for the last edge is numBins+1 (if pos is exactly the max edge).
        % Set any points falling outside the main range (0 to 140) to 0 or remove them.
        % Here we enforce valid bins (1 to numBins).
        validIndices = (binID > 0) & (binID <= numBins);
        
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

response.lapRunningSpeed = lapSpeedFinal;
response.lapPositionRunningSpeed = binnedSpeed;

disp('Running speed data has been successfully processed into session, lap, and binned structures.');
disp('Saving response with updated lapRunningSpeed and lapPositionRunningSpeed...');
save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');


%% Visualization Block
if nLaps > 0 && ~isempty(binnedSpeed)
    
    figure('Name', 'Running Speed Analysis');
    subplot(121);
    
    %mean speed profile across all laps
    meanSpeedProfile = mean(binnedSpeed, 1, 'omitnan');
    
    if any(~isnan(meanSpeedProfile))
        plot(meanSpeedProfile, 'LineWidth', 2);
        
        % Set X-axis to match position (1-140)
        xlim([1 numBins]);
        % set(gca, 'XTick', xTickPositions);
        % set(gca, 'XTickLabel', xTickLabels);
        
        title('Mean Speed');
        ylabel('Mean Speed (cm/s)');
        xlabel('Position on Track (cm)');
        xticks([0 50 70 90 110 140]);
        xticklabels({'0', '50', '70', '90', '110', '140'});
        xline(50, 'k--', 'LineWidth', 2.5);
        xline(70, 'k--', 'LineWidth', 2.5);
        xline(90, 'k--', 'LineWidth', 2.5);
        xline(110, 'k--', 'LineWidth', 2.5);
        grid on;
    else
        text(0.5, 0.5, 'Mean speed profile is empty.', 'HorizontalAlignment', 'center');
    end
    
     %  Heatmap of Binned Speed (Lap vs. Position)
    subplot(122);
    imagesc(binnedSpeed);
    caxis([1 60]); 
    % xTickPositions = 1:20:numBins;
    % xTickLabels = posBinEdges(xTickPositions);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    % set(gca, 'XTick', xTickPositions);
    % set(gca, 'XTickLabel', xTickLabels);

    yTickPositions = 1:5:nLaps; % Show every 5th lap
    set(gca, 'YTick', yTickPositions);

    % colormap('hot'); 
    colorbar;
    xline(50, 'k--', 'LineWidth', 2.5);
    xline(70, 'k--', 'LineWidth', 2.5);
    xline(90, 'k--', 'LineWidth', 2.5);
    xline(110, 'k--', 'LineWidth', 2.5);
    title('Binned Running Speed Across Laps');
    xlabel('Position on Track (cm)');
    ylabel('Lap Number');
    axis tight;
    
    % subplot(133)
    % lapActivity = response.lapPositionActivity.dFFNeuropilCorrected;
    % 
    % % Optional plot-time smoothing
    % applySmoothing = true; 
    % if applySmoothing
    %     w = gausswin(9); w = w / sum(w);
    %     for iCell = 1:size(lapActivity, 1)
    %         for iLap = 1:size(lapActivity, 2)
    %             trace = squeeze(lapActivity(iCell, iLap, :));
    %             if all(isnan(trace)), continue; end
    %             nanMask = isnan(trace);
    %             trace(nanMask) = 0;
    %             smoothed = filtfilt(w, 1, trace);
    %             smoothed(nanMask) = NaN;
    %             lapActivity(iCell, iLap, :) = smoothed;
    %         end
    %     end
    % end
    % 
    % % Find the mean running speed across laps 
    % meanRunningPerLap = mean(response.lapPositionRunningSpeed, 2);
    % [~, sortIdx] = sort(meanRunningPerLap);
    % 
    % meanAcrossAllLaps = squeeze(mean(lapActivity, 2, 'omitnan'));
    % % Normalize across position bins
    % normAll = normalize(meanAcrossAllLaps, 2, 'range');
    % 
    % imagesc(normAll(sortIdx, :));
    % caxis([0 1]); 
    % set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 18, 'YDir', 'normal');
    % xline(50, 'k--', 'LineWidth', 2.5);
    % xline(70, 'k--', 'LineWidth', 2.5);
    % xline(90, 'k--', 'LineWidth', 2.5);
    % xline(110, 'k--', 'LineWidth', 2.5);
    % xticks([0 50 70 90 110 140]);
    % xticklabels({'0', '50', '70', '90', '110', '140'});
    % xlabel('Position (cm)');
    % ylabel('ROIs');
    % title('Laps sorted by mean running speed');
    % colorbar; ylabel(colorbar, 'Activity (normalised)');


% %% Combine all laps
% meanAcrossAllLaps = squeeze(mean(lapActivity, 2, 'omitnan'));
% 
% % Normalize across position bins
% normAll = normalize(meanAcrossAllLaps, 2, 'range');


end

end

