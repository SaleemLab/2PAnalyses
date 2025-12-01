function plotRunningSpeedDistributionsVR_withSummary(MouseID)
% Structure to store metrics across sessions for later summary plots
sessionSummaryData = struct('sessionName', {}, 'percentageRun', {}, 'meanRunningSpeed', {}, 'numLaps', {}, 'timeRun', {});
sessionIdx = 1; % Counter for sessions
filteredTable = filterMasterTable('Exclude', 0, 'Suite2PPreprocessing', 1, 'MouseID', MouseID);
mouseInfo = sessionsToProcess(filteredTable);
sessionNames = mouseInfo{1, 2};
% Pre-calculate total plots and collect necessary data for plotting
allVrStimNames = {};
% Structure to hold data needed for the plotting loops, linking stims to file info
sessionPlotData = struct('sessionName', {}, 'stimNames', {}, 'fileInfo', {});
plotDataIdx = 1;
for thisSession = 1:length(sessionNames)
    sessionName = sessionNames{thisSession};
    sessionFileInfoFilePath = getSessionFileInfoFilePath(MouseID, sessionName);
    disp(['Loading sessionFileInfo for: ' sessionName])
    % Load sessionFileInfo once per session
    load(sessionFileInfoFilePath, 'sessionFileInfo');
    try
        stimList = {sessionFileInfo.stimFiles.name};
        % Filter for all relevant VR stimuli
        combinedStimNames = stimList(contains(stimList, 'CombinedRuns', 'IgnoreCase', true));
        vrCorrStimNames = stimList(contains(stimList, 'VRCorr', 'IgnoreCase', true));
        % Apply the selection logic: CombinedRuns OR VRCorr
        if ~isempty(combinedStimNames)
            vrStimNames = combinedStimNames(1); 
        elseif ~isempty(vrCorrStimNames)
            vrStimNames = vrCorrStimNames;
        else
            vrStimNames = {};
        end
        if ~isempty(vrStimNames)
            allVrStimNames = [allVrStimNames, vrStimNames];
            % Store file info and selected stim names for the main plotting loop
            sessionPlotData(plotDataIdx).sessionName = sessionName;
            sessionPlotData(plotDataIdx).stimNames = vrStimNames;
            sessionPlotData(plotDataIdx).fileInfo = sessionFileInfo;
            plotDataIdx = plotDataIdx + 1;
        end
    catch
        continue;
    end
end
% Running Speed Distributions (Histograms)
numSpeedPlots = length(allVrStimNames);
if numSpeedPlots == 0
    fprintf('No VR stimuli found for mouse %s.\n', MouseID);
    return;
end
numCols = ceil(sqrt(numSpeedPlots));
numRows = ceil(numSpeedPlots / numCols);
hFig1 = figure('Name', sprintf('1. Running Speed Distributions for Mouse: %s', MouseID));
tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');
fprintf('Processing Mouse: %s\n', MouseID);
% Lap Times Per Session (Line Plots)
% Use the same dimensions as the speed histograms for Figure 2
hFig2 = figure('Name', sprintf('2. Lap Time Performance Per Session for Mouse: %s', MouseID));
tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');
currentLapTile = 1;
% Plotting Loop 
for thisPlot = 1:length(sessionPlotData)
    sessionName = sessionPlotData(thisPlot).sessionName;
    vrStimNames = sessionPlotData(thisPlot).stimNames;
    sessionFileInfo = sessionPlotData(thisPlot).fileInfo; % Reuse loaded file info
    fprintf('\n-- Processing Session: %s --\n', sessionName);
    try
        fprintf('Found %d VR stimulus file(s).\n', length(vrStimNames));
        for thisVRStim = 1:length(vrStimNames)
            vrStimName = vrStimNames{thisVRStim};
            fprintf('Processing VR Stim: %s\n', vrStimName);
            stimIdx = find(strcmp(vrStimName, {sessionFileInfo.stimFiles.name}), 1); 
            if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
            % Load only the 'response' data
            load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');
            sessionWheelSpeed = response.wheelSpeed;
            % Define 'running' as speed > 1 
            runningSpeeds = sessionWheelSpeed(sessionWheelSpeed > 1);
            % Calculate Metrics
            % Percentage run 
            percentageRun = round(100 * length(runningSpeeds) / length(sessionWheelSpeed));
            % Mean Speed when Running
            meanSpeed = mean(runningSpeeds);
           %  Number of Laps 
            numLaps = length(response.completedLaps); 
            % Time Run (s)
            lapDurations = response.completedEndTimes - response.completedStartTimes;
            timeRun = mean(lapDurations, 'omitnan');
            if thisVRStim == 1 % Store metrics once per session
                sessionSummaryData(sessionIdx).sessionName = sessionName;
                sessionSummaryData(sessionIdx).percentageRun = percentageRun;
                sessionSummaryData(sessionIdx).meanRunningSpeed = meanSpeed;
                sessionSummaryData(sessionIdx).numLaps = numLaps;
                sessionSummaryData(sessionIdx).timeRun = timeRun; 
                sessionIdx = sessionIdx + 1;
            end
            % Speed distribution 
            figure(hFig1); % switch to fig 1
            nexttile; 
            histogram(runningSpeeds, 'Normalization', 'count', 'DisplayName', vrStimName, 'BinEdges', 1:5:60); 
            correctedVRStimName = replace(vrStimName, '_', '\_');
            title([correctedVRStimName ' (' num2str(percentageRun) '% Run)']);
            ylabel('Count');
            xlabel('Speed (cm/s)'); 
            % Time Spent per Lap in fig 2
            figure(hFig2); % Switch to fig 2
            nexttile(currentLapTile);
            lapNumbers = 1:length(lapDurations);
            plot(lapNumbers, lapDurations, '-o', 'MarkerSize', 4, 'LineWidth', 1.5);
            title(['Lap Times: ' correctedVRStimName]);
            xlabel('Lap Number');
            ylabel('Time per Lap (s)');
            grid on;
            currentLapTile = currentLapTile + 1;
        end
    catch ME
        fprintf('Error processing session %s: %s\n', sessionName, ME.message);
    end
end
figure(hFig1); sgtitle(sprintf('1. Running Speed Distributions for Mouse: %s', MouseID));
figure(hFig2); sgtitle(sprintf('2. Lap Time Performance Per Session for Mouse: %s', MouseID));
% Summary plots across sessions 
if ~isempty(sessionSummaryData)
    % Extract data into arrays 
    days = 1:length(sessionSummaryData);
    percentages = [sessionSummaryData.percentageRun];
    meanSpeeds = [sessionSummaryData.meanRunningSpeed];
    laps = [sessionSummaryData.numLaps];
    timeRuns = [sessionSummaryData.timeRun]; 
    sessionLabels = {sessionSummaryData.sessionName};
    
    % Create a 4x1 tiled layout for the summary plots in Figure 3
    hFig3 = figure('Name', sprintf('3. Summary Metrics Across Days for Mouse: %s', MouseID));
    tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    %Percentage Run Across Days
    nexttile;
    plot(days, percentages, '-o', 'MarkerSize',10, 'LineWidth',2);
    title('Percentage Run Across Days');
    ylabel('Percentage Run (%)');
    xticks(days);
    xticklabels(replace(sessionLabels, '_', '\_'));
    xtickangle(45);
    %Time Run Across Days 
    nexttile;
    plot(days, timeRuns, '-o', 'MarkerSize',10, 'LineWidth',2);
    title('Time Run Across Days (mean across laps)');
    ylabel('Time Run (s)');
    xticks(days);
    xticklabels(replace(sessionLabels, '_', '\_'));
    xtickangle(45);
    %Mean Speed When Running Across Days
    nexttile;
    plot(days, meanSpeeds, '-o', 'MarkerSize',10, 'LineWidth',2);
    title('Mean Speed When Running Across Days');
    ylabel('Mean Speed (cm/s)');
    xticks(days);
    xticklabels(replace(sessionLabels, '_', '\_'));
    xtickangle(45);
    %Number of Laps Run Across Days
    nexttile;
    plot(days, laps, '-o', 'MarkerSize',10,'LineWidth',2);
    title('Number of Laps Run Across Days');
    ylabel('Laps Run');
    xlabel('Session / Day');
    xticks(days);
    xticklabels(replace(sessionLabels, '_', '\_'));
    xtickangle(45);
    sgtitle(sprintf('3. Summary Running Metrics for Mouse: %s', MouseID));
end
end