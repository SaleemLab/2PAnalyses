% %% Fixed threshold speed modulation analysis
% % Low = <15 cm/s, High = >25 cm/s, discard 15-25 cm/s
% 
% lowThresh  = 18;
% highThresh = 25;
% 
% allSpeedMod   = [];
% allAnimalID   = [];
% animalCounter = 0;
% 
% for thisMouse = 1:size(mouseInfo, 1)
%     mousenumber  = mouseInfo{thisMouse, 1};
%     sessionNames = mouseInfo{thisMouse, 2};
%     animalCounter = animalCounter + 1;
%     fprintf('Processing Mouse: %s\n', mousenumber);
% 
%     for thisSession = 1:length(sessionNames)
%         sessionName = sessionNames{thisSession};
%         fprintf('  Processing Session: %s\n', sessionName);
% 
%         infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
%         if ~isfile(infoPath), fprintf('  Info Missing — skipping\n'); continue; end
%         loadedInfo      = load(infoPath, 'sessionFileInfo');
%         sessionFileInfo = loadedInfo.sessionFileInfo;
%         stimNames       = {sessionFileInfo.stimFiles.name};
% 
%         targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "CombinedRuns"), 1);
%         if isempty(targetIdx)
%             allCorridorIdx = find(contains(stimNames, "Corridor"));
%             if isscalar(allCorridorIdx)
%                 targetIdx = allCorridorIdx;
%             elseif length(allCorridorIdx) > 1
%                 targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "00002"), 1);
%             end
%         end
%         if isempty(targetIdx), fprintf('  No valid corridor — skipping\n'); continue; end
% 
%         filePath = sessionFileInfo.stimFiles(targetIdx).Response;
% 
%         %% Load raw data
%         rawData        = load(filePath, 'lapPositionActivity', 'lapPositionRunningSpeed', 'trialIndicesByCondition', 'SMI_Metrics');
%         rawROIData     = rawData.lapPositionActivity.dFFNeuropilCorrected; % [numROIs x numTrials x numPosBins]
%         speedData      = rawData.lapPositionRunningSpeed;                  % [numTrials x numPosBins]
%         baselineTrials = rawData.trialIndicesByCondition.Baseline;
% 
%         rawROIData = rawROIData(:, baselineTrials, :);
%         speedData  = speedData(baselineTrials, :);
% 
%         [numROIs, ~, numPosBins] = size(rawROIData);
% 
%         %% Check frame counts at landmarks
%         % landmarks = [40, 80, 120, 160];
%         % fprintf('  Frame counts per landmark:\n');
%         % for lm = 1:length(landmarks)
%         %     window = landmarks(lm)-15 : landmarks(lm)+15;
%         %     speedAtLandmark = mean(speedData(:, window), 2);
%         %     nLow  = sum(speedAtLandmark < lowThresh);
%         %     nHigh = sum(speedAtLandmark > highThresh);
%         %     fprintf('    Landmark %d cm: nLow = %d, nHigh = %d\n', landmarks(lm), nLow, nHigh);
%         % end
% 
%         allSpeeds   = speedData(:);
%         nLow  = sum(allSpeeds < lowThresh  & ~isnan(allSpeeds));
%         nHigh = sum(allSpeeds > highThresh & ~isnan(allSpeeds));
%         fprintf('  Total frames: nLow = %d, nHigh = %d\n', nLow, nHigh);
% 
%         % exclude sessions with unbalanced sampling
%         maxRatio = 3;
%         frameRatio = max(nLow, nHigh) / min(nLow, nHigh);
%         if nLow < 2 || nHigh < 2 || frameRatio > maxRatio
%             fprintf('  Skipping — unbalanced sampling: ratio = %.1f\n', frameRatio);
%             continue
%         end
% 
%         %% Load crossValExpVar and filter ROIs
%         rawData2       = load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'crossValExpVar');
%         crossValExpVar = rawData2.crossValExpVar;
% 
%         FilteredROIS_idx = find(crossValExpVar.dFFNeuropilCorrected.pValues <= 0.01 & ...
%             crossValExpVar.dFFNeuropilCorrected.medianExpVar > 0.1 & ...
%             ~rawData.SMI_Metrics.dFFNeuropilCorrected.ExcludeEdgePeakCells);
% 
%         if isempty(FilteredROIS_idx)
%             fprintf('  No spatially modulated ROIs — skipping\n'); continue;
%         end
%         fprintf('  nROIs: %d\n', length(FilteredROIS_idx));
% 
%         %% Build low and high matrices
%         matrixLow  = nan(numROIs, numPosBins);
%         matrixHigh = nan(numROIs, numPosBins);
% 
%         for b = 1:numPosBins
%             binSpeeds = speedData(:, b);
% 
%             lowIdx  = binSpeeds < lowThresh  & ~isnan(binSpeeds);
%             highIdx = binSpeeds > highThresh & ~isnan(binSpeeds);
% 
%             if sum(lowIdx) >= 2
%                 matrixLow(:, b)  = mean(squeeze(rawROIData(:, lowIdx, b)), 2, 'omitnan');
%             end
%             if sum(highIdx) >= 2
%                 matrixHigh(:, b) = mean(squeeze(rawROIData(:, highIdx, b)), 2, 'omitnan');
%             end
%         end
% 
%         %% Compute speed modulation per ROI
%         filteredLow  = matrixLow(FilteredROIS_idx, :);
%         filteredHigh = matrixHigh(FilteredROIS_idx, :);
% 
% 
%         % restrict to landmark windows ±10cm
%         landmarkBins = [30:50, 70:90, 110:130, 150:170];
% 
%         filteredLow_lm  = filteredLow(:,  landmarkBins);
%         filteredHigh_lm = filteredHigh(:, landmarkBins);
% 
%         speedModPerROI = mean(filteredLow_lm - filteredHigh_lm, 2, 'omitnan');
% 
%         % speedModPerROI = mean(filteredLow - filteredHigh, 2, 'omitnan');
% 
%         allSpeedMod   = [allSpeedMod;   speedModPerROI];
%         allAnimalID   = [allAnimalID;   repmat(animalCounter, length(speedModPerROI), 1)];
% 
%         fprintf('  Mean modulation: %.4f\n', mean(speedModPerROI, 'omitnan'));
%     end
% end
% 
% %% Plot per animal and pooled
% figure('Position', [100 100 800 400], 'Color', 'w');
% 
% subplot(1,3,1)
% histogram(allSpeedMod(allAnimalID == 1), 30, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none')
% xline(0, '--r', 'HandleVisibility', 'off')
% xline(mean(allSpeedMod(allAnimalID == 1), 'omitnan'), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
% [~, p1] = ttest(allSpeedMod(allAnimalID == 1));
% title(sprintf('M25132\np = %.4f', p1))
% xlabel('Mean (Low - High) \DeltaF/F')
% ylabel('Number of ROIs')
% set(gca, 'Box', 'off', 'TickDir', 'out')
% 
% subplot(1,3,2)
% histogram(allSpeedMod(allAnimalID == 2), 30, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none')
% xline(0, '--r', 'HandleVisibility', 'off')
% xline(mean(allSpeedMod(allAnimalID == 2), 'omitnan'), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
% [~, p2] = ttest(allSpeedMod(allAnimalID == 2));
% title(sprintf('M26003\np = %.4f', p2))
% xlabel('Mean (Low - High) \DeltaF/F')
% set(gca, 'Box', 'off', 'TickDir', 'out')
% 
% subplot(1,3,3)
% histogram(allSpeedMod, 30, 'FaceColor', [0.4 0.6 0.4], 'EdgeColor', 'none')
% xline(0, '--r', 'HandleVisibility', 'off')
% xline(mean(allSpeedMod, 'omitnan'), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
% [~, p] = ttest(allSpeedMod);
% title(sprintf('Pooled\np = %.4f', p))
% xlabel('Mean (Low - High) \DeltaF/F')
% set(gca, 'Box', 'off', 'TickDir', 'out')
% 
% sgtitle(sprintf('Speed modulation: <%.0f cm/s vs >%.0f cm/s', lowThresh, highThresh))



%% Fixed threshold speed modulation analysis — landmark positions only
% Low = <18 cm/s, High = >25 cm/s
% Only landmark position bins included throughout

lowThresh    = 18;
highThresh   = 25;
maxRatio     = 3;
landmarkBins = [30:50, 70:90, 110:130, 150:170];

allSpeedMod  = [];
allAnimalID  = [];
animalCounter = 0;

for thisMouse = 1:size(mouseInfo, 1)
    mousenumber  = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};
    animalCounter = animalCounter + 1;
    fprintf('Processing Mouse: %s\n', mousenumber);

    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        fprintf('  Processing Session: %s\n', sessionName);

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if ~isfile(infoPath), fprintf('  Info Missing — skipping\n'); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};

        targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "CombinedRuns"), 1);
        if isempty(targetIdx)
            allCorridorIdx = find(contains(stimNames, "Corridor"));
            if isscalar(allCorridorIdx)
                targetIdx = allCorridorIdx;
            elseif length(allCorridorIdx) > 1
                targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "00002"), 1);
            end
        end
        if isempty(targetIdx), fprintf('  No valid corridor — skipping\n'); continue; end

        filePath = sessionFileInfo.stimFiles(targetIdx).Response;

        %% Load raw data
        rawData        = load(filePath, 'lapPositionActivity', 'lapPositionRunningSpeed', 'trialIndicesByCondition', 'SMI_Metrics');
        rawROIData     = rawData.lapPositionActivity.dFFNeuropilCorrected; % [numROIs x numTrials x numPosBins]
        speedData      = rawData.lapPositionRunningSpeed;                  % [numTrials x numPosBins]
        baselineTrials = rawData.trialIndicesByCondition.Baseline;

        rawROIData = rawROIData(:, baselineTrials, :);
        speedData  = speedData(baselineTrials, :);

        [numROIs, ~, ~] = size(rawROIData);

        %% Check frame counts at landmark positions only
        landmarkSpeeds    = speedData(:, landmarkBins);
        allLandmarkSpeeds = landmarkSpeeds(:);

        nLow  = sum(allLandmarkSpeeds < lowThresh  & ~isnan(allLandmarkSpeeds));
        nHigh = sum(allLandmarkSpeeds > highThresh & ~isnan(allLandmarkSpeeds));
        fprintf('  Landmark frames: nLow = %d, nHigh = %d\n', nLow, nHigh);

        frameRatio = max(nLow, nHigh) / min(nLow, nHigh);
        if nLow < 2 || nHigh < 2 || frameRatio > maxRatio
            fprintf('  Skipping — unbalanced sampling: ratio = %.1f\n', frameRatio);
            continue
        end

        %% Load crossValExpVar and filter ROIs
        rawData2       = load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'crossValExpVar');
        crossValExpVar = rawData2.crossValExpVar;

        FilteredROIS_idx = find( ...
            crossValExpVar.dFFNeuropilCorrected.pValues    <= 0.01 & ...
            crossValExpVar.dFFNeuropilCorrected.medianExpVar > 0.1 & ...
            ~rawData.SMI_Metrics.dFFNeuropilCorrected.ExcludeEdgePeakCells);

        if isempty(FilteredROIS_idx)
            fprintf('  No spatially modulated ROIs — skipping\n'); continue;
        end
        fprintf('  nROIs: %d\n', length(FilteredROIS_idx));

        %% Build low and high matrices at landmark positions only
        matrixLow  = nan(numROIs, length(landmarkBins));
        matrixHigh = nan(numROIs, length(landmarkBins));

        for bi = 1:length(landmarkBins)
            b         = landmarkBins(bi);
            binSpeeds = speedData(:, b);

            lowIdx  = binSpeeds < lowThresh  & ~isnan(binSpeeds);
            highIdx = binSpeeds > highThresh & ~isnan(binSpeeds);

            if sum(lowIdx) >= 2
                matrixLow(:, bi)  = mean(squeeze(rawROIData(:, lowIdx, b)), 2, 'omitnan');
            end
            if sum(highIdx) >= 2
                matrixHigh(:, bi) = mean(squeeze(rawROIData(:, highIdx, b)), 2, 'omitnan');
            end
        end

        %% Compute speed modulation per ROI
        filteredLow  = matrixLow(FilteredROIS_idx, :);
        filteredHigh = matrixHigh(FilteredROIS_idx, :);

        speedModPerROI = mean(filteredLow - filteredHigh, 2, 'omitnan');

        allSpeedMod  = [allSpeedMod;  speedModPerROI];
        allAnimalID  = [allAnimalID;  repmat(animalCounter, length(speedModPerROI), 1)];

        fprintf('  Mean modulation: %.4f\n', mean(speedModPerROI, 'omitnan'));
    end
end

%% Plot per animal and pooled
figure('Position', [100 100 800 400], 'Color', 'w');

subplot(1,3,1)
histogram(allSpeedMod(allAnimalID == 1), 30, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none')
xline(0, '--r', 'HandleVisibility', 'off')
xline(mean(allSpeedMod(allAnimalID == 1), 'omitnan'), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
[~, p1] = ttest(allSpeedMod(allAnimalID == 1));
title(sprintf('M25132\np = %.4f', p1))
xlabel('Mean (Low - High) \DeltaF/F')
ylabel('Number of ROIs')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(1,3,2)
histogram(allSpeedMod(allAnimalID == 2), 30, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none')
xline(0, '--r', 'HandleVisibility', 'off')
xline(mean(allSpeedMod(allAnimalID == 2), 'omitnan'), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
[~, p2] = ttest(allSpeedMod(allAnimalID == 2));
title(sprintf('M26003\np = %.4f', p2))
xlabel('Mean (Low - High) \DeltaF/F')
set(gca, 'Box', 'off', 'TickDir', 'out')

subplot(1,3,3)
histogram(allSpeedMod, 30, 'FaceColor', [0.4 0.6 0.4], 'EdgeColor', 'none')
xline(0, '--r', 'HandleVisibility', 'off')
xline(mean(allSpeedMod, 'omitnan'), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
[~, p] = ttest(allSpeedMod);
title(sprintf('Pooled\np = %.4f', p))
xlabel('Mean (Low - High) \DeltaF/F')
set(gca, 'Box', 'off', 'TickDir', 'out')

sgtitle(sprintf('Landmark speed modulation: <%.0f cm/s vs >%.0f cm/s', lowThresh, highThresh))

%% Print summary
fprintf('\nSummary:\n')
fprintf('  M25132: n = %d ROIs, mean = %.4f, p = %.4f\n', ...
        sum(allAnimalID==1), mean(allSpeedMod(allAnimalID==1), 'omitnan'), p1)
fprintf('  M26003: n = %d ROIs, mean = %.4f, p = %.4f\n', ...
        sum(allAnimalID==2), mean(allSpeedMod(allAnimalID==2), 'omitnan'), p2)
fprintf('  Pooled: n = %d ROIs, mean = %.4f, p = %.4f\n', ...
        length(allSpeedMod), mean(allSpeedMod, 'omitnan'), p)