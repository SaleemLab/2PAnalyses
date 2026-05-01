function plotROIDifferencesAcrossConditions_HalvesStableOnly(sessionFileInfo, response, signalToUse, applySmoothing)
if nargin < 3; signalToUse = 'spks'; end
if nargin < 4; applySmoothing = true; end
%% Parameters
landmarks = [40, 80, 120, 160]; 
landNames = {'L1', 'L2', 'L3', 'L4', 'ALL_Landmarks', 'Non_Landmark'};
windowSize = 7; 
cMapDiff = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); 
            ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];
standardTicks = [40 80 120 160];
stableThresh = 0.5; % change to shuffle TODO
%% Load and pick stable rois
try
    vars = load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'lapCorr_Halves');
    stableIdx = find(vars.lapCorr_Halves.rho >= stableThresh);
catch
    stableIdx = 1:size(response.lapPositionActivity.(signalToUse), 1);
end
lapActivity = response.lapPositionActivity.(signalToUse)(stableIdx, :, :);
[nROIs, ~, nBins] = size(lapActivity);
conds = fieldnames(response.trialIndicesByCondition);
if applySmoothing
    w = gausswin(15); w = w / sum(w);
    for iCell = 1:nROIs
        for iLap = 1:size(lapActivity, 2)
            trace = squeeze(lapActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace); trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
            lapActivity(iCell, iLap, :) = smoothed;
        end
    end
end
%% Identify baseline and manipulation trials
baseIdx = find(contains(lower(conds), 'baseline') | contains(lower(conds), 'norm'), 1);
if isempty(baseIdx), baseIdx = 1; end
baseLaps = response.trialIndicesByCondition.(conds{baseIdx});
baseOddLaps = baseLaps(1:2:end);
baseEvenLaps = baseLaps(2:2:end);
otherCondNames = conds(setdiff(1:length(conds), baseIdx));
validManipNames = {};
for i = 1:length(otherCondNames)
    if ~isempty(response.trialIndicesByCondition.(otherCondNames{i}))
        validManipNames{end+1} = otherCondNames{i};
    end
end
nValidManip = length(validManipNames);
% Compute templates
meanOddAll = squeeze(mean(lapActivity(:, baseOddLaps, :), 2, 'omitnan'));
normOddAll = normalize(meanOddAll, 2, 'range');
meanEvenAll = squeeze(mean(lapActivity(:, baseEvenLaps, :), 2, 'omitnan'));
normEvenAll = normalize(meanEvenAll, 2, 'range');
allPeakPos = zeros(nROIs, 1);
for i = 1:nROIs
    [~, p] = max(normOddAll(i, :));
    allPeakPos(i) = p;
end
%% Iterate through groups
for thisLandmark = 1:6 
    if thisLandmark <= 4
        L_pos = landmarks(thisLandmark);
        currName = landNames{thisLandmark};
        batchIdx = find(allPeakPos >= (L_pos - windowSize) & allPeakPos <= (L_pos + windowSize));
    elseif thisLandmark == 5
        currName = landNames{5}; 
        batchIdx = [];
        for m = 1:4
            batchIdx = [batchIdx; find(allPeakPos >= (landmarks(m)-windowSize) & allPeakPos <= (landmarks(m)+windowSize))];
        end
        batchIdx = unique(batchIdx);
    else
        % Non-Landmark ROIs
        currName = landNames{6};
        landmarkROIs = [];
        for m = 1:4
            landmarkROIs = [landmarkROIs; find(allPeakPos >= (landmarks(m)-windowSize) & allPeakPos <= (landmarks(m)+windowSize))];
        end
        batchIdx = setdiff(1:nROIs, unique(landmarkROIs));
    end
    
    if isempty(batchIdx), continue; end
    
    [~, localSort] = sort(allPeakPos(batchIdx), 'ascend');
    targetIdx = batchIdx(localSort); 
    numInGroup = length(targetIdx);
    nCols = nValidManip + 2;
    fig = figure('Color', 'w', 'Position', [50 50 250*nCols 900]);
    t = tiledlayout(4, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('Batch: %s (n=%d)', currName, numInGroup), 'FontWeight', 'bold');
    %% ROW 1: HEATMAPS
    % Base odd
    nexttile(1);
    imagesc(1:nBins, 1:numInGroup, normOddAll(targetIdx, :));
    colormap(gca, flipud(gray)); clim([0 1]);
    title('Base Odd'); ylabel('ROIs'); set(gca, 'XTick', []); 
    % Base even
    nexttile(2);
    imagesc(1:nBins, 1:numInGroup, normEvenAll(targetIdx, :));
    colormap(gca, flipud(gray)); clim([0 1]);
    title('Base Even'); set(gca, 'XTick', [], 'YTick', []);
    % Manipulations
    for iM = 1:nValidManip
        nexttile(2+iM);
        laps = response.trialIndicesByCondition.(validManipNames{iM});
        mMean = squeeze(mean(lapActivity(targetIdx, laps, :), 2, 'omitnan'));
        imagesc(1:nBins, 1:numInGroup, normalize(mMean, 2, 'range'));
        colormap(gca, flipud(gray)); clim([0 1]);
        title(strrep(validManipNames{iM}, '_', ' ')); set(gca, 'XTick', [], 'YTick', []);
    end
    %% ROW 2: MEANS
    % Base odd
    nexttile(nCols + 1); 
    popMean = mean(normOddAll(targetIdx, :), 1, 'omitnan');
    popSEM = std(normOddAll(targetIdx, :), 0, 1, 'omitnan') ./ sqrt(numInGroup);
    xBins = 1:nBins;
    hold on;
    fill([xBins, fliplr(xBins)], [popMean + popSEM, fliplr(popMean - popSEM)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(xBins, popMean, 'k', 'LineWidth', 1.5);
    ylabel('Norm. \DeltaF/F'); xticks(standardTicks); box off;
    % Base even
    nexttile(nCols + 2); 
    popMean = mean(normEvenAll(targetIdx, :), 1, 'omitnan');
    popSEM = std(normEvenAll(targetIdx, :), 0, 1, 'omitnan') ./ sqrt(numInGroup);
    hold on;
    fill([xBins, fliplr(xBins)], [popMean + popSEM, fliplr(popMean - popSEM)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(xBins, popMean, 'k', 'LineWidth', 1.5);
    xticks(standardTicks); box off;
    
    % Manip Mean
    for iM = 1:nValidManip
        nexttile(nCols + 2 + iM);
        laps = response.trialIndicesByCondition.(validManipNames{iM});
        mMeanAll = squeeze(mean(lapActivity(targetIdx, laps, :), 2, 'omitnan'));
        normData = normalize(mMeanAll, 2, 'range');
        popMean = mean(normData, 1, 'omitnan');
        popSEM = std(normData, 0, 1, 'omitnan') ./ sqrt(numInGroup);
        hold on;
        fill([xBins, fliplr(xBins)], [popMean + popSEM, fliplr(popMean - popSEM)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        plot(xBins, popMean, 'k', 'LineWidth', 1.5); 
        ylim([0.1 0.9]); xticks(standardTicks); box off;
    end
    %% ROW 3: DIFFERENCE MAPS
    % placeholder 
    nexttile(2*nCols + 1); axis off; 
    % Base even repeat (ref)
    nexttile(2*nCols + 2); 
    imagesc(1:nBins, 1:numInGroup, normEvenAll(targetIdx, :));
    colormap(gca, flipud(gray)); clim([0 1]);
    title('Ref: Base Even'); set(gca, 'XTick', [], 'YTick', []);
    % difference
    diffStore = {};
    for iM = 1:nValidManip
        nexttile(2*nCols + 2 + iM);
        laps = response.trialIndicesByCondition.(validManipNames{iM});
        mMean = squeeze(mean(lapActivity(targetIdx, laps, :), 2, 'omitnan'));
        mNorm = normalize(mMean, 2, 'range');
        diffMap = mNorm - normEvenAll(targetIdx, :);
        diffStore{iM} = diffMap;
        
        imagesc(1:nBins, 1:numInGroup, diffMap);
        colormap(gca, cMapDiff); clim([-0.4 0.4]); 
        title('Difference'); set(gca, 'XTick', [], 'YTick', []);
    end
    %% ROW 4: MEAN DIFFERENCES
    nexttile(3*nCols + 1); axis off;
    
    % ref even mean repeat
    nexttile(3*nCols + 2); 
    popMean = mean(normEvenAll(targetIdx, :), 1, 'omitnan');
    popSEM = std(normEvenAll(targetIdx, :), 0, 1, 'omitnan') ./ sqrt(numInGroup);
    hold on;
    fill([xBins, fliplr(xBins)], [popMean + popSEM, fliplr(popMean - popSEM)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(xBins, popMean, 'k', 'LineWidth', 1.5);
    xticks(standardTicks); box off; ylabel('Ref Mean');
    for iM = 1:nValidManip
        nexttile(3*nCols + 2 + iM);
        dMean = mean(diffStore{iM}, 1, 'omitnan');
        dSEM = std(diffStore{iM}, 0, 1, 'omitnan') ./ sqrt(numInGroup);
        hold on;
        fill([xBins, fliplr(xBins)], [dMean + dSEM, fliplr(dMean - dSEM)], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        plot(xBins, dMean, 'k', 'LineWidth', 1.2);
        yline(0, '--', 'Color', [0.5 0.5 0.5]);
        ylim([-0.4 0.4]); ylabel('Diff'); xlabel('Pos (cm)');
        xticks(standardTicks); box off;
    end
    saveName = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', ...
        [sessionFileInfo.animal_name '_' currName '_Labeled_Diffs.png']);
    exportgraphics(fig, saveName, 'Resolution', 300);
    close(fig);
end 
end