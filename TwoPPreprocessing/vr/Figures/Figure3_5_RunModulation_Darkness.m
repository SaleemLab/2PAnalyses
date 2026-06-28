%% RSP speed tuning in darkness

%% Parameters
pairs = struct();
pairs.M25132 = {'20260219','20260223','20260226','20260228','20260303','20260313','20260306'};
pairs.M25133 = {'20260219','20260223','20260221'};
pairs.M26003 = {'20260316','20260322','20260324','20260325'};

targetStruct          = 'tuningCurve';
useField              = 'dFFNeuropilCorrected';
minR2Threshold        = 0.5;
pval_shuffleThreshold = 0.01;
saveFolder            = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.5\';

colors.highpass        = [0.89, 0.10, 0.11];   % red
colors.bandpass        = [0.99, 0.60, 0.17];   % orange
colors.lowpass         = [0.12, 0.47, 0.71];   % blue
colors.trough_inverted = [0.20, 0.63, 0.17];   % green

%% Filter and load session info
filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'Darkness', 'GrayScreen'});
mouseInfo     = sessionsToProcess(filteredTable);

%% accumulators — must be before session loop
tuningTypesPlot    = {'highpass', 'bandpass', 'lowpass', 'trough_inverted', 'unclassified'};
faceColors_fig     = {colors.highpass, colors.bandpass, colors.lowpass, colors.trough_inverted, [0.75 0.75 0.75]};
labelNames         = {'High-pass', 'Band-pass', 'Low-pass', 'Trough-inv', 'Unclassified'};
allSessionPcts_fig = [];

allTypes      = {};
allPrefSpeeds = [];
allR2s        = [];
allOccupancy = zeros(0, 4); 

sessionsProcessedCount = 0;
totalSessionsToProcess = sum(cellfun(@length, mouseInfo(:, 2)));
fprintf('Parsing %d mouse cohorts for tuned ROIs...\n', size(mouseInfo, 1));

%% Load data across sessions
for thisMouse = 1:size(mouseInfo, 1)
    mousenumber  = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};
    fprintf('Processing Mouse: %s\n', mousenumber);

    for thisSession = 1:length(sessionNames)
        sessionName            = sessionNames{thisSession};
        sessionsProcessedCount = sessionsProcessedCount + 1;

        fprintf('  Processing Session [%d/%d]: %s  %s\n', ...
            sessionsProcessedCount, totalSessionsToProcess, mousenumber, sessionName);

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if isempty(infoPath) || ~isfile(infoPath)
            fprintf('    WARNING: File path missing or invalid for %s - %s. Skipping.\n', mousenumber, sessionName);
            continue;
        end

        matchRow = find(strcmp(filteredTable.MouseID, mousenumber) & strcmp(filteredTable.Session, sessionName), 1);
        if isempty(matchRow)
            fprintf('    WARNING: Failed to find matched row: %s - %s. Skipping.\n', mousenumber, sessionName);
            continue;
        end

        fprintf('    Imaged structure: %s\n', filteredTable.TypeImaged{matchRow});

        try
            loadedInfo      = load(infoPath, 'sessionFileInfo');
            sessionFileInfo = loadedInfo.sessionFileInfo;
            stimNames       = {sessionFileInfo.stimFiles.name};
            targetIdx       = find(contains(stimNames, 'Darkness'));

            if ~any(targetIdx), continue; end
            fprintf('    Found %d target stimulus file(s).\n', length(targetIdx));

            for thisStim = 1:length(targetIdx)
                thisStimName = stimNames{targetIdx(thisStim)};
                stimFileName = sprintf('%s_%s_Response_%s.mat', mousenumber, sessionName, thisStimName);
                fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

                if ~isfile(fileFullPath)
                    warning('Could not find response file: %s. Skipping.', stimFileName);
                    continue;
                end

                loadedData = load(fileFullPath, 'response');
                resp       = loadedData.response;

                if ~isfield(resp, targetStruct) || ...
                   ~isfield(resp.(targetStruct), useField) || ...
                   ~isfield(resp.(targetStruct).(useField), 'classification')
                    continue;
                end

                cls        = resp.(targetStruct).(useField).classification;
                pvalMoving = resp.(targetStruct).(useField).pValMoving;

                classifiedIdx   = find(cls.R2 >= minR2Threshold & pvalMoving <= pval_shuffleThreshold);
                unclassifiedIdx = find(cls.R2 <  minR2Threshold & pvalMoving <= pval_shuffleThreshold);
                untunedIdx      = find(pvalMoving > pval_shuffleThreshold);

                fprintf('    ROIs — classified: %d | unclassified: %d | untuned: %d\n', ...
                    numel(classifiedIdx), numel(unclassifiedIdx), numel(untunedIdx));

                % --- Per-session percentages for figure ---
                totalSig = numel(classifiedIdx) + numel(unclassifiedIdx);
                if totalSig > 0
                    sessionCounts = zeros(1, numel(tuningTypesPlot));
                    for k = 1:numel(tuningTypesPlot)-1
                        sessionCounts(k) = sum(strcmp(cls.tuningType(classifiedIdx), tuningTypesPlot{k}));
                    end
                    sessionCounts(end) = numel(unclassifiedIdx);
                    allSessionPcts_fig = [allSessionPcts_fig; sessionCounts / totalSig * 100]; %#ok<AGROW>
                end

                % Check occupancy
                occ = resp.(targetStruct).occupancy.moving;
                nBins = numel(occ);
                fprintf('    Occupancy — min=%.1f  max=%.1f  mean=%.1f seconds | nBins=%d\n', ...
                    min(occ), max(occ), mean(occ), nBins);
                allOccupancy = [allOccupancy; min(occ), max(occ), mean(occ), nBins];
                if any(occ == 0)
                    fprintf('    WARNING: session has empty speed bins — tuning curves unreliable!\n');
                end

                if nBins < 5
                    warning('    WARNING: only %d speed bin(s) — skipping session.\n', nBins);
                    continue;
                end
                % pool info 
                if ~isempty(classifiedIdx)
                    allTypes      = [allTypes;      cls.tuningType(classifiedIdx)];            
                    allPrefSpeeds = [allPrefSpeeds; cls.preferredSpeed(classifiedIdx)];        
                    allR2s        = [allR2s;        cls.R2(classifiedIdx)];                    
                end
                if ~isempty(unclassifiedIdx)
                    allTypes      = [allTypes;      repmat({'unclassified'}, numel(unclassifiedIdx), 1)]; %#ok<AGROW>
                    allPrefSpeeds = [allPrefSpeeds; cls.preferredSpeed(unclassifiedIdx)];                 %#ok<AGROW>
                    allR2s        = [allR2s;        cls.R2(unclassifiedIdx)];                            %#ok<AGROW>
                end
                if ~isempty(untunedIdx)
                    allTypes      = [allTypes;      repmat({'untuned'}, numel(untunedIdx), 1)]; %#ok<AGROW>
                    allPrefSpeeds = [allPrefSpeeds; nan(numel(untunedIdx), 1)];                 %#ok<AGROW>
                    allR2s        = [allR2s;        cls.R2(untunedIdx)];                        %#ok<AGROW>
                end

            end % thisStim
        catch ME
            fprintf('    ERROR on %s-%s: %s\n', mousenumber, sessionName, ME.message);
        end
    end % thisSession
end % thisMouse

%% Summary
tuningTypes      = {'lowpass', 'trough_inverted', 'bandpass', 'highpass'};
tuningLabels     = {'Low-pass', 'Trough-inverted', 'Band-pass', 'High-pass'};
totalROIsCounted = numel(allTypes);
nClassified      = sum(~strcmp(allTypes, 'untuned') & ~strcmp(allTypes, 'unclassified'));
nUnclassified    = sum(strcmp(allTypes, 'unclassified'));
nUntuned         = sum(strcmp(allTypes, 'untuned'));

fprintf('\nTotal ROIs: %d | Classified: %d | Sessions: %d\n', ...
    totalROIsCounted, nClassified, size(allSessionPcts_fig,1));
fprintf('Classified (sig + R2 pass):  %d  (%.1f%%)\n', nClassified,   100*nClassified/totalROIsCounted);
fprintf('Unclassified (sig, low R2):  %d  (%.1f%%)\n', nUnclassified, 100*nUnclassified/totalROIsCounted);
fprintf('Untuned (not sig):           %d  (%.1f%%)\n', nUntuned,      100*nUntuned/totalROIsCounted);
fprintf('Classified breakdown:\n');
for k = 1:numel(tuningTypes)
    n = sum(strcmp(allTypes, tuningTypes{k}));
    fprintf('  %-20s %d  (%.1f%% of classified)\n', tuningLabels{k}, n, 100*n/nClassified);
end
fprintf('\nPreferred speeds:\n');
for k = 1:numel(tuningTypes)
    speeds = allPrefSpeeds(strcmp(allTypes, tuningTypes{k}));
    fprintf('  %-20s min=%.1f  max=%.1f  median=%.1f\n', tuningLabels{k}, min(speeds), max(speeds), median(speeds));
end


fprintf('\nOccupancy across all sessions:\n');
fprintf('  Mean of means: %.1f seconds\n', mean(allOccupancy(:,3)));
fprintf('  Overall min:   %.1f seconds\n', min(allOccupancy(:,1)));
fprintf('  Overall max:   %.1f seconds\n', max(allOccupancy(:,2)));
fprintf('  Bin counts:    ');
fprintf('%d ', allOccupancy(:,4)');
fprintf('\n');
%% Figure A: Running speed and dF/F example
proc2PData   = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_processed2PData_M26003_Darkness_20260322_00001.mat");
peripheraleg = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_PeripheralData_M26003_Darkness_20260322_00001.mat");
respeg       = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_Darkness_20260322_00001.mat");

isSigShuffle = find(respeg.response.tuningCurve.dFFNeuropilCorrected.pValFull < 0.001);
dFF          = proc2PData.processedSignals.dFFNeuropilCorrected(isSigShuffle, :);
timeVec      = proc2PData.TwoPFrameTime;

tickToCm   = 3.1416 * 20 / 1024;
wheelTimes = peripheraleg.peripheralData.Wheel.sampleTimes;
runSpeed   = [0; diff(peripheraleg.peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(wheelTimes)];
runSpeed(abs(runSpeed) > 150) = NaN;

tStart = 47; tEnd = 360;

keep_dff   = timeVec >= tStart & timeVec <= tEnd;
timeVec    = timeVec(keep_dff);
dFF        = dFF(:, keep_dff);

keep_wheel = wheelTimes >= tStart & wheelTimes <= tEnd;
wheelTimes = wheelTimes(keep_wheel);
runSpeed   = runSpeed(keep_wheel);

[~, nearestIdx]  = min(abs(wheelTimes - timeVec'), [], 1);
runSpeed_matched = runSpeed(nearestIdx);
runSpeed_matched = runSpeed_matched(:);

corrs = corr(dFF', runSpeed_matched, 'rows', 'complete');
[~, sortOrder] = sort(corrs, 'ascend');
dFF = dFF(sortOrder, :);

dFF_norm = dFF - min(dFF, [], 2);
dFF_norm = dFF_norm ./ max(dFF_norm, [], 2);

figA = figure('Color', 'w', 'Position', [100 100 900 400]);
ax1  = subplot(2,1,1);
plot(wheelTimes, runSpeed, 'k', 'LineWidth', 0.8);
ylabel('Speed cm/s', 'FontSize', 9);
title('RSC bouton activity in darkness', 'FontSize', 10);
set(ax1, 'Box', 'off', 'TickDir', 'out');

ax2 = subplot(2,1,2);
imagesc(timeVec, 1:size(dFF_norm,1), dFF_norm);
colormap(ax2, flipud(gray));
ylabel('Boutons', 'FontSize', 9);
set(ax2, 'Box', 'off', 'TickDir', 'out', 'YDir', 'normal');

cb = colorbar(ax2, 'eastoutside');
cb.Label.String = '\DeltaF/F (norm.)';
lo = prctile(dFF(:), 10);
hi = prctile(dFF(:), 95);
clim([lo hi]);
cb.Ticks = [lo hi]; cb.TickLabels = {'0', '>0.7'};

linkaxes([ax1 ax2], 'x');
drawnow;
pos1 = get(ax1, 'Position'); pos2 = get(ax2, 'Position');
set(ax2, 'Position', [pos1(1), pos2(2), pos1(3), pos2(4)]);
set(cb,  'Position', [pos1(1)+pos1(3)+0.01, pos2(2), 0.02, pos2(4)]);

outputDir = fullfile(saveFolder, 'exampleWheel_dFF\');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figA, fullfile(outputDir, 'wheelSpeed_dff_m26003_day5'));

%% Figure B: Gaussian fit examples
outputDir = fullfile(saveFolder, 'exampleFits\');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fits1 = load('Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260223\M25132_20260223_Response_M25132_Darkness_20260223_00001.mat');
figB1 = plotAllSpeedTuningCategories(fits1.response, useField);
saveFigureFormats(figB1, fullfile(outputDir, 'M25132_20260223_takebandpass'));

fits2 = load('Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260220\M25132_20260220_Response_M25132_Darkness_20260220_00001.mat');
figB2 = plotAllSpeedTuningCategories(fits2.response, useField);
saveFigureFormats(figB2, fullfile(outputDir, 'M25132_20260220_taketrough'));

fits3 = load('Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_Darkness_20260322_00001.mat');
figB3 = plotAllSpeedTuningCategories(fits3.response, useField);
saveFigureFormats(figB3, fullfile(outputDir, 'M26003_20260322_takelowpass'));

fits4 = load('Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260316\M26003_20260316_Response_M26003_Darkness_20260316_00001.mat');
figB4 = plotAllSpeedTuningCategories(fits4.response, useField);
saveFigureFormats(figB4, fullfile(outputDir, 'M26003_20260316_takehighpass'));

%% Figure C: Category counts with per-session woth dots
meanPct = mean(allSessionPcts_fig, 1);
semPct  = std(allSessionPcts_fig, 0, 1) / sqrt(size(allSessionPcts_fig, 1));

figC1 = figure('Position', [100, 100, 380, 360], 'Color', 'w');
ax_c1 = gca;
hold on;

for k = 1:numel(tuningTypesPlot)
    bar(k, meanPct(k), 'FaceColor', faceColors_fig{k}, 'EdgeColor', 'none', 'BarWidth', 0.6);
end

errorbar(1:numel(tuningTypesPlot), meanPct, semPct, 'k.', 'LineWidth', 1, 'CapSize', 4);

for k = 1:numel(tuningTypesPlot)
    xJitter = k + (rand(size(allSessionPcts_fig,1), 1) - 0.5) * 0.25;
    scatter(xJitter, allSessionPcts_fig(:,k), 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5);
end

ax_c1.XTick              = 1:numel(tuningTypesPlot);
ax_c1.XTickLabel         = labelNames;
ax_c1.XTickLabelRotation = 30;
ax_c1.FontSize           = 8;
ax_c1.FontName           = 'Arial';
ax_c1.TickDir            = 'out';
ax_c1.LineWidth          = 0.75;
ax_c1.Position           = [0.18, 0.28, 0.75, 0.62];
box off;

ylabel('% of sig. boutons', 'FontSize', 9, 'FontName', 'Arial');
title('RSP',                'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('%d rois out of %d (sig + R2 pass) | n=%d sessions', ...
    nClassified, totalROIsCounted, size(allSessionPcts_fig,1)), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');

defaultAxesProperties(ax_c1, 0);
offsetAxes(ax_c1);

saveFigureFormats(figC1, fullfile(saveFolder, 'speedModulated_histogram_darkness\category_counts_darkness'));

%% Figure D: Band-pass preferred speed histogram
prefSpeedBins = [2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5];
binCentres    = [5, 10, 15, 20, 25, 30];
bp_speeds     = allPrefSpeeds(strcmp(allTypes, 'bandpass'));
histBP        = histcounts(bp_speeds, prefSpeedBins);

figC2 = figure('Position', [100, 100, 380, 360], 'Color', 'w');
ax_c2 = gca;
hold on;

bar(binCentres, histBP, 'FaceColor', colors.bandpass, 'EdgeColor', 'none', 'BarWidth', 0.5);

xlim([prefSpeedBins(1), prefSpeedBins(end)]);
ax_c2.XTick      = binCentres;
ax_c2.XTickLabel = {'5','10','15','20','25','30'};
ylim([0, max(histBP) * 1.2]);
ax_c2.FontSize   = 8;
ax_c2.FontName   = 'Arial';
ax_c2.TickDir    = 'out';
ax_c2.LineWidth  = 0.75;
ax_c2.Layer      = 'top';
ax_c2.Position   = [0.18, 0.22, 0.75, 0.65];
box off;

xlabel('Preferred speed (cm s^{-1})', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Number of boutons',           'FontSize', 9, 'FontName', 'Arial');
title('RSP — Band-pass',              'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('n = %d boutons', numel(bp_speeds)), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');

defaultAxesProperties(ax_c2, 0);
offsetAxes(ax_c2);

saveFigureFormats(figC2, fullfile(saveFolder, 'speedModulated_histogram_darkness\preferred_speed_bandpass_darkness'));

%% Figure E: Preferred speed histogram — stacked (all categories)
prefSpeedBins = [2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5];
binCentres    = [5, 10, 15, 20, 25, 30];

histMatrix = zeros(numel(prefSpeedBins)-1, numel(tuningTypes));
for t = 1:numel(tuningTypes)
    speeds = allPrefSpeeds(strcmp(allTypes, tuningTypes{t}));
    histMatrix(:, t) = histcounts(speeds, prefSpeedBins);
end

figE = figure('Position', [100, 100, 480, 360], 'Color', 'w');
ax_e = gca;
hold on;

faceColors_stacked = {colors.lowpass, colors.trough_inverted, colors.bandpass, colors.highpass};
bHandle = bar(binCentres, histMatrix, 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.5);
for t = 1:numel(tuningTypes)
    bHandle(t).FaceColor = faceColors_stacked{t};
end

xlim([prefSpeedBins(1), prefSpeedBins(end)]);
ax_e.XTick      = binCentres;
ax_e.XTickLabel = {'5','10','15','20','25','30'};
yMax = max(sum(histMatrix, 2)) * 1.18;
ylim([0, yMax]);
ax_e.FontSize   = 8;
ax_e.FontName   = 'Arial';
ax_e.TickDir    = 'out';
ax_e.LineWidth  = 0.75;
ax_e.Layer      = 'top';
ax_e.Position   = [0.13, 0.22, 0.83, 0.65];
box off;

xlabel('Speed at peak tuning (cm s^{-1})', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Number of boutons',                'FontSize', 9, 'FontName', 'Arial');
title('RSP',                               'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('%d rois out of %d (sig + R2 pass)', nClassified, totalROIsCounted), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');

dummyHandles = [
    patch(NaN, NaN, colors.lowpass,         'EdgeColor', 'none'), ...
    patch(NaN, NaN, colors.trough_inverted, 'EdgeColor', 'none'), ...
    patch(NaN, NaN, colors.bandpass,        'EdgeColor', 'none'), ...
    patch(NaN, NaN, colors.highpass,        'EdgeColor', 'none')
];
leg = legend(dummyHandles, {'Low-pass','Trough-inverted','Band-pass','High-pass'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 7, 'FontName', 'Arial');
leg.ItemTokenSize = [10, 8];

counts_stacked = zeros(1, numel(tuningTypes));
for k = 1:numel(tuningTypes)
    counts_stacked(k) = sum(strcmp(allTypes, tuningTypes{k}));
end
annoStr = sprintf('Low-pass:    %d (%.1f%%)\nTrough-inv:  %d (%.1f%%)\nBand-pass:   %d (%.1f%%)\nHigh-pass:   %d (%.1f%%)', ...
    counts_stacked(1), 100*counts_stacked(1)/nClassified, ...
    counts_stacked(2), 100*counts_stacked(2)/nClassified, ...
    counts_stacked(3), 100*counts_stacked(3)/nClassified, ...
    counts_stacked(4), 100*counts_stacked(4)/nClassified);
text(ax_e, prefSpeedBins(end-2)*0.98, yMax*0.95, annoStr, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 7, 'FontName', 'Arial', 'Color', [0.3 0.3 0.3]);

defaultAxesProperties(ax_e, 0);
offsetAxes(ax_e);

saveFigureFormats(figE, fullfile(saveFolder, 'speedModulated_histogram_darkness\preferred_speed_histogram_darkness'));
%% --- plot tuning curves for all classified ri
% pdfSavePath = "Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.5\exampleFits\all_classified_tuningCurves.pdf";
% 
% % Delete existing file first or exportgraphics will append to old content
% if isfile(pdfSavePath)
%     delete(pdfSavePath);
%     fprintf('Deleted existing QC PDF.\n');
% end
% 
% pageCount = 0;
% 
% for thisMouse = 1:size(mouseInfo, 1)
%     mousenumber  = mouseInfo{thisMouse, 1};
%     sessionNames = mouseInfo{thisMouse, 2};
% 
%     for thisSession = 1:length(sessionNames)
%         sessionName = sessionNames{thisSession};
% 
%         infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
%         if isempty(infoPath) || ~isfile(infoPath), continue; end
% 
%         try
%             loadedInfo      = load(infoPath, 'sessionFileInfo');
%             sessionFileInfo = loadedInfo.sessionFileInfo;
%             stimNames       = {sessionFileInfo.stimFiles.name};
%             targetIdx       = find(contains(stimNames, {'Darkness'}));
% 
%             for thisStim = 1:length(targetIdx)
%                 thisStimName = stimNames{targetIdx(thisStim)};
%                 stimFileName = sprintf('%s_%s_Response_%s.mat', mousenumber, sessionName, thisStimName);
%                 fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
%                 if ~isfile(fileFullPath), continue; end
% 
%                 loadedData = load(fileFullPath, 'response');
%                 resp       = loadedData.response;
% 
%                 if ~isfield(resp, targetStruct) || ...
%                    ~isfield(resp.(targetStruct), useField) || ...
%                    ~isfield(resp.(targetStruct).(useField), 'classification')
%                     continue;
%                 end
% 
%                 cls        = resp.(targetStruct).(useField).classification;
%                 pvalMoving = resp.(targetStruct).(useField).pValMoving;
% 
%                 classifiedIdx = find(cls.R2 >= minR2Threshold & pvalMoving <= pval_shuffleThreshold);
%                 if isempty(classifiedIdx), continue; end
% 
%                 edges         = resp.(targetStruct).speedBins;
%                 movingCenters = (edges(1:end-1) + diff(edges)/2)';
%                 yMean         = resp.(targetStruct).(useField).moveMean;
% 
%                 typeColorMap = containers.Map(...
%                     {'lowpass','trough_inverted','bandpass','highpass'}, ...
%                     {colors.lowpass, colors.trough_inverted, colors.bandpass, colors.highpass});
% 
%                 nPlot = numel(classifiedIdx);
%                 nCols = 8;
%                 nRows = ceil(nPlot / nCols);
% 
%                 fig = figure('Color', 'w', ...
%                              'Position', [100 100 1600 max(nRows*180, 300)], ...
%                              'Visible', 'off');
% 
%                 for k = 1:nPlot
%                     r   = classifiedIdx(k);
%                     typ = cls.tuningType{r};
% 
%                     if isKey(typeColorMap, typ)
%                         c = typeColorMap(typ);
%                     else
%                         c = [0.5 0.5 0.5];
%                     end
% 
%                     ax_k = subplot(nRows, nCols, k);
%                     plot(movingCenters, yMean(r,:), 'o-', ...
%                         'Color', c, 'MarkerSize', 3, ...
%                         'MarkerFaceColor', c, 'LineWidth', 1);
%                     title(sprintf('%s %d | R²=%.2f', typ, r, cls.R2(r)), ...
%                         'FontSize', 5.5, 'Color', c, 'Interpreter', 'none');
%                     box off;
%                     set(ax_k, 'TickDir', 'out', 'FontSize', 5);
%                 end
% 
%                 sgtitle(sprintf('%s  %s  —  %d / %d ROIs classified', ...
%                     mousenumber, sessionName, nPlot, numel(pvalMoving)), ...
%                     'FontSize', 10, 'FontWeight', 'bold');
% 
%                 exportgraphics(fig, pdfSavePath, 'Append', true, 'Resolution', 120);
%                 close(fig);
%                 pageCount = pageCount + 1;
%                 fprintf('  Page %d written: %s %s (%d ROIs)\n', pageCount, mousenumber, sessionName, nPlot);
%             end
% 
%         catch ME
%             fprintf('ERROR on %s-%s: %s\n', mousenumber, sessionName, ME.message);
%         end
%     end
% end
% 
% fprintf('\nDone. QC PDF saved (%d pages):\n  %s\n', pageCount, pdfSavePath);

%% other response functions 

% targetStruct = 'tuningCurve';              % or 'tuningCurveFixedBins'
% useField = 'dFFNeuropilCorrected';         % change if needed
% 
% plotOpts = struct();
% plotOpts.maxPlots = 24;                    % how many examples to show
% plotOpts.nCols = 4;
% plotOpts.showSEM = true;
% plotOpts.onlySignificantMoving = true;    % set true if you want shuffle-significant only
% plotOpts.sortBy = 'R2';                    % 'R2', 'preferredSpeed', 'dynamicRange'
% 
% plotSpeedTuningCategoryExamples(response, targetStruct, useField, 'lowpass', plotOpts);
% plotSpeedTuningCategoryExamples(response, targetStruct, useField, 'highpass', plotOpts);
% plotSpeedTuningCategoryExamples(response, targetStruct, useField, 'bandpass', plotOpts);
% plotSpeedTuningCategoryExamples(response, targetStruct, useField, 'trough_inverted', plotOpts);


plotBandpassTroughCheck(pairs, targetStruct, useField, minR2Threshold, pval_shuffleThreshold)
    
        

