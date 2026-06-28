%% V1 Population Speed Tuning Landscape
pairs = struct();
pairs.M26005 = {'20260305', '20260306', '20260311', '20260318', '20260321', '20260322'};
pairs.M26004 = {'20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322'};
pairs.M25131 = {'20260312', '20260313', '20260314', '20260318', '20260321', '20260322'};
pairs.M25126 = {'20260311', '20260312', '20260313'};

targetStruct          = 'tuningCurve';
useField              = 'dFFNeuropilCorrected';
minR2Threshold        = 0.5;
pval_shuffleThreshold = 0.01;
saveFolder_v1         = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.5\V1\';
if ~exist(saveFolder_v1, 'dir'), mkdir(saveFolder_v1); end

colors.lowpass         = [0.20, 0.53, 0.74];
colors.trough_inverted = [0.85, 0.37, 0.01];
colors.bandpass        = [0.13, 0.63, 0.36];
colors.highpass        = [0.58, 0.18, 0.55];

%% Filter and load session info
filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'Darkness', 'GrayScreen'});
mouseInfo     = sessionsToProcess(filteredTable);

%% Initialise accumulators — before loop
tuningTypesPlot_v1 = {'highpass', 'bandpass', 'lowpass', 'trough_inverted', 'unclassified'};
faceColors_v1      = {colors.highpass, colors.bandpass, colors.lowpass, colors.trough_inverted, [0.75 0.75 0.75]};
labelNames_v1      = {'High-pass', 'Band-pass', 'Low-pass', 'Trough-inv', 'Unclassified'};
allSessionPcts_v1  = [];

allTypes      = {};
allPrefSpeeds = [];
allR2s        = [];

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
            targetIdx       = find(contains(stimNames, 'GrayScreen'));

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

                % Check occupancy
                if isfield(resp.(targetStruct), 'occupancy')
                    occ   = resp.(targetStruct).occupancy.moving;
                    nBins = numel(occ);
                    fprintf('    Occupancy — min=%.1f  max=%.1f  mean=%.1f | nBins=%d\n', ...
                        min(occ), max(occ), mean(occ), nBins);
                    if nBins < 5
                        fprintf('    WARNING: only %d bin(s) — skipping session.\n', nBins);
                        continue;
                    end
                end

                classifiedIdx   = find(cls.R2 >= minR2Threshold & pvalMoving <= pval_shuffleThreshold);
                unclassifiedIdx = find(cls.R2 <  minR2Threshold & pvalMoving <= pval_shuffleThreshold);
                untunedIdx      = find(pvalMoving > pval_shuffleThreshold);

                fprintf('    ROIs — classified: %d | unclassified: %d | untuned: %d\n', ...
                    numel(classifiedIdx), numel(unclassifiedIdx), numel(untunedIdx));

                % --- Per-session percentages ---
                totalSig = numel(classifiedIdx) + numel(unclassifiedIdx);
                if totalSig > 0
                    sessionCounts = zeros(1, numel(tuningTypesPlot_v1));
                    for k = 1:numel(tuningTypesPlot_v1)-1
                        sessionCounts(k) = sum(strcmp(cls.tuningType(classifiedIdx), tuningTypesPlot_v1{k}));
                    end
                    sessionCounts(end) = numel(unclassifiedIdx);
                    allSessionPcts_v1  = [allSessionPcts_v1; sessionCounts / totalSig * 100]; 
                end

                % --- Pooled accumulation ---
                if ~isempty(classifiedIdx)
                    allTypes      = [allTypes;      cls.tuningType(classifiedIdx)];           
                    allPrefSpeeds = [allPrefSpeeds; cls.preferredSpeed(classifiedIdx)];       
                    allR2s        = [allR2s;        cls.R2(classifiedIdx)];                  
                end
                if ~isempty(unclassifiedIdx)
                    allTypes      = [allTypes;      repmat({'unclassified'}, numel(unclassifiedIdx), 1)];
                    allPrefSpeeds = [allPrefSpeeds; cls.preferredSpeed(unclassifiedIdx)];                
                    allR2s        = [allR2s;        cls.R2(unclassifiedIdx)];                          
                end
                if ~isempty(untunedIdx)
                    allTypes      = [allTypes;      repmat({'untuned'}, numel(untunedIdx), 1)]; 
                    allPrefSpeeds = [allPrefSpeeds; nan(numel(untunedIdx), 1)];                
                    allR2s        = [allR2s;        cls.R2(untunedIdx)];                        
                end

            end % thisStim
        catch ME
            fprintf('    ERROR on %s-%s: %s\n', mousenumber, sessionName, ME.message);
            fprintf('    Line: %d\n', ME.stack(1).line);  % ADD THIS
        end
    end % thisSession
end % thisMouse

%% Summary
tuningTypes_v1   = {'lowpass', 'trough_inverted', 'bandpass', 'highpass'};
tuningLabels_v1  = {'Low-pass', 'Trough-inverted', 'Band-pass', 'High-pass'};
totalROIsCounted = numel(allTypes);
nClassified      = sum(~strcmp(allTypes, 'untuned') & ~strcmp(allTypes, 'unclassified'));
nUnclassified    = sum(strcmp(allTypes, 'unclassified'));
nUntuned         = sum(strcmp(allTypes, 'untuned'));

fprintf('\nTotal ROIs: %d | Classified: %d | Sessions: %d\n', ...
    totalROIsCounted, nClassified, size(allSessionPcts_v1,1));
fprintf('Classified (sig + R2 pass):  %d  (%.1f%%)\n', nClassified,   100*nClassified/totalROIsCounted);
fprintf('Unclassified (sig, low R2):  %d  (%.1f%%)\n', nUnclassified, 100*nUnclassified/totalROIsCounted);
fprintf('Untuned (not sig):           %d  (%.1f%%)\n', nUntuned,      100*nUntuned/totalROIsCounted);
fprintf('Classified breakdown:\n');
for k = 1:numel(tuningTypes_v1)
    n = sum(strcmp(allTypes, tuningTypes_v1{k}));
    fprintf('  %-20s %d  (%.1f%% of classified)\n', tuningLabels_v1{k}, n, 100*n/nClassified);
end
fprintf('\nPreferred speeds:\n');
for k = 1:numel(tuningTypes_v1)
    speeds = allPrefSpeeds(strcmp(allTypes, tuningTypes_v1{k}));
    fprintf('  %-20s min=%.1f  max=%.1f  median=%.1f\n', tuningLabels_v1{k}, min(speeds), max(speeds), median(speeds));
end

%% Figure 1 — Category counts with per-session dots
meanPct = mean(allSessionPcts_v1, 1);
semPct  = std(allSessionPcts_v1, 0, 1) / sqrt(size(allSessionPcts_v1, 1));

figV1_C1 = figure('Position', [100, 100, 380, 360], 'Color', 'w');
ax_v1c1  = gca;
hold on;

for k = 1:numel(tuningTypesPlot_v1)
    bar(k, meanPct(k), 'FaceColor', faceColors_v1{k}, 'EdgeColor', 'none', 'BarWidth', 0.6);
end

errorbar(1:numel(tuningTypesPlot_v1), meanPct, semPct, 'k.', 'LineWidth', 1, 'CapSize', 4);

for k = 1:numel(tuningTypesPlot_v1)
    xJitter = k + (rand(size(allSessionPcts_v1,1), 1) - 0.5) * 0.25;
    scatter(xJitter, allSessionPcts_v1(:,k), 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5);
end

ax_v1c1.XTick              = 1:numel(tuningTypesPlot_v1);
ax_v1c1.XTickLabel         = labelNames_v1;
ax_v1c1.XTickLabelRotation = 30;
ax_v1c1.FontSize           = 8;
ax_v1c1.FontName           = 'Arial';
ax_v1c1.TickDir            = 'out';
ax_v1c1.LineWidth          = 0.75;
ax_v1c1.Position           = [0.18, 0.28, 0.75, 0.62];
box off;

ylabel('% of sig. boutons', 'FontSize', 9, 'FontName', 'Arial');
title('VISp',               'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('%d rois out of %d (sig + R2 pass) | n=%d sessions', ...
    nClassified, totalROIsCounted, size(allSessionPcts_v1,1)), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');

defaultAxesProperties(ax_v1c1, 0);
offsetAxes(ax_v1c1);
saveFigureFormats(figV1_C1, fullfile(saveFolder_v1, 'category_counts_v1'));

%% Figure 2 — Band-pass preferred speed histogram
prefSpeedBins_v1 = [0, 10, 20, 30, 40, 50, 60];
binCentres_v1    = [5, 15, 25, 35, 45, 55];
bp_speeds_v1     = allPrefSpeeds(strcmp(allTypes, 'bandpass'));
histBP_v1        = histcounts(bp_speeds_v1, prefSpeedBins_v1);

figV1_C2 = figure('Position', [100, 100, 380, 360], 'Color', 'w');
ax_v1c2  = gca;
hold on;

bar(binCentres_v1, histBP_v1, 'FaceColor', colors.bandpass, 'EdgeColor', 'none', 'BarWidth', 0.5);

xlim([prefSpeedBins_v1(1), prefSpeedBins_v1(end)]);
ax_v1c2.XTick      = binCentres_v1;
ax_v1c2.XTickLabel = {'5','15','25','35','45','55'};
ylim([0, max(histBP_v1) * 1.2]);
ax_v1c2.FontSize   = 8;
ax_v1c2.FontName   = 'Arial';
ax_v1c2.TickDir    = 'out';
ax_v1c2.LineWidth  = 0.75;
ax_v1c2.Layer      = 'top';
ax_v1c2.Position   = [0.18, 0.22, 0.75, 0.65];
box off;

xlabel('Preferred speed bin centre (cm s^{-1})', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Number of boutons',                       'FontSize', 9, 'FontName', 'Arial');
title('VISp — Band-pass',                         'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('n = %d boutons', numel(bp_speeds_v1)), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');

defaultAxesProperties(ax_v1c2, 0);
offsetAxes(ax_v1c2);
saveFigureFormats(figV1_C2, fullfile(saveFolder_v1, 'preferred_speed_bandpass_v1'));

%% Figure 3 — Stacked bar all categories
histMatrix_v1 = zeros(numel(prefSpeedBins_v1)-1, numel(tuningTypes_v1));
for t = 1:numel(tuningTypes_v1)
    speeds = allPrefSpeeds(strcmp(allTypes, tuningTypes_v1{t}));
    histMatrix_v1(:, t) = histcounts(speeds, prefSpeedBins_v1);
end

figV1_C3 = figure('Position', [100, 100, 480, 360], 'Color', 'w');
ax_v1c3  = gca;
hold on;

faceColors_stacked = {colors.lowpass, colors.trough_inverted, colors.bandpass, colors.highpass};
bHandle_v1 = bar(binCentres_v1, histMatrix_v1, 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.5);
for t = 1:numel(tuningTypes_v1)
    bHandle_v1(t).FaceColor = faceColors_stacked{t};
end

xlim([prefSpeedBins_v1(1), prefSpeedBins_v1(end)]);
ax_v1c3.XTick      = binCentres_v1;
ax_v1c3.XTickLabel = {'5','15','25','35','45','55'};
yMax_v1 = max(sum(histMatrix_v1, 2)) * 1.18;
ylim([0, yMax_v1]);
ax_v1c3.FontSize   = 8;
ax_v1c3.FontName   = 'Arial';
ax_v1c3.TickDir    = 'out';
ax_v1c3.LineWidth  = 0.75;
ax_v1c3.Layer      = 'top';
ax_v1c3.Position   = [0.13, 0.22, 0.83, 0.65];
box off;

xlabel('Speed at peak tuning (cm s^{-1})', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Number of boutons',                'FontSize', 9, 'FontName', 'Arial');
title('VISp',                              'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('%d rois out of %d (sig + R2 pass)', nClassified, totalROIsCounted), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');

dummyHandles_v1 = [
    patch(NaN, NaN, colors.lowpass,         'EdgeColor', 'none'), ...
    patch(NaN, NaN, colors.trough_inverted, 'EdgeColor', 'none'), ...
    patch(NaN, NaN, colors.bandpass,        'EdgeColor', 'none'), ...
    patch(NaN, NaN, colors.highpass,        'EdgeColor', 'none')
];
leg_v1 = legend(dummyHandles_v1, {'Low-pass','Trough-inv','Band-pass','High-pass'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 7, 'FontName', 'Arial');
leg_v1.ItemTokenSize = [10, 8];

counts_v1 = zeros(1, numel(tuningTypes_v1));
for k = 1:numel(tuningTypes_v1)
    counts_v1(k) = sum(strcmp(allTypes, tuningTypes_v1{k}));
end
annoStr_v1 = sprintf('Low-pass:    %d (%.1f%%)\nTrough-inv:  %d (%.1f%%)\nBand-pass:   %d (%.1f%%)\nHigh-pass:   %d (%.1f%%)', ...
    counts_v1(1), 100*counts_v1(1)/nClassified, ...
    counts_v1(2), 100*counts_v1(2)/nClassified, ...
    counts_v1(3), 100*counts_v1(3)/nClassified, ...
    counts_v1(4), 100*counts_v1(4)/nClassified);
text(ax_v1c3, prefSpeedBins_v1(end-1)*0.98, yMax_v1*0.95, annoStr_v1, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 7, 'FontName', 'Arial', 'Color', [0.3 0.3 0.3]);

defaultAxesProperties(ax_v1c3, 0);
offsetAxes(ax_v1c3);
saveFigureFormats(figV1_C3, fullfile(saveFolder_v1, 'preferred_speed_histogram_v1'));
%%
pdfSavePath = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.5\example_session\all_classified_tuningCurves_somadff.pdf';

% Delete existing file first or exportgraphics will append to old content
if isfile(pdfSavePath)
    delete(pdfSavePath);
    fprintf('Deleted existing QC PDF.\n');
end

pageCount = 0;

for thisMouse = 1:size(mouseInfo, 1)
    mousenumber  = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};

    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if isempty(infoPath) || ~isfile(infoPath), continue; end

        try
            loadedInfo      = load(infoPath, 'sessionFileInfo');
            sessionFileInfo = loadedInfo.sessionFileInfo;
            stimNames       = {sessionFileInfo.stimFiles.name};
            targetIdx       = find(contains(stimNames, {'GrayScreen'}));

            for thisStim = 1:length(targetIdx)
                thisStimName = stimNames{targetIdx(thisStim)};
                stimFileName = sprintf('%s_%s_Response_%s.mat', mousenumber, sessionName, thisStimName);
                fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
                if ~isfile(fileFullPath), continue; end

                loadedData = load(fileFullPath, 'response');
                resp       = loadedData.response;

                if ~isfield(resp, targetStruct) || ...
                   ~isfield(resp.(targetStruct), useField) || ...
                   ~isfield(resp.(targetStruct).(useField), 'classification')
                    continue;
                end

                cls        = resp.(targetStruct).(useField).classification;
                pvalMoving = resp.(targetStruct).(useField).pValMoving;

                classifiedIdx = find(cls.R2 >= minR2Threshold & pvalMoving <= pval_shuffleThreshold);
                if isempty(classifiedIdx), continue; end

                edges         = resp.(targetStruct).speedBins;
                movingCenters = (edges(1:end-1) + diff(edges)/2)';
                yMean         = resp.(targetStruct).(useField).moveMean;

                typeColorMap = containers.Map(...
                    {'lowpass','trough_inverted','bandpass','highpass'}, ...
                    {colors.lowpass, colors.trough_inverted, colors.bandpass, colors.highpass});

                nPlot = numel(classifiedIdx);
                nCols = 8;
                nRows = ceil(nPlot / nCols);

                fig = figure('Color', 'w', ...
                             'Position', [100 100 1600 max(nRows*180, 300)], ...
                             'Visible', 'off');

                for k = 1:nPlot
                    r   = classifiedIdx(k);
                    typ = cls.tuningType{r};

                    if isKey(typeColorMap, typ)
                        c = typeColorMap(typ);
                    else
                        c = [0.5 0.5 0.5];
                    end

                    ax_k = subplot(nRows, nCols, k);
                    plot(movingCenters, yMean(r,:), 'o-', ...
                        'Color', c, 'MarkerSize', 3, ...
                        'MarkerFaceColor', c, 'LineWidth', 1);
                    title(sprintf('%s %d | R²=%.2f', typ, r, cls.R2(r)), ...
                        'FontSize', 5.5, 'Color', c, 'Interpreter', 'none');
                    box off;
                    set(ax_k, 'TickDir', 'out', 'FontSize', 5);
                end

                sgtitle(sprintf('%s  %s  —  %d / %d ROIs classified', ...
                    mousenumber, sessionName, nPlot, numel(pvalMoving)), ...
                    'FontSize', 10, 'FontWeight', 'bold');

                exportgraphics(fig, pdfSavePath, 'Append', true, 'Resolution', 120);
                close(fig);
                pageCount = pageCount + 1;
                fprintf('  Page %d written: %s %s (%d ROIs)\n', pageCount, mousenumber, sessionName, nPlot);
            end

        catch ME
            fprintf('ERROR on %s-%s: %s\n', mousenumber, sessionName, ME.message);
        end
    end
end

fprintf('\nDone. QC PDF saved (%d pages):\n  %s\n', pageCount, pdfSavePath);