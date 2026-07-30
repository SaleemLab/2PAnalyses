pairs=struct;
pairs.M25132 = ['20260219','20260220','20260221','20260223', '20260226','20260228', '20260313']; %
pairs.M25133 = ['20260220','20260221','20260223''20260219','20260224'];
pairs.M26003 = ['20260316','20260317', '20260320','20260321','20260322', '20260324', '20260325']; %


RSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
RSPDataAcrossDays = getTuningDataByCondition(RSPSessions);

RSPData = appendFilteredROIs(RSPDataAcrossDays,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);
daysToPlot      = [1 2 3 4 5 200]; 
daysToPlot(daysToPlot == 200) = 5;
daysToPlot = unique(daysToPlot, 'stable');
nDays           = numel(daysToPlot);
trueEVThreshold = 0.1;   % same as cvExpvarThreshold you used in appendFilteredROIs
nullThreshold   = 0.01;  % same as ExpVarSigThreshold you used in appendFilteredROIs

%% colours
dayColorPalette = {'r','m', [0.4 0.7 0.2],'k', 'b'};
dayColors = cell(1, nDays);
for d = 1:nDays
    dayColors{d} = dayColorPalette{daysToPlot(d)};
end
%% Behaviour 
for s = 1:length(RSPDataAcrossDays)
    nLaps = RSPDataAcrossDays(s).ConditionData.Baseline.NumLaps;  % or whatever your lap field is
    fprintf('Day %d: %d laps | Mouse %s \n', RSPDataAcrossDays(s).Day, nLaps, RSPDataAcrossDays(s).MouseID);
end

%% Pool NumLaps by Day, then plot laps across days 
lapsByDay      = cell(1, nDays);  % per-session lap counts for each day
meanLapsByDay  = nan(1, nDays);
for d = 1:nDays
    thisDay = daysToPlot(d);
    if thisDay == 5
        dayIdx = find([RSPDataAcrossDays.Day] == 5 | [RSPDataAcrossDays.Day] == 200);
    else
        dayIdx = find([RSPDataAcrossDays.Day] == thisDay);
    end

    theseLaps = [];
    for s = dayIdx
        nLaps = RSPDataAcrossDays(s).ConditionData.Baseline.NumLaps;
        theseLaps(end+1) = nLaps;
        fprintf('Day %d: %d laps | Mouse %s \n', thisDay, nLaps, RSPDataAcrossDays(s).MouseID);
    end
    lapsByDay{d}     = theseLaps;
    meanLapsByDay(d) = mean(theseLaps);
end

% Plot: laps across days ---
fig5 = figure('Color', 'w', 'Position', [200, 200, 500, 400]);
hold on;

% individual mouse/session points, jittered, colored by day
for d = 1:nDays
    if isempty(lapsByDay{d}); continue; end
    jitter = (rand(1, numel(lapsByDay{d})) - 0.5) * 0.15;
    scatter(daysToPlot(d) + jitter, lapsByDay{d}, 30, dayColors{d}, 'filled', ...
        'MarkerFaceAlpha', 0.6);
end

plot(daysToPlot, meanLapsByDay, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
for d = 1:nDays
    if isnan(meanLapsByDay(d)); continue; end
    plot(daysToPlot(d), meanLapsByDay(d), 'o', 'Color', dayColors{d}, ...
        'MarkerFaceColor', dayColors{d}, 'MarkerSize', 8, 'LineWidth', 1.2);
end

xlabel('Day of Experience');
ylabel('Number of laps');
xticks(daysToPlot);
xlim([min(daysToPlot) - 0.5, max(daysToPlot) + 0.5]);
defaultAxesProperties(gca, true);
hold off;

saveFigureFormats(fig5, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.1_Learning\allVersions\lapsAcrossDays_3mice_12345200');

%% plot running for the 3 days 
respDay1 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260219\M25132_20260219_Response_M25132_BaselineCorridor_20260219_00002.mat", 'lapRunningSpeed', 'speedPositionActivity', 'completedLaps_AbsoluteIdx');
respDay3 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260221\M25132_20260221_Response_M25132_BaselineCorridor_20260221_CombinedRuns.mat", 'lapRunningSpeed', 'speedPositionActivity');
respDay5 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_BaselineCorridor_20260226_CombinedRuns.mat", 'lapRunningSpeed', 'speedPositionActivity');
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.1_Learning\running_trejectories\';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end


figA1 = plotSpeedTrajectoriesSimple(respDay1);
set(figA1, 'Visible', 'on');
saveFigureFormats(figA1, fullfile(outputDir, 'runBeh_m25132_day1'));


figA2 = plotSpeedTrajectoriesSimple(respDay3);
set(figA2, 'Visible', 'on');
saveFigureFormats(figA2, fullfile(outputDir, 'runBeh_m25132_day3'));

figA3 = plotSpeedTrajectoriesSimple(respDay5);
set(figA3, 'Visible', 'on');
saveFigureFormats(figA3, fullfile(outputDir, 'runBeh_m25132_day5'));



%%
% %downsampleCVR2_ForLearning
% 
% 
% for s = 1:nSessions
%     RSPDataAcrossDays(s).cvExpVar_downsampled.meanExpVar     = downsampledResults(s).meanExpVar;
%     RSPDataAcrossDays(s).cvExpVar_downsampled.pValues        = downsampledResults(s).pValues;
%     RSPDataAcrossDays(s).cvExpVar_downsampled.targetNumLaps  = targetNumLaps;
%     RSPDataAcrossDays(s).cvExpVar_downsampled.nFolds         = nFolds;
%     RSPDataAcrossDays(s).cvExpVar_downsampled.nShuffles      = nShuffles;
%     RSPDataAcrossDays(s).cvExpVar_downsampled.nRepeats       = nRepeats;
% end
% 
% savePath = 'Z:\ibn-vision\USERS\Sonali\Figures\downsampled_cvev\RSPDataAcrossDays_withDownsampledCV.mat';
% save(savePath, 'RSPDataAcrossDays', '-v7.3');
% fprintf('Saved RSPDataAcrossDays (with downsampled cvEV) to %s\n', savePath);

%% Pool DOWNSAMPLED true EV by day (same day grouping as before)
% pooledDownsampledEV_byDay = cell(1, nDays);
% 
% for d = 1:nDays
%     thisDay = daysToPlot(d);
%     dayIdx  = find([RSPDataAcrossDays.Day] == thisDay);
% 
%     theseVals = [];
%     for s = dayIdx
%         if isfield(RSPDataAcrossDays(s), 'cvExpVar_downsampled') && ...
%                 ~isempty(RSPDataAcrossDays(s).cvExpVar_downsampled.meanExpVar)
%             theseVals = [theseVals; RSPDataAcrossDays(s).cvExpVar_downsampled.meanExpVar(:)]; %#ok<AGROW>
%         end
%     end
%     pooledDownsampledEV_byDay{d} = theseVals;
% end
%% Pool cvExpVar data by Day, then plot % ROIs passing threshold across days,

pctPassingByDay      = nan(1, nDays);   % pooled (all ROIs, all mice) % passing per day
mousePctByDay        = cell(1, nDays);  % per-session (per mouse) % passing, for scatter/error bars
pooledTrueEV_byDay    = cell(1, nDays);
pooledNullEV_byDay    = cell(1, nDays);
for d = 1:nDays
    thisDay = daysToPlot(d);

    % Day 200 is an alias for Day 5 — pool its sessions in with day 5's.
    if thisDay == 5
        dayIdx = find([RSPData.Day] == 5 | [RSPData.Day] == 200);
    else
        dayIdx = find([RSPData.Day] == thisDay);
    end

    pooled_true_EV          = [];
    pooled_null_EV           = [];
    pooled_null_max_per_ROI  = [];
    mousePct                 = [];
    for iSess = dayIdx
        current_cv = RSPData(iSess).cvExpVar;
        if isempty(current_cv) || isempty(current_cv.meanExpVar)
            continue;
        end
        true_means  = current_cv.meanExpVar;
        null_matrix = current_cv.cvExpVarNull;
        if isempty(true_means) || isempty(null_matrix)
            continue;
        end
        nROIs    = length(true_means);
        null_dims = size(null_matrix);
        roi_dim   = find(null_dims == nROIs, 1, 'first');
        if isempty(roi_dim)
            error('Session %d: Could not find a dimension matching %d ROIs.', iSess, nROIs);
        end
        dims_to_collapse = setdiff(1:length(null_dims), roi_dim);
        if isempty(dims_to_collapse)
            max_null_per_ROI = null_matrix;
        else
            max_null_per_ROI = max(null_matrix, [], dims_to_collapse);
        end
        pooled_true_EV         = [pooled_true_EV;  true_means(:)];
        pooled_null_EV          = [pooled_null_EV;  null_matrix(:)];
        pooled_null_max_per_ROI = [pooled_null_max_per_ROI; max_null_per_ROI(:)];
        % per-session (mouse) passing percentage for this day
        sessPassIdx = (true_means(:) > trueEVThreshold) & (max_null_per_ROI(:) > nullThreshold);
        mousePct(end+1) = 100 * sum(sessPassIdx) / numel(true_means);
    end
    if isempty(pooled_true_EV)
        continue;
    end
    passing_idx = (pooled_true_EV > trueEVThreshold) & (pooled_null_max_per_ROI > nullThreshold);
    pctPassingByDay(d) = 100 * sum(passing_idx) / numel(pooled_true_EV);
    mousePctByDay{d}   = mousePct;
    pooledTrueEV_byDay{d} = pooled_true_EV;
    pooledNullEV_byDay{d}  = pooled_null_EV;
end

%% EV: Plot 1: % ROIs passing threshold across days 
fig2 = figure('Color', 'w', 'Position', [200, 200, 500, 400]);
hold on;

% individual mouse/session points, jittered, colored by day, behind the mean line
for d = 1:nDays
    if isempty(mousePctByDay{d}); continue; end
    jitter = (rand(1, numel(mousePctByDay{d})) - 0.5) * 0.15;
    scatter(daysToPlot(d) + jitter, mousePctByDay{d}, 30, dayColors{d}, 'filled', ...
        'MarkerFaceAlpha', 0.6);
end

plot(daysToPlot, pctPassingByDay, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
for d = 1:nDays
    if isnan(pctPassingByDay(d)); continue; end
    plot(daysToPlot(d), pctPassingByDay(d), 'o', 'Color', dayColors{d}, ...
        'MarkerFaceColor', dayColors{d}, 'MarkerSize', 8, 'LineWidth', 1.2);
end

xlabel('Day of experience');
ylabel('ROIs passing EV threshold (%)');
xticks(daysToPlot);
xlim([min(daysToPlot) - 0.5, max(daysToPlot) + 0.5]);
defaultAxesProperties(gca, true);
hold off;

saveFigureFormats(fig2, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.1_Learning\allVersions\pctPassingAcrossDays_3mice_12345200');

%% EV Plot 2: Superimposed CDFs of true EV across days ---
fig3 = figure('Color', 'w', 'Position', [200, 200, 600, 450]);
hold on;

hCDF      = gobjects(1, nDays);
dayLabels = cell(1, nDays);

for d = 1:nDays
    if isempty(pooledTrueEV_byDay{d}); continue; end
    [f, x]   = ecdf(pooledTrueEV_byDay{d});
    hCDF(d)  = plot(x, f, 'LineWidth', 1.8, 'Color', dayColors{d});
    dayLabels{d} = sprintf('Day %d', daysToPlot(d));
end

xline(trueEVThreshold, 'k--', 'LineWidth', 1.2);
xlabel('Explained variance (R^2)');
ylabel('Cumulative probability');
xlim([-0.8, 1.0]);
validDays = isgraphics(hCDF);
legend(hCDF(validDays), dayLabels(validDays), 'Location', 'SouthEast');
legend boxoff;
defaultAxesProperties(gca, true);
hold off;

saveFigureFormats(fig3, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.1_Learning\allVersions\pctPassingAcrossDays_3mice_12345200_cdf');

%% EV Plot 3: Superimposed EV distributions (histograms) across days
fig4 = figure('Color', 'w', 'Position', [200, 200, 600, 450]);
hold on;

bin_width = 0.02;
bin_edges = -0.8:bin_width:1.0;

hDist = gobjects(1, nDays);

for d = 1:nDays
    if isempty(pooledTrueEV_byDay{d}); continue; end
    hDist(d) = histogram(pooledTrueEV_byDay{d}, 'BinEdges', bin_edges, ...
        'Normalization', 'probability', 'DisplayStyle', 'bar', ...
        'FaceColor', dayColors{d}, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
end

xline(trueEVThreshold, 'k--', 'LineWidth', 1.2);
xlabel('Explained variance (R^2)');
ylabel('Probability');
xlim([-0.8, 1.0]);
xticks([-0.8, -0.5, 0, 0.5, 1.0]);
validDist = isgraphics(hDist);
legend(hDist(validDist), dayLabels(validDist), 'Location', 'NorthEast');
legend boxoff;
defaultAxesProperties(gca, true);
hold off;

% saveFigureFormats(fig4, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.2\EVHistogram\evDistAcrossDays');


%%
%% Compute mean EV per day: original vs downsampled, then plot as lines across days

meanEV_orig_byDay = nan(1, nDays);
meanEV_down_byDay = nan(1, nDays);

for d = 1:nDays
    if ~isempty(pooledTrueEV_byDay{d})
        meanEV_orig_byDay(d) = mean(pooledTrueEV_byDay{d}, 'omitnan');
    end
    if ~isempty(pooledDownsampledEV_byDay{d})
        meanEV_down_byDay(d) = mean(pooledDownsampledEV_byDay{d}, 'omitnan');
    end
end

%% --- Plot: original vs downsampled EV across days ---
fig7 = figure('Color', 'w', 'Position', [200, 200, 500, 400]);
hold on;

% connecting lines (neutral, different style per condition)
plot(daysToPlot, meanEV_orig_byDay, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
plot(daysToPlot, meanEV_down_byDay, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

% markers colored per day: filled circle = original, open square = downsampled
for d = 1:nDays
    if ~isnan(meanEV_orig_byDay(d))
        plot(daysToPlot(d), meanEV_orig_byDay(d), 'o', 'Color', dayColors{d}, ...
            'MarkerFaceColor', dayColors{d}, 'MarkerSize', 8, 'LineWidth', 1.2);
    end
    if ~isnan(meanEV_down_byDay(d))
        plot(daysToPlot(d), meanEV_down_byDay(d), 's', 'Color', dayColors{d}, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
end

xlabel('Day of learning');
ylabel('Mean explained variance (R^2)');
xticks(daysToPlot);
xlim([min(daysToPlot) - 0.5, max(daysToPlot) + 0.5]);

% custom legend entries for line style (condition), separate from day colors
h1 = plot(nan, nan, '-ok', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
h2 = plot(nan, nan, '--sk', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
legend([h1, h2], {'Original', 'Downsampled'}, 'Location', 'best');
legend boxoff;

defaultAxesProperties(gca, true);
hold off;

% saveFigureFormats(fig7, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.2\EVHistogram\evOrigVsDownsampledLine');

%% tuning curve
plotPooledPopulation_DayWise_Even(RSPDataAcrossDays, ...
    'RSP', ...
    'DaysToPlot',[1 2 3 4 5 200], ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.1_Learning\allVersions\tuningcurves_acrossdays_3mice_12345200');



%% smi
% RSPDataAcrossDays = getTuningDataByCondition(RSPSessions);
pooledRSP = poolSMIAcrossSessions(RSPDataAcrossDays);
fig5=plotPooledSMI_CDF(pooledRSP, 'RSP', 'Days', [1 3 5]);
defaultAxesProperties(gca, true);

saveFigureFormats(fig5, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.1_Learning\allVersions\smi_3mice_135200_cdf');
