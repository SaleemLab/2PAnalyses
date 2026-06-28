function fig = plotSpatialAutocorrelation_Pooled(RSPData, V1Data)
% Computes and plots population mean spatial autocorrelation for RSP and V1
% overlaid on the same figure
%
% Input:
%   RSPData - struct array with fields ConditionData.Baseline.LapActivity
%             and FilteredROIs
%   V1Data  - same format as RSPData

fprintf('Computing spatial autocorrelation for RSP and V1...\n');

binRange = 30:170;  % exclude unreliable edges
[~, lagsFull] = xcorr(zeros(1, length(binRange)), 'normalized');

%% --- Helper function to compute ACF for a region ---
function allACF = computeACF(RegionData)
    allACF = [];
    w_space = gausswin(15); w_space = w_space / sum(w_space);
    
    for s = 1:length(RegionData)
        sess = RegionData(s);
        
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
           ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end
        
        lapActivity = sess.ConditionData.Baseline.LapActivity;
        [numROIsTotal, numLaps, ~] = size(lapActivity);
        
        % smooth
        smoothedActivity = lapActivity;
        for iCell = 1:numROIsTotal
            for iLap = 1:numLaps
                trace = squeeze(lapActivity(iCell, iLap, :));
                if all(isnan(trace)), continue; end
                nanMask = isnan(trace); trace(nanMask) = 0;
                smoothed = filtfilt(w_space, 1, trace);
                smoothed(nanMask) = NaN;
                smoothedActivity(iCell, iLap, :) = smoothed;
            end
        end
        
        roisToAnalyze = sess.FilteredROIs;
        roiActivity   = smoothedActivity(roisToAnalyze, :, :);
        numROIs       = length(roisToAnalyze);
        if numROIs == 0, continue; end
        
        % mean tuning curve across all laps
        meanTuning = squeeze(mean(roiActivity, 2, 'omitnan'));  % numROIs x numBins
        if numROIs == 1, meanTuning = meanTuning'; end
        
        % compute ACF for each cell
        sessACF = NaN(numROIs, length(lagsFull));
        for i = 1:numROIs
            trace = meanTuning(i, binRange);
            if all(isnan(trace)), continue; end
            demeaned = trace - mean(trace, 'omitnan');
            sessACF(i, :) = xcorr(demeaned, 'normalized');
        end
        
        allACF = [allACF; sessACF];
    end
end

%% --- Compute for both regions ---
ACF_RSP = computeACF(RSPData);
ACF_V1  = computeACF(V1Data);

fprintf('RSP: %d cells\n', size(ACF_RSP, 1));
fprintf('V1:  %d cells\n', size(ACF_V1, 1));

%% --- Compute population mean and SEM ---
meanACF_RSP = mean(ACF_RSP, 1, 'omitnan');
semACF_RSP  = std(ACF_RSP, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(ACF_RSP), 1));

meanACF_V1  = mean(ACF_V1, 1, 'omitnan');
semACF_V1   = std(ACF_V1, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(ACF_V1), 1));

% clip negative values to 0
% meanACF_RSP = max(0, meanACF_RSP);
% meanACF_V1  = max(0, meanACF_V1);

%% --- Statistical tests at landmark lags ---
lags = lagsFull;
lag40idx  = lags == 40;
lag80idx  = lags == 80;
lag120idx = lags == 120;

[p40_RSP,  ~] = signrank(ACF_RSP(:, lag40idx));
[p80_RSP,  ~] = signrank(ACF_RSP(:, lag80idx));
[p40_V1,   ~] = signrank(ACF_V1(:, lag40idx));
[p80_V1,   ~] = signrank(ACF_V1(:, lag80idx));

fprintf('\n=== RSP ===\n');
fprintf('ACF lag 40: median=%.3f, p=%.4f\n', median(ACF_RSP(:,lag40idx),'omitnan'), p40_RSP);
fprintf('ACF lag 80: median=%.3f, p=%.4f\n', median(ACF_RSP(:,lag80idx),'omitnan'), p80_RSP);
fprintf('\n=== V1 ===\n');
fprintf('ACF lag 40: median=%.3f, p=%.4f\n', median(ACF_V1(:,lag40idx),'omitnan'), p40_V1);
fprintf('ACF lag 80: median=%.3f, p=%.4f\n', median(ACF_V1(:,lag80idx),'omitnan'), p80_V1);

%% --- Plot ---
fig = figure('Color', 'w', 'Position', [100 100 500 350]);
ax = axes('Position', [0.15 0.18 0.75 0.72]);
hold(ax, 'on');

% SEM shading

% fix legend - name the fill patches as empty string to hide from legend
fill(ax, [lags, fliplr(lags)], ...
    [meanACF_RSP + semACF_RSP, fliplr(meanACF_RSP - semACF_RSP)], ...
    [0 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

fill(ax, [lags, fliplr(lags)], ...
    [meanACF_V1 + semACF_V1, fliplr(meanACF_V1 - semACF_V1)], ...
    [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% mean lines
plot(ax, lags, meanACF_RSP, 'k-', 'LineWidth', 2, 'DisplayName', sprintf('RSP (n=%d)', size(ACF_RSP,1)));
plot(ax, lags, meanACF_V1, 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'DisplayName', sprintf('V1 (n=%d)', size(ACF_V1,1)));

% landmark lines
% xline(ax, 40,  'r--', 'LineWidth', 1.2);
% xline(ax, 80,  'b--', 'LineWidth', 1.2);
% xline(ax, 120, 'r--', 'LineWidth', 1.2);
% xline(ax, 0,   'k-',  'LineWidth', 1);

% only show positive lags - clip at 0
xlim(ax, [0 140]);
ylim(ax, [-0.5 max([meanACF_RSP, meanACF_V1])*1.2]);

xlabel(ax, 'Lag (cm)', 'FontName', 'Arial', 'FontSize', 12);
ylabel(ax, 'Spatial autocorrelation', 'FontName', 'Arial', 'FontSize', 12);
title(ax, 'Population spatial autocorrelation', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');

legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 10);

set(ax, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 11);

%% --- Save ---
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4\SpatialAutocorrelation';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(fig, fullfile(outputDir, 'SpatialAutocorrelation_RSP_V1'));

end