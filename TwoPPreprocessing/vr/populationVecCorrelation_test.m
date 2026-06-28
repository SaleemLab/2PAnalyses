%% RSP Population Vector Correlation - Full Track
PVcorr_allSess = {};

for s = 1:length(VISpData)
    sess = VISpData(s);
    
    if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
       ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
        continue;
    end
    
    lapActivity = sess.ConditionData.Baseline.LapActivity;
    [numROIsTotal, numLaps, numPosBins] = size(lapActivity);
    
    if numPosBins < 200
        fprintf('Session %d skipped - only %d bins\n', s, numPosBins);
        continue;
    end
    
    % smooth per lap
    w = gausswin(15); w = w/sum(w);
    smoothed = lapActivity;
    for iCell = 1:numROIsTotal
        for iLap = 1:numLaps
            tr = squeeze(lapActivity(iCell, iLap, :));
            if all(isnan(tr)), continue; end
            nm = isnan(tr); tr(nm) = 0;
            sm = filtfilt(w, 1, tr);
            sm(nm) = NaN;
            smoothed(iCell, iLap, :) = sm;
        end
    end
    
    % mean across all laps for filtered ROIs
    roiActivity = smoothed(sess.FilteredROIs, :, :);
    meanT = squeeze(mean(roiActivity, 2, 'omitnan'));  % nROIs x 200
    if length(sess.FilteredROIs) == 1, meanT = meanT'; end
    
    if size(meanT, 1) < 5
        fprintf('Session %d skipped - too few cells (%d)\n', s, size(meanT,1));
        continue;
    end
    
    % smooth each cell's tuning curve
    w2 = gausswin(21); w2 = w2/sum(w2);
    meanT_smooth = zeros(size(meanT));
    for iCell = 1:size(meanT, 1)
        tr = meanT(iCell, :);
        if all(isnan(tr)), continue; end
        nm = isnan(tr); tr(nm) = 0;
        sm = filtfilt(w2, 1, tr);
        sm(nm) = NaN;
        meanT_smooth(iCell, :) = sm;
    end
    
    % compute PV correlation for this session (200 x 200)
    PVcorr_sess = corr(meanT_smooth, 'rows', 'pairwise');
    PVcorr_allSess{end+1} = PVcorr_sess;
    
    fprintf('Session %d: %d cells included\n', s, size(meanT,1));
end

fprintf('\nComputed PV correlation for %d sessions\n', length(PVcorr_allSess));

%% --- Average across sessions ---
PVcorr_stack = cat(3, PVcorr_allSess{:});  % 200 x 200 x nSessions
meanPVcorr = mean(PVcorr_stack, 3, 'omitnan');
meanPVcorr_smooth = imgaussfilt(meanPVcorr, 3);

%% --- Plot ---
figure('Color', 'w', 'Position', [100 100 500 450]);
ax = axes('Position', [0.15 0.15 0.65 0.75]);
imagesc(1:200, 1:200, meanPVcorr_smooth);
colormap(ax, flipud(gray));
caxis([0 1]);

cb = colorbar;
cb.Label.String = 'Population vector correlation';
cb.TickDirection = 'out';
cb.Box = 'off';
cb.FontName = 'Arial';
cb.FontSize = 10;

hold on;
xline(40,  'r--', 'LineWidth', 1.5);
xline(80,  'b--', 'LineWidth', 1.5);
xline(120, 'r--', 'LineWidth', 1.5);
xline(160, 'b--', 'LineWidth', 1.5);
yline(40,  'r--', 'LineWidth', 1.5);
yline(80,  'b--', 'LineWidth', 1.5);
yline(120, 'r--', 'LineWidth', 1.5);
yline(160, 'b--', 'LineWidth', 1.5);

xlabel('Position (cm)', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Position (cm)', 'FontName', 'Arial', 'FontSize', 12);
title(sprintf('RSP mean PV correlation (%d sessions)', length(PVcorr_allSess)), ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');
axis square;
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 11);

%% --- Save ---
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4\PVcorrelation';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(gcf, fullfile(outputDir, 'RSP_PVcorrelation_fulltrack'));