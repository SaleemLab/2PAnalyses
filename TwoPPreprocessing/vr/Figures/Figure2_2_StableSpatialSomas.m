% Figure3.2: Responsive boutons from day 5 onwards
% use all sessions even ones where the model didnt analyse? all spikes 
pairs=struct;
% pairs.M26005 = [ '20260318', '20260321', '20260322']; % unique fovs  '20260305', '20260306',
% pairs.M26004 = [ '20260318', '20260321', '20260322']; % unique fovs '20260305', '20260307', '20260312', '20260313', '20260314',
% pairs.M25131 = ['20260318', '20260321', '20260322']; % unique fovs  these sessions were excluded from juline's analyses['20260312', '20260313', '20260314', 
% pairs.M25126 = ['20260311', '20260312', '20260313']; % unique fovs 

pairs.M26005 = ['20260305', '20260306', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M26004 = ['20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs
pairs.M25131 = ['20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M25126 = ['20260311', '20260312', '20260313']; % unique fovs 
 


VISpSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
% VISpData = getTuningDataByCondition(VISpSessions);
VISpSpksData = getTuningDataByCondition(VISpSessions, 'signalToUse', 'spks'); 
% bin data by condition and load inclusion critera
trueEVThreshold = 0.1;
dataToUse = VISpSpksData;
%% ev

pooled_true_EV = [];
pooled_null_EV = [];
pooled_null_max_per_ROI = []; % Track the highest null value an ROI achieved across shuffles

num_sessions = length(dataToUse);

% Loop through each session and pool the ROIs 
for iSess = 1:num_sessions
    current_cv = dataToUse(iSess).cvExpVar;
    true_means = current_cv.meanExpVar; % Vector of length = nROIs
    null_matrix = current_cv.cvExpVarNull; 
    
    if isempty(true_means) || isempty(null_matrix)
        continue; 
    end
    
    nROIs = length(true_means);
    null_dims = size(null_matrix);
    roi_dim = find(null_dims == nROIs, 1, 'first');
    
    if isempty(roi_dim)
        error('Session %d: Could not find a dimension matching %d ROIs.', iSess, nROIs);
    end
    
    dims_to_collapse = setdiff(1:length(null_dims), roi_dim);
    
    % Instead of the MEAN, let's find the MAX null value this ROI hit during shuffles
    if isempty(dims_to_collapse)
        max_null_per_ROI = null_matrix; 
    else
        max_null_per_ROI = max(null_matrix, [], dims_to_collapse);
    end
    
    % Force to column vectors and pool
    pooled_true_EV = [pooled_true_EV; true_means(:)];
    pooled_null_EV = [pooled_null_EV; null_matrix(:)]; 
    pooled_null_max_per_ROI = [pooled_null_max_per_ROI; max_null_per_ROI(:)];
end

% Keep the threshold calculation logic exactly the same for annotations
passing_idx = (pooled_true_EV > 0.1) & (pooled_null_max_per_ROI > 0.01);

% --- Plotting ---
fig1 = figure('Color', 'w', 'Position', [200, 200, 600, 450]);
hold on;

bin_width = 0.02; 
bin_edges = -0.8:bin_width:1.0; 

% ev null
h_null = histogram(pooled_null_EV, 'BinEdges', bin_edges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
    'EdgeColor', [0.55, 0.55, 0.55], 'LineWidth', 1.5);

xlabel('Explained variance (R^2)');
ylabel('Probability');

% ev true
h_true = histogram(pooled_true_EV, 'BinEdges', bin_edges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
    'EdgeColor', '#008080', 'LineWidth', 1.5);


pooled_99th = prctile(pooled_null_EV, 99);
xline(pooled_99th, 'r--', 'LineWidth', 1.5);
xline(trueEVThreshold, 'k--', 'LineWidth', 1.5);

y_limits = ylim;
text_y_position = y_limits(2) * 0.9; 

text(pooled_99th + 0.03, text_y_position, 'Null 99%', ...
     'Color', 'r', 'FontSize', 9, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% Calculate passing percentages
num_passing = sum(passing_idx);
total_boutons = length(pooled_true_EV);
pct_passing = (num_passing / total_boutons) * 100;

text(pooled_99th + 0.05, text_y_position * 0.8, sprintf('Passing: %.1f%%', pct_passing), ...
     'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% Set up Legend for the two outline distributions
lgd = legend([h_null, h_true], {'Shuffled', 'True EV'}, ...
    'Location', 'NorthEast');
defaultAxesProperties(gca, true);
legend boxoff;

xticks(-0.8:0.2:1.0);  
xticks([-0.8, -0.5, 0, 0.5, 1.0]);  
xlim([-0.8, 1.0]);
hold off;

% Save Figure
% saveFigureFormats(fig1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\EVHistogram\evDistributions_spks');
saveFigureFormats(fig1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\EVHistogram\evDistributions_spks');
%% Find  middle and high sig rois from pool that pass the shuffle test as well 
target_R2 = -0.15; % Change this to find an ROI near any value (e.g., -0.1, 0.02, 0.45)

all_ev_values = [];
all_session_idx = [];
all_roi_idx = [];

for iSess = 1:length(dataToUse)
    current_cv = dataToUse(iSess).cvExpVar;
    true_means = current_cv.meanExpVar;
    if isempty(true_means), continue; end
    
    num_rois_in_session = length(true_means);
    
    % Identify which ROIs passed the shuffle test in this session
    passed_shuffle_mask = current_cv.pValues < 0.01;
    
    % Only accumulate data points that are statistically significant
    if any(passed_shuffle_mask)
        sig_indices = find(passed_shuffle_mask);
        
        all_ev_values   = [all_ev_values; true_means(sig_indices)];
        all_session_idx = [all_session_idx; repmat(iSess, length(sig_indices), 1)];
        all_roi_idx     = [all_roi_idx; sig_indices(:)];
    end
end

if isempty(all_ev_values)
    error('No ROIs passed the shuffle criteria across the provided sessions.');
end

[~, match_global_idx] = min(abs(all_ev_values - target_R2));

matched_val     = all_ev_values(match_global_idx);
matched_session = all_session_idx(match_global_idx);
matched_roi     = all_roi_idx(match_global_idx);




fprintf('  Actual R^2 Value: %.4f\n', matched_val);
fprintf('  Session Index:    RSPData(%d) (%s, Day %d)\n', ...
    matched_session, dataToUse(matched_session).MouseID, dataToUse(matched_session).Day);
fprintf('  ROI Number:       %d\n', matched_roi);

%% find target that also fails bootstrapping 
target_R2 = 0.85;
num_matches = 20;

all_ev_values   = [];
all_session_idx = [];
all_roi_idx     = [];

for iSess = 1:length(dataToUse)
    current_cv = dataToUse(iSess).cvExpVar;
    true_means = current_cv.meanExpVar;
    if isempty(true_means), continue; end

    % Mask for ROIs that are strictly less than your target threshold
    target_mask = true_means < target_R2;
    if any(target_mask)
        match_indices = find(target_mask);
        all_ev_values   = [all_ev_values; true_means(match_indices)];
        all_session_idx = [all_session_idx; repmat(iSess, length(match_indices), 1)];
        all_roi_idx     = [all_roi_idx; match_indices(:)];
    end
end

fprintf('=== TOP %d ROIs WITH R^2 LESS THAN %.2f (closest to %.2f) ===\n', num_matches, target_R2, target_R2);

if isempty(all_ev_values)
    fprintf('  No ROIs found with an R^2 value less than %.2f\n', target_R2);
else
    % Sort by closeness to target_R2, ascending distance
    [sorted_dist, sort_idx] = sort(abs(all_ev_values - target_R2));

    n_available = length(sorted_dist);
    n_to_show = min(num_matches, n_available);

    if n_available < num_matches
        fprintf('  Only %d matches found (fewer than requested %d)\n', n_available, num_matches);
    end

    for k = 1:n_to_show
        idx = sort_idx(k);
        matched_val     = all_ev_values(idx);
        matched_session = all_session_idx(idx);
        matched_roi     = all_roi_idx(idx);

        fprintf('%2d. R^2 = %.4f | Session: VISp(%d) (%s, Day %d) | ROI: %d\n', ...
            k, matched_val, matched_session, dataToUse(matched_session).MouseID, ...
            dataToUse(matched_session).Day, matched_roi);
    end
end
%% Plot as heatmaps and save (SELECT EXAMPLES AT THE VERY END)
% plot raw? 
figLow = plotRoiSpatialTuning(dataToUse, 10,344, 'Low_Significant');
saveFigureFormats(figLow, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\EV_HighMidLowROIs\ExampleTuning_Low');

figLowSplot = plotRoiSpatialTuning_oddEvenSplit(dataToUse, 22,53, 'Low_Significant_split')

figMid = plotRoiSpatialTuning(dataToUse, 17, 438, 'Middle_Significant'); %17, 438
saveFigureFormats(figMid, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\EV_HighMidLowROIs\ExampleTuning_Middle');

figHigh = plotRoiSpatialTuning(dataToUse,15 ,11, 'High_Significant');%14, 109
saveFigureFormats(figHigh, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\EV_HighMidLowROIs\ExampleTuning_High');

% single peak 
figSingle = plotRoiSpatialTuning(dataToUse,22 ,53, 'singlePeak');%14, 109
saveFigureFormats(figSingle, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\typesOfResponses\ExampleTuning_singlePeak');

% 

%%
% Loop plotRoiSpatialTuning across all ROIs for a given session/condition
dataToUse=VISpSpksData
sessIdx = 21;
condName = 'Baseline';

% Determine number of ROIs from the LapActivity array (dim 1 = ROIs)
numROIs = size(dataToUse(sessIdx).ConditionData.(condName).LapActivity, 1);

outDir = fullfile('Z:\ibn-vision\DATA\SUBJECTS\M26004\Analysis\20260318\Figures', sprintf('RoiSpatialTuning_Sess%d_%s', sessIdx, condName));
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for roiIdx = 1:numROIs
    try
        figLow = plotRoiSpatialTuning(dataToUse, sessIdx, roiIdx, 'spks');

        % Save each figure, then close it to avoid piling up open windows
        figName = fullfile(outDir, sprintf('ROI_%03d.png', roiIdx));
        exportgraphics(figLow, figName, 'Resolution', 150);
        close(figLow);

    catch ME
        fprintf('Skipping ROI %d (session %d): %s\n', roiIdx, sessIdx, ME.message);
    end
end


%% Laps per session, %tuned, and laps-vs-R2 relationship

combinedEVThreshold = 0.1;
nullMaxThreshold    = 0.01;
pValThreshold        = 0.01;

num_sessions = length(dataToUse);

session_laps       = nan(num_sessions,1);
session_meanR2      = nan(num_sessions,1);
session_medianR2     = nan(num_sessions,1);
session_pctTuned     = nan(num_sessions,1);
session_nROIs        = nan(num_sessions,1);

% for the supplementary ROI-level (pseudoreplicated) scatter
roi_R2_all    = [];
roi_laps_all  = [];
roi_sess_all  = [];
roi_tuned_all = [];

for iSess = 1:num_sessions

    % --- laps for this session ---
    try
        nLaps = dataToUse(iSess).ConditionData.Baseline.NumLaps;
    catch
        warning('Session %d: no ConditionData.Baseline.NumLaps found, skipping.', iSess);
        continue;
    end
    if isempty(nLaps), continue; end
    session_laps(iSess) = nLaps;

    % --- R2 / tuning data for this session ---
    current_cv = dataToUse(iSess).cvExpVar;
    true_means = current_cv.meanExpVar;
    if isempty(true_means), continue; end

    nROIs = length(true_means);
    null_matrix = current_cv.cvExpVarNull;
    null_dims = size(null_matrix);
    roi_dim = find(null_dims == nROIs, 1, 'first');
    if isempty(roi_dim)
        error('Session %d: could not match null matrix dim to nROIs.', iSess);
    end
    dims_to_collapse = setdiff(1:length(null_dims), roi_dim);
    if isempty(dims_to_collapse)
        max_null_per_ROI = null_matrix;
    else
        max_null_per_ROI = max(null_matrix, [], dims_to_collapse);
    end
    max_null_per_ROI = max_null_per_ROI(:);
    true_means = true_means(:);

    pVals = current_cv.pValues(:);

    % Combined "tuned" criterion: EV threshold + null-max check + shuffle p-value
    tuned_mask = (true_means > combinedEVThreshold) & ...
                 (max_null_per_ROI > nullMaxThreshold) & ...
                 (pVals < pValThreshold);

    session_nROIs(iSess)     = nROIs;
    session_meanR2(iSess)    = mean(true_means);
    session_medianR2(iSess)  = median(true_means);
    session_pctTuned(iSess)  = 100 * sum(tuned_mask) / nROIs;

    % accumulate ROI-level (pseudoreplicated) data
    roi_R2_all    = [roi_R2_all; true_means];
    roi_laps_all  = [roi_laps_all; repmat(nLaps, nROIs, 1)];
    roi_sess_all  = [roi_sess_all; repmat(iSess, nROIs, 1)];
    roi_tuned_all = [roi_tuned_all; tuned_mask];
end

valid = ~isnan(session_laps) & ~isnan(session_meanR2);

%% --- Report text-ready stats ---

fprintf('\n=== LAPS SUMMARY (n = %d sessions) ===\n', sum(valid));
fprintf('Mean +/- SD laps: %.1f +/- %.1f\n', mean(session_laps(valid)), std(session_laps(valid)));
fprintf('Range: %d - %d laps\n', min(session_laps(valid)), max(session_laps(valid)));

fprintf('\n=== %% TUNED SUMMARY ===\n');
overall_pct_tuned = 100 * sum(roi_tuned_all) / length(roi_tuned_all);
fprintf('Overall %% tuned (pooled across all ROIs): %.1f%%\n', overall_pct_tuned);
fprintf('Per-session %% tuned range: %.1f - %.1f%%\n', ...
    min(session_pctTuned(valid)), max(session_pctTuned(valid)));

%% --- (a) Session-level: laps vs mean R2, and laps vs %tuned ---

[r_meanR2, p_meanR2]   = corr(session_laps(valid), session_meanR2(valid), 'Type', 'Pearson');
[rho_meanR2, p_rho_meanR2] = corr(session_laps(valid), session_meanR2(valid), 'Type', 'Spearman');

[r_pctTuned, p_pctTuned] = corr(session_laps(valid), session_pctTuned(valid), 'Type', 'Pearson');
[rho_pctTuned, p_rho_pctTuned] = corr(session_laps(valid), session_pctTuned(valid), 'Type', 'Spearman');

fprintf('\n=== SESSION-LEVEL: Laps vs mean R2 (n=%d sessions) ===\n', sum(valid));
fprintf('Pearson r = %.3f, p = %.4f\n', r_meanR2, p_meanR2);
fprintf('Spearman rho = %.3f, p = %.4f\n', rho_meanR2, p_rho_meanR2);

fprintf('\n=== SESSION-LEVEL: Laps vs %% tuned (n=%d sessions) ===\n', sum(valid));
fprintf('Pearson r = %.3f, p = %.4f\n', r_pctTuned, p_pctTuned);
fprintf('Spearman rho = %.3f, p = %.4f\n', rho_pctTuned, p_rho_pctTuned);

% Plot: laps vs mean R2
figLapsR2 = figure('Color','w','Position',[200,200,500,400]);
scatter(session_laps(valid), session_meanR2(valid), 60, 'filled', ...
    'MarkerFaceColor', '#008080', 'MarkerFaceAlpha', 0.8);
lsline;
xlabel('Number of laps (session)');
ylabel('Mean R^2 across ROIs');
title(sprintf('r = %.2f, p = %.3f', r_meanR2, p_meanR2));
defaultAxesProperties(gca, true);
saveFigureFormats(figLapsR2, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\LapsVsR2\LapsVsMeanR2_sessionLevel');

% Plot: laps vs %tuned
figLapsTuned = figure('Color','w','Position',[200,200,500,400]);
scatter(session_laps(valid), session_pctTuned(valid), 60, 'filled', ...
    'MarkerFaceColor', '#008080', 'MarkerFaceAlpha', 0.8);
lsline;
xlabel('Number of laps (session)');
ylabel('% tuned ROIs');
title(sprintf('r = %.2f, p = %.3f', r_pctTuned, p_pctTuned));
defaultAxesProperties(gca, true);
saveFigureFormats(figLapsTuned, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\LapsVsR2\LapsVsPctTuned_sessionLevel');

%% --- (b) Supplementary ROI-level scatter (PSEUDOREPLICATED - descriptive only) ---
% NOTE: each session contributes many ROIs sharing one lap count, so a naive
% correlation here treats ROIs as independent when they are not - do not
% report a p-value from this as if it were a valid inferential test.
% Shown only to visualize whether the session-level trend holds up at the
% single-ROI level, colored by session.

figROI = figure('Color','w','Position',[250,200,550,450]);
hold on;
cmap = parula(num_sessions);
for iSess = 1:num_sessions
    idx = roi_sess_all == iSess;
    if any(idx)
        scatter(roi_laps_all(idx), roi_R2_all(idx), 12, cmap(iSess,:), 'filled', ...
            'MarkerFaceAlpha', 0.35);
    end
end
xlabel('Number of laps (session)');
ylabel('R^2 (single ROI)');
title('ROI-level (pseudoreplicated - descriptive only)');
defaultAxesProperties(gca, true);
hold off;
saveFigureFormats(figROI, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_2_Section1\LapsVsR2\LapsVsR2_ROIlevel_descriptive');
fprintf('\n=== OVERALL CELL COUNT ===\n');
total_cells_pooled = length(roi_R2_all);
fprintf('Total ROIs/neurons pooled across all sessions: %d\n', total_cells_pooled);
fprintf('Per-session ROI count range: %d - %d\n', ...
    min(session_nROIs(valid)), max(session_nROIs(valid)));
fprintf('Mean +/- SD ROIs per session: %.1f +/- %.1f\n', ...
    mean(session_nROIs(valid)), std(session_nROIs(valid)));