% Figure3.2: Responsive boutons from day 5 onwards

pairs=struct;

pairs.M25132 = ['20260313', '20260228']; 
pairs.M26003 = ['20260324', '20260325'];

% compute across all sessions because this chapter now includes the
% learning days as well. 
% pairs.M25132 = ['20260219','20260220','20260221','20260223', '20260226','20260228', '20260313']; %
% pairs.M25133 = ['20260220','20260221','20260223''20260219','20260224'];
% pairs.M26003 = ['20260316','20260317', '20260320','20260321','20260322', '20260324', '20260325']; %


ExpRSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
RSPData = getTuningDataByCondition(ExpRSPSessions);
RSPDataSpks = getTuningDataByCondition(ExpRSPSessions, 'signalToUse', 'spks'); 
% bin data by condition and load inclusion critera

%% ev
trueEVThreshold = 0.1;

% --- Pool ΔF/F ---
pooled_true_EV = [];
pooled_null_EV = [];
pooled_null_max_per_ROI = []; % Track the highest null value an ROI achieved across shuffles

% --- Pool spks ---
pooled_true_EV_spks = [];
pooled_null_EV_spks = [];
pooled_null_max_per_ROI_spks = [];

num_sessions = length(RSPData);

% Loop through each session and pool the ROIs (ΔF/F)
for iSess = 1:num_sessions
    current_cv = RSPData(iSess).cvExpVar;
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

% Loop through each session and pool the ROIs (spks)
num_sessions_spks = length(RSPDataSpks);
for iSess = 1:num_sessions_spks
    current_cv = RSPDataSpks(iSess).cvExpVar;
    true_means = current_cv.meanExpVar;
    null_matrix = current_cv.cvExpVarNull;
    
    if isempty(true_means) || isempty(null_matrix)
        continue;
    end
    
    nROIs = length(true_means);
    null_dims = size(null_matrix);
    roi_dim = find(null_dims == nROIs, 1, 'first');
    
    if isempty(roi_dim)
        error('Session %d (spks): Could not find a dimension matching %d ROIs.', iSess, nROIs);
    end
    
    dims_to_collapse = setdiff(1:length(null_dims), roi_dim);
    
    if isempty(dims_to_collapse)
        max_null_per_ROI = null_matrix;
    else
        max_null_per_ROI = max(null_matrix, [], dims_to_collapse);
    end
    
    pooled_true_EV_spks = [pooled_true_EV_spks; true_means(:)];
    pooled_null_EV_spks = [pooled_null_EV_spks; null_matrix(:)];
    pooled_null_max_per_ROI_spks = [pooled_null_max_per_ROI_spks; max_null_per_ROI(:)];
end

% Keep the threshold calculation logic exactly the same for annotations
passing_idx = (pooled_true_EV > 0.1) & (pooled_null_max_per_ROI > 0.01);
passing_idx_spks = (pooled_true_EV_spks > 0.1) & (pooled_null_max_per_ROI_spks > 0.01);

% --- Plotting ---
fig1 = figure('Color', 'w', 'Position', [200, 200, 600, 450]);
hold on;

bin_width = 0.02; 
bin_edges = -0.8:bin_width:1.0; 

% Shuffled Baseline Profile Line (ΔF/F)
h_null = histogram(pooled_null_EV, 'BinEdges', bin_edges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
    'EdgeColor', [0.55, 0.55, 0.55], 'LineWidth', 1.5);

xlabel('Explained variance (R^2)');
ylabel('Probability');

% True EV Profile Line (ΔF/F)
h_true = histogram(pooled_true_EV, 'BinEdges', bin_edges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
    'EdgeColor', '#008080', 'LineWidth', 1.5);

% Shuffled Baseline Profile Line (spks)
h_null_spks = histogram(pooled_null_EV_spks, 'BinEdges', bin_edges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
    'EdgeColor', [0.75, 0.75, 0.75], 'LineWidth', 1.5, 'LineStyle', ':');

% True EV Profile Line (spks)
h_true_spks = histogram(pooled_true_EV_spks, 'BinEdges', bin_edges, ...
    'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
    'EdgeColor', '#D95F02', 'LineWidth', 1.5);

% --- Annotations & Markers ---
pooled_99th = prctile(pooled_null_EV, 99);
xline(pooled_99th, 'r--', 'LineWidth', 1.5);
xline(trueEVThreshold, 'k--', 'LineWidth', 1.5);

y_limits = ylim;
text_y_position = y_limits(2) * 0.9; 

text(pooled_99th + 0.03, text_y_position, 'Null 99%', ...
     'Color', 'r', 'FontSize', 9, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% Calculate passing percentages (dF/F)
num_passing = sum(passing_idx);
total_boutons = length(pooled_true_EV);
pct_passing = (num_passing / total_boutons) * 100;

text(pooled_99th + 0.05, text_y_position * 0.8, sprintf('%s Passing: %.1f%%', 'dF/F', pct_passing), ...
     'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% Calculate passing percentages (spks)
num_passing_spks = sum(passing_idx_spks);
total_boutons_spks = length(pooled_true_EV_spks);
pct_passing_spks = (num_passing_spks / total_boutons_spks) * 100;

text(pooled_99th + 0.05, text_y_position * 0.65, sprintf('%s Passing: %.1f%%', 'spks', pct_passing_spks), ...
     'Color', '#D95F02', 'FontSize', 10, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% Set up Legend for the four outline distributions
lgd = legend([h_null, h_true, h_null_spks, h_true_spks], ...
    {'Shuffled (dF/F)', 'True EV (dF/F)', 'Shuffled (spks)', 'True EV (spks)'}, ...
    'Location', 'NorthEast');
defaultAxesProperties(gca, true);
legend boxoff;

xticks(-0.8:0.2:1.0);  
xticks([-0.8, -0.5, 0, 0.5, 1.0]);  
xlim([-0.8, 1.0]);
hold off;

% Save Figure
saveFigureFormats(fig1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.2\EVHistogram_spks_dff_manipulationonly\evDistributions_spks_dff');


%% Scatter of true EV, dF/F vs spikes
fig2 = figure('Color', 'w', 'Position', [200, 200, 400, 450]);
hold on;
jitter_width = 0.15;

% Only boutons that pass the R^2 > 0.1 threshold and the shuffled null (p < 0.01)
tuned_dff = pooled_true_EV(passing_idx);
tuned_spks = pooled_true_EV_spks(passing_idx_spks);

n_dff = length(tuned_dff);
n_spks = length(tuned_spks);
x_dff = 1 + (rand(n_dff, 1) - 0.5) * jitter_width * 2;
x_spks = 2 + (rand(n_spks, 1) - 0.5) * jitter_width * 2;
scatter(x_dff, tuned_dff, 15, 'MarkerFaceColor', '#74B448', ...
'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
scatter(x_spks, tuned_spks, 15, 'MarkerFaceColor', 'k', ...
'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
% Overlay median markers for each group
plot([1 - jitter_width, 1 + jitter_width], [median(tuned_dff), median(tuned_dff)], ...
'Color', '#8896A5', 'LineWidth', 2.5);
plot([2 - jitter_width, 2 + jitter_width], [median(tuned_spks), median(tuned_spks)], ...
'Color', '#8896A5', 'LineWidth', 2.5);
xlim([0.5, 2.5]);
xticks([1, 2]);
xticklabels({'dF/F', 'spks'});
ylabel('Explained variance (R^2)');
defaultAxesProperties(gca, true);
hold off;

[p, h, stats] = signrank(pooled_true_EV, pooled_true_EV_spks);
fprintf('Wilcoxon signed-rank test: p = %.3g, n = %d\n', p, length(pooled_true_EV));
fprintf('Median dF/F EV = %.3f, Median spks EV = %.3f\n', ...
    median(pooled_true_EV), median(pooled_true_EV_spks));
% That result says explained variance is significantly higher for ΔF/F than for deconvolved spikes in this same set of boutons (p = 7.42×10⁻¹⁷, Wilcoxon signed-rank), with the median dropping from 0.103 to 0.051, roughly half. That's consistent with the direction you already reported (51.1% vs 36% tuned) and supports the point that deconvolution reduces reliability for bouton signals.

% Save Figure
saveFigureFormats(fig2, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.2\EVHistogram_spks_dff_manipulationonly\evScatter_spks_dff');