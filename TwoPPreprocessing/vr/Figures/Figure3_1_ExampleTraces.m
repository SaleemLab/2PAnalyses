% Example traces Fig 2.1 [Bouton traces and data streams]
% Load example sessions 
sessionFileInfo = get2PsessionFilePaths('M25132', '20260226');
VRIdx = find(contains({sessionFileInfo.stimFiles.name}, 'Corridor',  'IgnoreCase', true));
if length(VRIdx) > 1
    VRIdx = VRIdx(1); 
end 
response = load(sessionFileInfo.stimFiles(VRIdx).Response, "wheelSpeed", "lapPositionActivity", "lapPositionRunningSpeed", "pupilArea", 'mouseVirtualPosition'); 
proc2PData = load(sessionFileInfo.stimFiles(VRIdx).processedMergedBonsaiSuite2pData, 'zScoredProcessedSignals', 'TwoPFrameTime', 'ops');
load(sessionFileInfo.stimFiles(VRIdx).BonsaiData);
load(sessionFileInfo.stimFiles(VRIdx).processedPeripheralData);
 
%% Define your custom time window here (in seconds)
plot_start_time = 1000; 
plot_end_time   = 1436; 
%% Define data metrics
roi_idx = [6 7 9 1 30];
selected = proc2PData.zScoredProcessedSignals.dFFNeuropilCorrected(roi_idx, :);
num_rois = length(roi_idx);
time_vector = proc2PData.TwoPFrameTime;
%% Scale Bar Constants
t_scale = plot_end_time; 
t_pad   = (plot_end_time - plot_start_time) * 0.005; 
tick_w  = (plot_end_time - plot_start_time) * 0.007; 
%% Colour maps 
color_neural = [0.25, 0.25, 0.25];   % Dark Gray
color_vr     = "#ff4500";             % Orange
color_wheel  = '#87cefa';             % Sky Blue
color_pupil  = '#9370db';             % Purple
%% Set up Figure 1 Window (Data Streams Tracking)
fig = figure('Color', 'w', 'Position', [100, 100, 800, 650]);
line_width = 1;
%% Neural  
ax_neural = axes('Position', [0.05, 0.45, 0.85, 0.50]); 
hold on;
spacing = max(selected(:)) * 1.3; 
for i = 1:num_rois
    plot(time_vector, selected(i, :) + (i-1)*spacing, 'Color', color_neural, 'LineWidth', line_width);
end
scale_value = 2; 
scale_bar_y_start = (num_rois-1) * spacing; 
scale_bar_y_end   = scale_bar_y_start + scale_value; 
plot([t_scale, t_scale], [scale_bar_y_start, scale_bar_y_end], 'k', 'LineWidth', 1.5);
plot([t_scale, t_scale + tick_w], [scale_bar_y_start, scale_bar_y_start], 'k', 'LineWidth', 1.2);
plot([t_scale, t_scale + tick_w], [scale_bar_y_end, scale_bar_y_end], 'k', 'LineWidth', 1.2);
text(t_scale + t_pad + tick_w, (scale_bar_y_start + scale_bar_y_end)/2, '2 \DeltaF/F', 'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', 'k');
axis off;
%% VR Position
ax_vr = axes('Position', [0.05, 0.31, 0.85, 0.11]);
hold on;
plot(time_vector, response.mouseVirtualPosition, 'Color', color_vr, 'LineWidth', line_width);
plot([t_scale, t_scale], [1, 200], 'k', 'LineWidth', 1.5);
plot([t_scale, t_scale + tick_w], [1, 1], 'k', 'LineWidth', 1.2);       
plot([t_scale, t_scale + tick_w], [200, 200], 'k', 'LineWidth', 1.2); 
text(t_scale + t_pad + tick_w, 200, '200 cm', 'VerticalAlignment', 'middle', 'FontSize', 9);
text(t_scale + t_pad + tick_w, 1, '1 cm', 'VerticalAlignment', 'middle', 'FontSize', 9);
axis off;
%% Run Speed
ax_wheel = axes('Position', [0.05, 0.17, 0.85, 0.11]);
hold on;
plot(time_vector, response.wheelSpeed, 'Color', color_wheel, 'LineWidth', line_width);
plot([t_scale, t_scale], [0, 30], 'k', 'LineWidth', 1.5);
plot([t_scale, t_scale + tick_w], [0, 0], 'k', 'LineWidth', 1.2);   
plot([t_scale, t_scale + tick_w], [30, 30], 'k', 'LineWidth', 1.2); 
text(t_scale + t_pad + tick_w, 30, '30 cm/s', 'VerticalAlignment', 'middle', 'FontSize', 9);
text(t_scale + t_pad + tick_w, 0, '0 cm/s', 'VerticalAlignment', 'middle', 'FontSize', 9);
axis off;
%% Pupil Area (Z-Scored)
ax_pupil = axes('Position', [0.05, 0.03, 0.85, 0.11]); 
hold on;

% Z-score pupil data across time
raw_pupil_data = peripheralData.Pupil.Value.Area;
mu_pupil = mean(raw_pupil_data, 'omitnan');
sig_pupil = std(raw_pupil_data, 'omitnan');
if sig_pupil == 0, sig_pupil = 1; end
zscore_pupil = (raw_pupil_data - mu_pupil) / sig_pupil;

% Pull metrics inside window bounds for plotting scale offsets
visible_indices = (time_vector >= plot_start_time) & (time_vector <= plot_end_time);
visible_zscore = zscore_pupil(visible_indices);
p_min = min(visible_zscore);
p_max = max(visible_zscore);

plot(time_vector, zscore_pupil, 'Color', color_pupil, 'LineWidth', line_width);

% Scale bar showing 2 standard deviations (\sigma) relative to z-score metric
plot([t_scale, t_scale], [0, 2], 'k', 'LineWidth', 1.5);
plot([t_scale, t_scale + tick_w], [0, 0], 'k', 'LineWidth', 1.2);
plot([t_scale, t_scale + tick_w], [2, 2], 'k', 'LineWidth', 1.2);
text(t_scale + t_pad + tick_w, 1, '2 \sigma (\Delta pupil)', 'VerticalAlignment', 'middle', 'FontSize', 9);

% Horizontal 30-second scale bar positioned precisely below the window minimum
time_bar_length = 30; 
time_bar_start = plot_start_time + (plot_end_time - plot_start_time) * 0.45; 
time_bar_y = p_min - 1.5; 

plot([time_bar_start, time_bar_start + time_bar_length], [time_bar_y, time_bar_y], 'k', 'LineWidth', 2);
text(time_bar_start + (time_bar_length/2), time_bar_y, '30 s', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);

% Bound the window coordinates so elements remain within focus
ylim([time_bar_y * 1.3, p_max * 1.3]);
axis off;
%% Sync across all axes and lock view to limits 
all_axes = [ax_neural, ax_vr, ax_wheel, ax_pupil];
linkaxes(all_axes, 'x');
set(all_axes, 'XLim', [plot_start_time, plot_end_time + (t_pad * 14)]); 
saveFigureFormats(fig, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section1_Fig3.1\Fig3.1_ExpSetup\traces\ExampleDataTraces');

%% Set up Figure 2 Window (FOV Reference Image)
sessionFileInfo = get2PsessionFilePaths('M26003', '20260322');
VRIdx = find(contains({sessionFileInfo.stimFiles.name}, 'Corridor',  'IgnoreCase', true));
if length(VRIdx) > 1
    VRIdx = VRIdx(1); 
end 

% Scale bar representin
img_pixels = 256;       % Width/height in pixels
img_microns = 81;       % Width/height in microns/ At zoom 9; 95 microns at zoom 8 
desired_bar_um = 20;
% Calculate how many pixels represent 1 micron
pixel_per_micron = img_pixels / img_microns;
bar_width_px = desired_bar_um * pixel_per_micron; 
bar_height_px = 6; % Thickness of the scale bar in pixels
edge_padding = 15; 

% Calculate the starting coordinates (Top-Left corner of the scale bar rectangle)
x_start = img_pixels - edge_padding - bar_width_px;
y_start = img_pixels - edge_padding - bar_height_px;

% 
proc2PData = load(sessionFileInfo.stimFiles(VRIdx).processedMergedBonsaiSuite2pData, 'ops');
fig2 = figure('Color', 'w', 'Position', [100, 100, 400, 400]);
imagesc(proc2PData.ops{1}.refImg);
colormap gray;
axis image off;
hold on 

rectangle('Position', [x_start, y_start, bar_width_px, bar_height_px], ...
          'FaceColor', 'w', 'EdgeColor', 'none')
hold off 


fig3 = figure('Color', 'w', 'Position', [100, 100, 400, 400]);
imagesc(proc2PData.ops{1}.meanImg_chan2);
colormap gray;
axis image off;
hold on; 

rectangle('Position', [x_start, y_start, bar_width_px, bar_height_px], ...
          'FaceColor', 'w', 'EdgeColor', 'none')
hold off; 
saveFigureFormats(fig2, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section1_Fig3.1\Fig3.1_ExpSetup\FOVs\ExampleFOV_GreenChan');
saveFigureFormats(fig3, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section1_Fig3.1\Fig3.1_ExpSetup\FOVs\ExampleFOV_RedChan');




%% multi plane 
suite2p_dir = 'Z:\ibn-vision\DATA\SUBJECTS\M25133\Processed\20260219\suite2p';
num_planes = 8;

plane_images = cell(1, num_planes);
sharpness_scores = nan(1, num_planes);
valid_planes = false(1, num_planes);

%% Load images and calculate Laplacian variance sharpness
fprintf('Loading Fall.mat files across planes...\n');

for p = 0:(num_planes - 1)
    % Construct folder and file paths (note: backup0, backup1, ...)
    fall_path = fullfile(suite2p_dir, sprintf('backup%d', p), 'Fall.mat');
    
    if ~exist(fall_path, 'file')
        warning('Plane %d: Fall.mat not found at %s - skipping', p, fall_path);
        continue;
    end
    
    % Load Fall.mat (ops is stored inside)
    data = load(fall_path, 'ops');
    if ~isfield(data, 'ops')
        warning('Plane %d: "ops" struct not found in %s', p, fall_path);
        continue;
    end
    
    ops = data.ops;
    
    % Extract meanImg (ensure double precision)
    mean_img = double(ops.meanImg);
    
    % Calculate Laplacian variance sharpness score
    % (fspecial/imfilter equivalent to scipy.ndimage.laplace)
    lap_filter = [0 1 0; 1 -4 1; 0 1 0];
    lap_img = imfilter(mean_img, lap_filter, 'replicate');
    sharpness = var(lap_img(:), 'omitnan');
    
    % MATLAB cell indices are 1-based (Plane 0 -> Index 1)
    plane_idx = p + 1; 
    plane_images{plane_idx} = mean_img;
    sharpness_scores(plane_idx) = sharpness;
    valid_planes(plane_idx) = true;
    
    fprintf('Plane %d: sharpness = %.2f\n', p, sharpness);
end

%% Print Sharpness Summary
fprintf('\n--- Sharpness Summary (ascending = potential fly-back plane) ---\n');
[sorted_scores, sort_idx] = sort(sharpness_scores(valid_planes), 'ascend');
valid_plane_nums = find(valid_planes) - 1; % Convert back to 0-indexed plane numbers

for k = 1:length(sorted_scores)
    fprintf('Plane %d: %.2f\n', valid_plane_nums(sort_idx(k)), sorted_scores(k));
end

%% Plot Grid Diagnostic Figure
num_found = sum(valid_planes);
n_cols = 4;
n_rows = ceil(num_found / n_cols);

fig_grid = figure('Color', 'w', 'Position', [100, 100, 1200, 300 * n_rows]);
t = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

found_indices = find(valid_planes);

for i = 1:num_found
    idx = found_indices(i);
    plane_num = idx - 1;
    img = plane_images{idx};
    
    nexttile(t);
    
    % 1st and 99th percentile contrast stretch (similar to vmin/vmax in pyplot)
    p_limits = prctile(img(:), [1 99]);
    
    imagesc(img, p_limits);
    colormap gray;
    axis image off;
    title(sprintf('Plane %d\nsharpness = %.1f', plane_num, sharpness_scores(idx)), ...
        'FontSize', 10);
end

title(t, 'Mean image per plane — look for visibly blurred/smeared planes (fly-back)', ...
    'FontSize', 12, 'FontWeight', 'bold');

% Save grid plot
saveFigureFormats(fig_grid, fullfile(suite2p_dir, 'all_plane_mean_images'));

%% Pass collected images into your stack function
% Filter to valid images only
valid_images = plane_images(valid_planes);

% Call isometric stack plot (adjust z_spacing and rot_deg as needed)
fig_stack = plot_isometric_stack(valid_images, 0.25, 180);



%% supp figure:
%% Set up Figure 2 Window (FOV Reference Image)
sessionFileInfo = get2PsessionFilePaths('MI268', '20260721B');
RF = find(contains({sessionFileInfo.stimFiles.name}, 'RFMapping', 'IgnoreCase', true));
if length(RF) > 1
    RF = RF(1); 
end 

% Scale bar parameters (Zoom 9)
img_pixels = 256;       
img_microns = 81;       
desired_bar_um = 20;

pixel_per_micron = img_pixels / img_microns;
bar_width_px = desired_bar_um * pixel_per_micron; 
bar_height_px = 6; 
edge_padding = 15; 

x_start = img_pixels - edge_padding - bar_width_px;
y_start = img_pixels - edge_padding - bar_height_px;

% Load ops and stat
proc2PData = load(sessionFileInfo.stimFiles(RF).processedMergedBonsaiSuite2pData, 'ops', 'stat');

% Extract ops
if iscell(proc2PData.ops) && iscell(proc2PData.ops{1})
    ops_p1 = proc2PData.ops{1}{1};
elseif iscell(proc2PData.ops)
    ops_p1 = proc2PData.ops{1};
else
    ops_p1 = proc2PData.ops;
end

% FIX: Extract all ROIs without discarding elements by taking {1} prematurely
stat_all = proc2PData.stat;
if iscell(stat_all) && length(stat_all) == 1 && iscell(stat_all{1})
    stat_all = stat_all{1}; % Unwrap plane layer if multi-plane cell
end

num_rois = length(stat_all);

fig2 = figure('Color', 'w', 'Position', [100, 100, 400, 400]);
imagesc(ops_p1.meanImgE);
colormap gray;
axis image off;
hold on;

% Concatenate all ROI outline points into single vectors for fast plotting
all_x = [];
all_y = [];

for roi_idx = 1:num_rois
    if iscell(stat_all)
        st = stat_all{roi_idx};
        if iscell(st), st = st{1}; end
    else
        st = stat_all(roi_idx);
    end
    
    x = double(st.xpix);
    y = double(st.ypix);
    
    % Adjust 0-based Python indices to 1-based MATLAB indices
    if min(x) == 0 || min(y) == 0
        x = x + 1;
        y = y + 1;
    end
    
    try
        k = convhull(x, y);
        bx = x(k);
        by = y(k);
    catch
        bx = x;
        by = y;
    end
    
    % Append polygon coordinates separated by NaN
    all_x = [all_x, bx(:)', NaN];
    all_y = [all_y, by(:)', NaN];
end

% Plot ALL 435 masks in one shot
plot(all_x, all_y, 'Color', [1, 0.2, 0.2], 'LineWidth', 0.8);

% Scale bar
rectangle('Position', [x_start, y_start, bar_width_px, bar_height_px], ...
          'FaceColor', 'w', 'EdgeColor', 'none');
hold off;

%% Export Figure
save_path_fig2 = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Supp_Section1_Fig4_1_VISp\fovs\ExampleFOV_Somas';
saveFigureFormats(fig2, save_path_fig2);