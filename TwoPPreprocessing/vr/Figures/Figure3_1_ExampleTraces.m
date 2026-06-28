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