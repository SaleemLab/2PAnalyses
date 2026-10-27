%Example traces Fig 2.1 [Bouton traces and data streams]
%Load example sessions 
sessionFileInfo = get2PsessionFilePaths('M25132', '20260226');
VRIdx = find(contains({sessionFileInfo.stimFiles.name}, 'Corridor',  'IgnoreCase', true));
if length(VRIdx) > 1
    VRIdx = VRIdx(1);
end
response = load(sessionFileInfo.stimFiles(VRIdx).Response, "wheelSpeed", "lapPositionActivity", "lapPositionRunningSpeed", "pupilArea", 'mouseVirtualPosition');
proc2PData = load(sessionFileInfo.stimFiles(VRIdx).processedMergedBonsaiSuite2pData, 'zScoredProcessedSignals', 'TwoPFrameTime', 'ops', 'spks');
load(sessionFileInfo.stimFiles(VRIdx).BonsaiData);
load(sessionFileInfo.stimFiles(VRIdx).processedPeripheralData);

%% Define your custom time window here (in seconds)
plot_start_time = 1000; 
plot_end_time   = 1436; 

%% Define data metrics — just 3 example ROIs
roi_idx = [6 7 9 1 30];
selected = proc2PData.zScoredProcessedSignals.dFFNeuropilCorrected(roi_idx, :);
%num_rois = length(roi_idx);
%time_vector = proc2PData.TwoPFrameTime;
dff_data  = proc2PData.zScoredProcessedSignals.dFFNeuropilCorrected(roi_idx, :); % dF/F
spks_data = proc2PData.spks(roi_idx, :);                       % spikes
num_rois = length(roi_idx);
time_vector = proc2PData.TwoPFrameTime;

%% Scale Bar Constants
t_scale = plot_end_time;
t_pad   = (plot_end_time - plot_start_time) * 0.005;
tick_w  = (plot_end_time - plot_start_time) * 0.007;

%% Colour maps
color_neural = [0.25, 0.25, 0.25];   % Dark Gray (dF/F)
color_spikes = [0.85, 0.10, 0.10];   % Red (spike events)
color_vr     = "#ff4500";
color_wheel  = '#87cefa';
color_pupil  = '#9370db';

%% Set up Figure 1 Window (Data Streams Tracking)
fig = figure('Color', 'w', 'Position', [100, 100, 800, 650]);
line_width = 1;

%% Neural: dF/F traces with spike events overlaid
ax_neural = axes('Position', [0.05, 0.45, 0.85, 0.50]);
hold on;

spacing = max(dff_data(:)) * 1.3;

for i = 1:num_rois
    offset = (i-1) * spacing;

    % dF/F trace (already z-scored)
    plot(time_vector, dff_data(i, :) + offset, 'Color', color_neural, 'LineWidth', line_width);

    % Z-score spikes across time, same as dFF, then plot the same way
    mu_spk  = mean(spks_data(i, :), 'omitnan');
    sig_spk = std(spks_data(i, :), 'omitnan');
    if sig_spk == 0, sig_spk = 1; end
    zscore_spk = (spks_data(i, :) - mu_spk) / sig_spk;

    plot(time_vector, zscore_spk + offset, 'Color', color_spikes, 'LineWidth', line_width);
end

% dF/F scale bar (kept as before, still useful as a visual reference)
scale_value = 2;
scale_bar_y_start = (num_rois-1) * spacing;
scale_bar_y_end   = scale_bar_y_start + scale_value;
plot([t_scale, t_scale], [scale_bar_y_start, scale_bar_y_end], 'k', 'LineWidth', 1.5);
plot([t_scale, t_scale + tick_w], [scale_bar_y_start, scale_bar_y_start], 'k', 'LineWidth', 1.2);
plot([t_scale, t_scale + tick_w], [scale_bar_y_end, scale_bar_y_end], 'k', 'LineWidth', 1.2);
text(t_scale + t_pad + tick_w, (scale_bar_y_start + scale_bar_y_end)/2, '2 \sigma', ...
    'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', 'k');

% --- Real y-axis for ROI 1 only, showing actual z-score units ---
roi_ref = 1;                      % which ROI to use as the scale reference
offset_ref = (roi_ref - 1) * spacing;

box off;
ax_neural.XColor = 'none';        % hide x-axis (time handled by scale bar elsewhere)
ax_neural.YColor = 'k';           % show y-axis
ax_neural.TickDir = 'out';
ax_neural.YTick = offset_ref + [-2 0 2];        % z-score reference ticks for ROI 1
ax_neural.YTickLabel = {'-2', '0', '2'};
ylabel(ax_neural, 'z-score (\sigma)', 'FontSize', 9);
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
% ax_pupil = axes('Position', [0.05, 0.03, 0.85, 0.11]); 
% hold on;
% 
% % Z-score pupil data across time
% raw_pupil_data = peripheralData.Pupil.Value.Area;
% mu_pupil = mean(raw_pupil_data, 'omitnan');
% sig_pupil = std(raw_pupil_data, 'omitnan');
% if sig_pupil == 0, sig_pupil = 1; end
% zscore_pupil = (raw_pupil_data - mu_pupil) / sig_pupil;
% 
% % Pull metrics inside window bounds for plotting scale offsets
% visible_indices = (time_vector >= plot_start_time) & (time_vector <= plot_end_time);
% visible_zscore = zscore_pupil(visible_indices);
% p_min = min(visible_zscore);
% p_max = max(visible_zscore);
% 
% plot(time_vector, zscore_pupil, 'Color', color_pupil, 'LineWidth', line_width);
% 
% % Scale bar showing 2 standard deviations (\sigma) relative to z-score metric
% plot([t_scale, t_scale], [0, 2], 'k', 'LineWidth', 1.5);
% plot([t_scale, t_scale + tick_w], [0, 0], 'k', 'LineWidth', 1.2);
% plot([t_scale, t_scale + tick_w], [2, 2], 'k', 'LineWidth', 1.2);
% text(t_scale + t_pad + tick_w, 1, '2 \sigma (\Delta pupil)', 'VerticalAlignment', 'middle', 'FontSize', 9);
% 
% % Horizontal 30-second scale bar positioned precisely below the window minimum
% time_bar_length = 30; 
% time_bar_start = plot_start_time + (plot_end_time - plot_start_time) * 0.45; 
% time_bar_y = p_min - 1.5; 
% 
% plot([time_bar_start, time_bar_start + time_bar_length], [time_bar_y, time_bar_y], 'k', 'LineWidth', 2);
% text(time_bar_start + (time_bar_length/2), time_bar_y, '30 s', ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);
% 
% % Bound the window coordinates so elements remain within focus
% ylim([time_bar_y * 1.3, p_max * 1.3]);
% axis off;
ax_pupil = axes('Position', [0.05, 0.03, 0.85, 0.11]);
hold on;

% --- Clean pupil data before z-scoring ---
raw_pupil_data = peripheralData.Pupil.Value.Area;

clean_pupil = raw_pupil_data;
clean_pupil(clean_pupil <= 0) = NaN;

% Detect blink dips directly: find sharp downward deflections using
% peak-prominence on the inverted signal (targets dip shape, not just
% magnitude, so it won't miss dips that MAD/isoutlier under-flag)
inv_pupil = -clean_pupil;
inv_pupil(isnan(inv_pupil)) = min(inv_pupil, [], 'omitnan'); % temporarily fill NaNs so findpeaks doesn't choke
prominence_thresh = 0.5 * std(clean_pupil, 'omitnan'); % tune: lower = catches smaller dips too
[~, dip_locs, dip_widths] = findpeaks(inv_pupil, 'MinPeakProminence', prominence_thresh);

% Remove each detected dip plus a padding margin on both sides (blinks
% have ramp-up/ramp-down edges beyond just the peak sample)
pad_samples = 5; % samples of padding each side; increase if edges still show
for k = 1:length(dip_locs)
    half_width = ceil(dip_widths(k)/2) + pad_samples;
    lo = max(1, dip_locs(k) - half_width);
    hi = min(length(clean_pupil), dip_locs(k) + half_width);
    clean_pupil(lo:hi) = NaN;
end

clean_pupil = fillmissing(clean_pupil, 'pchip', 'EndValues', 'nearest');

% --- Z-score pupil data across time ---
mu_pupil = mean(clean_pupil, 'omitnan');
sig_pupil = std(clean_pupil, 'omitnan');
if sig_pupil == 0, sig_pupil = 1; end
zscore_pupil = (clean_pupil - mu_pupil) / sig_pupil;

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
time_bar_length = 10;
time_bar_start = plot_start_time + (plot_end_time - plot_start_time) * 0.45;
time_bar_y = p_min - 1.5;
plot([time_bar_start, time_bar_start + time_bar_length], [time_bar_y, time_bar_y], 'k', 'LineWidth', 2);
text(time_bar_start + (time_bar_length/2), time_bar_y, '10 s', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);

% Bound the window coordinates so elements remain within focus
ylim([time_bar_y * 1.3, p_max * 1.3]);
axis off;

%% Sync across all axes and lock view to limits 
all_axes = [ax_neural, ax_vr, ax_wheel, ax_pupil];
linkaxes(all_axes, 'x');
set(all_axes, 'XLim', [plot_start_time, plot_end_time + (t_pad * 14)]); 
saveFigureFormats(fig, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section1_Fig3.1\Fig3.1_ExpSetup\traces\ExampleDataTraces_withSpikes');

