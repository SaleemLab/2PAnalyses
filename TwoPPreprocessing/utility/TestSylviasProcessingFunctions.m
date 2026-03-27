load(sessionFileInfo.stimFiles(1).processedMergedBonsaiSuite2pData, 'F', 'Fneu', 'ops');

%% 1. Setup & ROI Selection
roi_idx = 200; 
fs = 60; 
absZero = 23; 

F_raw = F(roi_idx, :)';
N_raw = Fneu(roi_idx, :)';
t = (0:length(F_raw)-1)' / fs;

%% changed smoothning to 15 (250ms) 
% Offset subtraction and Gaussian smoothing (gausswin 15 ~ 250ms at 60Hz)
F_smooth = filtfilt(gausswin(15)/sum(gausswin(15)), 1, F_raw - absZero);
N_smooth = filtfilt(gausswin(15)/sum(gausswin(15)), 1, N_raw - absZero);

%% 
% with neuropil correction
[signal_corr, regPars, F_bins, N_bins] = correct_neuropil(F_smooth, N_smooth, fs);
F0_corr = get_F0(signal_corr, 8, 60);
dFF_corr = (signal_corr - F0_corr) ./ max(absZero, F0_corr);

% without neuropil correction 
signal_no_np = F_smooth; 
F0_no_np = get_F0(signal_no_np, 8, 60);
% Use max(absZero, F0) to prevent noise amplification in dim boutons.
% This  (23) prevents tiny fluctuations in dark pixels from looking like huge % changes.
dFF = (signal - F0) ./ max(absZero, F0);
dFF_no_np = (signal_no_np - F0_no_np) ./ max(absZero, F0_no_np); 

%% plots 
figure('Color', 'w', 'Name', ['Full Pipeline Comparison ROI ' num2str(roi_idx)], 'Position', [50 50 1100 950]);

% smoothed f and n 
subplot(4,1,1);
plot(t, F_smooth, 'r', 'DisplayName', 'Signal (Smoothed)'); hold on;
plot(t, N_smooth, 'Color', [0.4 0.4 0.4], 'DisplayName', 'Neuropil (Smoothed)');
plot(t, signal_corr, 'k', 'LineWidth', 1, 'DisplayName', 'Corrected (F - rN)');
ylabel('Intensity'); title(['Raw Intensity Traces (ROI ' num2str(roi_idx) ')']);
legend('Location', 'northeastoutside'); grid on;

% regression
subplot(4,1,2);
scatter(N_bins, F_bins, 15, 'filled', 'MarkerFaceAlpha', 0.5); hold on;
% Plot the regression line
x_vals = [min(N_bins) max(N_bins)];
y_vals = regPars(1) + regPars(2)*x_vals;
plot(x_vals, y_vals, 'r', 'LineWidth', 2);
xlabel('Neuropil Intensity'); ylabel('Signal Intensity');
title(['Neuropil Regression (Slope r = ', num2str(regPars(2), '%.2f'), ')']);
grid on;

% dF/F comparison with and without np correction 
subplot(4,1,3);
plot(t, dFF_no_np * 100, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'No Correction'); hold on;
plot(t, dFF_corr * 100, 'b', 'LineWidth', 1, 'DisplayName', 'Neuropil Corrected');
ylabel('\DeltaF/F (%)'); title('Final Activity Comparison');
legend('Location', 'northeastoutside'); 

% difference
subplot(4,1,4);
dFF_diff = (dFF_no_np - dFF_corr) * 100;
area(t, dFF_diff, 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'r', 'DisplayName', 'Signal Subtracted');
ylabel('Diff (\DeltaF/F %)'); xlabel('Time (s)');
title('Difference (Amount of signal removed by correction)');
legend('Location', 'northeastoutside');