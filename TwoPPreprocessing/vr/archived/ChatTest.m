% =========================================================================
% SCRIPT: Verify Neuropil Correction Pipelines (Old vs. New)
% =========================================================================
% This script generates dummy data and runs it through both the slow,
% loop-based pipeline and the fast, vectorized pipeline to verify that
% their outputs are functionally identical.

clear;
close all;
clc;

%% 1. GENERATE REALISTIC DUMMY DATA
fprintf('Generating realistic dummy data...\n\n');

% Parameters
fs = 30;             % Frame rate (Hz)
n_rois = 10;         % Number of ROIs
n_seconds = 800;     % Duration of signal in seconds
nt = n_seconds * fs; % Total number of time points

% Create a shared slow, drifting baseline
baseline_drift = 30 * sin(linspace(0, 5*pi, nt))' + 150;
baseline_drift = baseline_drift + randn(nt, 1) * 3; % Add noise

% Create fluorescence signal (F) with sharp calcium events
F = repmat(baseline_drift, 1, n_rois);
for i = 1:n_rois
    num_spikes = randi([15, 25]);
    for j = 1:num_spikes
        spike_frame = randi([1, nt]);
        spike_amp = randi([30, 100]);
        spike_decay = exp(-(0:10*fs)' / (fs*1.5)); % 1.5-second decay
        end_frame = min(nt, spike_frame + length(spike_decay) - 1);
        F(spike_frame:end_frame, i) = F(spike_frame:end_frame, i) + spike_amp * spike_decay(1:end_frame-spike_frame+1);
    end
end

% Create a correlated neuropil signal (N)
% It should be a scaled, slightly noisier version of the true signal
N = F * 0.4 + randn(nt, n_rois) * 5 + repmat(baseline_drift*0.2, 1, n_rois);


%% 2. RUN BOTH PIPELINES AND TIME THEM

% --- Run OLD, slow pipeline ---
fprintf('Running OLD pipeline (slow loops)...\n');
tic;
[signal_old, regPars_old] = correct_neuropil_OLD(F, N, fs);
time_old = toc;
fprintf('-> Old pipeline finished in %.2f seconds.\n\n', time_old);

% --- Run NEW, fast pipeline ---
fprintf('Running NEW pipeline (vectorized)...\n');
tic;
[signal_new, regPars_new] = correct_neuropil_NEW(F, N, fs);
time_new = toc;
fprintf('-> New pipeline finished in %.2f seconds.\n\n', time_new);


%% 3. QUANTITATIVELY COMPARE THE OUTPUTS
fprintf('Comparing outputs...\n');

% Compare the final corrected signals
difference_signal = signal_new - signal_old;
max_abs_diff_signal = max(abs(difference_signal(:)));
fprintf('Maximum absolute difference in corrected signal: %e\n', max_abs_diff_signal);

% Compare the regression parameters (the slope 'b' is the most important)
difference_params = regPars_new - regPars_old;
max_abs_diff_params = max(abs(difference_params(:)));
fprintf('Maximum absolute difference in regression parameters: %e\n', max_abs_diff_params);
fprintf('--------------------------------------------------\n');
fprintf('Speedup factor: %.1fx faster\n', time_old / time_new);
fprintf('--------------------------------------------------\n');

%% 4. VISUALIZE THE COMPARISON
fprintf('Generating comparison plot...\n');
roi_to_plot = 1; % Which ROI to visualize
time_axis = (1:nt) / fs;

figure('Position', [50, 50, 1400, 900], 'Color', 'w');
sgtitle('Verification: Old vs. New Neuropil Correction Pipeline', 'FontSize', 18, 'FontWeight', 'bold');

% --- Top Plot: Corrected Signals Overlay ---
ax1 = subplot(2, 2, [1 2]);
plot(time_axis, F(:, roi_to_plot), 'Color', [0.8 0.8 0.8], 'DisplayName', 'Original F');
hold on;
plot(time_axis, signal_old(:, roi_to_plot), 'b-', 'LineWidth', 2.5, 'DisplayName', 'Corrected (OLD)');
plot(time_axis, signal_new(:, roi_to_plot), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Corrected (NEW)');
hold off;
title('Overlay of Final Corrected Signals (Should Overlap Perfectly)', 'FontSize', 14);
xlabel('Time (s)');
ylabel('Fluorescence (a.u.)');
legend('show', 'Location', 'northwest');
grid on; axis tight;

% --- Bottom-Left Plot: Difference Signal ---
ax2 = subplot(2, 2, 3);
plot(time_axis, difference_signal(:, roi_to_plot), 'k-', 'LineWidth', 1);
title('Difference (New Signal - Old Signal)', 'FontSize', 14);
xlabel('Time (s)');
ylabel('Difference (a.u.)');
grid on; axis tight;
ylim([-1e-9, 1e-9]); % Zoom in to show the difference is negligible

% --- Bottom-Right Plot: Regression Parameter Comparison ---
ax3 = subplot(2, 2, 4);
scatter(regPars_old(2, :), regPars_new(2, :), 50, 'filled');
hold on;
plot([0 2], [0 2], 'r--', 'LineWidth', 2, 'DisplayName', 'y=x line');
hold off;
title('Comparison of Neuropil Correction Coefficient (b)', 'FontSize', 14);
xlabel('Coefficient from OLD pipeline');
ylabel('Coefficient from NEW pipeline');
grid on; axis square;
legend('show', 'Location', 'northwest');


% #########################################################################
% #                 FUNCTION DEFINITIONS (PASTED BELOW)                   #
% #########################################################################

% --- NEW, FAST PIPELINE ---
function F0 = get_F0_NEW(Fc, fs, prctl_F, window_size)
    if nargin < 3, prctl_F = 8; end
    if nargin < 4, window_size = 180; end
    window_frames = round(fs * window_size);
    if mod(window_frames, 2) == 0, window_frames = window_frames + 1; end
    order = round((prctl_F / 100) * window_frames);
    order = max(order, 1);
    domain = ones(window_frames, 1);
    F0 = ordfilt2(Fc, order, domain, 'symmetric');
end

function [signal, regPars] = correct_neuropil_NEW(F, N, fs)
    numN = 20; minNp = 10; maxNp = 90; prctl_F = 5; prctl_F0 = 5; Npil_window_F0 = 180;
    [nt, nROIs] = size(F);
    F_binValues = NaN(numN, nROIs);
    regPars = NaN(2, nROIs);
    signal = NaN(nt, nROIs);
    F0 = get_F0_NEW(F, fs, prctl_F0, Npil_window_F0);
    N0 = get_F0_NEW(N, fs, prctl_F0, Npil_window_F0);
    Fc = F - F0;
    Nc = N - N0;
    for iROI = 1:nROIs
        iN = Nc(:, iROI); iF = Fc(:, iROI);
        N_prct = prctile(iN, [minNp, maxNp]);
        binSize = (N_prct(2) - N_prct(1)) / numN;
        if binSize <= 0, continue; end
        N_ind = floor((iN - N_prct(1)) / binSize) + 1;
        valid_mask = N_ind >= 1 & N_ind <= numN & ~isnan(iF);
        F_binValues(:, iROI) = accumarray(N_ind(valid_mask), iF(valid_mask), [numN 1], @(x) prctile(x, prctl_F), NaN);
        N_binValues = N_prct(1) + ((1:numN)' - 0.5) * binSize;
        validIdx = ~isnan(F_binValues(:, iROI));
        if any(validIdx)
            [a, b, ~] = linear_analytical_solution(N_binValues(validIdx), F_binValues(validIdx, iROI), false);
            b = min(max(b, 0), 2);
            regPars(:, iROI) = [a; b];
            signal(:, iROI) = F(:, iROI) - b * N(:, iROI) - a;
        end
    end
end

% --- OLD, SLOW PIPELINE ---
function F0 = get_F0_OLD(Fc, fs, prctl_F, window_size)
    if nargin < 3, prctl_F = 8; end
    if nargin < 4, window_size = 180; end
    [nt, nROIs] = size(Fc);
    F0 = zeros(nt, nROIs);
    window_frames = round(fs * window_size);
    for roi = 1:nROIs
        for t = 1:nt
            win_start = max(1, t - floor(window_frames / 2));
            win_end = min(nt, t + floor(window_frames / 2));
            F0(t, roi) = prctile(Fc(win_start:win_end, roi), prctl_F);
        end
    end
end

function [signal, regPars] = correct_neuropil_OLD(F, N, fs)
    numN = 20; minNp = 10; maxNp = 90; prctl_F = 5; prctl_F0 = 5; Npil_window_F0 = 180;
    [nt, nROIs] = size(F);
    F_binValues = NaN(numN, nROIs);
    regPars = NaN(2, nROIs);
    signal = NaN(nt, nROIs);
    F0 = get_F0_OLD(F, fs, prctl_F0, Npil_window_F0);
    N0 = get_F0_OLD(N, fs, prctl_F0, Npil_window_F0);
    Fc = F - F0;
    Nc = N - N0;
    for iROI = 1:nROIs
        iN = Nc(:, iROI); iF = Fc(:, iROI);
        N_prct = prctile(iN, [minNp, maxNp]);
        binSize = (N_prct(2) - N_prct(1)) / numN;
        if binSize <= 0, continue; end
        N_binValues = N_prct(1) + (0:numN-1)' * binSize;
        N_ind = floor((iN - N_prct(1)) / binSize);
        for Ni = 1:numN
            tmp = NaN(size(iF));
            tmp(N_ind == Ni-1) = iF(N_ind == Ni-1);
            F_binValues(Ni, iROI) = prctile(tmp(~isnan(tmp)), prctl_F);
        end
        validIdx = ~isnan(F_binValues(:, iROI));
        if any(validIdx)
            [a, b, ~] = linear_analytical_solution(N_binValues(validIdx), F_binValues(validIdx, iROI), false);
            b = min(max(b, 0), 2);
            regPars(:, iROI) = [a; b];
            signal(:, iROI) = F(:, iROI) - b * N(:, iROI) - a;
        end
    end
end

% --- SHARED DEPENDENCY ---
function [a, b, mse] = linear_analytical_solution(x, y, noIntercept)
    if nargin < 3, noIntercept = false; end
    x = x(:); y = y(:); n = length(x);
    if noIntercept
        b = sum(x .* y) / sum(x .^ 2); a = 0;
    else
        sum_x = sum(x); sum_y = sum(y); sum_x2 = sum(x .^ 2); sum_xy = sum(x .* y);
        denom = n * sum_x2 - sum_x^2;
        if denom == 0, a = NaN; b = NaN; mse = NaN; return; end
        a = (sum_y * sum_x2 - sum_x * sum_xy) / denom;
        b = (n * sum_xy - sum_x * sum_y) / denom;
    end
    y_fit = a + b * x;
    mse = mean((y - y_fit) .^ 2);
end