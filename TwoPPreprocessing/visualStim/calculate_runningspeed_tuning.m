%%
function [speed_response, speed_response_std, speed_pos_FR, speed_tuning_type, speed_tuning_code] = calculate_runningspeed_tuning(temporal_activity, tvec, speed, position, track_ID_all, start_time_all, end_time_all)
% CALCULATE_SPEED_TUNING_ANALYSIS Bins firing rates, fits a speed tuning
% curve averaged across all conditions, and provides detailed classification.
%
% SYNTAX:
% [speed_response, speed_response_std, speed_pos_FR, speed_tuning_type, speed_tuning_code] = calculate_speed_tuning_analysis(...)
%
% OUTPUTS:
%   speed_response:     An N x S x (K+1) matrix of the mean firing rates binned by speed.
%   speed_response_std: An N x S x (K+1) matrix of the standard deviation of firing rates binned by speed.
%   speed_pos_FR:       An N x S x P x K matrix of firing rates binned by speed and position.
%   speed_tuning_type:  An N x 1 cell array with a string classification for each neuron
%                       ('lowpass', 'bandpass_5', 'highpass', 'untuned', etc.).
%   speed_tuning_code:  An N x 1 numeric matrix with a code for each classification.
%       - 0: untuned
%       - 2.5: lowpass (s_max < 2.5)
%       - 5, 10, 15, 20, 25, 30: bandpass (center of preferred speed bin)
%       - 40: highpass (s_max > 30)

% --- Input Validation and Setup ---
assert(ismatrix(temporal_activity), '`temporal_activity` must be an N x T matrix.');

tvec = tvec(:)';
speed = speed(:)';
position = position(:)';

num_neurons = size(temporal_activity, 1);
no_lap = length(start_time_all);
no_track = length(unique(track_ID_all(~isnan(track_ID_all))));

% --- Smoothing ---
w = gausswin(11);
w = w / sum(w);
speed(isnan(speed)) = 0;
speed_smoothed = filtfilt(w, 1, speed);
firing_rate = filtfilt(w, 1, temporal_activity')';

% --- Binning Setup ---
speed_bins = 0:1:40; % Updated speed bins
pos_bins = 0:1:ceil(max(position));
speed_indices = discretize(speed_smoothed, speed_bins);
pos_indices = discretize(position, pos_bins);
num_speed_bins = length(speed_bins) - 1;
num_pos_bins = length(pos_bins) - 1;

% --- 1. Calculate speed_response (1D Speed Tuning) ---
speed_response = nan(num_neurons, num_speed_bins, no_track + 1);
speed_response_std = nan(num_neurons, num_speed_bins, no_track + 1); % Initialize STD matrix
all_tracks_mask = false(size(tvec));

for iT = 1:no_track
    track_mask = false(size(tvec));
    for iLap = 1:no_lap
        if track_ID_all(iLap) == iT
            track_mask = track_mask | (tvec >= start_time_all(iLap) & tvec <= end_time_all(iLap));
        end
    end
    all_tracks_mask = all_tracks_mask | track_mask;

    for iSpeed = 1:num_speed_bins
        speed_mask = speed_indices == iSpeed;
        combined_mask = track_mask & speed_mask;
        time_in_bin = sum(combined_mask);
        if time_in_bin > 0
            % Calculate mean firing rate
            speed_response(:, iSpeed, iT) = mean(firing_rate(:, combined_mask), 2);
            % Calculate standard deviation of firing rate
            speed_response_std(:, iSpeed, iT) = std(firing_rate(:, combined_mask), 0, 2)/sqrt(time_in_bin);
        end
    end
end

non_track_mask = ~all_tracks_mask;
for iSpeed = 1:num_speed_bins
    speed_mask = speed_indices == iSpeed;
    combined_mask = non_track_mask & speed_mask;
    time_in_bin = sum(combined_mask);
    if time_in_bin > 0
        % Calculate mean firing rate for non-track time
        speed_response(:, iSpeed, no_track + 1) = mean(firing_rate(:, combined_mask), 2);
        % Calculate standard deviation for non-track time
        speed_response_std(:, iSpeed, no_track + 1) = std(firing_rate(:, combined_mask), 0, 2)/sqrt(time_in_bin);
    end
end

% --- 2. Calculate speed_pos_FR (2D Speed-Position Tuning) ---
speed_pos_FR = nan(num_neurons, num_speed_bins, num_pos_bins, no_track);
for iT = 1:no_track
    track_mask_pos = false(size(tvec));
    for iLap = 1:no_lap
        if track_ID_all(iLap) == iT
            if iLap < no_lap
                track_mask_pos = track_mask_pos | (tvec >= start_time_all(iLap) & tvec < start_time_all(iLap+1));
            else
                track_mask_pos = track_mask_pos | (tvec >= start_time_all(iLap));
            end
        end
    end

    for iSpeed = 1:num_speed_bins
        speed_mask = speed_indices == iSpeed;
        for iPos = 1:num_pos_bins
            pos_mask = pos_indices == iPos;
            combined_mask_2D = track_mask_pos & speed_mask & pos_mask;
            time_in_bin_2D = sum(combined_mask_2D);
            if time_in_bin_2D > 0
                % Note: std is not calculated for the 2D case as per the request
                speed_pos_FR(:, iSpeed, iPos, iT) = mean(firing_rate(:, combined_mask_2D), 2);
            end
        end
    end
end

% --- 3. Fit Descriptive Function and Classify ---
% Exclude the first near-stationary speed bin from the descriptive fit
% so the skewed-Gaussian classification reflects moving-speed tuning
% rather than the transition out of a stationary state.
avg_speed_response = nanmean(speed_response, 3);

speed_tuning_type = cell(num_neurons, 1);
speed_tuning_code = nan(num_neurons, 1);
speed_bin_centers = speed_bins(1:end-1) + diff(speed_bins)/2;

skewed_gaussian = @(params, s_vec) ...
    params(1) * exp(-(s_vec - params(2)).^2 ./ ...
    ((s_vec < params(2)) .* params(3) + (s_vec >= params(2)) .* params(4)));

options = optimoptions('lsqcurvefit', 'Display', 'off', 'Algorithm', 'trust-region-reflective');

for iNeuron = 1:num_neurons
    y_data = avg_speed_response(iNeuron, :);
    x_data = speed_bin_centers;

    if numel(y_data) >= 2
        y_data = y_data(2:end);
        x_data = x_data(2:end);
    end

    valid_indices = ~isnan(y_data) & y_data >= 0;
    y_data = y_data(valid_indices);
    x_data = x_data(valid_indices);

    if length(y_data) < 4
        speed_tuning_type{iNeuron, 1} = 'untuned';
        speed_tuning_code(iNeuron, 1) = 0;
        continue;
    end

    y_max_guess = max(y_data);
    if isempty(y_max_guess) || y_max_guess == 0; y_max_guess = 1; end
    [~, max_idx] = max(y_data);
    s_max_guess = x_data(max_idx);

    p0 = [y_max_guess, s_max_guess, 10, 10];
    lb = [0, 0, 1e-3, 1e-3];
    ub = [y_max_guess*2.5, 50, 1000, 1000];

    [params, resnorm] = lsqcurvefit(skewed_gaussian, p0, x_data, y_data, lb, ub, options);

    ss_total = sum((y_data - mean(y_data)).^2);
    if ss_total < 1e-9; ss_total = 1; end
    R2 = 1 - (resnorm / ss_total);

    s_max = params(2);

    if R2 < 0.15
        type = 'untuned';
        code = 0;
    elseif s_max < 2.5
        type = 'lowpass';
        code = 2.5;
    elseif s_max > 30
        type = 'highpass';
        code = 40;
    else % Bandpass case: 2.5 <= s_max <= 30
        bandpass_edges = 2.5:5:32.5; % Edges: 2.5, 7.5, 12.5, ..., 32.5
        bandpass_centers = 5:5:30;   % Centers: 5, 10, 15, ..., 30

        % Find which bin the s_max falls into
        bin_idx = discretize(s_max, bandpass_edges);

        if isnan(bin_idx) % Should not happen given the logic, but as a safeguard
            type = 'untuned';
            code = 0;
        else
            code = bandpass_centers(bin_idx);
            type = sprintf('bandpass_%d', code);
        end
    end

    speed_tuning_type{iNeuron, 1} = type;
    speed_tuning_code(iNeuron, 1) = code;
end
end