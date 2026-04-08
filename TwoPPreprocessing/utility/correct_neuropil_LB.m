%% Older version (which was apparently liad's version! and not what was described in the paper)
function [signal, regPars, F_binValues, N_binValues] = correct_neuropil_LB(F, N, fs, numN, minNp, maxNp, prctl_F, prctl_F0, Npil_window_F0, verbose)
    % Estimates the correction factor r for neuropil correction.
    % This version is OPTIMIZED to use `accumarray` to avoid a slow inner loop.
    %
    % Parameters
    % ----------
    % F : np.ndarray [t x nROIs]
    %     Calcium traces (measured signal) of ROIs.
    % N : np.ndarray [t x nROIs]
    %     Neuropil traces of ROIs.
    % ... (rest of parameters) ...
    %
    % Returns
    % -------
    % signal : np.ndarray [t x nROIs]
    %     Neuropil corrected calcium traces.
    % regPars : np.ndarray [2 x nROIs]
    %     Intercept and slope of linear fits.
    % F_binValues : np.array [numN, nROIs]
    %     Low percentile values for each calcium trace bin.
    % N_binValues : np.array [numN, nROIs]
    %     Values for each neuropil bin.

    if nargin < 4, numN = 20; end
    if nargin < 5, minNp = 10; end
    if nargin < 6, maxNp = 90; end
    if nargin < 7, prctl_F = 5; end
    if nargin < 8, prctl_F0 = 5; end
    if nargin < 9, Npil_window_F0 = 180; end
    if nargin < 10, verbose = true; end

    [nt, nROIs] = size(F);
    N_binValues = NaN(numN, nROIs);
    F_binValues = NaN(numN, nROIs);
    regPars = NaN(2, nROIs);
    signal = NaN(nt, nROIs);

    % This now calls your new, fast version of get_F0
    F0 = get_F0(F, fs, prctl_F0, Npil_window_F0);
    N0 = get_F0(N, fs, prctl_F0, Npil_window_F0);

    % Correct for slow drift
    Fc = F - F0;
    Nc = N - N0;

    for iROI = 1:nROIs
        iN = Nc(:, iROI);
        iF = Fc(:, iROI);

        % Get low and high percentile of neuropil trace
        N_prct = prctile(iN, [minNp, maxNp]);
        binSize = (N_prct(2) - N_prct(1)) / numN;

        if binSize <= 0 % Avoid division by zero if neuropil is flat
            continue; 
        end

        % Vectorized calculation of bin indices (+1 to make them 1-based)
        N_ind = floor((iN - N_prct(1)) / binSize) + 1;

        % Create a mask for valid data points (within bin range and not NaN)
        valid_mask = N_ind >= 1 & N_ind <= numN & ~isnan(iF);

        % Use accumarray to calculate percentile for all bins at once.
        % This replaces the entire inner 'for Ni = 1:numN' loop.
        F_binValues(:, iROI) = accumarray(N_ind(valid_mask), iF(valid_mask), [numN 1], ...
                                          @(x) prctile(x, prctl_F), NaN);

        % Use bin centers for a more accurate regression
        N_binValues(:, iROI) = N_prct(1) + ((1:numN)' - 0.5) * binSize; 

        % Fit non-NaN values
        validIdx = ~isnan(F_binValues(:, iROI));
        if any(validIdx)
            [a, b, ~] = linear_analytical_solution(N_binValues(validIdx, iROI), F_binValues(validIdx, iROI), false);
            b = min(max(b, 0), 2); % Constrain the slope
            regPars(:, iROI) = [a; b];

            % Apply correction to the original, non-F0-subtracted traces
            corrected_sig = F(:, iROI) - b * N(:, iROI) - a;
            signal(:, iROI) = corrected_sig;
        end
    end
end

function [a, b, mse] = linear_analytical_solution(x, y, noIntercept)
    % Fits a robust line to data using least squares.
    %
    % Inputs:
    %   x - [n x 1] or [1 x n] array of x-values
    %   y - [n x 1] or [1 x n] array of y-values
    %   noIntercept - logical (optional), if true, fit without intercept
    %
    % Outputs:
    %   a - intercept
    %   b - slope
    %   mse - mean squared error of fit

    if nargin < 3
        noIntercept = false;
    end

    x = x(:);  % Ensure column vectors
    y = y(:);
    n = length(x);

    if noIntercept
        b = sum(x .* y) / sum(x .^ 2);
        a = 0;
    else
        sum_x = sum(x);
        sum_y = sum(y);
        sum_x2 = sum(x .^ 2);
        sum_xy = sum(x .* y);
        denom = n * sum_x2 - sum_x^2;

        if denom == 0
            a = NaN;
            b = NaN;
            mse = NaN;
            return;
        end

        a = (sum_y * sum_x2 - sum_x * sum_xy) / denom;
        b = (n * sum_xy - sum_x * sum_y) / denom;
    end

    y_fit = a + b * x;
    mse = mean((y - y_fit) .^ 2);
end