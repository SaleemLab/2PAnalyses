function [signal, regPars, F_binValues, N_binValues] = correct_neuropil(F, N, fs, numN, minNp, maxNp, prctl_F, prctl_F0, Npil_window_F0, verbose)
    % Estimates the correction factor r for neuropil correction, so that:
    %     C = S - rN
    %     with C: actual signal from the ROI, S: measured signal, N: neuropil 
    % Estimates the correction factor r for neuropil correction, so that:
    %     Parameters
    %     ----------
    %     F : np.ndarray [t x nROIs]
    %         Calcium traces (measured signal) of ROIs.
    %     N : np.ndarray [t x nROIs]
    %         Neuropil traces of ROIs.
    %     numN : int, optional
    %         Number of bins used to partition the distribution of neuropil values.
    %         Each bin will be associated with a mean neuropil value and a mean
    %         signal value. The default is 20.
    %     minNp : int, optional
    %         Minimum values of neuropil considered, expressed in percentile.
    %         0 < minNp < 100. The default is 10.
    %     maxNp : int, optional
    %         Maximum values of neuropil considered, expressed in percentile.
    %         0 < maxNp < 100, minNp < maxNp. The
    %         default is 90.
    %     prctl_F : int, optional
    %         Percentile of the measured signal that will be matched to neuropil.
    %         The default is 5.
    %     prctl_F0 : int, optional
    %         Percentile of the measured signal that will be taken as F0.
    %         The default is 8
    %     window_F0 : int, optional
    %         The window size for the calculation of F0 for both signal and neuropil.
    %         The default is 60.
    %     verbose : boolean, optional
    %         Feedback on fitting. The default is True.
    % 
    %     Returns
    %     -------
    %     signal : np.ndarray [t x nROIs]
    %         Neuropil corrected calcium traces.
    %     regPars : np.ndarray [2 x nROIs], each row: [intercept, slope]
    %         Intercept and slope of linear fits of neuropil (N) to measured calcium
    %         traces (F)
    %     F_binValues : np.array [numN, nROIs]
    %         Low percentile (prctl_F) values for each calcium trace bin. These
    %         values were used for linear regression.
    %     N_binValues : np.array [numN, nROIs]
    %         Values for each neuropil bin. These values were used for linear
    %         regression.
    % 
    %     Based on Matlab function estimateNeuropil (in +preproc) written by Mario
    %     Dipoppa and Sylvia Schroeder
        %%%

    if nargin < 4, numN = 20; end
    if nargin < 5, minNp = 10; end
    if nargin < 6, maxNp = 90; end
    if nargin < 7, prctl_F = 5; end
    if nargin < 8, prctl_F0 = 5; end
    if nargin < 9, Npil_window_F0 = 60; end
    if nargin < 10, verbose = true; end
    
    [nt, nROIs] = size(F);
    N_binValues = NaN(numN, nROIs);
    F_binValues = NaN(numN, nROIs);
    regPars = NaN(2, nROIs);
    signal = NaN(nt, nROIs);

    % Compute F0 traces for Calcium and Neuropil traces
    F0 = get_F0(F, fs, Npil_window_F0);
    N0 = get_F0(N, fs, Npil_window_F0);
    
    % Correct for slow drift
    Fc = F - F0;
    Nc = N - N0;

    % Determine where the minimum normalized difference between F0 and N0 occurs
    [~, ti] = min((F0 - N0) ./ N0, [], 1, 'omitnan');
    
    for iROI = 1:nROIs
        iN = Nc(:, iROI);
        iF = Fc(:, iROI);

        % Get low and high percentile of neuropil trace
        N_prct = prctile(iN, [minNp, maxNp]);
        binSize = (N_prct(2) - N_prct(1)) / numN;
        N_binValues(:, iROI) = N_prct(1) + (0:numN-1)' * binSize;

        % Discretize neuropil values into bins
        N_ind = floor((iN - N_prct(1)) / binSize);
        
        % Find the matching low percentile value from F trace for each neuropil bin
        for Ni = 1:numN
            tmp = NaN(size(iF));
            tmp(N_ind == Ni-1) = iF(N_ind == Ni-1);
            F_binValues(Ni, iROI) = prctile(tmp(~isnan(tmp)), prctl_F);

        end
        
        % Fit non-NaN values
        validIdx = ~isnan(F_binValues(:, iROI)) & ~isnan(N_binValues(:, iROI));
        if any(validIdx)
            [a, b, ~] = linear_analytical_solution(N_binValues(validIdx, iROI), F_binValues(validIdx, iROI), false);
            b = min(max(b, 0), 2);
            regPars(:, iROI) = [a; b];
            corrected_sig = iF - (b * iN + a) + F0(:, iROI);
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
