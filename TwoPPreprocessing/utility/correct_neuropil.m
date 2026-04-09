function [signal, regPars, F_binValues, N_binValues] = correct_neuropil(F, F0, N, plane_rate, prctl_F0, prctl_F, Npil_window_F0, numN, minNp, maxNp)
% Based on Sylvia Schroeder's Python function [correct_neuropil] in repository 'depth-for-2p'
% [Ref: correct_neuropil() in preprocess_traces.py] 
% Estimates the correction factor r for neuropil correction, so that:
%    C = S - rN
%    with C: actual signal from the ROI, S: measured signal, N: neuropil
%  RESTORE BASELINE: Add f0_F back to the corrected centered signal 


%Parameters
%----------
%F :  [t x nROIs]
%    Calcium traces (measured signal) of ROIs.
%F0: [t x nROIs]
%    Baseline F of rois 
%N :  [t x nROIs]
%    Neuropil traces of ROIs.
%prctl_F0 : int, optional
%    Percentile of the measured signal that will be taken as F0.
%    The default is 8
%window_F0 : int, optional
%    The window size for the calculation of F0 for both signal and neuropil.
%    The default is 60.
%prctl_F : int, optional
%    Percentile of the measured signal that will be matched to neuropil.
%    The default is 5.
%numN : int, optional
%    Number of bins used to partition the distribution of neuropil values.
%    Each bin will be associated with a mean neuropil value and a mean
%    signal value. The default is 20.
%minNp : int, optional
%    Minimum values of neuropil considered, expressed in percentile.
%    0 < minNp < 100. The default is 10.
%maxNp : int, optional
%    Maximum values of neuropil considered, expressed in percentile.
%    0 < maxNp < 100, minNp < maxNp. The
%    default is 90.


%Returns
%-------
%signal :  [t x nROIs]
%    Neuropil corrected calcium traces.
%regPars :  [2 x nROIs], each row: [intercept, slope]
%    Intercept and slope of linear fits of neuropil (N) to measured calcium
%    traces (F)
%F_binValues : [numN, nROIs]
%    Low percentile (prctl_F) values for each calcium trace bin. These
%    values were used for linear regression.
%N_binValues : [numN, nROIs]
%    Values for each neuropil bin. These values were used for linear
%    regression.

% Originally based on Matlab function estimateNeuropil (in +preproc) written by Mario Dipoppa and Sylvia Schroeder

if nargin < 5 || isempty(prctl_F0), prctl_F0 = 8; end 
if nargin < 6 || isempty(prctl_F), prctl_F = 5; end % why is the prctile f lower than for f0? @Sylvia  
if nargin < 7 || isempty(Npil_window_F0), Npil_window_F0 = 60; end % the same window size as F0_window size [60s] 
if nargin < 8 || isempty(numN), numN = 20; end
if nargin < 9 || isempty(minNp), minNp = 10; end
if nargin < 10 || isempty(maxNp), maxNp = 90; end

% Initialise outputs
[nt, nROIs] = size(F);

N_binValues = NaN(numN, nROIs);
F_binValues = NaN(numN, nROIs);
regPars = NaN(2, nROIs);
signal = NaN(nt, nROIs);

% Correct for slow drift in rois and neuropil traces separately
% Compute N0 inside this function; 
% F0 = get_F0(F, prctl_F0, Npil_window_F0, plane_rate); % using the same percentile for F and FNeu [defulat 8th percentile]
N0 = get_F0(N, prctl_F0, Npil_window_F0, plane_rate);
Fc = F - F0;
Nc = N - N0;

for iROI = 1:nROIs
    iN = Nc(:, iROI);
    iF = Fc(:, iROI);

    % Get range of neuropil values (default: between 10th and 90th percentile) + divide range into numN
    % (default: 20) binsn of equal width.
    % NaNs can be present for boutons sessions after z-motion correction pipeline: 'bad_frames' are converted to nans.  
    iN_NoNaN = iN(~isnan(iN)); 
    
    N_prct = prctile(iN_NoNaN, [minNp, maxNp],1); %first dimension
    binSize = (N_prct(2) - N_prct(1)) / numN;

    % Exact bin values (start of bin)
    N_binValues(:, iROI) = N_prct(1) + (0:numN-1)' * binSize;

    % Associate each neuropil value with a bin number.
    N_ind = floor((iN - N_prct(1)) / binSize) + 1;

    % Bin ROI signal values the same way as neuropil values. For each bin, find the prctl_F value (default: 5). The
    % idea is to match ROI signal values without spiking activity to the corresponding neuropil values.
   for Ni = 1:numN
        mask = (N_ind == Ni);
        if any(mask)
            % Clean the bin before calculating percentile
            temp_F = iF(mask);
            temp_F = temp_F(~isnan(temp_F));
            if ~isempty(temp_F)
                F_binValues(Ni, iROI) = prctile(temp_F, prctl_F, 1); % dim 1 
            end
        end
    end

    % Determine relation between neuropil and ROI signal using linear regression.
    noNan = ~isnan(F_binValues(:, iROI));
    if sum(noNan) > 1
        p = polyfit(N_binValues(noNan, iROI), F_binValues(noNan, iROI), 1); % Degree of polynomial fit=1
        b = p(1); 
        a = p(2); 

        % Restrict slope: allows up to 2.0
        % (@Aman - the allen white paper restricts slope to 0.7 to prevent over-correcting)
        b = min(max(b, 0), 2.0);
        regPars(:, iROI) = [a; b];

        % signal = (iF - (b * iN + a)) + F0
        % Correct ROI signal by subtracting prediction based on neuropil.
        signal(:, iROI) = (iF - (b * iN + a)) + F0(:, iROI);
    end
end
end
