function F0 = get_F0(signal, F0_percentile, window_size_sec, plane_rate)
% % Based on Sylvia Schroeder's Python function [get_F0] in repository
% 'depth-for-2p'; This version is adapted to use with ordfilt2 for speed.
% The python version usues pd.rolling() [pandas dataframe]
% Determines the baseline fluorescence (F0) using optimised rank filtering.


% signal: [frames x ROIs]
% F0_percentile: 8 [Sylvia's defult: 'depth-for-2p']
% window_size_sec: usually 60 [Sylvia's defult: 'depth-for-2p']
if nargin < 2 || isempty(F0_percentile), F0_percentile = 8; end
if nargin < 3 || isempty(window_size_sec), window_size_sec = 60; end
% Here we pass the interpolatd rate [60hz] 


nWindow = round(window_size_sec * plane_rate);

% Determine the rank (index) for the required percentile
% ordfilt2 uses 1-based indexing for the rank
rank_idx = max(1, fix(nWindow * F0_percentile / 100));

% Create a domain for the filter (1D window applied across time)
domain = true(nWindow, 1);


% Apply ordfilt2 to each ROI
% ordfilt2 expects 2D input. We process columns to maintain [frames x ROIs]
% We pad the signal to handle edge effects similar to sliding windows
% 'symmetric' padding 
Fc_q = ordfilt2(signal, rank_idx, domain, 'symmetric');


% Smoothing step
% if using plane_rate at 7.5 * 60 = 450 would be ok; 
% Since we're running this function on 60hz data I have changed this line of
% code to smooth slightly otherwise the baseline looked like a flat line @Aman 
smoothWin = round(5 * plane_rate); %60*5=300
F0 = smoothdata(Fc_q, 1, 'gaussian', smoothWin);

end

% THIS VERSION TAKE TOO LONG: Iterates through each frame to compute 8th
% percentile in a window of 60s 
% function F0 = get_F0(signal, F0_percentile, window_size_sec, fs)
% % Determines the baseline fluorescence (F0) for computing deltaF/F.
% % signal: [frames x ROIs]
% % F0_percentile: usually 5 [sylvia's defult]
% % window_size_sec: usually 60 [sylvia's defult]
% % From cortex lab: removeSlowDrift from +preproc [removed the filetered
% % trace component and smoothed F0 similar to Sylvia's function in Python from 'depth-for-2p'; currently private repo] 
% % [https://github.com/sylviaschroeder/CortexLab/blob/master/%2Bpreproc/removeSlowDrift.m]
% % 
% 
% if nargin < 2 || isempty(F0_percentile), F0_percentile = 5; end
% if nargin < 3 || isempty(window_size_sec), window_size_sec = 60; end
% 
% [nFrames, nROIs] = size(signal);
% 
% % extract centered window 
% n = round(window_size_sec * fs);
% if mod(n,2) == 0
%     n = n + 1; % Ensure window is odd for perfect centering
% end
% nBefore = floor((n-1)/2);
% nAfter = n - nBefore - 1;
% 
% Fc_q = zeros(nFrames, nROIs);
% 
% % Check for parallel pool
% poolObj = gcp('nocreate');
% 
% if isempty(poolObj)
%     fprintf('Calculating percentile frame-by-frame (Serial)...');
%     for k = 1:nFrames
%         % Exact slice logic from removeSlowDrift
%         tmpTraces = signal(max(1, k-nBefore) : min(nFrames, k+nAfter), :);
%         Fc_q(k, :) = prctile(tmpTraces, F0_percentile, 1);
%     end
% else
%     fprintf('Calculating percentile frame-by-frame (Parallel)...');
%     parfor k = 1:nFrames
%         % Exact slice logic from removeSlowDrift
%         tmpTraces = signal(max(1, k-nBefore) : min(nFrames, k+nAfter), :);
%         Fc_q(k, :) = prctile(tmpTraces, F0_percentile, 1);
%     end
% end
% 
% % sigma_val = 450 (Matches original 7.5Hz logic for a 60s window)
% sigma_val = 450; 
% % Using 'smoothdata' with a window of 8*sigma to match SciPy truncation
% F0 = smoothdata(Fc_q, 1, 'gaussian', 8 * sigma_val);
% 
% end
% 
