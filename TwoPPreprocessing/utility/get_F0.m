function F0 = get_F0(Fc, fs, prctl_F, window_size)
    % Determines the baseline fluorescence to use for computing deltaF/F.
    % Function converted to Matlab from Schroeder Lab 
    %
    % Parameters:
    % ----------
    % Fc : matrix [t x nROIs]
    %     Calcium traces (measured signal) of ROIs.
    % fs : float
    %     The frame rate (frames/second/plane).
    % prctl_F : int, optional
    %     The percentile from which to take F0. The default is 8.
    % window_size : int, optional
    %     The rolling window over which to calculate F0. The default is 60.
    % framesPerFolder : array, optional
    %     An array with the number of frames in each experiment. If not empty,
    %     then gets individual F0 for each experiment. Default is empty.
    % verbose : bool, optional
    %     Whether or not to provide detailed processing information.
    %
    % Returns:
    % -------
    % F0 : matrix [t x nROIs]
    %     The baseline fluorescence (F0) traces for each ROI.
    if nargin < 3 || isempty(prctl_F)
        prctl_F = 8;
    end
    if nargin < 4 || isempty(window_size)
        window_size = 60;
    end
  
    % Initialize F0 array
    F0 = zeros(size(Fc));

    % Translate window size from seconds into frames
    window_size = round(fs * window_size);
    
    % Compute rolling percentile for baseline fluorescence
    for roi = 1:size(Fc, 2)
        for t = 1:size(Fc, 1)
            win_start = max(1, t - floor(window_size / 2));
            win_end = min(size(Fc, 1), t + floor(window_size / 2));
            F0(t, roi) = prctile(Fc(win_start:win_end, roi), prctl_F);
        end
    end
end
