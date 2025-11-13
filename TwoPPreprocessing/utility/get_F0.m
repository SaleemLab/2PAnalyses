function F0 = get_F0(Fc, fs, prctl_F, window_size)
    % Determines the baseline fluorescence (F0) for computing deltaF/F.
    %
    % This is the FASTEST version, using the `ordfilt2` function from the
    % Image Processing Toolbox to perform a highly optimized rolling percentile.
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
    %     The rolling window over which to calculate F0, in seconds. Default is 180.
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
    
    % Check if the required toolbox is available
    if ~license('test', 'image_toolbox')
        error('This fast version of get_F0 requires the Image Processing Toolbox.');
    end

    % 1. Translate window size from seconds into an odd number of frames
    window_frames = round(fs * window_size);
    if mod(window_frames, 2) == 0
        window_frames = window_frames + 1;
    end
    
    % 2. Determine the order (k-th smallest value) needed for the percentile
    % For a percentile 'p' in a window of size 'N', the order is (p/100)*N.
    order = round((prctl_F / 100) * window_frames);
    order = max(order, 1); % The order must be at least 1.
    
    % 3. Define the filtering neighborhood (a 1D vertical window)
    domain = ones(window_frames, 1);
    
    % 4. Apply the 2D order-statistic filter. It's designed for images but
    % works perfectly on our [time x ROIs] matrix with a 1D domain.
    % The 'symmetric' option handles the edges of the data correctly.
    F0 = ordfilt2(Fc, order, domain, 'symmetric');
end