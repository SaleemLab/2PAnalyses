function [tuningCurves, positionAxis] = calculateTuningFromBinned(lapPositionActivity, occupancyMatrix)
    % calculateTuningFromBinned Creates tuning curves from pre-binned data.
    %
    % INPUTS:
    %   lapPositionActivity    - 3D array (neurons x laps x position bins) of mean activity.
    %   occupancyMatrix        - 2D array (laps x position bins) of time spent in each bin.
    %
    % OUTPUTS:
    %   tuningCurves           - 2D array (neurons x position bins) of the final
    %                            smoothed and normalized tuning curves.
    %   positionAxis           - A 1D vector for plotting the x-axis.

    %% Sum Activity and Occupancy Across Laps 
    totalActivity = squeeze(sum(lapPositionActivity, 2));
    totalOccupancy = sum(occupancyMatrix, 1);

    %%  Perform Occupancy Normalization 
    rawTuningCurves = totalActivity ./ totalOccupancy;
    
    unvisitedBins = (totalOccupancy == 0);
    rawTuningCurves(:, unvisitedBins) = 0;

    %% Apply Spatial Smoothing ---
    binSize = 1; 
    smoothingWindow = 5; % cm, as per Diamanti paper
    smoothingBins = round(smoothingWindow / binSize);
    
    tuningCurves = smoothdata(rawTuningCurves, 2, 'gaussian', smoothingBins);

    %% Prepare Axis for Plotting ---
    numBins = size(lapPositionActivity, 3);
    positionAxis = (1:numBins) - (1-binSize)/2; % Bin centers
end