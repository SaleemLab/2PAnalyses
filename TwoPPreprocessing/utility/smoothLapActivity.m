function smoothedActivity = smoothLapActivity(lapActivity, win)
% smoothLapActivity: Reusable function that applies zero-phase spatial 
% smoothing across position bins for each cell and lap individually, 
% cleanly preserving and isolating NaN missing data masks.
%
% Input:
%   lapActivity - 3D matrix [nCells x nLaps x nPositionBins]
% Output:
%   smoothedActivity - 3D matrix of identical dimensions, spatially smoothed

    % Allocate output matrix to match dimensions
    if nargin < 2, win=15; end 
    smoothedActivity = lapActivity;
    
    nROIs = size(lapActivity, 1);
    nLaps = size(lapActivity, 2);
    
    % Initialize your exact 10-point Gaussian window profile
    w = gausswin(win); 
    w = w / sum(w);
    
    % --- NESTED CELL & LAP TRACE FILTERING LOOP ---
    for iCell = 1:nROIs
        for iLap = 1:nLaps
            % Extract raw spatial position bin trace (Dimension 3)
            trace = squeeze(lapActivity(iCell, iLap, :));
            
            % If the entire lap is empty/unvisited, skip it
            if all(isnan(trace)), continue; end
            
            % Isolate and swap out NaNs to safeguard against filter bleeding
            nanMask = isnan(trace); 
            trace(nanMask) = 0;
            
            % Apply zero-phase filter along the spatial track bins
            smoothed = filtfilt(w, 1, trace); 
            
            % Restore original NaN values back into their unvisited coordinate slots
            smoothed(nanMask) = NaN;
            
            % Reassign to output matrix
            smoothedActivity(iCell, iLap, :) = smoothed;
        end
    end
end