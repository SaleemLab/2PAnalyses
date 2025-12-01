function plot_discarded_groups(roisToDiscard, groups, ...
                               planeROIIndices, planeYpixAll, planeXpixAll, imgH, imgW)
% PLOT_DISCARDED_GROUPS plots the boundaries of ONLY the discarded correlated ROIs, 
% color-coded by their group ID, to visually confirm spatial overlap.

% Define a diverse color map (Subset of MATLAB's colors)
colorMap = [
    0.8500, 0.3250, 0.0980; % Red-Orange
    0.0000, 0.4470, 0.7410; % Blue
    0.9290, 0.6940, 0.1250; % Yellow
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330; % Cyan
    0.6350, 0.0780, 0.1840; % Dark Red
    0.75, 0.75, 0;          % Olive
];

% 1. Identify ONLY the discarded ROIs in the current plane
planeROIsToDiscard = intersect(roisToDiscard, planeROIIndices);

fprintf('Plotting %d discarded ROIs (color-coded by group ID)...\n', length(planeROIsToDiscard));

if ~iscell(planeYpixAll) || ~iscell(planeXpixAll)
    error('Internal error: Pixel arrays must be cell arrays.');
end

% 2. Iterate through all discarded ROIs in the plane
for i = 1:length(planeROIsToDiscard)
    globalIdx = planeROIsToDiscard(i);
    
    % --- Determine Color ---
    groupID = groups(globalIdx);
    colorIndex = mod(groupID - 1, size(colorMap, 1)) + 1;
    groupColor = colorMap(colorIndex, :);
    
    % --- Retrieve Position and Plot Boundary ---
    localIdx = find(planeROIIndices == globalIdx, 1); 

    if ~isempty(localIdx)
        ypix_list_int = planeYpixAll{localIdx};
        xpix_list_int = planeXpixAll{localIdx};
        
        % Convert coordinates to 'double' for poly2mask
        ypix_list = double(ypix_list_int);
        xpix_list = double(xpix_list_int);

        % Create binary mask and find boundaries
        mask = poly2mask(xpix_list, ypix_list, imgH, imgW);
        B = bwboundaries(mask);
        
        for k = 1:length(B)
            boundary = B{k};
            % Plot the boundary (Line width increased for visibility)
            plot(boundary(:, 2), boundary(:, 1), 'Color', groupColor, 'LineWidth', 2);
        end
    end
end

end % End of function