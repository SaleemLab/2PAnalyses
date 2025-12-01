function plot_correlated_groups(roisToKeep, roisToDiscard, groups, ...
                                planeROIIndices, planeYpixAll, planeXpixAll, imgH, imgW)
% PLOT_CORRELATED_GROUPS plots the boundaries of correlated ROIs 
% (kept and discarded) using a unique color for each group ID.

% Define a diverse color map 
colorMap = [
    0.8500, 0.3250, 0.0980; 0.0000, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 
    0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 
    0.6350, 0.0780, 0.1840; 0.75, 0.75, 0; 0, 0.5, 0; 0, 0.75, 0.75;
];

% 1. Identify ALL ROIs in the plane that belong to a non-singleton group
correlatedROIsGlobal = union(roisToKeep, roisToDiscard);
planeCorrelatedROIs = intersect(correlatedROIsGlobal, planeROIIndices);

% 2. Get all ROI positions for the target plane (used for indexing below)
if ~iscell(planeYpixAll) || ~iscell(planeXpixAll)
    error('Internal error: Pixel arrays must be cell arrays.');
end

% 3. Iterate through all correlated ROIs in the plane
for i = 1:length(planeCorrelatedROIs)
    globalIdx = planeCorrelatedROIs(i);
    
    % --- Determine Color ---
    groupID = groups(globalIdx);
    colorIndex = mod(groupID - 1, size(colorMap, 1)) + 1;
    groupColor = colorMap(colorIndex, :);
    
    % --- Retrieve Position and Plot Boundary ---
    localIdx = find(planeROIIndices == globalIdx, 1); 

    if ~isempty(localIdx)
        ypix_list_int = planeYpixAll{localIdx};
        xpix_list_int = planeXpixAll{localIdx};
        
        % CRUCIAL FIX: Convert coordinates to 'double' for poly2mask
        ypix_list = double(ypix_list_int);
        xpix_list = double(xpix_list_int);

        % Create binary mask and find boundaries
        mask = poly2mask(xpix_list, ypix_list, imgH, imgW);
        B = bwboundaries(mask);
        
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:, 2), boundary(:, 1), 'Color', groupColor, 'LineWidth', 1.5);
        end
        
        % --- Mark Winner ROI with a thick circle (Visual Confirmation) ---
        if ismember(globalIdx, roisToKeep)
            % Find the centroid of the mask for marking
            stats = regionprops(mask, 'Centroid');
            if ~isempty(stats)
                centroid = stats(1).Centroid;
                % Plot a magenta circle over the centroid
                plot(centroid(1), centroid(2), 'mo', 'MarkerSize', 8, 'LineWidth', 2);
            end
        end
    end
end

end % End of function