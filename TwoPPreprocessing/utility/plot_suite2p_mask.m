function plot_suite2p_mask(roi_idx, stat, ops)
    % Get image dimensions
    Ly = double(ops.Ly);
    Lx = double(ops.Lx);

    % ROI pixels (Convert to double to avoid class errors)
    x_roi = double(stat{roi_idx}.xpix) + 1;
    y_roi = double(stat{roi_idx}.ypix) + 1;

    %
    % Convert to double before performing division/mod
    neu_idx_py = double(stat{roi_idx}.neuropil_mask); 
    
    % [y, x]
    y_neu = floor(neu_idx_py / Ly) + 1;
    x_neu = mod(neu_idx_py, Lx) + 1;

    % 4. Create the Figure
    figure('Color', 'w', 'Name', ['ROI ' num2str(roi_idx) 'F & Neu Masks']);
    
    % Background
    if isfield(ops, 'meanImgE')
        imagesc(ops.meanImgE); 
    elseif isfield(ops, 'max_proj') 
        imagesc(ops.max_proj);
    end
    colormap('gray'); hold on;

    % neuropil in cyan 
    scatter(x_neu, y_neu, 3, 'c', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Neuropil Mask');
    
    % ori mask in red 
    scatter(x_roi, y_roi, 6, 'r', 'filled', 'DisplayName', 'ROI Mask');

    %  median
    center_y = double(stat{roi_idx}.med(1)) + 1;
    center_x = double(stat{roi_idx}.med(2)) + 1;
    
    padding = 60; 
    xlim([center_x - padding, center_x + padding]);
    ylim([center_y - padding, center_y + padding]);

    axis image; 
    set(gca, 'YDir', 'reverse'); 
    title(['ROI ' num2str(roi_idx) ':Mask Alignment']);
    xlabel('Pixels (X)'); ylabel('Pixels (Y)');
    
    legend('Location', 'northeast');
end