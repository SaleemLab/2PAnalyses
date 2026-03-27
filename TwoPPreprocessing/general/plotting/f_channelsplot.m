fall = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Processed\20260313\suite2p\plane_z\Fall.mat")
if iscell(fall.ops)
    ops_p = fall.ops{1}; 
    stat_p = fall.stat{1};
else
    ops_p = fall.ops;
    stat_p = fall.stat;
end

% 2. Define the ROI index you want to plot (e.g., ROI 1)
roi_idx = 300; 

figure('Name', 'ROI Comparison', 'Color', 'w', 'Position', [100 100 1000 450]);

% Define channel fields and matching colors from your image
channels = {'meanImg', 'meanImg_chan2'};
titles = {['ROI ' num2str(roi_idx) ' on Ch1 (Mean)'], ...
          ['ROI ' num2str(roi_idx) ' on Ch2 (Structural)']};
outline_colors = {'g', 'r'}; % Green for Ch1, Red for Ch2

for c = 1:2
    fieldName = channels{c};
    
    if isfield(ops_p, fieldName)
        subplot(1, 2, c);
        
        % Plot Background Image
        img = ops_p.(fieldName);
        imagesc(img); 
        colormap(gca, 'gray');
        hold on;
        
        % 3. Plot the specific ROI outline
        % Get coordinates and adjust for MATLAB 1-based indexing (+1)
        ypix = double(stat_p{roi_idx}.ypix) + 1;
        xpix = double(stat_p{roi_idx}.xpix) + 1;
        
        % Use convhull to create the clean boundary line seen in your image
        if ~isempty(ypix)
            k = convhull(xpix, ypix);
            plot(xpix(k), ypix(k), outline_colors{c}, 'LineWidth', 2);
        end
        
        title(titles{c});
        axis image; % Keep aspect ratio square like your image
        set(gca, 'FontSize', 12);
    end
end

sgtitle('Comparing ROI across Channels', 'FontSize', 16);