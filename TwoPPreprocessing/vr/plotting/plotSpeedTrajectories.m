function figA = plotSpeedTrajectories(response, strat, metrics, smoothSigma)
    % Replicates Extended Data Figure 5a with 3 subplots for clarity
    if nargin < 4, smoothSigma = 2; end

    figA = figure('Name', 'Panel A: Stratified Speed Trajectories', 'Position', [100, 100, 600, 800]);
    
    numLaps = length(response.lapRunningSpeed);
    groupColors = [0.8 0 0.8; 0.2 0.2 0.2; 0 0.8 0.8]; % High (Mag), Med (Blk), Low (Cyan)
    titles = {'High Speed (>20 cm/s)', 'Medium Speed (10-20 cm/s)', 'Low Speed (1-10 cm/s)'};
    
    % Prepare the categorical indices in an array for looping
    masks = {strat.highIdx, strat.medIdx, strat.lowIdx};
    thresholds = [20, 50; 10, 20; 1, 10]; % [min, max] for shading/y-axis

    for i = 1:3
        subplot(3, 1, i); hold on;
        currentMask = masks{i};
        lapIndices = find(currentMask);
        
        for l = lapIndices'
            rawV = response.lapRunningSpeed{l};
            % Resample to 200 bins
            binnedV = interp1(linspace(1, 200, length(rawV)), rawV, 1:200, 'linear', 'extrap');
            % Smooth
            smoothedV = smoothdata(binnedV, 'gaussian', smoothSigma * 3);
            
            % Plot with specific group color
            plot(1:200, smoothedV, 'Color', [groupColors(i,:), 0.3], 'LineWidth', 1.5);
            
            % Plot the robust trimmean dot at the end for verification
            plot(205, metrics.robustSpeedsPerLap(l), 'o', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', groupColors(i,:), 'MarkerSize', 4);
        end
        
        % Formatting each subplot
        yline(thresholds(i,1), '--r');
        if i < 3, yline(thresholds(i,2), '--r'); end
        
        ylabel('Speed (cm/s)');
        title([titles{i}, ' - n=', num2str(length(lapIndices))]);
        grid on; box off;
        xlim([0 210]); 
        ylim([0 60]); % Keep Y-axis same across all to see relative difference
    end
    xlabel('Position Bin (1-200)');
end

% function figA = plotSpeedTrajectories(response, strat, metrics)
%     % Replicates Extended Data Figure 5a with Gaussian smoothing
%     % smoothSigma: standard deviation for smoothing (default = 2 bins)
% 
%     if nargin < 4, smoothSigma = 2; end
% 
%     figA = figure('Name', 'Panel A: Smoothed Speed Trajectories');
%     hold on;
% 
%     numLaps = length(response.lapRunningSpeed);
% 
%     % Colors matching paper: Cyan (Low), Black (Med), Magenta (High)
%     groupColors = [0 0.8 0.8; 0.2 0.2 0.2; 0.8 0 0.8; 0.7 0.7 0.7]; 
% 
%     for l = 1:numLaps
%         % 1. Determine which speed group the lap belongs to
%         if strat.lowIdx(l), colorID = 1;
%         elseif strat.medIdx(l), colorID = 2;
%         elseif strat.highIdx(l), colorID = 3;
%         else, colorID = 4; % Excluded (< 1 cm/s)
%         end
% 
%         % 2. Process the Speed Vector
%         rawV = response.lapRunningSpeed{l};
% 
%         % Interpolate raw treadmill samples to the 200 position bins
%         % This ensures the x-axis (1:200) matches your dFF matrix
%         binnedV = interp1(linspace(1, 200, length(rawV)), rawV, 1:200, 'linear', 'extrap');
% 
%         % 3. Apply Gaussian Smoothing to the line
%         % Window size is typically 3*sigma to capture the full curve
%         smoothedV = smoothdata(binnedV, 'gaussian', smoothSigma * 3);
% 
%         % 4. Plot the Trajectory with transparency (Alpha = 0.25)
%         % This allows you to see the "density" of the lines
%         plot(1:200, smoothedV, 'Color', [groupColors(colorID,:), 0.25], 'LineWidth', 2);
%     end
% 
%     % --- Shading Bands (Thresholds) ---
%     % Add the 1, 10, and 30 cm/s lines
%     yline([1, 10, 30], '--r', {'1 cm/s', '10 cm/s', '20 cm/s'}, ...
%         'LabelHorizontalAlignment', 'left', 'Alpha', 0.6, 'LineWidth', 1);
% 
%     % --- Verification Dots ---
%     % Plot the 'Robust Speed' (trimmean) value at the very end
%     for l = 1:numLaps
%         if strat.lowIdx(l), c = groupColors(1,:);
%         elseif strat.medIdx(l), c = groupColors(2,:);
%         elseif strat.highIdx(l), c = groupColors(3,:);
%         else, c = groupColors(4,:);
%         end
%         % We place these at 205 to keep them clear of the corridor plot
%         plot(205, metrics.robustSpeedsPerLap(l), 'o', ...
%             'MarkerEdgeColor', 'none', 'MarkerFaceColor', c, 'MarkerSize', 5);
%     end
% 
%     % --- Formatting ---
%     xlabel('Position Bin (1-200)');
%     ylabel('Running Speed (cm/s)');
%     title(sprintf('Stratification: %d Low, %d Med, %d High', ...
%         sum(strat.lowIdx), sum(strat.medIdx), sum(strat.highIdx)));
%     grid on; box off;
%     xlim([0 215]); ylim([-2 60]);
% end