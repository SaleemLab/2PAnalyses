function figA = plotSpeedTrajectories(response, smoothSigma)
% Loads lap running speed data and plots traces across position bins.
% Colour patches use the fixed speed edges saved by getLowMedHigh_SpeedPosActivityMatrix
% so the shading is consistent with the analysis (not recomputed from this session's quantiles).
% Y-axis ticks: 1, 15, 27, and max.

if nargin < 2, smoothSigma = 2; end

%% Interpolate and smooth lap speeds
numLaps    = length(response.lapRunningSpeed);
numPosBins = 200;
allInterpSpeeds = nan(numLaps, numPosBins);

for l = 1:numLaps
    rawV = response.lapRunningSpeed{l};
    if length(rawV) < 5, continue; end
    binnedV = interp1(linspace(1, numPosBins, length(rawV)), rawV, 1:numPosBins, 'linear', 'extrap');
    allInterpSpeeds(l, :) = smoothdata(binnedV, 'gaussian', smoothSigma * 3);
end

%% Get speed edges from saved analysis
% Use the edges saved by getLowMedHigh_SpeedPosActivityMatrix so patches
% match exactly what was used to split the data
if isfield(response, 'speedPositionActivity') && isfield(response.speedPositionActivity, 'speedEdges')
    speedEdges = response.speedPositionActivity.speedEdges;
    lowThresh  = speedEdges(2);   % e.g. 15 cm/s
    highThresh = speedEdges(3);   % e.g. 27 cm/s
    fprintf('Using saved speed edges: low < %.1f | med %.1f-%.1f | high > %.1f\n', ...
        lowThresh, lowThresh, highThresh, highThresh);
else
    % fallback: recompute from session quantiles with a warning
    warning('No saved speedEdges found in response — falling back to session quantiles. Run getLowMedHigh_SpeedPosActivityMatrix first.');
    runningSpeeds  = allInterpSpeeds(allInterpSpeeds > 1 & ~isnan(allInterpSpeeds));
    sessionQuantiles = quantile(runningSpeeds, linspace(0, 1, 4));
    lowThresh  = sessionQuantiles(2);
    highThresh = sessionQuantiles(3);
end

%% Axis limits
minY = 0;
maxY = max(allInterpSpeeds(:)) + 1;
if isnan(maxY), maxY = 50; end

%% Plot
figA = figure('Name', 'Speed Trajectories', ...
    'Position', [100, 100, 650, 520], 'Color', 'w', 'Visible', 'off');
ax = axes('Position', [0.15, 0.15, 0.75, 0.75]);
hold on;

% Background patches: low (blue), med (grey), high (pink)
patch([1 numPosBins numPosBins 1], [0         0         lowThresh  lowThresh],  [0.85 0.95 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
patch([1 numPosBins numPosBins 1], [lowThresh  lowThresh  highThresh highThresh], [0.92 0.92 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
patch([1 numPosBins numPosBins 1], [highThresh highThresh maxY       maxY],       [0.98 0.85 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Speed traces
for l = 1:numLaps
    if all(isnan(allInterpSpeeds(l,:))), continue; end
    plot(ax, 1:numPosBins, allInterpSpeeds(l,:), 'Color', [0.2, 0.2, 0.2, 0.25], 'LineWidth', 1.5);
end

% Add horizontal dashed lines at thresholds for clarity
yline(ax, lowThresh,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'Alpha', 0.8);
yline(ax, highThresh, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'Alpha', 0.8);

%% Formatting
ylabel(ax, 'Speed (cm/s)', 'FontSize', 14);
xlabel(ax, 'Position (cm)',  'FontSize', 14);
box off;
set(ax, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 12);
xlim([1, numPosBins]);
xticks([1, 40, 80, 120, 160, 200]);
xticklabels({'1', '40', '80', '120', '160', '200'});
ylim([minY, maxY]);

% Y ticks match the fixed thresholds used in the analysis
cleanTicks = unique([1, lowThresh, highThresh, round(maxY)]);
yticks(cleanTicks);

if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax, 0); end
if exist('offsetAxes',            'file') == 2, offsetAxes(ax); end

set(figA, 'Visible', 'on');
end
% function figA = plotSpeedTrajectories(response, smoothSigma)
%     
%     % this function loads the lap data for a session and plots the traces across the position bins 
%     % The lap data is interpolated using the num of bins because they have different lengths 
%     % the colours patches shoould use the speed edges for low medium and
%     % high from the getLowMedHigh_SpeedPosActivity function; 
%     % Tick 1, 15, 30 and max can be labelled on the y-axis 
%     
%     if nargin < 2, smoothSigma = 2; end
%     
%     numLaps = length(response.lapRunningSpeed);
%     numPosBins = 200;
%     allInterpSpeeds = nan(numLaps, numPosBins);
%     
%     for l = 1:numLaps
%         rawV = response.lapRunningSpeed{l};
%         if length(rawV) < 5, continue; end
%         binnedV = interp1(linspace(1, numPosBins, length(rawV)), rawV, 1:numPosBins, 'linear', 'extrap');
%         allInterpSpeeds(l, :) = smoothdata(binnedV, 'gaussian', smoothSigma * 3);
%     end
%     
%     runningSpeeds = allInterpSpeeds(allInterpSpeeds > 1 & ~isnan(allInterpSpeeds));
%     sessionQuantiles = quantile(runningSpeeds, linspace(0, 1, 4));
%     lowThresh  = sessionQuantiles(2); 
%     highThresh = sessionQuantiles(3); 
%     
%     minY = 0;
%     maxY = max(allInterpSpeeds(:)) + 1; if isnan(maxY), maxY = 50; end
%     
%     figA = figure('Name', 'Quantile Speed Separation Zones', ...
%                   'Position', [100, 100, 650, 520], 'Color', 'w', 'Visible', 'off');
%     ax = axes('Position', [0.15, 0.15, 0.75, 0.75]);
%     hold on;
%     
%     patch([1 numPosBins numPosBins 1], [0 0 lowThresh lowThresh], ...
%           [0.85 0.95 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
%     patch([1 numPosBins numPosBins 1], [lowThresh lowThresh highThresh highThresh], ...
%           [0.92 0.92 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
%     patch([1 numPosBins numPosBins 1], [highThresh highThresh maxY maxY], ...
%           [0.98 0.85 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
%       
%     for l = 1:numLaps
%         if all(isnan(allInterpSpeeds(l,:))), continue; end
%         plot(ax, 1:numPosBins, allInterpSpeeds(l,:), 'Color', [0.2, 0.2, 0.2, 0.25], 'LineWidth', 1.5);
%     end
%     
%     ylabel(ax, 'Speed (cm/s)', 'FontSize', 14);
%     xlabel(ax, 'Position (cm)', 'FontSize', 14); 
%     box off; 
%     set(ax, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 12);
%     
%     xlim([1, numPosBins]); 
%     xticks([1, 40, 80, 120, 160, 200]);
%     xticklabels({'1', '40', '80', '120', '160', '200'});
%     
%     ylim([minY, maxY]);
%     
%     cleanTicks = [1, 15, 30, round(maxY)];
%     yticks(unique(cleanTicks));
%     
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax, 0); end 
%     if exist('offsetAxes', 'file') == 2, offsetAxes(ax); end
%     
%     set(figA, 'Visible', 'on');
% end

% function figA = plotSpeedTrajectories(response, smoothSigma)
%     % Replicates runBeh_example1.png format using QUANTILE split boundaries
%     % but visually color-codes and highlights absolute physical speed ranges.
%     
%     if nargin < 2, smoothSigma = 2; end
%     
%     numLaps = length(response.lapRunningSpeed);
%     numPosBins = 200;
%     allInterpSpeeds = nan(numLaps, numPosBins);
%     medianSpeedsPerLap = zeros(numLaps, 1);
%     
%     % Process and spatially interpolate every single lap
%     for l = 1:numLaps
%         rawV = response.lapRunningSpeed{l};
%         if length(rawV) < 5, continue; end
%         binnedV = interp1(linspace(1, numPosBins, length(rawV)), rawV, 1:numPosBins, 'linear', 'extrap');
%         allInterpSpeeds(l, :) = smoothdata(binnedV, 'gaussian', smoothSigma * 3);
%         
%         % Lap overall speed metric
%         medianSpeedsPerLap(l) = median(allInterpSpeeds(l, :), 'omitnan');
%     end
%     
%     % Filter out the laps where median speed is 0 before calculating quantiles
%     validLapSpeeds = medianSpeedsPerLap(medianSpeedsPerLap > 1 & ~isnan(medianSpeedsPerLap));
%     
%     % Split into 3 equal data groups (33.3% / 66.6%); Quantile based binning
%     lapEdges = quantile(validLapSpeeds, linspace(0, 1, 4));
%     
%     lowThresh  = lapEdges(2); % Upper bound of Low / Lower bound of Med
%     highThresh = lapEdges(3); % Upper bound of Med / Lower bound of High
%     
%     % -------------------------------------------------------------------------
%     % RIGID PHYSICAL BOUNDARIES FOR COLOR-CODING
%     % -------------------------------------------------------------------------
%     COLOR_LOW_LIMIT  = 15; % cm/s
%     COLOR_HIGH_LIMIT = 27; % cm/s
%     
%     minY = 0;
%     maxY = max(allInterpSpeeds(:)) + 1; if isnan(maxY), maxY = 50; end
%     
%     % Figure 
%     figA = figure('Name', 'Quantile Stratified Speed Trajectories', ...
%                   'Position', [100, 100, 650, 520], 'Color', 'w', 'Visible', 'off');
%     ax = axes('Position', [0.15, 0.15, 0.75, 0.75]);
%     hold on;
%     
%     % 4. Create Background Shading Bands using RIGID physical milestones
%     % Blue patch (Low Speed Range: 0 to 15 cm/s)
%     patch([1 numPosBins numPosBins 1], [0 0 COLOR_LOW_LIMIT COLOR_LOW_LIMIT], ...
%           [0.85 0.95 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
%     % Grey patch (Medium Speed Range: 15 to 27 cm/s)
%     patch([1 numPosBins numPosBins 1], [COLOR_LOW_LIMIT COLOR_LOW_LIMIT COLOR_HIGH_LIMIT COLOR_HIGH_LIMIT], ...
%           [0.92 0.92 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
%     % Pink patch (High Speed Range: 27 to maxY cm/s)
%     patch([1 numPosBins numPosBins 1], [COLOR_HIGH_LIMIT COLOR_HIGH_LIMIT maxY maxY], ...
%           [0.98 0.85 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
%       
%     % 5. Plot Individual Lap Lines
%     for l = 1:numLaps
%         if all(isnan(allInterpSpeeds(l,:))), continue; end
%         plot(ax, 1:numPosBins, allInterpSpeeds(l,:), 'Color', [0 0 0 0.3], 'LineWidth', 2);
%     end
%     
%     % 6. Assign Laps via Quantile Masks and Plot Means
%     lowIdx  = medianSpeedsPerLap >= 1  & medianSpeedsPerLap < lowThresh;
%     medIdx  = medianSpeedsPerLap >= lowThresh & medianSpeedsPerLap <= highThresh;
%     highIdx = medianSpeedsPerLap > highThresh;
%     
%     % Updated labels so they cleanly display your actual, current quantile split numbers
%     groupProps = { ...
%         {lowIdx,  [0, 0.6, 0.8],  sprintf('Quantile Low Mean (<%.1f cm/s)', lowThresh)}, ...
%         {medIdx,  [0.2, 0.2, 0.2], sprintf('Quantile Med Mean (%.1f-%.1f cm/s)', lowThresh, highThresh)}, ...
%         {highIdx, [0.8, 0, 0.6],   sprintf('Quantile High Mean (>%.1f cm/s)', highThresh)} ...
%     };
%     
% %     for g = 1:3
% %         currentMask = groupProps{g}{1};
% %         groupColor  = groupProps{g}{2};
% %         groupLabel  = groupProps{g}{3};
% %     
% %         if sum(currentMask) > 0
% %             stratifiedMean = mean(allInterpSpeeds(currentMask, :), 1, 'omitnan');
% %             plot(ax, 1:numPosBins, stratifiedMean, 'Color', groupColor, ...
% %                 'LineWidth', 4.0, 'DisplayName', groupLabel);
% %         end
% %     end
% %     
%     % 7. Formatting
%     ylabel(ax, 'Speed (cm/s)', 'FontSize', 14);
%     xlabel(ax, 'Position (cm)', 'FontSize', 14);
%     box off; set(ax, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 12);
%     
%     xlim([1, numPosBins]);
%     xticks([1, 40, 80, 120, 160, 200]);
%     xticklabels({'1', '40', '80', '120', '160', '200'});
%     
%     ylim([minY, maxY]);
%     
%     % Forces your y-axis ticks to match your physical color-coding lines perfectly
%     yticks([1, COLOR_LOW_LIMIT, 30, 45]); 
%     yticklabels({'1', '15', '30', '45'});
%     
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax, 0); end
%     if exist('offsetAxes', 'file') == 2, offsetAxes(ax); end
%     
% %     legend('show', 'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
% end

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


% function figA = plotSpeedTrajectories(response, smoothSigma)
% % Loads lap running speed data and plots traces across position bins.
% % Shows the per-position median threshold line from getLowHighSpeedPositionMatrix
% % instead of fixed horizontal patches — so the shading reflects the actual
% % occupancy-based split used in the analysis.
% 
% if nargin < 2, smoothSigma = 2; end
% 
% %% Interpolate and smooth lap speeds
% numLaps    = length(response.lapRunningSpeed);
% numPosBins = 200;
% allInterpSpeeds = nan(numLaps, numPosBins);
% 
% for l = 1:numLaps
%     rawV = response.lapRunningSpeed{l};
%     if length(rawV) < 5, continue; end
%     binnedV = interp1(linspace(1, numPosBins, length(rawV)), rawV, 1:numPosBins, 'linear', 'extrap');
%     allInterpSpeeds(l, :) = smoothdata(binnedV, 'gaussian', smoothSigma * 3);
% end
% 
% %% Get per-position median threshold from saved analysis
% if isfield(response, 'speedPositionActivity') && ...
%    isfield(response.speedPositionActivity, 'lowHigh') && ...
%    isfield(response.speedPositionActivity.lowHigh, 'medianThreshLine')
%     medianThreshLine = response.speedPositionActivity.lowHigh.medianThreshLine;
%     fprintf('Using saved per-position median threshold line\n');
% else
%     error('No medianThreshLine found. Run getLowHighSpeedPositionMatrix first.');
% end
% 
% %% Axis limits
% minY = 0;
% maxY = max(allInterpSpeeds(:)) + 1;
% if isnan(maxY), maxY = 50; end
% 
% %% Plot
% figA = figure('Name', 'Speed Trajectories', ...
%     'Position', [100, 100, 650, 520], 'Color', 'w', 'Visible', 'off');
% ax = axes('Position', [0.15, 0.15, 0.75, 0.75]);
% hold on;
% 
% % Speed traces
% for l = 1:numLaps
%     if all(isnan(allInterpSpeeds(l,:))), continue; end
%     plot(ax, 1:numPosBins, allInterpSpeeds(l,:), 'Color', [0.2, 0.2, 0.2, 0.25], 'LineWidth', 1.5);
% end
% 
% % Per-position median threshold line
% plot(ax, 1:numPosBins, medianThreshLine, '-', 'Color', [0.8 0.0 0.6], 'LineWidth', 2.5);
% 
% %% Formatting
% ylabel(ax, 'Speed (cm/s)', 'FontSize', 14);
% xlabel(ax, 'Position (cm)', 'FontSize', 14);
% box off;
% set(ax, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 12);
% xlim([1, numPosBins]);
% xticks([1, 40, 80, 120, 160, 200]);
% xticklabels({'1', '40', '80', '120', '160', '200'});
% ylim([minY, maxY]);
% 
% % Y ticks: min, session median summary, max
% sessionMedian = response.speedPositionActivity.lowHigh.sessionMedian;
% cleanTicks    = unique([1, round(sessionMedian), round(maxY)]);
% yticks(cleanTicks);
% 
% if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax, 0); end
% if exist('offsetAxes',            'file') == 2, offsetAxes(ax); end
% 
% set(figA, 'Visible', 'on');
% end


% function figA = plotSpeedTrajectories(response, smoothSigma)
% % Loads lap running speed data and plots traces across position bins.
% % Shows the per-position median threshold line from getLowHighSpeedPositionMatrix
% % instead of fixed horizontal patches — so the shading reflects the actual
% % occupancy-based split used in the analysis.
% 
% if nargin < 2, smoothSigma = 2; end
% 
% %% Interpolate and smooth lap speeds
% numLaps    = length(response.lapRunningSpeed);
% numPosBins = 200;
% allInterpSpeeds = nan(numLaps, numPosBins);
% 
% for l = 1:numLaps
%     rawV = response.lapRunningSpeed{l};
%     if length(rawV) < 5, continue; end
%     binnedV = interp1(linspace(1, numPosBins, length(rawV)), rawV, 1:numPosBins, 'linear', 'extrap');
%     allInterpSpeeds(l, :) = smoothdata(binnedV, 'gaussian', smoothSigma * 3);
% end
% 
% %% Get per-position median threshold from saved analysis
% if isfield(response, 'speedPositionActivity') && ...
%    isfield(response.speedPositionActivity, 'lowHigh') && ...
%    isfield(response.speedPositionActivity.lowHigh, 'medianThreshLine')
%     medianThreshLine = response.speedPositionActivity.lowHigh.medianThreshLine;
%     fprintf('Using saved per-position median threshold line\n');
% else
%     error('No medianThreshLine found. Run getLowHighSpeedPositionMatrix first.');
% end
% 
% %% Axis limits
% minY = 0;
% maxY = max(allInterpSpeeds(:)) + 1;
% if isnan(maxY), maxY = 50; end
% 
% %% Plot
% figA = figure('Name', 'Speed Trajectories', ...
%     'Position', [100, 100, 650, 520], 'Color', 'w', 'Visible', 'off');
% ax = axes('Position', [0.15, 0.15, 0.75, 0.75]);
% hold on;
% 
% % Speed traces
% for l = 1:numLaps
%     if all(isnan(allInterpSpeeds(l,:))), continue; end
%     plot(ax, 1:numPosBins, allInterpSpeeds(l,:), 'Color', [0.2, 0.2, 0.2, 0.25], 'LineWidth', 1.5);
% end
% 
% % Per-position median threshold line
% plot(ax, 1:numPosBins, medianThreshLine, '-', 'Color', [0.8 0.0 0.6], 'LineWidth', 2.5);
% 
% %% Formatting
% ylabel(ax, 'Speed (cm/s)', 'FontSize', 14);
% xlabel(ax, 'Position (cm)', 'FontSize', 14);
% box off;
% set(ax, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 12);
% xlim([1, numPosBins]);
% xticks([1, 40, 80, 120, 160, 200]);
% xticklabels({'1', '40', '80', '120', '160', '200'});
% ylim([minY, maxY]);
% 
% % Y ticks: min, session median summary, max
% sessionMedian = response.speedPositionActivity.lowHigh.sessionMedian;
% cleanTicks    = unique([1, round(sessionMedian), round(maxY)]);
% yticks(cleanTicks);
% 
% if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(ax, 0); end
% if exist('offsetAxes',            'file') == 2, offsetAxes(ax); end
% 
% set(figA, 'Visible', 'on');
% end
