function offsetAxes(ax)
% thanks to Pierre Morel, undocumented Matlab
% and https://stackoverflow.com/questions/38255048/separating-axes-from-plot-area-in-matlab
% this only seems to work for older matlab versions 
%
% by Anne Urai, 2016

set(ax, 'TickDir','out', 'TickLength', [0.01 0.001]...
    ,  'color', 'none', 'box','off', 'XColor', 'k', 'YColor', 'k',...
    'FontName', 'Calibri', 'LineWidth', 0.5)

if ~exist('ax', 'var'), ax = gca; end
% modify the x and y limits to below the data (by a small amount)
ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/4;
ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/4;
% this will keep the changes constant even when resizing axes
addlistener (ax, 'MarkedClean', @(obj,event)resetVertex(ax));
end

function resetVertex(ax)

set(ax, 'TickDir','out', 'TickLength', [0.01 0.001]...
    ,  'color', 'none', 'box','off', 'XColor', 'k', 'YColor', 'k',...
    'FontName', 'Calibri', 'LineWidth', 0.5)

% extract the x axis vertext data
% X, Y and Z row of the start and end of the individual axle.
ax.XRuler.Axle.VertexData(1,1) = min(get(ax, 'Xtick'));
% repeat for Y (set 2nd row)
ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
end


% function offsetAxes(ax)
%     % Designed for MATLAB 2025a/2026
%     if nargin < 1 || isempty(ax), ax = gca; end
% 
%     % 1. Force a draw to finalize tick positions
%     drawnow; 
%     hold(ax, 'on');
%     box(ax, 'off');
% 
%     % 2. Get the current Ticks
%     xTicks = ax.XTick;
%     yTicks = ax.YTick;
%     if isempty(xTicks) || isempty(yTicks), return; end
% 
%     % 3. THE FIX: Make the Spine invisible without losing Ticks
%     % We set the 'Interactions' to none so the axis doesn't "repair" itself
%     ax.Interactions = []; 
% 
%     % Set the Spine (Axle) to transparent
%     % In 2026, setting the color to [0 0 0 0] (RGBA) hides the line 
%     % while keeping the tick marks intact.
%     ax.XAxis.Color = 'none'; 
%     ax.YAxis.Color = 'none';
% 
%     % 4. Adjust Limits to create the physical gap
%     % This moves the labels away from the corner
%     xGap = (xTicks(2) - xTicks(1)) / 4;
%     yGap = (yTicks(2) - yTicks(1)) / 4;
%     ax.XLim(1) = xTicks(1) - xGap;
%     ax.YLim(1) = yTicks(1) - yGap;
% 
%     % 5. Draw the "Manual" Spines
%     % These lines start at the first tick and end at the last tick.
%     % Because they are separate, they cannot touch.
%     line(ax, [min(xTicks), max(xTicks)], [ax.YLim(1), ax.YLim(1)], ...
%         'Color', 'k', 'LineWidth', ax.LineWidth, 'Clipping', 'off', 'HandleVisibility', 'off');
% 
%     line(ax, [ax.XLim(1), ax.XLim(1)], [min(yTicks), max(yTicks)], ...
%         'Color', 'k', 'LineWidth', ax.LineWidth, 'Clipping', 'off', 'HandleVisibility', 'off');
% 
%     % 6. Restore Text Visibility
%     % Setting XColor to 'none' often hides labels; we force them back to black.
%     ax.XAxis.Label.Color = [0 0 0];
%     ax.YAxis.Label.Color = [0 0 0];
%     ax.XAxis.TickLabelColor = [0 0 0];
%     ax.YAxis.TickLabelColor = [0 0 0];
% 
%     % Re-force tick marks to appear outward
%     ax.XAxis.TickDirection = 'out';
%     ax.YAxis.TickDirection = 'out';
% end

% function offsetAxes(ax)
%     if nargin < 1 || isempty(ax), ax = gca; end
% 
%     % 1. Finalize and check scales
%     drawnow; 
%     hold(ax, 'on');
%     box(ax, 'off');
%     
%     xTicks = ax.XTick;
%     yTicks = ax.YTick;
%     if isempty(xTicks) || isempty(yTicks), return; end
% 
%     % Hide the default connecting lines
%     ax.XAxis.Color = 'none'; 
%     ax.YAxis.Color = 'none';
%     ax.XAxis.Label.Color = [0 0 0];
%     ax.YAxis.Label.Color = [0 0 0];
%     ax.XAxis.TickLabelColor = [0 0 0];
%     ax.YAxis.TickLabelColor = [0 0 0];
% 
%     % 2. Calculate the Gap (shifter)
%     xGap = diff(ax.XLim) * 0.02;
%     yGap = diff(ax.YLim) * 0.02;
%     
%     ax.XLim(1) = xTicks(1) - xGap;
%     ax.YLim(1) = yTicks(1) - yGap;
% 
%     % 3. Draw the "Manual" Spines
%     line(ax, [min(xTicks), max(xTicks)], [ax.YLim(1), ax.YLim(1)], ...
%         'Color', 'k', 'LineWidth', ax.LineWidth, 'Clipping', 'off', 'HandleVisibility', 'off');
%     line(ax, [ax.XLim(1), ax.XLim(1)], [min(yTicks), max(yTicks)], ...
%         'Color', 'k', 'LineWidth', ax.LineWidth, 'Clipping', 'off', 'HandleVisibility', 'off');
% 
%     % 4. Draw Manual Ticks
%     tickLen = 0.02 * diff(ax.YLim); % Slightly longer ticks for a "Nature" look
%     for i = 1:length(xTicks)
%         line(ax, [xTicks(i), xTicks(i)], [ax.YLim(1), ax.YLim(1)-tickLen], ...
%             'Color', 'k', 'LineWidth', ax.LineWidth, 'Clipping', 'off', 'HandleVisibility', 'off');
%     end
%     
%     tickLenX = 0.02 * diff(ax.XLim);
%     for i = 1:length(yTicks)
%         line(ax, [ax.XLim(1), ax.XLim(1)-tickLenX], [yTicks(i), yTicks(i)], ...
%             'Color', 'k', 'LineWidth', ax.LineWidth, 'Clipping', 'off', 'HandleVisibility', 'off');
%     end
% 
%     % 5. ADJUST SPACING (The "More Space" Fix)
%     % Increase the gap between the tick marks and the numbers (0.5, 1, 10, etc.)
%     ax.XAxis.TickLabelGapOffset = 6; % Doubled from previous
%     ax.YAxis.TickLabelGapOffset = 6;
%     
%     % Move the main Axis Labels (e.g., "Speed (cm/s)") further out
%     % We multiply the tick length to ensure they don't collide with the numbers
%     ax.XAxis.Label.VerticalAlignment = 'top';
%     ax.XAxis.Label.Position(2) = ax.YLim(1) - (tickLen * 4); % Increased multiplier
%     
%     ax.YAxis.Label.HorizontalAlignment = 'right';
%     ax.YAxis.Label.Position(1) = ax.XLim(1) - (tickLenX * 4); % Increased multiplier
% end