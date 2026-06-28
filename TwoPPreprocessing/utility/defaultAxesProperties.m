function defaultAxesProperties(ax, offsetFlag)
set(ax, 'TickDir','out', 'Layer', 'top', 'TickLength', [0.02 0.001]...
    ,  'color', 'none', 'box','off', 'XColor', 'k', 'YColor', 'k',...
    'FontName', 'Arial', 'LineWidth', 1, 'FontSize', 12)
if nargin>1
    if offsetFlag
        offsetAxes(ax);
    end
end

end


% function defaultAxesProperties(ax, offsetFlag)
%     % Standard styling
%     set(ax, 'TickDir','out', 'TickLength', [0.01 0.001], ...
%             'Color', 'none', 'Box','off', 'XColor', 'k', 'YColor', 'k', ...
%             'FontName', 'Calibri', 'LineWidth', 0.5);
% 
%     if nargin > 1 && offsetFlag
%         offsetAxes(ax);
%     end
% end

% function defaultAxesProperties(ax, offsetFlag)
%     % 
%     if ~ishandle(ax)
%         ax = gca; 
%     end
% 
%     % 2. Apply your core styling
%     set(ax, 'TickDir','out', 'TickLength', [0.01 0.001], ...
%             'Color', 'none', 'Box','off', 'XColor', 'k', 'YColor', 'k', ...
%             'FontName', 'Ariel', 'FontSize', 12, 'LineWidth', 0.5);
%             
%     % 3. Set specific Label sizes
%     ax.XAxis.Label.FontSize = 11;
%     ax.YAxis.Label.FontSize = 11;
% 
%     % 4. Run the offset logic
%     if nargin > 1 && offsetFlag
%         offsetAxes(ax);
%     end
% end