function standardizeGridAxes(axHandles)
% standardizeGridAxes Ensures all subplots in a grid have the same width.
% Useful when some columns lack colorbars or y-labels.
    drawnow;
    validAx = axHandles(isgraphics(axHandles));
    if isempty(validAx), return; end
    
    posAll = get(validAx, 'Position');
    if iscell(posAll), posMat = vertcat(posAll{:}); else, posMat = posAll; end
    
    % Find minimum width among all valid axes
    minW = min(posMat(:,3));
    
    for i = 1:numel(validAx)
        p = get(validAx(i), 'Position');
        p(3) = minW; 
        set(validAx(i), 'Position', p);
    end
end