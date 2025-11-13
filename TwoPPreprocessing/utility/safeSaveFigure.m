function safeSaveFigure(figH, fullPath)
    if isempty(fullPath), return; end
    try
        folder = fileparts(fullPath);
        if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
        saveas(figH, fullPath);
        fprintf('Figure saved to: %s\n', fullPath);
    catch ME
        warning('Failed to save figure to %s: %s', fullPath, ME.message);
    end
end