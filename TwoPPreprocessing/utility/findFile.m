function filePath = findFile(fileArray, keyword)
    % findFile Searches for a file path containing the specified keyword.
    %
    % fileArray - A cell array of file paths.
    % keyword   - A string or pattern to search for in the file paths
    % filePath  - The file path that matches the keyword. Returns empty if not found.

    % Initialize output
    filePath = '';
    
    % Lowercase for case-insensitive search
    keyword = lower(keyword);
    for i = 1:length(fileArray)
        currentFilePath = fileArray{i};
        if contains(lower(currentFilePath), keyword)
            filePath = currentFilePath;
            return; % 
        end
    end
    
    if isempty(filePath)
        warning('No file matching the keyword "%s" was found.', keyword);
    end
end
