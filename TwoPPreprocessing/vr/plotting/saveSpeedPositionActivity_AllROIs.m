function saveSpeedPositionActivity_AllROIs(sessionFileInfo, response, options)
% saveallroistopdf - loops through all rois and saves them into a single pdf.
% every roi gets its own page.

    % 1. setup defaults and directory
    if nargin < 3, options = struct(); end
    if ~isfield(options, 'applySmoothing'), options.applySmoothing = false; end
    if ~isfield(options, 'smoothSigma'),    options.smoothSigma = [1.1, 1.5]; end

    figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'SpeedPositionSummaries');
    if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end

    % 2. define the master pdf filename
    pdfName = sprintf('%s_%s_AllROIs_SpeedPosition.pdf', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name);
    fullPDFPath = fullfile(figSaveDir, pdfName);

    % 3. delete existing pdf if it exists (so you don't keep appending to old files)
    if exist(fullPDFPath, 'file'), delete(fullPDFPath); end

    % 4. get total number of rois from the response matrix
    numROIs = size(response.speedPositionActivity.matrix, 3);
    fprintf('starting pdf generation for %d rois...\n', numROIs);

    % 5. the loop
    for targetROI = 1:numROIs
        % call your plotting function
        % note: we use 'Visible', 'off' inside a loop to speed things up
        % if your function creates the figure, we capture it here
        plotSpeedPositionActivity_ForROI(sessionFileInfo, response, targetROI, ...
            options.applySmoothing, options.smoothSigma);
        
        % find the figure that was just created
        figHandle = gcf; 
        
        % append to pdf
        % the first time this runs, it creates the file; after that, it appends pages.
        exportgraphics(figHandle, fullPDFPath, 'ContentType', 'vector', 'Append', true);
        
        % close figure immediately to save ram
        close(figHandle);
        
        % console progress update
        if mod(targetROI, 10) == 0 || targetROI == numROIs
            fprintf('processed %d/%d rois...\n', targetROI, numROIs);
        end
    end

    fprintf('success! all rois saved to: %s\n', fullPDFPath);
end