function plotDotFieldRFs(sessionFileInfo, doSmooth)
% plotDotFieldRFs: Automatically loads TL/TR/BL/BR responses and plots a 
% compact 2x2 comparison grid with shared scaling (0-5s stimulus window).

if nargin < 2; doSmooth = true; end

%% --- 1. Automated Data Loading ---
locNames = {'TL', 'TR', 'BL', 'BR'};
dataStructs = cell(1, 4); % To store response.psthData for each location

fprintf('Finding and loading quadrant response files...\n');
for iLoc = 1:length(locNames)
    % Search stimFiles for names containing the specific quadrant tag
    % Format assumes: ..._RFMapping_TL_... 
    matchIdx = find(contains({sessionFileInfo.stimFiles.name}, ['_' locNames{iLoc} '_']), 1);
    
    if ~isempty(matchIdx)
        resPath = sessionFileInfo.stimFiles(matchIdx).Response;
        fprintf('  Loading %s: %s\n', locNames{iLoc}, resPath);
        temp = load(resPath);
        % Standardize the field access 
        dataStructs{iLoc} = temp.response; 
    else
        error('Could not find response file for location: %s', locNames{iLoc});
    end
end

%% --- 2. Setup Directory and PDF Path ---
saveFolder = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(saveFolder, 'dir'), mkdir(saveFolder); end
pdfPath = fullfile(saveFolder, [sessionFileInfo.animal_name, '_' ...
    sessionFileInfo.session_name '_DotFieldRFs.pdf']);
if exist(pdfPath, 'file'), delete(pdfPath); end 

%% --- 3. Parameter Setup ---
tWin = 5; 
nROI = size(dataStructs{1}.psthData(1).alignedResponses, 1);
timeVec = dataStructs{1}.alignedTimes;

% Requested comparison window (Stimulus period)
plotWin = [0 5]; 
plotIdx = timeVec >= plotWin(1) & timeVec <= plotWin(2);

%% --- 4. Figure Generation ---
% Compact figure size for better PDF presentation
hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 850 700]);

for iROI = 1:nROI
    clf(hFig);
    
    % A. Calculate Global Y-Limits based ONLY on the 0-5s window traces
    % This ensures consistent amplitude comparison across positions 
    allMeans = [];
    for iLoc = 1:4
        pData = dataStructs{iLoc}.psthData;
        for s = 1:length(pData)
            % Use only the plotting window to determine scale
            trials = squeeze(pData(s).alignedResponses(iROI, plotIdx, :));
            allMeans = [allMeans; mean(trials, 2, 'omitnan')];
        end
    end
    
    % Buffer scaling: 15% head room, ensure 0 is visible for context
    yMax = max(allMeans) * 1.15; 
    yMin = min(allMeans) * 1.15;
    if yMin > 0, yMin = -0.05; end 
    yRange = [yMin, yMax];
    
    % B. Plot quadrants
    for iLoc = 1:4
        % Precise positions for touching panels [cite: 1252-1255]
        if iLoc == 1, pos = [0.12, 0.52, 0.38, 0.38]; % TL
        elseif iLoc == 2, pos = [0.50, 0.52, 0.38, 0.38]; % TR
        elseif iLoc == 3, pos = [0.12, 0.10, 0.38, 0.38]; % BL
        elseif iLoc == 4, pos = [0.50, 0.10, 0.38, 0.38]; % BR
        end
        
        ax = axes('Position', pos); hold on;
        currPsth = dataStructs{iLoc}.psthData;
        
        % Sort speeds numerically for consistent colormapping 
        rawNums = arrayfun(@(x) str2double(string(x.stimValue)), currPsth);
        [~, sIdx] = sort(rawNums);
        sPSTH = currPsth(sIdx);
        
        nS = length(sPSTH);
        cmap = parula(nS);
        pHTemp = gobjects(nS, 1);
        
        for s = 1:nS
            trials = squeeze(sPSTH(s).alignedResponses(iROI, :, :));
            mT = mean(trials, 2, 'omitnan');
            sT = std(trials, 0, 2, 'omitnan') ./ sqrt(size(trials, 2));
            if doSmooth, mT = smoothdata(mT, 'gaussian', tWin); end
            
            % Force Red for speed +64 [cite: 638-659]
            lCol = cmap(s,:);
            if strcmp(string(sPSTH(s).stimValue), "64"), lCol = [1 0 0]; end
            
            fill([timeVec(:); flipud(timeVec(:))], [mT(:)-sT(:); flipud(mT(:)+sT(:))], ...
                lCol, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            pHTemp(s) = plot(timeVec, mT, 'Color', lCol, 'LineWidth', 2.5);
        end
        
        % Formatting
        ylim(yRange); 
        xlim(plotWin); 
        set(ax, 'TickDir', 'out', 'Box', 'on', 'FontSize', 8);
        
        % Location label in top-left of quadrant 
        text(0.05, 0.90, locNames{iLoc}, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 10);
        
        % Only show outer labels to keep inner touching edges clean 
        if iLoc == 1 || iLoc == 3, ylabel('\DeltaF/F'); else set(ax, 'YTickLabel', []); end
        if iLoc == 3 || iLoc == 4, xlabel('Time (s)'); else set(ax, 'XTickLabel', []); end
        
        % 
        if iLoc == 2
            lgd = legend(pHTemp, arrayfun(@(x) string(x.stimValue), sPSTH), ...
                'Position', [0.89, 0.70, 0.08, 0.1], 'FontSize', 7);
            title(lgd, 'Speed');
        end
    end
    
    sgtitle(sprintf('%s | ROI %d | Dot field Rfs', sessionFileInfo.session_name, iROI), ...
        'FontSize', 11, 'Interpreter', 'none', 'FontWeight', 'bold');
    
    % Compact export [cite: 1364]
    exportgraphics(hFig, pdfPath, 'Append', true, 'Resolution', 150);
end

close(hFig);
fprintf('RF Comparison Complete. PDF saved: %s\n', pdfPath);
end