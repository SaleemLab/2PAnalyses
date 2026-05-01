function plotSpeedTuning(sessionFileInfo, response, doSmooth)
% plotSpeedTuning: PSTHs, Disconnected Tuning Curves, and Trial Rasters.

if nargin < 3; doSmooth = true; end

%% --- Setup Directory and PDF Path ---
saveFolder = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(saveFolder, 'dir'), mkdir(saveFolder); end
pdfPath = fullfile(saveFolder, [sessionFileInfo.animal_name, '_' ...
    sessionFileInfo.session_name '_SpeedTuning.pdf']);

if exist(pdfPath, 'file'), delete(pdfPath); end 

%% --- Data Preparation ---
psthData  = response.psthData; 
timeVec   = psthData(1).timeVector;
nROI      = size(psthData(1).alignedResponses, 1);

% Clean Speed Names (Remove '-' signs and convert to numeric)
rawNames = {psthData.stimValue};
cleanSpeeds = zeros(1, length(rawNames));
for i = 1:length(rawNames)
    % Remove '-' and convert to number
    valStr = strrep(string(rawNames{i}), '-', '');
    cleanSpeeds(i) = str2double(valStr);
end

% Sort logic: 0 first, then 1, then others descending
is0 = (cleanSpeeds == 0);
is1 = (cleanSpeeds == 1);
othersIdx = find(~is0 & ~is1);

[~, sortOrder] = sort(cleanSpeeds(othersIdx), 'ascend');
sortedIndices = [find(is0), find(is1), othersIdx(sortOrder)];

uSpeedsSorted = cleanSpeeds(sortedIndices);
psthDataSorted = psthData(sortedIndices);

% Analysis Window
respWin = [0.5 2]; 
respIdx = timeVec >= respWin(1) & timeVec <= respWin(2);
tWin = 5; 

%% --- Figure Generation ---
hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 800]);

for iROI = 1:nROI
    clf(hFig);
    means = zeros(length(uSpeedsSorted), 1);
    sems  = zeros(length(uSpeedsSorted), 1);
    
    % --- Panel A: PSTH (Top Left) ---
    subplot(2, 4, [1, 2]); hold on;
    cmap = parula(length(uSpeedsSorted));
    pHandles = gobjects(length(uSpeedsSorted), 1);
    
    for s = 1:length(uSpeedsSorted)
        trials = squeeze(psthDataSorted(s).alignedResponses(iROI, :, :)); 
        meanTrace = mean(trials, 2, 'omitnan');
        semTrace  = std(trials, 0, 2, 'omitnan') ./ sqrt(size(trials, 2));
        
        if doSmooth
            meanTrace = smoothdata(meanTrace, 'gaussian', tWin);
        end
        
        trialAverages = mean(trials(respIdx, :), 1, 'omitnan');
        means(s) = mean(trialAverages);
        sems(s)  = std(trialAverages) / sqrt(length(trialAverages));
        
        fill([timeVec(:); flipud(timeVec(:))], [meanTrace(:)-semTrace(:); flipud(meanTrace(:)+semTrace(:))], ...
            cmap(s,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        pHandles(s) = plot(timeVec, meanTrace, 'Color', cmap(s,:), 'LineWidth', 2);
    end
    
    xline(0, 'k--', 'LineWidth', 1); 
    xlabel('Time (s)'); ylabel('\DeltaF/F');
    
    % FIXED LEGEND CALL
    lgd = legend(pHandles, string(uSpeedsSorted), 'Location', 'bestoutside');
    title(lgd, 'Speed'); 
    title('PSTH');

    % --- Panel B: Disconnected Tuning Curve (Top Right) ---
    subplot(2, 4, [3, 4]); hold on;
    
    idx0 = find(uSpeedsSorted == 0);
    idx1 = find(uSpeedsSorted == 1);
    idxRest = find(uSpeedsSorted ~= 0 & uSpeedsSorted ~= 1);
    
    % Line for the rest
    if ~isempty(idxRest)
        errorbar(idxRest, means(idxRest), sems(idxRest), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    end
    % Points for 0 and 1
    if ~isempty(idx0)
        errorbar(idx0, means(idx0), sems(idx0), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
    end
    if ~isempty(idx1)
        errorbar(idx1, means(idx1), sems(idx1), 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
    end
    
    set(gca, 'XTick', 1:length(uSpeedsSorted), 'XTickLabel', string(uSpeedsSorted));
    xlabel('Visual Speeds (deg/s)'); ylabel('Avg \DeltaF/F');
    title('Speed Tuning Curve');

    rasterTimeWin = [-1, 4]; 
    tIdx = timeVec >= rasterTimeWin(1) & timeVec <= rasterTimeWin(2);
    nSpeeds = length(uSpeedsSorted);
    
    for s = 1:nSpeeds
        % Create subplot in the bottom row (row 2 of a 2xN grid)
        ax = subplot(2, nSpeeds, nSpeeds + s);
        
        currData = squeeze(psthDataSorted(s).alignedResponses(iROI, :, :));
        if ismatrix(currData)
            % Transpose so trials are on Y-axis
            currData = currData'; 
            
            imagesc(timeVec(tIdx), 1:size(currData,1), currData(:, tIdx));
            colormap(ax, parula); hold on;
            xline(0, 'w:', 'LineWidth', 1.2);
            
            % Robust scaling
            cLim = [quantile(currData(:), 0.05), quantile(currData(:), 0.98) + 1e-6];
            if cLim(1) < cLim(2), set(ax, 'CLim', cLim); end
        end
        
        set(ax, 'YDir', 'reverse', 'FontSize', 8, 'XLim', rasterTimeWin);
        title(sprintf('Sp: %g', uSpeedsSorted(s)));
        
        if s > 1, set(ax, 'YTickLabel', []); end
        if s == 1, ylabel('Trials'); end
        xlabel('Time (s)');
    end

    sgtitle(sprintf('%s | %s | ROI %d', sessionFileInfo.animal_name, ...
        sessionFileInfo.session_name, iROI), 'Interpreter', 'none', 'FontWeight', 'bold');
    
    exportgraphics(hFig, pdfPath, 'Append', true);
end

close(hFig);
fprintf('PDF Export Finished: %s\n', pdfPath);
end