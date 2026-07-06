function plotSpeedTuning_TemporNasalSpeedInc(sessionFileInfo, response, doSmooth)
% plotSpeedTuning: Custom order, negative connectivity, Red +64, and shared raster colorbar.
%
% PATCHED: alignedResponses is RAW (not baseline-subtracted) upstream.
% This version now baseline-subtracts each trial (using baseWin, defined
% below) before using it for EITHER the PSTH panel or the rasters, and
% derives the raster color scale (globalCLim) from that baseline-
% subtracted data -- not from raw absolute fluorescence pooled across
% baseline+response+all-speeds, which was swamping the color scale with
% baseline-offset variability rather than actual response magnitude.

if nargin < 3; doSmooth = true; end

%% --- Setup Directory and PDF Path ---
saveFolder = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(saveFolder, 'dir'), mkdir(saveFolder); end
pdfPath = fullfile(saveFolder, [sessionFileInfo.animal_name, '_' ...
    sessionFileInfo.session_name '_SpeedTuningWithTemporNasalTest.pdf']);
if exist(pdfPath, 'file'), delete(pdfPath); end 

%% --- Data Preparation ---
psthData = response.psthData; 
timeVec  = psthData(1).timeVector;
nROI     = size(psthData(1).alignedResponses, 1);

% Specific order requested [cite: 41-47, 133-139]
targetOrder = {'0', '1', '-16', '-32', '-64', '-128', '-256', '64'};

% Initialize sorted structure to avoid dissimilar structure errors
sortedPsthData = psthData; 
sortedPsthData(:) = []; 
finalLabels = {};
rawStimValues = arrayfun(@(x) string(x.stimValue), psthData);

for i = 1:length(targetOrder)
    matchIdx = find(rawStimValues == targetOrder{i}, 1);
    if ~isempty(matchIdx)
        sortedPsthData(end+1) = psthData(matchIdx);
        finalLabels{end+1} = targetOrder{i};
    end
end

nSpeeds = length(finalLabels);
baseWin = [-0.5 0];   % pre-stim baseline window, consistent with the DirTuning analysis
respWin = [0.5 5]; 
baseIdx = timeVec >= baseWin(1) & timeVec <= baseWin(2);
respIdx = timeVec >= respWin(1) & timeVec <= respWin(2);
tWin = 5; 

%% --- Figure Generation ---
hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1350 850]);

for iROI = 1:nROI
    clf(hFig);
    means = zeros(nSpeeds, 1);
    sems  = zeros(nSpeeds, 1);

    % Pre-compute baseline-subtracted trials for every speed ONCE per ROI
    % (used by both Panel A/B and Panel C, and by the shared CLim calc)
    subTrialsPerSpeed = cell(nSpeeds, 1);
    for s = 1:nSpeeds
        trials = squeeze(sortedPsthData(s).alignedResponses(iROI, :, :)); % [nTimepoints x nTrials]
        perTrialBaseline = mean(trials(baseIdx, :), 1, 'omitnan');        % [1 x nTrials]
        subTrialsPerSpeed{s} = trials - perTrialBaseline;                 % baseline-subtracted, same shape
    end

    % --- Panel A: PSTH ---
    subplot(2, 2, 1); hold on;
    cmap = parula(nSpeeds);
    pHandles = gobjects(nSpeeds, 1);
    
    for s = 1:nSpeeds
        trials = subTrialsPerSpeed{s}; % now baseline-subtracted
        mTrace = mean(trials, 2, 'omitnan');
        sTrace = std(trials, 0, 2, 'omitnan') ./ sqrt(size(trials, 2));
        if doSmooth, mTrace = smoothdata(mTrace, 'gaussian', tWin); end
        
        trialAverages = mean(trials(respIdx, :), 1, 'omitnan');
        means(s) = mean(trialAverages, 'omitnan');
        sems(s)  = std(trialAverages, 'omitnan') / sqrt(sum(~isnan(trialAverages)));
        
        % Force +64 to Red [cite: 35, 122, 214]
        if strcmp(finalLabels{s}, '64'), lCol = [1 0 0]; else, lCol = cmap(s,:); end 
        
        fill([timeVec(:); flipud(timeVec(:))], [mTrace(:)-sTrace(:); flipud(mTrace(:)+sTrace(:))], ...
            lCol, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        pHandles(s) = plot(timeVec, mTrace, 'Color', lCol, 'LineWidth', 1.5);
    end
    xline(0, 'k--', 'LineWidth', 1); ylabel('\DeltaF/F (baseline-subtracted)'); xlabel('Time (s)');
    lgd = legend(pHandles, finalLabels, 'Location', 'bestoutside');
    title(lgd, 'Speeds'); title('PSTH');

    % --- Panel B: Tuning Curve ---
    subplot(2, 2, 2); hold on;
    xPts = 1:nSpeeds;
    negIdx = find(startsWith(finalLabels, '-'));
    
    % 1. Plot connected line for all negatives (including -64) 
    if ~isempty(negIdx)
        plot(xPts(negIdx), means(negIdx), 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
        errorbar(xPts(negIdx), means(negIdx), sems(negIdx), 'ko', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off','LineWidth', 1.5);
    end
    
    % 2. Plot standalone points (0, 1, and +64) [cite: 27, 35, 41]
    standaloneIdx = find(ismember(finalLabels, {'0', '1', '64'}));
    for idx = standaloneIdx
        mCol = 'k'; if strcmp(finalLabels{idx}, '64'), mCol = 'r'; end
        errorbar(xPts(idx), means(idx), sems(idx), 'o', 'Color', mCol, ...
            'MarkerFaceColor', mCol, 'MarkerSize', 6, 'LineStyle', 'none', 'LineWidth', 1.5);
    end
    
    set(gca, 'XTick', xPts, 'XTickLabel', finalLabels);
    ylabel('Avg \DeltaF/F (baseline-subtracted)'); title('Speed Tuning Curve');
    xlabel('Visual Speeds (deg/s)');

    % --- Panel C: Trial Rasters ---
    startX = 0.05; 
    widthPerRaster = 0.85 / nSpeeds; 
    rasterHeight = 0.25;
    bottomY = 0.1;

    % Global color limits, now computed from BASELINE-SUBTRACTED data
    % (previously this was raw absolute fluorescence pooled across
    % baseline+response+all-speeds, which swamped the scale with
    % baseline-offset variability rather than actual response magnitude).
    allDataAcrossSpeeds = [];
    for s = 1:nSpeeds
        allDataAcrossSpeeds = [allDataAcrossSpeeds; subTrialsPerSpeed{s}(:)];
    end
    globalCLim = [quantile(allDataAcrossSpeeds, 0.05), quantile(allDataAcrossSpeeds, 0.98) + 1e-6];

    for s = 1:nSpeeds
        axPos = [startX + (s-1)*widthPerRaster, bottomY, widthPerRaster*0.85, rasterHeight];
        axR = axes('Position', axPos);
        
        currData = subTrialsPerSpeed{s}'; % [nTrials x nTimepoints], baseline-subtracted
        
        if ~isempty(currData)
            imagesc(timeVec, 1:size(currData,1), currData);
            colormap(axR, redblue_local()); hold on; xline(0, 'k:', 'LineWidth', 1.2);
            
            if globalCLim(1) < globalCLim(2)
                set(axR, 'CLim', globalCLim); 
            end
        end
        
        set(axR, 'YDir', 'reverse', 'FontSize', 7, 'XLim', [-1 6]);
        tCol = 'k'; if strcmp(finalLabels{s}, '64'), tCol = 'r'; end % Red title for +64 [cite: 55, 147]
        title(sprintf('Sp: %s', finalLabels{s}), 'Color', tCol, 'FontSize', 8);
        
        if s == 1; ylabel('Trials'); else set(axR, 'YTickLabel', []); end
        xlabel('Time (s)');
        % Shared Colorbar for all rasters
        if s == nSpeeds
            cbPos = [startX + s*widthPerRaster, bottomY, 0.012, rasterHeight];
            cb = colorbar(axR, 'Position', cbPos);
            ylabel(cb, '\DeltaF/F (baseline-subtracted)', 'FontSize', 7);
        end
    end
    
    sgtitle(sprintf('%s | %s | ROI %d', sessionFileInfo.animal_name, sessionFileInfo.session_name, iROI), ...
        'Interpreter', 'none', 'FontWeight', 'bold'); 
    exportgraphics(hFig, pdfPath, 'Append', true);
end
close(hFig);
fprintf('Finished processing %d ROIs. PDF: %s\n', nROI, pdfPath);
end

function cmap = redblue_local()
% Simple diverging blue-white-red colormap, since MATLAB has no built-in
% 'redblue'. Blue = below baseline, white = at baseline, red = above.
n = 128;
r = [linspace(0, 1, n), ones(1, n)];
g = [linspace(0, 1, n), linspace(1, 0, n)];
b = [ones(1, n), linspace(1, 0, n)];
cmap = [r(:), g(:), b(:)];
end