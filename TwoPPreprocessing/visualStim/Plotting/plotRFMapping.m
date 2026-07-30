function plotRFMapping(sessionFileInfo, RFMapping, RFMappingMetadata, doSmooth, onlyTuned)
% plotRFMapping: Summary (Left) + 16 Stim Rasters + 1 Blank Raster (Right)
% Overlays mean blank PSTH (white) and includes a colorbar.

if nargin < 4; doSmooth = false; end
if nargin < 5; onlyTuned = true; end 


rasterTimeWin = [-1 2.5]; 


saveFolder = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~isfolder(saveFolder); mkdir(saveFolder); end
pdfPath = fullfile(saveFolder, [sessionFileInfo.animal_name, '_' sessionFileInfo.session_name '_RFMapping.pdf']);
if exist(pdfPath, 'file'), delete(pdfPath); end

uAz = RFMappingMetadata.uAz;
uEl_plot = RFMappingMetadata.uEl;
timeVector = RFMappingMetadata.timeVector;
nROI = numel(RFMapping);
nAz = length(uAz); nEl = length(uEl_plot);
sigma_smooth = 0.8;
if doSmooth; tWin = 5; else; tWin = 1; end

dAz = abs(uAz(2) - uAz(1)); dEl = abs(uEl_plot(1) - uEl_plot(2));
azLimFull = [min(uAz) - dAz/2, max(uAz) + dAz/2];
elLimFull = [min(uEl_plot) - dEl/2, max(uEl_plot) + dEl/2];

for iROI = 1:nROI
    isTuned = RFMapping(iROI).isResponsive;
    if onlyTuned && ~isTuned; continue; end

    fig = figure('Color', 'w', 'Position', [50 50 1500 800], 'Visible', 'off');

    %% --- LEFT PANEL: Heatmap + Red Stim PSTH + White Blank PSTH ---
    ax1 = axes('Position', [0.05 0.15 0.35 0.7]); 
    imagesc(uAz, uEl_plot, imgaussfilt(RFMapping(iROI).meanGridResponse, sigma_smooth)); 
    hold on; colormap(ax1, gray);
    set(ax1, 'YDir', 'normal', 'CLim', [0, max(RFMapping(iROI).meanGridResponse(:)) + 1e-6]);

    if isTuned
        rectangle('Position', [RFMapping(iROI).centerAz-dAz/2, RFMapping(iROI).centerEl-dEl/2, dAz, dEl], ...
                  'EdgeColor', 'y', 'LineWidth', 2.5);
    end

    % Overlay PSTHs
    vS_base = dEl * 0.4; hS = dAz * 0.85;     
    for r = 1:nEl
        for c = 1:nAz
            % Blank PSTH (White)
            trB = RFMapping(iROI).meanBlankResponse;
            % Stim PSTH (Red)
            trS = RFMapping(iROI).meanTemporalResponse(:, r, c);

            if doSmooth
                trB = smoothdata(trB, 'gaussian', tWin);
                trS = smoothdata(trS, 'gaussian', tWin);
            end

            tNorm = (timeVector - min(timeVector)) / (max(timeVector) - min(timeVector));
            sX = (tNorm - 0.5) * hS + uAz(c);
            sY_B = (trB / RFMapping(iROI).peakAmplitude * vS_base) + uEl_plot(r);
            sY_S = (trS / RFMapping(iROI).peakAmplitude * vS_base) + uEl_plot(r);

            plot(ax1, sX, sY_B, 'w', 'LineWidth', 0.8, 'LineStyle', ':'); % Blank in White
            plot(ax1, sX, sY_S, 'r', 'LineWidth', 1.2);                   % Stim in Red
        end
    end
    title(['ROI ' num2str(iROI)], 'FontSize', 14);

    %% --- RIGHT PANEL: 4x4 Stim Rasters + 1 Blank Raster ---
    rStartX = 0.45; rStartY = 0.2;
    rW = 0.08; rH = 0.12; 

    % Main 4x4 Grid
    for r = 1:nEl
        for c = 1:nAz
            axR = axes('Position', [rStartX + (c-1)*(rW+0.02), rStartY + (nEl-r)*(rH+0.05), rW, rH]);
            currData = RFMapping(iROI).baselineSubtracted{r,c};

            if ~isempty(currData)
                tIdx = timeVector(1:size(currData,2)) >= rasterTimeWin(1) & timeVector(1:size(currData,2)) <= rasterTimeWin(2);
                imagesc(timeVector(tIdx), 1:size(currData,1), currData(:,tIdx));
                colormap(axR, parula); hold on;
                xline(0, 'w:', 'LineWidth', 1.2);
                % Consistent color scaling across the grid
                set(axR, 'CLim', [quantile(currData(:), 0.05), quantile(currData(:), 0.98) + 1e-6]);
            end
            set(axR, 'YDir', 'reverse', 'FontSize', 7, 'XLim', rasterTimeWin);
            if r < nEl; set(axR, 'XTickLabel', []); end
            if c > 1; set(axR, 'YTickLabel', []); end
            title(sprintf('%d, %d', uAz(c), uEl_plot(r)), 'FontSize', 7);
        end
    end

    % --- ADDITIONAL: Blank Trial Raster ---
    % Positioned below the main grid
    axBlank = axes('Position', [rStartX, 0.05, rW, rH]);
    blankTrials = RFMapping(iROI).baselineSubtractedBlank;
    if ~isempty(blankTrials)
        tIdxB = timeVector(1:size(blankTrials,2)) >= rasterTimeWin(1) & timeVector(1:size(blankTrials,2)) <= rasterTimeWin(2);
        imagesc(timeVector(tIdxB), 1:size(blankTrials,1), blankTrials(:,tIdxB));
        colormap(axBlank, parula); hold on;
        set(axBlank, 'CLim', [quantile(blankTrials(:), 0.05), quantile(blankTrials(:), 0.98) + 1e-6]);
    end
    set(axBlank, 'YDir', 'reverse', 'FontSize', 7, 'XLim', rasterTimeWin);
    title('BLANK TRIALS', 'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.4 0.4 0.4]);

    % --- COLORBAR ---
    cb = colorbar(axR, 'Position', [rStartX + (nAz)*(rW+0.02), rStartY, 0.015, rH*nEl]);
    cb.Label.String = 'dF/F';

    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType', 'vector'); 
    close(fig);
end
fprintf('PDF saved to: %s\n', pdfPath);
end


