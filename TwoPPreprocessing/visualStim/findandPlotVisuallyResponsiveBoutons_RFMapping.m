function findandPlotVisuallyResponsiveBoutons_RFMapping(sessionFileInfo)
% findandPlotVisuallyResponsiveBoutons_RFMapping: Summary reporter pulling 
% directly from the pre-calculated RFMapping structure.

% 1. Load the pre-calculated RFMapping data
savePath = sessionFileInfo.otherSessFilePaths.sessionROIData;
if ~exist(savePath, 'file')
    error('RFMapping data not found. Please run analyseRFMapping first.');
end
load(savePath, 'RFMapping', 'RFMappingMetadata');

% Extract metadata for plotting
uAz = RFMappingMetadata.uAz;
uEl_plot = RFMappingMetadata.uEl;
timeVector = RFMappingMetadata.timeVector;
nAz = length(uAz); nEl = length(uEl_plot);
nROI = numel(RFMapping);

% Calculate grid offsets for proper scaling
dAz = 20; dEl = 20; 
if nAz > 1, dAz = abs(uAz(2) - uAz(1)); end
if nEl > 1, dEl = abs(uEl_plot(1) - uEl_plot(2)); end
azLimFull = [min(uAz) - dAz/2, max(uAz) + dAz/2];
elLimFull = [min(uEl_plot) - dEl/2, max(uEl_plot) + dEl/2];

% Setup PDF path
outDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end
pdfPath = fullfile(outDir, [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_Responsive_Bouton_Report.pdf']);
if exist(pdfPath, 'file'), delete(pdfPath); end

%% 2. IDENTIFY RESPONSIVE (TUNED) ROIs
% We use the 'isTuned' flag we already calculated 
responsiveIdx = find([RFMapping.isTuned]);
nTuned = length(responsiveIdx);

%% 3. GENERATE SUMMARY COVER PAGE
coverFig = figure('Color', 'w', 'Position', [100 100 600 800], 'Visible', 'off');
summaryText = { ...
    ['Session: ' sessionFileInfo.session_name], ...
    ['Total Boutons: ' num2str(nROI)], ...
    ['Responsive (Tuned) Boutons: ' num2str(nTuned)], ...
    ['% Responsive: ' num2str(round(nTuned/nROI*100, 1)) '%'], ...
    '', ...
    'Criteria: Peak > Blank + 2SD', ... % Consistent with your saved gating 
    ['Response Window: ' num2str(RFMappingMetadata.respWin(1)) ' to ' num2str(RFMappingMetadata.respWin(2)) 's']}; 

annotation('textbox', [0.1 0.1 0.8 0.8], 'String', summaryText, ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
exportgraphics(coverFig, pdfPath); 
close(coverFig);

%% 4. GENERATE INDIVIDUAL ROI PAGES
allCenters = [];
for p = 1:nTuned
    iROI = responsiveIdx(p);
    
    % Store centers for the population map
    allCenters = [allCenters; RFMapping(iROI).centerAz, RFMapping(iROI).centerEl];
    
    fig = figure('Color', 'w', 'Position', [50 50 750 650], 'Visible', 'off');
    ax1 = axes('Position', [0.15 0.15 0.7 0.7]);
    
    % Plot Heatmap
    imagesc(uAz, uEl_plot, imgaussfilt(RFMapping(iROI).meanGridResponse, 0.8)); 
    hold on; colormap(ax1, gray);
    set(ax1, 'YDir', 'normal', 'CLim', [0, max(RFMapping(iROI).meanGridResponse(:)) + 1e-6]);
    xlim(azLimFull); ylim(elLimFull); colorbar;
    
    % PSTH Overlay Positioning
    vS_base = dEl * 0.4; hS = dAz * 0.85;     
    
    for r = 1:nEl
        for c = 1:nAz
            % Blank (White Dotted)
            trB = RFMapping(iROI).meanBlankResponse;
            % Stim (Red Solid)
            trS = RFMapping(iROI).meanTemporalResponse(:, r, c);
            
            % Smooth for visualization
            trB = smoothdata(trB, 'gaussian', 3);
            trS = smoothdata(trS, 'gaussian', 3);
            
            tNorm = (timeVector - min(timeVector)) / (max(timeVector) - min(timeVector));
            sX = (tNorm - 0.5) * hS + uAz(c);
            sY_B = (trB / RFMapping(iROI).peakAmplitude * vS_base) + uEl_plot(r);
            sY_S = (trS / RFMapping(iROI).peakAmplitude * vS_base) + uEl_plot(r);
            
            plot(ax1, sX, sY_B, 'w:', 'LineWidth', 0.8); % Blank baseline
            plot(ax1, sX, sY_S, 'r', 'LineWidth', 1.2);  % Stim response
        end
    end
    
    title(['ROI ' num2str(iROI) ' Receptive Field Summary']);
    xlabel('Azimuth (°)'); ylabel('Elevation (°)');
    
    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType', 'vector');
    close(fig);
end

%% 5. GENERATE POPULATION DENSITY SUMMARY
if ~isempty(allCenters)
    popFig = figure('Color', 'w', 'Position', [100 100 750 650], 'Visible', 'off');
    
    % Use histcounts2 to create the density grid
    % We use uEl_plot (descending) but histcounts2 needs ascending edges
    uEl_asc = sort(uEl_plot, 'ascend');
    [counts, ~, ~] = histcounts2(allCenters(:,1), allCenters(:,2), ...
        [uAz; max(uAz)+dAz], [uEl_asc; max(uEl_asc)+dEl]);
    
    % Plot density map (using bone colormap for "raw" feel)
    imagesc(uAz, uEl_plot, flipud(counts')); 
    colormap(flipud(bone)); colorbar;
    set(gca, 'YDir', 'normal'); xlim(azLimFull); ylim(elLimFull);
    
    % Overlay numerical counts
    for r = 1:nEl
        for c = 1:nAz
            val = counts(c, nEl-r+1);
            if val > 0
                text(uAz(c), uEl_plot(r), num2str(val), 'Color', 'r', ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center'); 
            end
        end
    end
    
    title('Population RF Density (Total Tuned Boutons)');
    xlabel('Azimuth (°)'); ylabel('Elevation (°)');
    exportgraphics(popFig, pdfPath, 'Append', true); 
    close(popFig);
end

fprintf('Report generated: %s\n', pdfPath);
end