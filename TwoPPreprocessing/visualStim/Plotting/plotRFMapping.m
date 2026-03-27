function plotRFMapping(sessionFileInfo, RFMapping, RFMappingMetadata, doSmooth)
% plotRFMapping: Standalone plotting with optional temporal smoothing.
% doSmooth = true applies a 3-sample Gaussian filter

if nargin < 4; doSmooth = true; end

% load save path from sessionFileInfo
saveFolder = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~isfolder(saveFolder); mkdir(saveFolder); end
pdfPath = fullfile(saveFolder, [sessionFileInfo.animal_name, '_' sessionFileInfo.session_name '_RF_Summary.pdf']);
if exist(pdfPath, 'file'), delete(pdfPath); end


%% 
uAz = RFMappingMetadata.uAz;
uEl_plot = RFMappingMetadata.uEl;
timeVector = RFMappingMetadata.timeVector;
stimName = RFMappingMetadata.stimName;
nROI = numel(RFMapping);
nAz = length(uAz); nEl = length(uEl_plot);
sigma_smooth = 0.8;

% Set temporal smoothing window
% At 60Hz, 3 samples = 50ms (Very safe)
if doSmooth; tWin = 3; else; tWin = 1; end

% Grid spacing for plotting
dAz = 20; dEl = 20; 
if nAz > 1, dAz = abs(uAz(2) - uAz(1)); end
if nEl > 1, dEl = abs(uEl_plot(1) - uEl_plot(2)); end
azLimFull = [min(uAz) - dAz/2, max(uAz) + dAz/2];
elLimFull = [min(uEl_plot) - dEl/2, max(uEl_plot) + dEl/2];



for iROI = 1:nROI
    fig = figure('Color', 'w', 'Position', [50 50 850 650], 'Visible', 'off');
    ax1 = axes('Position', [0.15 0.15 0.7 0.7]);
    
    % Background Heatmap (The Mean Grid Response)
    imagesc(uAz, uEl_plot, imgaussfilt(RFMapping(iROI).meanGridResponse, sigma_smooth)); 
    hold on; colormap(ax1, gray);
    set(ax1, 'YDir', 'normal', 'CLim', [0, max(RFMapping(iROI).meanGridResponse(:)) + 1e-6]);
    xlim(azLimFull); ylim(elLimFull); 
    
    cb = colorbar;
    cb.Label.String = 'z-scored dF/F';
    
    vS_base = dEl * 0.4; hS = dAz * 0.85;     
    for r = 1:nEl
        for c = 1:nAz
            % Blank (Reference)
            trB = RFMapping(iROI).meanBlankResponse;
            % Stimulus
            trS = RFMapping(iROI).meanTemporalResponse(:, r, c);
            
            % Optional Smoothing
            if doSmooth
                trB = smoothdata(trB, 'gaussian', tWin);
                trS = smoothdata(trS, 'gaussian', tWin);
            end

            tNorm = (timeVector - min(timeVector)) / (max(timeVector) - min(timeVector));
            sX = (tNorm - 0.5) * hS + uAz(c);
            
            sY_B = (trB / RFMapping(iROI).peakAmplitude * vS_base) + uEl_plot(r);
            sY_S = (trS / RFMapping(iROI).peakAmplitude * vS_base) + uEl_plot(r);

            % Plot Blank (Gray Dotted)
            plot(ax1, sX, sY_B, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'LineStyle', ':'); 
            % Plot Stimulus (Red)
            plot(ax1, sX, sY_S, 'r', 'LineWidth', 1.2); 
        end
    end
    title(['ROI ' num2str(iROI) ' Heatmap + PSTHs']);
    xlabel('azimuth (°)'); ylabel('elevation (°)');

    exportgraphics(fig, pdfPath, 'Append', true); close(fig);
end
end