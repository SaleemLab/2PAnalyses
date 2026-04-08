filePath = "Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260214A\M25132_20260214A_Response_M25132_RFMapping_20260214_00001.mat";
load(filePath, 'response');

% 
psth = response.psthData;
nROI = size(psth(1).alignedResponses, 1);
stimVs = vertcat(psth.stimValue);


blankIdx = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
if isempty(blankIdx)
    error('Blank trials ([200, 0]) not found. Verify your getTrialGroups logic.');
end

gridMask = stimVs(:,1) ~= 200;
gridPSTH = psth(gridMask);
gridStim = stimVs(gridMask, :);

% Define Grid Dimensions
uAz = sort(unique(gridStim(:,1)), 'ascend');  
uEl = sort(unique(gridStim(:,2)), 'descend'); 
nAz = length(uAz);
nEl = length(uEl);

%% 
respWin = [0.5 2]; % Relative to stimulus onset (0s)
alphaThresh = 0.05;
responsiveROIs = [];

% Pre-calculate blank distribution mask
tVecBlank = psth(blankIdx).timeVector(:);
bMask = tVecBlank >= respWin(1) & tVecBlank <= respWin(2);

fprintf('Analyzing %d ROIs for responsiveness...\n', nROI);

for i = 1:nROI
    %
    rawBlank = psth(blankIdx).alignedResponses(i, :, :);
    blankData = reshape(rawBlank, size(rawBlank, 2), []); 
    blankDist = mean(blankData(bMask, :), 1, 'omitnan');
    
    pVals = ones(numel(gridPSTH), 1);
    for k = 1:numel(gridPSTH)
        rawPos = gridPSTH(k).alignedResponses(i, :, :);
        posData = reshape(rawPos, size(rawPos, 2), []);
        % Calculate response in the same window
        posDist = mean(posData(bMask, :), 1, 'omitnan');
        
        if numel(posDist) > 1 && numel(blankDist) > 1
            [~, pVals(k)] = ttest2(posDist, blankDist, 'Tail', 'right');
        end
    end
    
    % ROI is responsive if it beats the blank in at least one position
    if any(pVals < alphaThresh)
        responsiveROIs = [responsiveROIs; i]; 
    end
end

fprintf('Found %d responsive ROIs out of %d.\n', length(responsiveROIs), nROI);

%% PSTH Grids
for p = 1:min(15, length(responsiveROIs))
    roiIdx = responsiveROIs(p);
    
    figure('Color', 'w', 'Name', sprintf('ROI %d Receptive Field', roiIdx), 'Position', [50 50 1200 700]);
    tlo = tiledlayout(nEl, nAz, 'TileSpacing', 'none', 'Padding', 'compact');
    
    % Global scale for this ROI to make subplots comparable
    allMeans = [];
    
    for r = 1:nEl
        for c = 1:nAz
            targetAz = uAz(c);
            targetEl = uEl(r);
            k = find(gridStim(:,1) == targetAz & gridStim(:,2) == targetEl, 1);
            
            ax = nexttile(tlo);
            if ~isempty(k)
                tVec = gridPSTH(k).timeVector;
                rawP = gridPSTH(k).alignedResponses(roiIdx, :, :);
                pData = reshape(rawP, size(rawP, 2), []);
                mu = mean(pData, 2, 'omitnan');
                
                plot(ax, tVec, mu, 'k', 'LineWidth', 1.2);
                hold(ax, 'on');
                xline(ax, 0, 'r:', 'HandleVisibility', 'off'); 
                
                % Record for global scaling
                allMeans = [allMeans; mu];
                axis(ax, 'tight');
            else
                axis(ax, 'off');
            end
            
            % Labels
            if r == 1, title(ax, sprintf('%d° Az', targetAz), 'FontSize', 9); end
            if c == 1, ylabel(ax, sprintf('%.1f° El', targetEl), 'Visible', 'on', 'FontWeight', 'bold'); end
            set(ax, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 7);
        end
    end
    
    % Apply uniform Y-limits to all tiles for this figure
    set(findobj(gcf, 'type', 'axes'), 'YLim', [min(allMeans) max(allMeans)*1.1]);
    sgtitle(sprintf('ROI %d Spatial PSTH Grid (Time: %.1fs to %.1fs)', roiIdx, min(tVec), max(tVec)));
end




%% 
sigma = 0.8; 
respWin = [0.5 2];
baseWin = [-2 0]; 

for p = 10:min(20, length(responsiveROIs))
    roiIdx = responsiveROIs(p);
    
    % --- STEP 1: Calculate Heatmap & Global Max ---
    rfMatrix = nan(nEl, nAz);
    traceMax = 1e-6; % Initialize with tiny value to avoid div by zero
    
    for k = 1:numel(gridPSTH)
        % Ensure muTrace and tVec are column vectors
        muTrace = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
        muTrace = muTrace(:); 
        tVec = gridPSTH(k).timeVector(:);
        

        bVal = mean(muTrace(tVec >= baseWin(1) & tVec < baseWin(2)), 'omitnan');
        muTrace = muTrace - bVal;
        
        % global peak for scaling red traces
        traceMax = max(traceMax, max(muTrace));
        
        % 
        rCoord = find(uEl == gridStim(k,2), 1);
        cCoord = find(uAz == gridStim(k,1), 1);
        
        if ~isempty(rCoord) && ~isempty(cCoord)
            % Use the cleaned trace for the heatmap value
            wMask = tVec >= respWin(1) & tVec <= respWin(2);
            rfMatrix(rCoord, cCoord) = mean(muTrace(wMask), 'omitnan');
        end
    end
    
    % Smooth the Heatmap
    rfSmoothed = imgaussfilt(rfMatrix, sigma, 'Padding', 'replicate');
    

    figure('Color', 'w', 'Name', sprintf('ROI %d Standard RF', roiIdx), 'Position', [100 100 800 600]);
    

    imagesc(uAz, uEl, rfSmoothed); 
    hold on;
    colormap(gray); 
    if max(rfSmoothed(:)) > 0, clim([0, max(rfSmoothed(:))]); end
    
    cb = colorbar;
    ylabel(cb, '$\Delta F/F$ (Baseline Subtracted)', 'Interpreter', 'latex');
    

    vStep = abs(uEl(1)-uEl(2));
    hStep = abs(uAz(1)-uAz(2));
    vScale = vStep * 0.7; % Height of red trace relative to grid cell
    hScale = hStep * 0.9; % Width of red trace relative to grid cell


    for r = 1:nEl
        for c = 1:nAz
            k = find(gridStim(:,1) == uAz(c) & gridStim(:,2) == uEl(r), 1);
            if ~isempty(k)
                muTrace = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
                muTrace = muTrace(:);
                tVec = gridPSTH(k).timeVector(:);
                
     
                muTrace = muTrace - mean(muTrace(tVec >= baseWin(1) & tVec < baseWin(2)), 'omitnan');
                
                % Global scaling: Biggest peak in whole session = vScale height
                normTrace = muTrace / traceMax;
                
                % ANCHORING: Baseline sits exactly on the elevation coordinate line
                scaledY = (normTrace * vScale) + uEl(r); 
                
                % Horizontal centering on the azimuth
                normTime = (tVec - min(tVec)) / (max(tVec) - min(tVec));
                scaledX = (normTime - 0.5) * hScale + uAz(c);
    
                plot(scaledX, scaledY, 'r', 'LineWidth', 1.2); 
                
                % Stimulus onset (Time 0) reference marker
                [~, t0Idx] = min(abs(tVec));
                plot(scaledX(t0Idx), scaledY(t0Idx), 'k.', 'MarkerSize', 6);
            end
        end
    end
    

    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
    xlabel('Azimuth (deg)');
    ylabel('Elevation (deg)');
    title(sprintf('ROI %d: Receptive Field Mapping (Baseline Corrected)', roiIdx));
    axis tight;
    % grid on;
    set(gca, 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.4);
end


%% 
sigma_smooth = 0.8; 
respWin = [0.5 2]; 
baseWin = [-2 0];  

for p = 1:min(10, length(responsiveROIs))
    roiIdx = responsiveROIs(p);
    
    % --- STEP 1: Calculate Heatmap & Global Max ---
    rfMatrix = nan(nEl, nAz);
    traceMax = 1e-6; 
    
    for k = 1:numel(gridPSTH)
        muTrace = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
        tVec = gridPSTH(k).timeVector(:);
        
        % Baseline Subtraction
        bVal = mean(muTrace(tVec >= baseWin(1) & tVec < baseWin(2)), 'omitnan');
        muTrace = muTrace(:) - bVal;
        traceMax = max(traceMax, max(muTrace));
        
        % Corrected indexing for 7x3 grid 
        rIdx = find(uEl == gridStim(k,2), 1);
        cIdx = find(uAz == gridStim(k,1), 1);
        
        if ~isempty(rIdx) && ~isempty(cIdx)
            wMask = tVec >= respWin(1) & tVec <= respWin(2);
            rfMatrix(rIdx, cIdx) = mean(muTrace(wMask), 'omitnan');
        end
    end


    [X, Y] = meshgrid(uAz, uEl);
    xdata(:,:,1) = X; xdata(:,:,2) = Y;
    

    [maxAmp, maxIdx] = max(rfMatrix(:));
    x0_guess = [maxAmp, X(maxIdx), Y(maxIdx), 20, 20]; 
    gauss2D = @(x, xdata) x(1) * exp( -((xdata(:,:,1)-x(2)).^2/(2*x(4)^2) + (xdata(:,:,2)-x(3)).^2/(2*x(5)^2)) );
    
    opts = optimset('Display','off');
    fitP = lsqcurvefit(gauss2D, x0_guess, xdata, rfMatrix, [0 min(uAz) min(uEl) 2 2], [inf max(uAz) max(uEl) 60 60], opts);

    figure('Color', 'w', 'Position', [100 100 1100 500]);
    tlo = tiledlayout(1, 2, 'TileSpacing', 'compact');
    
    axLeft = nexttile(tlo);
    imagesc(uAz, uEl, imgaussfilt(rfMatrix, sigma_smooth)); hold on;
    colormap(axLeft, gray); clim([0, max(rfMatrix(:)) + 1e-6]);
    
    vScale = abs(uEl(1)-uEl(2)) * 0.7; hScale = abs(uAz(1)-uAz(2)) * 0.9;
    for k = 1:numel(gridPSTH)
        tr = mean(gridPSTH(k).alignedResponses(roiIdx, :, :), 3, 'omitnan');
        tr = tr(:) - mean(tr(tVec >= baseWin(1) & tVec < baseWin(2)), 'omitnan');
        sY = (tr/traceMax * vScale) + gridStim(k,2);
        sX = (normalize(tVec, 'range') - 0.5) * hScale + gridStim(k,1);
        plot(axLeft, sX, sY, 'r', 'LineWidth', 1);
    end
    title('Raw Data: Heatmap + Traces'); xlabel('azimuth (°)'); ylabel('elevation (°)');
    set(axLeft, 'YDir', 'normal'); axis tight;

    % 
    axRight = nexttile(tlo);
 
    [Xf, Yf] = meshgrid(linspace(min(uAz), max(uAz), 100), linspace(min(uEl), max(uEl), 100));
    xfdata(:,:,1) = Xf; xfdata(:,:,2) = Yf;
    Zfit = gauss2D(fitP, xfdata);
    
    imagesc(linspace(min(uAz), max(uAz), 100), linspace(min(uEl), max(uEl), 100), Zfit); hold on;
    colormap(axRight, parula); % The "Yellow-Blue" look from the paper
    
    % Draw the 1-sigma contour (RF Extent) 
    contour(Xf, Yf, Zfit, [fitP(1)*exp(-0.5) fitP(1)*exp(-0.5)], 'w--', 'LineWidth', 2);
    
    title(sprintf('Fit: \\sigma_{az}=%.1f, \\sigma_{el}=%.1f', fitP(4), fitP(5)));
    xlabel('azimuth (°)'); set(axRight, 'YDir', 'normal');
    axis tight; colorbar;
end