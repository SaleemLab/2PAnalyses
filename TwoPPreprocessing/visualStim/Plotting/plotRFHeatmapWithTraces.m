function plotRFHeatmapWithTraces(ax, boutonData, uAz, uEl_plot, timeVector, varargin)
% plotRFHeatmapWithTraces: Draws the spatial tuning heatmap for a single
% bouton, with the mean temporal trace at each Az/El position overlaid
% directly on top at its correct screen coordinate, plus the blank trace
% for comparison.
% 
% INPUTS:
%   ax          - target axes handle to plot into
%   boutonData  - single element of an RFMapping struct array, must contain:
%                   .meanGridResponse (nEl x nAz)
%                   .meanTemporalResponse (nTpts x nEl x nAz)
%                   .meanBlankResponse (nTpts x 1)
%                   .peakAmplitude
%                   .centerAz, .centerEl, .isResponsive
%   uAz         - vector of unique azimuth positions
%   uEl_plot    - vector of unique elevation positions
%   timeVector  - shared time vector (nTpts x 1)
% 
% NAME-VALUE OPTIONS:
%   'Colormap'     - colormap name or Nx3 matrix. Default: 'bone'.
%   'MapMethod'    - 'mean' (use meanGridResponse directly, default), 'svd'
%                    (denoise via computeSVDMap), or 'gaussianfit'
%   'RespWin'      - [start end] response window (s), used for filtering and SVD.
%   'PlotRespWin'  - true/false. If true, ONLY plots the segment of the traces 
%                    that falls within the RespWin range.
%   'Smooth'       - true/false.
%   'AxisMode'     - 'square' (default, matches prior behavior: square AXES BOX,
%                    independent of data range) or 'image' (1 data-unit-az = 1
%                    data-unit-el, box cropped to data extent -- use this to match
%                    panels built elsewhere with axis(ax,'image'), e.g. a Gaussian
%                    fit surface plotted via imagesc + axis image).
%   'LabelUnits'   - '\circ' (default, matches prior behavior: 'Azimuth (\circ)') or
%                    'deg' (produces 'Azimuth (deg)') -- use 'deg' to match panels
%                    elsewhere that label axes as "(deg)" instead of the degree symbol.
%   'StimDiameterDeg' - []  (default, matches prior behavior: each tile drawn at the
%                    GRID STEP size, i.e. tiles are contiguous and non-overlapping,
%                    regardless of how wide the actual stimulus was). If set to a
%                    number (e.g. 30), each tile is instead drawn at that width/height,
%                    centered on its grid position -- so tiles VISUALLY OVERLAP exactly
%                    as much as the real stimuli did on screen (e.g. 30 deg stimuli on
%                    a 20 deg grid overlap by 10 deg with each neighbor). Overlapping
%                    tiles are alpha-blended (FaceAlpha) so overlap is visible rather
%                    than one tile silently occluding another. Only applies to
%                    MapMethod='mean'/'svd' with Smooth=false; ignored for 'gaussianfit'
%                    or Smooth=true, which already use their own continuous rendering.
%                    NOTE: superseded by ScreenAzLimits/ScreenElLimits below if both are
%                    given -- ScreenAzLimits/ScreenElLimits is the recommended option.
%   'ScreenAzLimits'/'ScreenElLimits' - [] (default, matches prior behavior: tiles are
%                    CONTIGUOUS with a fixed half-grid-step border, leaving blank axis
%                    space beyond the outermost tile centers). If set to a 2-element
%                    range (e.g. [-80 30] / [-25 40], matching the true screen extent),
%                    tiles remain CONTIGUOUS (no overlap, no alpha blending, no striping)
%                    but the OUTER edge of the first/last tile on each axis is stretched
%                    to meet the given screen limit exactly, instead of stopping at
%                    half a grid step past the outermost center. This eliminates blank
%                    axis space entirely while staying geometrically honest: interior
%                    tile boundaries are still plain midpoints between neighboring grid
%                    centers (unchanged), only the two outermost edges move. Overlaid
%                    traces (including the blank-condition trace) are rescaled to fit
%                    each tile's own actual width/height, so edge-tile traces widen to
%                    fill their now-larger tile rather than staying pinned to the old,
%                    narrower grid-step size.
%   ... (other options inherit standard defaults)

    p = struct('Colormap', 'bone', 'MapMethod', 'mean', 'RespWin', [], ...
                'PlotRespWin', false, 'Smooth', false, 'ScaleFactor', 10, ...
                'Sigma', 3, 'GaussRes', 200, 'ShowFitCenter', true, ...
                'AxisMode', 'square', 'LabelUnits', '\circ', 'StimDiameterDeg', [], ...
                'ScreenAzLimits', [], 'ScreenElLimits', []);
    for k = 1:2:numel(varargin)
        p.(varargin{k}) = varargin{k+1};
    end
    nAz = length(uAz);
    nEl = length(uEl_plot);

    % check that the grid is evenly spaced
    azSteps = abs(diff(sort(uAz)));
    elSteps = abs(diff(sort(uEl_plot)));
    if numel(azSteps) > 1 && ~all(abs(azSteps - azSteps(1)) < 1e-6)
        error('plotRFHeatmapWithTraces:nonUniformGrid', ...
            'uAz is not evenly spaced — cell sizing/box placement will be wrong.');
    end
    if numel(elSteps) > 1 && ~all(abs(elSteps - elSteps(1)) < 1e-6)
        error('plotRFHeatmapWithTraces:nonUniformGrid', ...
            'uEl_plot is not evenly spaced — cell sizing/box placement will be wrong.');
    end
    
    dAz = abs(uAz(2) - uAz(1));
    dEl = abs(uEl_plot(1) - uEl_plot(2));
    tempStack  = boutonData.meanTemporalResponse;   % nTpts x nEl x nAz
    blankTrace = boutonData.meanBlankResponse;
    peakAmp    = max(boutonData.peakAmplitude, eps); % avoid divide-by-zero
    isGaussFitMap = strcmpi(p.MapMethod, 'gaussianfit');
    useScreenEdges = ~isempty(p.ScreenAzLimits) && ~isempty(p.ScreenElLimits) ...
        && ~isGaussFitMap && ~p.Smooth;

    if useScreenEdges
        % Non-uniform cell EDGES: interior boundaries are plain midpoints between
        % neighboring grid centers (same as contiguous imagesc behavior); only the two
        % outer edges are stretched out to the true screen limit instead of stopping at
        % half a grid step. Assumes uAz/uEl_plot are given in monotonic (ascending or
        % descending) order, consistent with how imagesc already interprets them elsewhere
        % in this function.
        azEdges = zeros(1, nAz+1);
        azEdges(2:end-1) = (uAz(1:end-1) + uAz(2:end)) / 2;
        azEdges([1 end]) = sort(p.ScreenAzLimits);
        if uAz(1) > uAz(end)
            azEdges([1 end]) = azEdges([end 1]);
        end
        elEdges = zeros(1, nEl+1);
        elEdges(2:end-1) = (uEl_plot(1:end-1) + uEl_plot(2:end)) / 2;
        elEdges([1 end]) = sort(p.ScreenElLimits);
        if uEl_plot(1) > uEl_plot(end)
            elEdges([1 end]) = elEdges([end 1]);
        end
    end

    switch lower(p.MapMethod)
        case 'mean'
            rawMap = boutonData.meanGridResponse;
        case 'svd'
            rawMap = computeSVDMap(tempStack, timeVector, p.RespWin);
        case 'gaussianfit'
            requiredFields = {'gaussFit_A', 'gaussFit_Az0', 'gaussFit_El0', 'gaussFit_sigmaAz', 'gaussFit_sigmaEl'};
            missingFields  = requiredFields(~isfield(boutonData, requiredFields));
            if ~isempty(missingFields) || any(cellfun(@(f) isempty(boutonData.(f)), requiredFields(isfield(boutonData, requiredFields))))
                error('plotRFHeatmapWithTraces:missingGaussFit', ...
                    'MapMethod=''gaussianfit'' requires gaussFit_A/Az0/El0/sigmaAz/sigmaEl fields.');
            end
            rawMap = []; 
        otherwise
            error('plotRFHeatmapWithTraces:badMapMethod', 'Unknown MapMethod.');
    end

    if isGaussFitMap
        azStep = uAz(2) - uAz(1);
        elStep = uEl_plot(2) - uEl_plot(1);
        xImg = [uAz(1) - azStep/2, uAz(end) + azStep/2];
        yImg = [uEl_plot(1) - elStep/2, uEl_plot(end) + elStep/2];
        xFine = linspace(min(xImg), max(xImg), p.GaussRes);
        yFine = linspace(min(yImg), max(yImg), p.GaussRes);
        [XFine, YFine] = meshgrid(xFine, yFine);
        displayMap = boutonData.gaussFit_A * exp( -( (XFine - boutonData.gaussFit_Az0).^2 ./ (2*boutonData.gaussFit_sigmaAz^2) + (YFine - boutonData.gaussFit_El0).^2 ./ (2*boutonData.gaussFit_sigmaEl^2) ) );
    elseif p.Smooth
        displayMap = smoothRFMapForDisplay(rawMap, p.ScaleFactor, p.Sigma);
        azStep = uAz(2) - uAz(1);
        elStep = uEl_plot(2) - uEl_plot(1);
        xImg = [uAz(1) - azStep/2, uAz(end) + azStep/2];
        yImg = [uEl_plot(1) - elStep/2, uEl_plot(end) + elStep/2];
    else
        displayMap = rawMap;
        xImg = uAz;
        yImg = uEl_plot;
    end

    axes(ax);
    if useScreenEdges
        % Contiguous tiles (no overlap, no alpha blending -- avoids the striping problem
        % from the earlier StimDiameterDeg overlap approach), but the outermost tile on
        % each side is stretched to the true screen edge instead of stopping at half a
        % grid step. One patch per cell, solid fill, colored the same way imagesc would.
        cmapData = colormap(ax, p.Colormap);
        climLo = min(displayMap(:));
        climHi = max(displayMap(:)) + 1e-6;
        hold(ax, 'on');
        for r = 1:nEl
            for c = 1:nAz
                val = displayMap(r, c);
                if isnan(val)
                    cellColor = cmapData(1, :);
                else
                    colorIdx = round( (val - climLo) / max(climHi - climLo, eps) * (size(cmapData,1)-1) ) + 1;
                    colorIdx = min(max(colorIdx, 1), size(cmapData,1));
                    cellColor = cmapData(colorIdx, :);
                end
                patch(ax, [azEdges(c), azEdges(c+1), azEdges(c+1), azEdges(c)], ...
                          [elEdges(r), elEdges(r), elEdges(r+1), elEdges(r+1)], ...
                          cellColor, 'EdgeColor', 'none');
            end
        end
        set(ax, 'YDir', 'normal', 'CLim', [climLo, climHi]);
    else
        imagesc(ax, xImg, yImg, displayMap);
        hold(ax, 'on');
        colormap(ax, p.Colormap);
        set(ax, 'YDir', 'normal', 'CLim', [min(displayMap(:)), max(displayMap(:)) + 1e-6]);
    end

    % ----------------------------------------------------
    % SELECTIVE WINDOW FILTERING
    % ----------------------------------------------------
    if p.PlotRespWin && ~isempty(p.RespWin)
        % Filter to only include time points inside the specified response window
        winIdx = (timeVector >= p.RespWin(1)) & (timeVector <= p.RespWin(2));
        subTime = timeVector(winIdx);
        subTempStack = tempStack(winIdx, :, :);
        subBlank = blankTrace(winIdx);
    else
        % Use the full trace
        subTime = timeVector;
        subTempStack = tempStack;
        subBlank = blankTrace;
    end

    % Trace overlay scaling and geometry setup
    tMin   = min(subTime);
    tMax   = max(subTime);
    tRange = max(tMax - tMin, eps);
    tNorm  = (subTime - tMin) / tRange;

    for r = 1:nEl
        for c = 1:nAz
            cellAz = uAz(c);
            cellEl = uEl_plot(r);

            if useScreenEdges
                % Use THIS cell's actual width/height (which may be wider at the outer
                % edges, now that they're stretched to the screen limit) instead of the
                % fixed global dAz/dEl -- so edge-tile traces fill their larger tile
                % rather than staying pinned to the old, narrower grid-step size.
                cellW = abs(azEdges(c+1) - azEdges(c));
                cellH = abs(elEdges(r+1) - elEdges(r));
            else
                cellW = dAz;
                cellH = dEl;
            end
            vScale = cellH * 0.4;
            hScale = cellW * 0.85;

            stimTrace = subTempStack(:, r, c);
            sX      = (tNorm - 0.5) * hScale + cellAz;
            sY_stim = (stimTrace  / peakAmp * vScale) + cellEl;

            % Blank trace comparison (cropped to same window) -- kept exactly as before,
            % just resized to match this cell's (possibly wider) footprint.
            sY_blank = (subBlank / peakAmp * vScale) + cellEl;
            plot(ax, sX, sY_blank, 'Color', [1 1 1 0.85], 'LineWidth', 0.8);
            
            % Active stimulus trace (cropped to same window)
            plot(ax, sX, sY_stim, 'Color', [1 0.3 0.2], 'LineWidth', 1.6);
        end
    end
    
    if isGaussFitMap && p.ShowFitCenter
        theta = linspace(0, 2*pi, 100);
        ellX  = boutonData.gaussFit_Az0 + boutonData.gaussFit_sigmaAz * cos(theta);
        ellY  = boutonData.gaussFit_El0 + boutonData.gaussFit_sigmaEl * sin(theta);
        plot(ax, ellX, ellY, '--', 'Color', [0.3 1 1], 'LineWidth', 1.2);
        plot(ax, boutonData.gaussFit_Az0, boutonData.gaussFit_El0, '+', ...
            'Color', [0.3 1 1], 'MarkerSize', 10, 'LineWidth', 1.5);
    end

    if useScreenEdges
        xlim(ax, sort(p.ScreenAzLimits));
        ylim(ax, sort(p.ScreenElLimits));
    else
        xlim(ax, [min(uAz)-dAz/2, max(uAz)+dAz/2]);
        ylim(ax, [min(uEl_plot)-dEl/2, max(uEl_plot)+dEl/2]);
    end
    if strcmpi(p.LabelUnits, 'deg')
        azLabelStr = 'Azimuth (deg)';
        elLabelStr = 'Elevation (deg)';
    else
        azLabelStr = 'Azimuth (\circ)';
        elLabelStr = 'Elevation (\circ)';
    end
    xlabel(ax, azLabelStr, 'FontName', 'Arial', 'FontSize', 9);
    ylabel(ax, elLabelStr, 'FontName', 'Arial', 'FontSize', 9);
    set(ax, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
    set(ax, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
    axis(ax, p.AxisMode);
    cb = colorbar(ax);
    cb.Label.String = '\Delta F/F';
    cb.FontSize = 7;
end

% function plotRFHeatmapWithTraces(ax, boutonData, uAz, uEl_plot, timeVector, varargin)
% % plotRFHeatmapWithTraces: Draws the spatial tuning heatmap for a single
% % bouton, with the mean temporal trace at each Az/El position overlaid
% % directly on top at its correct screen coordinate, plus the blank trace
% % for comparison.
% % 
% % INPUTS:
% %   ax          - target axes handle to plot into
% %   boutonData  - single element of an RFMapping struct array, must contain:
% %                   .meanGridResponse (nEl x nAz)
% %                   .meanTemporalResponse (nTpts x nEl x nAz)
% %                   .meanBlankResponse (nTpts x 1)
% %                   .peakAmplitude
% %                   .centerAz, .centerEl, .isResponsive
% %   uAz         - vector of unique azimuth positions
% %   uEl_plot    - vector of unique elevation positions
% %   timeVector  - shared time vector (nTpts x 1)
% % 
% % NAME-VALUE OPTIONS:
% %   'Colormap'    - colormap name or Nx3 matrix. Default: 'bone'.
% %   'MapMethod'   - 'mean' (use meanGridResponse directly, default), 'svd'
% %                   (denoise via computeSVDMap), or 'gaussianfit' (render
% %                   the fitted 2D Gaussian model from RFMapping_Gaussian2DFit.m
% %                   as a smooth analytic surface -- requires boutonData to
% %                   have gaussFit_A, gaussFit_Az0, gaussFit_El0,
% %                   gaussFit_sigmaAz, gaussFit_sigmaEl fields set).
% %   'RespWin'     - [start end] response window (s), passed to computeSVDMap
% %                   if MapMethod is 'svd'.
% %   'Smooth'      - true/false. If true (default), background is upsampled
% %                   + Gaussian-blurred for a continuous-looking display. If
% %                   false, background is plotted as the RAW grid — each
% %                   pixel is exactly one uniform-size Az/El cell (dAz x dEl),
% %                   with no interpolation, so cell boundaries are visible
% %                   and you can directly confirm each trace is centered in
% %                   its own cell. Recommended when checking alignment.
% %                   (Ignored when MapMethod is 'gaussianfit' -- that surface
% %                   is always evaluated smoothly/analytically, see GaussRes.)
% %   'ScaleFactor' - upsampling factor for display smoothing (only used if
% %                   Smooth=true). Default: 10.
% %   'Sigma'       - Gaussian blur sigma (only used if Smooth=true). Default: 3.
% %   'GaussRes'    - resolution (points per side) for the analytic Gaussian
% %                   surface when MapMethod='gaussianfit'. Default: 200.
% %   'ShowFitCenter' - true/false, only used when MapMethod='gaussianfit'.
% %                   If true (default), overlays the fitted center (x0,y0)
% %                   as a small marker and the 1-sigma ellipse, so you can
% %                   directly see whether the traces line up with where the
% %                   fit says the true peak is.
% % 
% % USAGE:
% %   plotRFHeatmapWithTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector);
% %   plotRFHeatmapWithTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector, ...
% %       'Smooth', false);   % raw uniform grid cells, for alignment checking
% %   plotRFHeatmapWithTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector, ...
% %       'MapMethod', 'gaussianfit');   % smooth fitted-model background instead of raw/blurred data
% % 
% % SGS code for SparseNoise analysis applied here - specifially the
% % smoothing and SVD 
% % 
%     p = struct('Colormap', 'bone', 'MapMethod', 'mean', 'RespWin', [], ...
%                 'Smooth', false, 'ScaleFactor', 10, 'Sigma', 3, ...
%                 'GaussRes', 200, 'ShowFitCenter', true);
%     for k = 1:2:numel(varargin)
%         p.(varargin{k}) = varargin{k+1};
%     end
% 
%     nAz = length(uAz);
%     nEl = length(uEl_plot);
% 
% %     check that the grid is evenly spaced; dAz/dEl below assume this ---
%     azSteps = abs(diff(sort(uAz)));
%     elSteps = abs(diff(sort(uEl_plot)));
%     if numel(azSteps) > 1 && ~all(abs(azSteps - azSteps(1)) < 1e-6)
%         error('plotRFHeatmapWithTraces:nonUniformGrid', ...
%             'uAz is not evenly spaced (steps: %s) — cell sizing/box placement will be wrong.', ...
%             mat2str(azSteps));
%     end
%     if numel(elSteps) > 1 && ~all(abs(elSteps - elSteps(1)) < 1e-6)
%         error('plotRFHeatmapWithTraces:nonUniformGrid', ...
%             'uEl_plot is not evenly spaced (steps: %s) — cell sizing/box placement will be wrong.', ...
%             mat2str(elSteps));
%     end
% 
%     dAz = abs(uAz(2) - uAz(1));
%     dEl = abs(uEl_plot(1) - uEl_plot(2));
% 
%     tempStack  = boutonData.meanTemporalResponse;   % nTpts x nEl x nAz
%     blankTrace = boutonData.meanBlankResponse;
%     peakAmp    = max(boutonData.peakAmplitude, eps); % avoid divide-by-zero
% 
%     isGaussFitMap = strcmpi(p.MapMethod, 'gaussianfit');
% 
%     switch lower(p.MapMethod)
%         case 'mean'
%             rawMap = boutonData.meanGridResponse;
%         case 'svd'
%             rawMap = computeSVDMap(tempStack, timeVector, p.RespWin);
%         case 'gaussianfit'
%             requiredFields = {'gaussFit_A', 'gaussFit_Az0', 'gaussFit_El0', 'gaussFit_sigmaAz', 'gaussFit_sigmaEl'};
%             missingFields  = requiredFields(~isfield(boutonData, requiredFields));
%             if ~isempty(missingFields) || any(cellfun(@(f) isempty(boutonData.(f)), requiredFields(isfield(boutonData, requiredFields))))
%                 error('plotRFHeatmapWithTraces:missingGaussFit', ...
%                     ['MapMethod=''gaussianfit'' requires gaussFit_A/Az0/El0/sigmaAz/sigmaEl on boutonData ' ...
%                      '(run RFMapping_Gaussian2DFit.m first). Missing/empty: %s'], strjoin(missingFields, ', '));
%             end
%             rawMap = []; % not used for this method -- built directly as displayMap below
%         otherwise
%             error('plotRFHeatmapWithTraces:badMapMethod', ...
%                 'MapMethod must be ''mean'', ''svd'', or ''gaussianfit'', got ''%s''.', p.MapMethod);
%     end
% % 
%     if isGaussFitMap
% %         Evaluate the fitted 2D Gaussian directly on a fine mesh -- this is an ANALYTIC
% %         surface, not a blurred version of discrete data, so it doesn't go through
% %         smoothRFMapForDisplay at all (Smooth/ScaleFactor/Sigma are ignored here).
%         azStep = uAz(2) - uAz(1);
%         elStep = uEl_plot(2) - uEl_plot(1);
%         xImg = [uAz(1) - azStep/2, uAz(end) + azStep/2];
%         yImg = [uEl_plot(1) - elStep/2, uEl_plot(end) + elStep/2];
% 
%         xFine = linspace(min(xImg), max(xImg), p.GaussRes);
%         yFine = linspace(min(yImg), max(yImg), p.GaussRes);
%         [XFine, YFine] = meshgrid(xFine, yFine);
% 
%         A0      = boutonData.gaussFit_A;
%         x0      = boutonData.gaussFit_Az0;
%         y0      = boutonData.gaussFit_El0;
%         sigmaX0 = boutonData.gaussFit_sigmaAz;
%         sigmaY0 = boutonData.gaussFit_sigmaEl;
% 
%         displayMap = A0 * exp( -( (XFine - x0).^2 ./ (2*sigmaX0^2) + (YFine - y0).^2 ./ (2*sigmaY0^2) ) );
%     elseif p.Smooth
%         displayMap = smoothRFMapForDisplay(rawMap, p.ScaleFactor, p.Sigma);
% %         imagesc only auto-extends half a cell beyond the outer edge when
% %         numel(x) matches size(C,2) (i.e. unsmoothed/raw grid). Once the
% %         map is upsampled, size(C) no longer matches uAz/uEl_plot, so we
% %         must give imagesc the true padded edges explicitly or it will
% %         only span from uAz(1) to uAz(end) (leaving a blank margin).
%         azStep = uAz(2) - uAz(1);
%         elStep = uEl_plot(2) - uEl_plot(1);
%         xImg = [uAz(1) - azStep/2, uAz(end) + azStep/2];
%         yImg = [uEl_plot(1) - elStep/2, uEl_plot(end) + elStep/2];
%     else
%         displayMap = rawMap; % raw grid: each pixel = exactly one uniform Az/El cell
%         xImg = uAz;
%         yImg = uEl_plot;
%     end
% 
%     axes(ax); 
%     imagesc(ax, xImg, yImg, displayMap);
%     hold(ax, 'on');
%     colormap(ax, p.Colormap);
%     set(ax, 'YDir', 'normal', 'CLim', [min(displayMap(:)), max(displayMap(:)) + 1e-6]);
% % 
% %     (yellow "responsive" box removed -- background + trace already agree)
% % 
% %     trace overlay geometry: each trace spans ~85% of a grid cell width,
% %     and ~40% of a cell height, centered on that cell's coordinate
%     vScale = dEl * 0.4;
%     hScale = dAz * 0.85;
%     tNorm  = (timeVector - min(timeVector)) / (max(timeVector) - min(timeVector));
% 
%     for r = 1:nEl
%         for c = 1:nAz
%             stimTrace = tempStack(:, r, c);
% 
%             sX      = (tNorm - 0.5) * hScale + uAz(c);
%             sY_stim = (stimTrace  / peakAmp * vScale) + uEl_plot(r);
% % 
% %             blank/baseline trace overlay removed for a cleaner look; to
% %             bring it back, uncompute sY_blank and re-add the dotted plot:
%               sY_blank = (blankTrace / peakAmp * vScale) + uEl_plot(r);
%               plot(ax, sX, sY_blank, 'Color', [1 1 1 0.85], 'LineWidth', 0.8, 'LineStyle', ':');
%             plot(ax, sX, sY_stim, 'Color', [1 0.3 0.2], 'LineWidth', 1.1);
%         end
%     end
% 
%     if isGaussFitMap && p.ShowFitCenter
% %         overlay the fitted center + 1-sigma ellipse so you can see directly whether the
% %         traces line up with where the model says the true peak/extent actually is
%         theta = linspace(0, 2*pi, 100);
%         ellX  = boutonData.gaussFit_Az0 + boutonData.gaussFit_sigmaAz * cos(theta);
%         ellY  = boutonData.gaussFit_El0 + boutonData.gaussFit_sigmaEl * sin(theta);
%         plot(ax, ellX, ellY, '--', 'Color', [0.3 1 1], 'LineWidth', 1.2);
%         plot(ax, boutonData.gaussFit_Az0, boutonData.gaussFit_El0, '+', ...
%             'Color', [0.3 1 1], 'MarkerSize', 10, 'LineWidth', 1.5);
%     end
% 
%     xlim(ax, [min(uAz)-dAz/2, max(uAz)+dAz/2]);
%     ylim(ax, [min(uEl_plot)-dEl/2, max(uEl_plot)+dEl/2]);
%     xlabel(ax, 'Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
%     ylabel(ax, 'Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
% %     ticks placed exactly at grid positions -- read the label directly
% %     under/beside a trace to confirm which Az/El it's centered on
%     set(ax, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
%     set(ax, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
%     axis(ax, 'square');
% 
%     cb = colorbar(ax);
%     cb.Label.String = '\Delta F/F';
%     cb.FontSize = 7;
% end