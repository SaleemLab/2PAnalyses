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
%   'Colormap'    - colormap name or Nx3 matrix. Default: 'bone'.
%   'MapMethod'   - 'mean' (use meanGridResponse directly, default) or
%                   'svd' (denoise via computeSVDMap).
%   'RespWin'     - [start end] response window (s), passed to computeSVDMap
%                   if MapMethod is 'svd'.
%   'Smooth'      - true/false. If true (default), background is upsampled
%                   + Gaussian-blurred for a continuous-looking display. If
%                   false, background is plotted as the RAW grid — each
%                   pixel is exactly one uniform-size Az/El cell (dAz x dEl),
%                   with no interpolation, so cell boundaries are visible
%                   and you can directly confirm each trace is centered in
%                   its own cell. Recommended when checking alignment.
%   'ScaleFactor' - upsampling factor for display smoothing (only used if
%                   Smooth=true). Default: 10.
%   'Sigma'       - Gaussian blur sigma (only used if Smooth=true). Default: 3.
%
% USAGE:
%   plotRFHeatmapWithTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector);
%   plotRFHeatmapWithTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector, ...
%       'Smooth', false);   % raw uniform grid cells, for alignment checking

% SGS code for SparseNoise analysis applied here - specifially the
% smoothing and SVD 

    p = struct('Colormap', 'bone', 'MapMethod', 'mean', 'RespWin', [], ...
                'Smooth', false, 'ScaleFactor', 10, 'Sigma', 3);
    for k = 1:2:numel(varargin)
        p.(varargin{k}) = varargin{k+1};
    end

    nAz = length(uAz);
    nEl = length(uEl_plot);

    % check that the grid is evenly spaced; dAz/dEl below assume this ---
    azSteps = abs(diff(sort(uAz)));
    elSteps = abs(diff(sort(uEl_plot)));
    if numel(azSteps) > 1 && ~all(abs(azSteps - azSteps(1)) < 1e-6)
        error('plotRFHeatmapWithTraces:nonUniformGrid', ...
            'uAz is not evenly spaced (steps: %s) — cell sizing/box placement will be wrong.', ...
            mat2str(azSteps));
    end
    if numel(elSteps) > 1 && ~all(abs(elSteps - elSteps(1)) < 1e-6)
        error('plotRFHeatmapWithTraces:nonUniformGrid', ...
            'uEl_plot is not evenly spaced (steps: %s) — cell sizing/box placement will be wrong.', ...
            mat2str(elSteps));
    end

    dAz = abs(uAz(2) - uAz(1));
    dEl = abs(uEl_plot(1) - uEl_plot(2));

    tempStack  = boutonData.meanTemporalResponse;   % nTpts x nEl x nAz
    blankTrace = boutonData.meanBlankResponse;
    peakAmp    = max(boutonData.peakAmplitude, eps); % avoid divide-by-zero

    switch lower(p.MapMethod)
        case 'mean'
            rawMap = boutonData.meanGridResponse;
        case 'svd'
            rawMap = computeSVDMap(tempStack, timeVector, p.RespWin);
        otherwise
            error('plotRFHeatmapWithTraces:badMapMethod', ...
                'MapMethod must be ''mean'' or ''svd'', got ''%s''.', p.MapMethod);
    end

    if p.Smooth
        displayMap = smoothRFMapForDisplay(rawMap, p.ScaleFactor, p.Sigma);
        % imagesc only auto-extends half a cell beyond the outer edge when
        % numel(x) matches size(C,2) (i.e. unsmoothed/raw grid). Once the
        % map is upsampled, size(C) no longer matches uAz/uEl_plot, so we
        % must give imagesc the true padded edges explicitly or it will
        % only span from uAz(1) to uAz(end) (leaving a blank margin).
        azStep = uAz(2) - uAz(1);
        elStep = uEl_plot(2) - uEl_plot(1);
        xImg = [uAz(1) - azStep/2, uAz(end) + azStep/2];
        yImg = [uEl_plot(1) - elStep/2, uEl_plot(end) + elStep/2];
    else
        displayMap = rawMap; % raw grid: each pixel = exactly one uniform Az/El cell
        xImg = uAz;
        yImg = uEl_plot;
    end

    axes(ax); 
    imagesc(ax, xImg, yImg, displayMap);
    hold(ax, 'on');
    colormap(ax, p.Colormap);
    set(ax, 'YDir', 'normal', 'CLim', [min(displayMap(:)), max(displayMap(:)) + 1e-6]);

    % (yellow "responsive" box removed -- background + trace already agree)

    % trace overlay geometry: each trace spans ~85% of a grid cell width,
    % and ~40% of a cell height, centered on that cell's coordinate
    vScale = dEl * 0.4;
    hScale = dAz * 0.85;
    tNorm  = (timeVector - min(timeVector)) / (max(timeVector) - min(timeVector));

    for r = 1:nEl
        for c = 1:nAz
            stimTrace = tempStack(:, r, c);

            sX      = (tNorm - 0.5) * hScale + uAz(c);
            sY_stim = (stimTrace  / peakAmp * vScale) + uEl_plot(r);

            % blank/baseline trace overlay removed for a cleaner look; to
            % bring it back, uncompute sY_blank and re-add the dotted plot:
            %   sY_blank = (blankTrace / peakAmp * vScale) + uEl_plot(r);
            %   plot(ax, sX, sY_blank, 'Color', [1 1 1 0.85], 'LineWidth', 0.8, 'LineStyle', ':');
            plot(ax, sX, sY_stim, 'Color', [1 0.3 0.2], 'LineWidth', 1.1);
        end
    end

    xlim(ax, [min(uAz)-dAz/2, max(uAz)+dAz/2]);
    ylim(ax, [min(uEl_plot)-dEl/2, max(uEl_plot)+dEl/2]);
    xlabel(ax, 'Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel(ax, 'Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
    % ticks placed exactly at grid positions -- read the label directly
    % under/beside a trace to confirm which Az/El it's centered on
    set(ax, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
    set(ax, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
    axis(ax, 'square');

    cb = colorbar(ax);
    cb.Label.String = '\Delta F/F';
    cb.FontSize = 7;
end

% function plotRFHeatmapWithTraces(ax, boutonData, uAz, uEl_plot, timeVector, varargin)
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
%   'Colormap'    - colormap name or Nx3 matrix. Default: 'bone'.
%   'MapMethod'   - 'mean' (use meanGridResponse directly, default), 'svd'
%                   (denoise via computeSVDMap), or 'gaussianfit' (render
%                   the fitted 2D Gaussian model from RFMapping_Gaussian2DFit.m
%                   as a smooth analytic surface -- requires boutonData to
%                   have gaussFit_A, gaussFit_Az0, gaussFit_El0,
%                   gaussFit_sigmaAz, gaussFit_sigmaEl fields set).
%   'RespWin'     - [start end] response window (s), passed to computeSVDMap
%                   if MapMethod is 'svd'.
%   'Smooth'      - true/false. If true (default), background is upsampled
%                   + Gaussian-blurred for a continuous-looking display. If
%                   false, background is plotted as the RAW grid — each
%                   pixel is exactly one uniform-size Az/El cell (dAz x dEl),
%                   with no interpolation, so cell boundaries are visible
%                   and you can directly confirm each trace is centered in
%                   its own cell. Recommended when checking alignment.
%                   (Ignored when MapMethod is 'gaussianfit' -- that surface
%                   is always evaluated smoothly/analytically, see GaussRes.)
%   'ScaleFactor' - upsampling factor for display smoothing (only used if
%                   Smooth=true). Default: 10.
%   'Sigma'       - Gaussian blur sigma (only used if Smooth=true). Default: 3.
%   'GaussRes'    - resolution (points per side) for the analytic Gaussian
%                   surface when MapMethod='gaussianfit'. Default: 200.
%   'ShowFitCenter' - true/false, only used when MapMethod='gaussianfit'.
%                   If true (default), overlays the fitted center (x0,y0)
%                   as a small marker and the 1-sigma ellipse, so you can
%                   directly see whether the traces line up with where the
%                   fit says the true peak is.
% 
% USAGE:
%   plotRFHeatmapWithTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector);
%   plotRFHeatmapWithTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector, ...
%       'Smooth', false);   % raw uniform grid cells, for alignment checking
%   plotRFHeatmapWithTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector, ...
%       'MapMethod', 'gaussianfit');   % smooth fitted-model background instead of raw/blurred data
% 
% SGS code for SparseNoise analysis applied here - specifially the
% smoothing and SVD 
% 
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
%     check that the grid is evenly spaced; dAz/dEl below assume this ---
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
% 
%     if isGaussFitMap
%         Evaluate the fitted 2D Gaussian directly on a fine mesh -- this is an ANALYTIC
%         surface, not a blurred version of discrete data, so it doesn't go through
%         smoothRFMapForDisplay at all (Smooth/ScaleFactor/Sigma are ignored here).
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
%         imagesc only auto-extends half a cell beyond the outer edge when
%         numel(x) matches size(C,2) (i.e. unsmoothed/raw grid). Once the
%         map is upsampled, size(C) no longer matches uAz/uEl_plot, so we
%         must give imagesc the true padded edges explicitly or it will
%         only span from uAz(1) to uAz(end) (leaving a blank margin).
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
% 
%     (yellow "responsive" box removed -- background + trace already agree)
% 
%     trace overlay geometry: each trace spans ~85% of a grid cell width,
%     and ~40% of a cell height, centered on that cell's coordinate
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
% 
%             blank/baseline trace overlay removed for a cleaner look; to
%             bring it back, uncompute sY_blank and re-add the dotted plot:
%               sY_blank = (blankTrace / peakAmp * vScale) + uEl_plot(r);
%               plot(ax, sX, sY_blank, 'Color', [1 1 1 0.85], 'LineWidth', 0.8, 'LineStyle', ':');
%             plot(ax, sX, sY_stim, 'Color', [1 0.3 0.2], 'LineWidth', 1.1);
%         end
%     end
% 
%     if isGaussFitMap && p.ShowFitCenter
%         overlay the fitted center + 1-sigma ellipse so you can see directly whether the
%         traces line up with where the model says the true peak/extent actually is
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
%     ticks placed exactly at grid positions -- read the label directly
%     under/beside a trace to confirm which Az/El it's centered on
%     set(ax, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
%     set(ax, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
%     axis(ax, 'square');
% 
%     cb = colorbar(ax);
%     cb.Label.String = '\Delta F/F';
%     cb.FontSize = 7;
% end