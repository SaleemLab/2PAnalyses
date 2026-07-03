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