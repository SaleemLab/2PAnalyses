function plotRFHeatmapSelectTraces(ax, boutonData, uAz, uEl_plot, timeVector, varargin)
% plotRFHeatmapSelectTraces: Full spatial heatmap (as in
% plotRFHeatmapWithTraces), but overlays traces only at the peak position
% and its 4-connected neighbors (up/down/left/right), reducing clutter
% while still showing real temporal data at the key spots.
%
% Falls back to the grid's amplitude-max position if the bouton is not
% responsive (no defined peak) — flagged as arbitrary in that case, since
% there's no principled "peak" for a non-responsive bouton.
%
% INPUTS:
%   ax          - target axes handle to plot into
%   boutonData  - single element of an RFMapping struct array (from
%                 analyseRFMapping.m), must contain:
%                   .meanGridResponse (nEl x nAz)
%                   .meanTemporalResponse (nTpts x nEl x nAz)
%                   .meanBlankResponse (nTpts x 1)
%                   .peakAmplitude
%                   .centerAz, .centerEl, .isResponsive
%   uAz         - vector of unique azimuth positions
%   uEl_plot    - vector of unique elevation positions
%   timeVector  - shared time vector (nTpts x 1)
%
% NAME-VALUE OPTIONS (same interface as plotRFHeatmapOnly):
%   'Colormap'    - colormap name or Nx3 matrix. Default: 'bone'.
%   'MapMethod'   - 'mean' (use meanGridResponse directly, default) or
%                   'svd' (denoise via computeSVDMap, using this
%                   function's timeVector automatically).
%   'RespWin'     - [start end] response window (s), passed to computeSVDMap
%                   if MapMethod is 'svd'. Optional even then.
%   'ScaleFactor' - upsampling factor for display smoothing. Default: 10.
%   'Sigma'       - Gaussian blur sigma (upsampled pixels). Default: 3.
%
% USAGE:
%   plotRFHeatmapSelectTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector);
%   plotRFHeatmapSelectTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector, 'Colormap', 'parula');
%   plotRFHeatmapSelectTraces(ax, allRFMapping(iBouton), uAz, uEl_plot, timeVector, ...
%       'MapMethod', 'svd', 'RespWin', respWin, 'Colormap', 'parula');

    p = struct('Colormap', 'bone', 'MapMethod', 'mean', 'RespWin', [], ...
                'ScaleFactor', 10, 'Sigma', 3);
    for k = 1:2:numel(varargin)
        p.(varargin{k}) = varargin{k+1};
    end

    nAz = length(uAz);
    nEl = length(uEl_plot);

    dAz = abs(uAz(2) - uAz(1));
    dEl = abs(uEl_plot(1) - uEl_plot(2));

    grid       = boutonData.meanGridResponse;  % used for peak-finding regardless of MapMethod
    tempStack  = boutonData.meanTemporalResponse;
    blankTrace = boutonData.meanBlankResponse;
    peakAmp    = max(boutonData.peakAmplitude, eps);

    switch lower(p.MapMethod)
        case 'mean'
            rawMap = grid;
        case 'svd'
            rawMap = computeSVDMap(tempStack, timeVector, p.RespWin);
        otherwise
            error('plotRFHeatmapSelectTraces:badMapMethod', ...
                'MapMethod must be ''mean'' or ''svd'', got ''%s''.', p.MapMethod);
    end

    smoothedMap = smoothRFMapForDisplay(rawMap, p.ScaleFactor, p.Sigma);

    axes(ax); %#ok<LAXES>
    imagesc(ax, uAz, uEl_plot, smoothedMap);
    hold(ax, 'on');
    colormap(ax, p.Colormap);
    set(ax, 'YDir', 'normal', 'CLim', [min(smoothedMap(:)), max(smoothedMap(:)) + 1e-6]);

    % determine peak position (use bouton's own peak if responsive, otherwise
    % just the grid's amplitude-max position so there's still something to show)
    if boutonData.isResponsive
        rPeak = find(uEl_plot == boutonData.centerEl, 1);
        cPeak = find(uAz == boutonData.centerAz, 1);
        rectangle(ax, 'Position', [boutonData.centerAz - dAz/2, boutonData.centerEl - dEl/2, dAz, dEl], ...
            'EdgeColor', [1 0.85 0], 'LineWidth', 2.5);
    else
        [~, mI] = max(grid(:), [], 'omitnan');
        [rPeak, cPeak] = ind2sub(size(grid), mI);
    end

    % peak + 4-connected neighbors (clipped to grid edges)
    selPositions = [rPeak, cPeak];
    neighborOffsets = [-1 0; 1 0; 0 -1; 0 1];
    for k = 1:size(neighborOffsets, 1)
        rN = rPeak + neighborOffsets(k, 1);
        cN = cPeak + neighborOffsets(k, 2);
        if rN >= 1 && rN <= nEl && cN >= 1 && cN <= nAz
            selPositions = [selPositions; rN, cN]; %#ok<AGROW>
        end
    end

    vScale = dEl * 0.4;
    hScale = dAz * 0.85;
    tNorm  = (timeVector - min(timeVector)) / (max(timeVector) - min(timeVector));

    for k = 1:size(selPositions, 1)
        r = selPositions(k, 1);
        c = selPositions(k, 2);
        stimTrace = tempStack(:, r, c);

        sX       = (tNorm - 0.5) * hScale + uAz(c);
        sY_stim  = (stimTrace  / peakAmp * vScale) + uEl_plot(r);
        sY_blank = (blankTrace / peakAmp * vScale) + uEl_plot(r);

        plot(ax, sX, sY_blank, 'Color', [1 1 1 0.85], 'LineWidth', 0.8, 'LineStyle', ':');
        plot(ax, sX, sY_stim,  'Color', [1 0.3 0.2], 'LineWidth', 1.3);
    end

    xlim(ax, [min(uAz)-dAz/2, max(uAz)+dAz/2]);
    ylim(ax, [min(uEl_plot)-dEl/2, max(uEl_plot)+dEl/2]);
    xlabel(ax, 'Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel(ax, 'Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
    set(ax, 'FontName', 'Arial', 'FontSize', 8, 'Box', 'off', 'TickDir', 'out');
    axis(ax, 'square');

    cb = colorbar(ax);
    cb.Label.String = 'dF/F';
    cb.FontSize = 7;
end
