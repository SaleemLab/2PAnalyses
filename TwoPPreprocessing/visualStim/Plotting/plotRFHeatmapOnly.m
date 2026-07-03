function plotRFHeatmapOnly(ax, boutonData, uAz, uEl_plot, varargin)
% plotRFHeatmapOnly: Draws just the spatial tuning heatmap for a single
% bouton, with no trace overlays — meant to be placed alongside
% plotRFHeatmapWithTraces so you get one clean, easy-to-read spatial map
% plus one panel with the underlying temporal traces.
%
% INPUTS:
%   ax          - target axes handle to plot into
%   boutonData  - single element of an RFMapping struct array, must contain:
%                   .meanGridResponse (nEl x nAz)
%                   .meanTemporalResponse (nTpts x nEl x nAz) — only needed
%                   if MapMethod is 'svd'
%                   .centerAz, .centerEl, .isResponsive
%   uAz         - vector of unique azimuth positions
%   uEl_plot    - vector of unique elevation positions
%
% NAME-VALUE OPTIONS:
%   'Colormap'    - colormap name or Nx3 matrix. Default: 'bone'.
%   'MapMethod'   - 'mean' (use meanGridResponse directly, default) or
%                   'svd' (denoise via computeSVDMap).
%   'TimeVector'  - required if MapMethod is 'svd'.
%   'RespWin'     - [start end] response window (s), passed to computeSVDMap
%                   if MapMethod is 'svd'.
%   'ScaleFactor' - upsampling factor for display smoothing. Default: 10.
%   'Sigma'       - Gaussian blur sigma. Default: 3.
%
% USAGE:
%   plotRFHeatmapOnly(ax, allRFMapping(iBouton), uAz, uEl_plot);
%   plotRFHeatmapOnly(ax, allRFMapping(iBouton), uAz, uEl_plot, ...
%       'MapMethod', 'svd', 'TimeVector', timeVector, 'RespWin', respWin, ...
%       'Colormap', 'parula');

    p = struct('Colormap', 'bone', 'MapMethod', 'mean', 'TimeVector', [], ...
                'RespWin', [], 'ScaleFactor', 10, 'Sigma', 3);
    for k = 1:2:numel(varargin)
        p.(varargin{k}) = varargin{k+1};
    end

    dAz = abs(uAz(2) - uAz(1));
    dEl = abs(uEl_plot(1) - uEl_plot(2));

    switch lower(p.MapMethod)
        case 'mean'
            rawMap = boutonData.meanGridResponse;
        case 'svd'
            if isempty(p.TimeVector)
                error('plotRFHeatmapOnly:missingTimeVector', ...
                    '''TimeVector'' is required when MapMethod is ''svd''.');
            end
            rawMap = computeSVDMap(boutonData.meanTemporalResponse, p.TimeVector, p.RespWin);
        otherwise
            error('plotRFHeatmapOnly:badMapMethod', ...
                'MapMethod must be ''mean'' or ''svd'', got ''%s''.', p.MapMethod);
    end

    smoothedMap = smoothRFMapForDisplay(rawMap, p.ScaleFactor, p.Sigma);

    % imagesc only auto-extends half a cell beyond the outer edge when
    % numel(x) matches size(C,2) (true for the raw grid). Once the map is
    % upsampled here, size(C) no longer matches uAz/uEl_plot, so we must
    % give imagesc the true padded edges explicitly, or it will only span
    % from uAz(1) to uAz(end) and leave a blank margin around the image.
    azStep = uAz(2) - uAz(1);
    elStep = uEl_plot(2) - uEl_plot(1);
    xImg = [uAz(1) - azStep/2, uAz(end) + azStep/2];
    yImg = [uEl_plot(1) - elStep/2, uEl_plot(end) + elStep/2];

    axes(ax);
    imagesc(ax, xImg, yImg, smoothedMap);
    hold(ax, 'on');
    colormap(ax, p.Colormap);
    set(ax, 'YDir', 'normal', 'CLim', [min(smoothedMap(:)), max(smoothedMap(:)) + 1e-6]);

    if boutonData.isResponsive
        rectangle(ax, 'Position', [boutonData.centerAz - dAz/2, boutonData.centerEl - dEl/2, dAz, dEl], ...
            'EdgeColor', [1 0.85 0], 'LineWidth', 2.5);
    end

    xlim(ax, [min(uAz)-dAz/2, max(uAz)+dAz/2]);
    ylim(ax, [min(uEl_plot)-dEl/2, max(uEl_plot)+dEl/2]);
    xlabel(ax, 'Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel(ax, 'Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
    set(ax, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
    set(ax, 'FontName', 'Arial', 'FontSize', 8, 'Box', 'off', 'TickDir', 'out');
    axis(ax, 'square');

    cb = colorbar(ax);
    cb.Label.String = 'dF/F';
    cb.FontSize = 7;
end