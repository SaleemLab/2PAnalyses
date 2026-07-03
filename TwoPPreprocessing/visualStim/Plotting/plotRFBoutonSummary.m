function plotRFBoutonSummary(fig, boutonData, uAz, uEl_plot, timeVector, varargin)
% plotRFBoutonSummary: Two-panel summary for a single bouton — left panel
% shows the actual per-position traces (ground truth), right panel shows a
% smoothed spatial map (mean or SVD-denoised) for a clean at-a-glance view
% of spatial selectivity. Combines plotRFHeatmapWithTraces (or
% plotRFHeatmapSelectTraces) with plotRFHeatmapOnly.
%
% INPUTS:
%   fig         - target figure handle to draw both panels into
%   boutonData  - single element of an RFMapping struct array (from
%                 analyseRFMapping.m)
%   uAz         - vector of unique azimuth positions
%   uEl_plot    - vector of unique elevation positions
%   timeVector  - shared time vector (nTpts x 1)
%
% NAME-VALUE OPTIONS:
%   'TracesMode'  - 'all' (plot every position's trace, via
%                   plotRFHeatmapWithTraces, default) or 'select' (only
%                   peak + neighbors, via plotRFHeatmapSelectTraces).
%   'MapMethod'   - 'mean' (default) or 'svd', for the RIGHT (smoothed) panel.
%                   The LEFT (traces) panel always shows the raw
%                   meanGridResponse as its background heatmap, since it's
%                   meant to be the "ground truth" reference.
%   'Colormap'    - colormap for both panels. Default: 'bone'.
%   'RespWin'     - [start end] response window (s); required if MapMethod
%                   is 'svd', unused otherwise.
%   'ScaleFactor' - upsampling factor for the smoothed panel. Default: 10.
%   'Sigma'       - Gaussian blur sigma for the smoothed panel. Default: 3.
%   'TitlePrefix' - string prepended to both panel titles, e.g. bouton ID.
%
% USAGE:
%   fig = figure('Color', 'w', 'Position', [50 50 1400 600]);
%   plotRFBoutonSummary(fig, allRFMapping(iBouton), uAz, uEl_plot, timeVector, ...
%       'MapMethod', 'svd', 'RespWin', respWin, 'Colormap', 'parula', ...
%       'TitlePrefix', sprintf('Bouton %d', iBouton));

    p = struct('TracesMode', 'all', 'MapMethod', 'mean', 'Colormap', 'bone', ...
                'RespWin', [], 'ScaleFactor', 10, 'Sigma', 3, 'TitlePrefix', '');
    for k = 1:2:numel(varargin)
        p.(varargin{k}) = varargin{k+1};
    end

    axTraces = subplot(1, 2, 1, 'Parent', fig);
    switch lower(p.TracesMode)
        case 'all'
            plotRFHeatmapWithTraces(axTraces, boutonData, uAz, uEl_plot, timeVector, ...
                'Colormap', p.Colormap, 'MapMethod', 'mean');  % traces panel always uses raw mean
        case 'select'
            plotRFHeatmapSelectTraces(axTraces, boutonData, uAz, uEl_plot, timeVector, ...
                'Colormap', p.Colormap, 'MapMethod', 'mean');
        otherwise
            error('plotRFBoutonSummary:badTracesMode', ...
                'TracesMode must be ''all'' or ''select'', got ''%s''.', p.TracesMode);
    end
    title(axTraces, [p.TitlePrefix ' — traces'], 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');

    axMap = subplot(1, 2, 2, 'Parent', fig);
    mapArgs = {'Colormap', p.Colormap, 'MapMethod', p.MapMethod, ...
               'ScaleFactor', p.ScaleFactor, 'Sigma', p.Sigma};
    if strcmpi(p.MapMethod, 'svd')
        mapArgs = [mapArgs, {'TimeVector', timeVector, 'RespWin', p.RespWin}];
    end
    plotRFHeatmapOnly(axMap, boutonData, uAz, uEl_plot, mapArgs{:});
    title(axMap, sprintf('%s — smoothed (%s)', p.TitlePrefix, upper(p.MapMethod)), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
end