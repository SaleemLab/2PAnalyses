% RFMapping_RawSliceInspection.m
%
% For a set of boutons, plots:
%   (1) the RAW (non-interpolated) grid response as a heatmap, with the
%       fitted center marked and dashed lines showing which row/column
%       the 1D slices below are taken from
%   (2) a 1D azimuth slice: raw response at each tested azimuth, for the
%       grid ROW closest to the fitted y0 (elevation center), with the
%       fitted Gaussian's azimuth profile overlaid as a curve
%   (3) a 1D elevation slice: same idea, along the grid COLUMN closest to
%       the fitted x0
%
% This directly shows the real data spread next to the fit, rather than
% relying on visual impression from an interpolated/smoothed heatmap.
%
% REQUIRES: allRFMapping, uAz, uEl_plot, gaussFitResults (from
% RFMapping_Gaussian2DFit.m -- run that first and keep it in the workspace).

if ~exist('gaussFitResults', 'var')
    error('gaussFitResults not found -- run RFMapping_Gaussian2DFit.m first and keep it in the workspace.');
end

%% params -- default to the same 8 pileup boutons from the multi-start check;
%% override boutonIDs with your own list of iROI values if you want different examples
boutonIDs = [5518, 5446, 1610, 2252, 3011, 1995, 5120, 2236];

azStep = mean(diff(sort(uAz)));
elStep = mean(diff(sort(uEl_plot)));
gaussFit1D_az = @(p, az) p(1) * exp( -(az - p(2)).^2 ./ (2*p(4)^2) ); % p = [A, x0, y0, sigmaX, sigmaY], evaluated at fixed el (via A scaling below)
gaussFit1D_el = @(p, el) p(1) * exp( -(el - p(3)).^2 ./ (2*p(5)^2) );

nBoutons = numel(boutonIDs);
figure('Color', 'w', 'Position', [50 50 1000 260*nBoutons], 'Name', 'Raw slice inspection');
tl = tiledlayout(nBoutons, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:nBoutons
    iROI = boutonIDs(k);
    gi = find([gaussFitResults.iROI] == iROI, 1);
    if isempty(gi)
        warning('Bouton %d not found in gaussFitResults -- skipping.', iROI);
        continue;
    end
    fitR = gaussFitResults(gi);
    b = allRFMapping(iROI);
    resp = double(b.meanGridResponse); % [nEl x nAz]

    % find the grid row/col closest to the fitted center
    [~, elIdxAtCenter] = min(abs(uEl_plot - fitR.y0));
    [~, azIdxAtCenter] = min(abs(uAz - fitR.x0));

    %% (1) raw heatmap
    nexttile;
    imagesc(uAz, uEl_plot, resp);
    set(gca, 'YDir', 'normal');
    colormap(gca, gray);
    hold on;
    yline(uEl_plot(elIdxAtCenter), 'c--', 'LineWidth', 1);
    xline(uAz(azIdxAtCenter), 'm--', 'LineWidth', 1);
    plot(fitR.x0, fitR.y0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
    xlabel('Azimuth (\circ)'); ylabel('Elevation (\circ)');
    title(sprintf('Bouton %d | R^2=%.2f | \\sigma_{az}=%.1f \\sigma_{el}=%.1f', ...
        iROI, fitR.R2, fitR.sigmaX, fitR.sigmaY), 'FontSize', 8);
    hold off;

    %% (2) azimuth slice at the row closest to fitted y0 (cyan dashed line above)
    nexttile;
    azSliceResp = resp(elIdxAtCenter, :);
    plot(uAz, azSliceResp, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    hold on;
    azFine = linspace(min(uAz)-azStep, max(uAz)+azStep, 200);
    % overlay fitted Gaussian's azimuth profile, scaled by the fit's predicted value at this
    % elevation slice (i.e. A * exp(-(el_slice - y0)^2 / (2*sigmaEl^2)) as the effective peak)
    elOffsetScale = exp( -(uEl_plot(elIdxAtCenter) - fitR.y0)^2 / (2*fitR.sigmaY^2) );
    azFitCurve = fitR.A * elOffsetScale * exp( -(azFine - fitR.x0).^2 / (2*fitR.sigmaX^2) );
    plot(azFine, azFitCurve, 'r-', 'LineWidth', 1.5);
    xline(fitR.x0, 'r:');
    xlabel('Azimuth (\circ)'); ylabel('\Delta F/F');
    title(sprintf('Azimuth slice @ el=%.0f\\circ (cyan line)', uEl_plot(elIdxAtCenter)), 'FontSize', 8);
    legend({'raw data', 'fitted Gaussian'}, 'FontSize', 6, 'Location', 'best');
    hold off;

    %% (3) elevation slice at the column closest to fitted x0 (magenta dashed line above)
    nexttile;
    elSliceResp = resp(:, azIdxAtCenter);
    plot(uEl_plot, elSliceResp, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    hold on;
    elFine = linspace(min(uEl_plot)-elStep, max(uEl_plot)+elStep, 200);
    azOffsetScale = exp( -(uAz(azIdxAtCenter) - fitR.x0)^2 / (2*fitR.sigmaX^2) );
    elFitCurve = fitR.A * azOffsetScale * exp( -(elFine - fitR.y0).^2 / (2*fitR.sigmaY^2) );
    plot(elFine, elFitCurve, 'r-', 'LineWidth', 1.5);
    xline(fitR.y0, 'r:');
    xlabel('Elevation (\circ)'); ylabel('\Delta F/F');
    title(sprintf('Elevation slice @ az=%.0f\\circ (magenta line)', uAz(azIdxAtCenter)), 'FontSize', 8);
    legend({'raw data', 'fitted Gaussian'}, 'FontSize', 6, 'Location', 'best');
    hold off;
end

title(tl, 'Raw grid data vs. fitted Gaussian: does the response actually spread to neighboring grid points, or drop off after one step?', ...
    'FontWeight', 'normal', 'FontSize', 10);
