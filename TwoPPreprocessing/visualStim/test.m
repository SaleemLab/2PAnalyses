% plotRFSummary.m
%
% Visualizes the trusted set from RFMapping_Gaussian2DFit.m (gaussFitResults with
% isTrusted computed). Run this AFTER the main fitting script's COMMIT block.
%
% Produces:
%   Figure 1: sigmaAz / sigmaEl distributions (histograms) with mean +/- SD marked,
%             plus median/IQR printed, since these distributions are typically
%             right-skewed (see earlier discussion).
%   Figure 2: preferred center positions (x0, y0) scatter across the grid, plus
%             marginal histograms of preferred azimuth and elevation separately.

%%
trustR2Threshold = 0.7;   % or whatever you want to test
ciLevel = 0.00;
lowerPct = 100 * (1 - ciLevel);

for gi = 1:numel(gaussFitResults)
    iROI = gaussFitResults(gi).iROI;
    validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));
    if numel(validBoot) >= 20
        lowerBound = prctile(validBoot, lowerPct);
    else
        lowerBound = NaN;
    end
    isRobust  = ~isnan(lowerBound) && lowerBound >= trustR2Threshold;
    isTrusted = (gaussFitResults(gi).R2 >= trustR2Threshold) && isRobust && ~gaussFitResults(gi).isDegenerateSigma;

    gaussFitResults(gi).bootR2LowerCI = lowerBound;
    gaussFitResults(gi).isRobust      = isRobust;
    gaussFitResults(gi).isTrusted     = isTrusted;

    allRFMapping(iROI).gaussFit_bootR2LowerCI = lowerBound;
    allRFMapping(iROI).gaussFit_isRobust      = isRobust;
    allRFMapping(iROI).gaussFit_isTrusted     = isTrusted;
end

isTrustedAll = [gaussFitResults.isTrusted];
fprintf('trustR2Threshold=%.2f -> %d / %d trusted\n', trustR2Threshold, sum(isTrustedAll), numel(gaussFitResults));

%% gather trusted set
trustR2Threshold = 0.1;
trustedIdx = [gaussFitResults.isTrusted];
if ~any(trustedIdx)
    error('No boutons pass isTrusted -- check trustR2Threshold/ciLevel in the commit block.');
end

azVals = [gaussFitResults(trustedIdx).sigmaX];
elVals = [gaussFitResults(trustedIdx).sigmaY];
x0Vals = [gaussFitResults(trustedIdx).x0];
y0Vals = [gaussFitResults(trustedIdx).y0];
n = sum(trustedIdx);

%% summary stats
mAz = mean(azVals); sdAz = std(azVals); semAz = sdAz/sqrt(n);
mEl = mean(elVals); sdEl = std(elVals); semEl = sdEl/sqrt(n);
medAz = median(azVals); iqrAz = prctile(azVals,[25 75]);
medEl = median(elVals); iqrEl = prctile(elVals,[25 75]);

fprintf('\n=== RF width summary (trusted set, n=%d) ===\n', n);
fprintf('sigAz = %.2f +/- %.2f (SD), SEM = %.2f | median = %.2f, IQR = [%.2f, %.2f]\n', ...
    mAz, sdAz, semAz, medAz, iqrAz(1), iqrAz(2));
fprintf('sigEl = %.2f +/- %.2f (SD), SEM = %.2f | median = %.2f, IQR = [%.2f, %.2f]\n', ...
    mEl, sdEl, semEl, medEl, iqrEl(1), iqrEl(2));

%% ===================== Figure 1: width distributions =====================
figure('Name', 'RF width distributions (trusted set)', 'Color', 'w');

subplot(1,2,1);
histogram(azVals, 'FaceColor', [0.8 0.4 0.3], 'EdgeColor', 'none');
hold on;
xline(mAz, '-k', 'LineWidth', 1.5);
xline(mAz - sdAz, '--k'); xline(mAz + sdAz, '--k');
xline(medAz, '-b', 'LineWidth', 1.2);
title(sprintf('sigma Azimuth (n=%d)\nmean=%.1f\\pm%.1f (SD), median=%.1f', n, mAz, sdAz, medAz));
xlabel('sigma Az (deg)'); ylabel('count');
legend({'', 'mean', 'mean\pm SD', '', 'median'}, 'Location', 'best');
hold off;

subplot(1,2,2);
histogram(elVals, 'FaceColor', [0.3 0.5 0.7], 'EdgeColor', 'none');
hold on;
xline(mEl, '-k', 'LineWidth', 1.5);
xline(mEl - sdEl, '--k'); xline(mEl + sdEl, '--k');
xline(medEl, '-b', 'LineWidth', 1.2);
title(sprintf('sigma Elevation (n=%d)\nmean=%.1f\\pm%.1f (SD), median=%.1f', n, mEl, sdEl, medEl));
xlabel('sigma El (deg)'); ylabel('count');
legend({'', 'mean', 'mean\pm SD', '', 'median'}, 'Location', 'best');
hold off;

%% ===================== Figure 2: preferred center positions =====================
figure('Name', 'Preferred RF center positions (trusted set)', 'Color', 'w');

% main scatter, with grid lines at actual tested positions for reference
subplot(2,2,[1 3]);
scatter(x0Vals, y0Vals, 25, [0.2 0.4 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
for a = uAz(:)', xline(a, 'Color', [0.85 0.85 0.85]); end
for e = uEl_plot(:)', yline(e, 'Color', [0.85 0.85 0.85]); end
xlabel('Preferred Azimuth (deg)'); ylabel('Preferred Elevation (deg)');
title(sprintf('Preferred center positions (n=%d)', n));
axis equal; box on;
hold off;

% marginal az histogram
subplot(2,2,2);
histogram(x0Vals, numel(uAz)*2, 'FaceColor', [0.2 0.4 0.7], 'EdgeColor', 'none');
hold on;
for a = uAz(:)', xline(a, 'Color', [0.85 0.85 0.85]); end
xlabel('Preferred Azimuth (deg)'); ylabel('count');
title('Marginal: preferred azimuth');
hold off;

% marginal el histogram
subplot(2,2,4);
histogram(y0Vals, numel(uEl_plot)*2, 'FaceColor', [0.2 0.4 0.7], 'EdgeColor', 'none');
hold on;
for e = uEl_plot(:)', xline(e, 'Color', [0.85 0.85 0.85]); end
xlabel('Preferred Elevation (deg)'); ylabel('count');
title('Marginal: preferred elevation');
hold off;


%% --------------------------------------
%% ===================== DIAGNOSTIC VISUALIZATION (MATCHED AXES) =====================
% Pick a few interesting ROIs to inspect
numROIsToPlot = 3;
passedIdx = find([gaussFitResults.isTrusted]);
failedIdx = find(~[gaussFitResults.isTrusted]);

% Pick a mix of indices from the gaussFitResults array itself
giToPlot = [];
if ~isempty(passedIdx), giToPlot = [giToPlot, passedIdx(randi(numel(passedIdx), 1, min(2, numel(passedIdx))))]; end
if ~isempty(failedIdx), giToPlot = [giToPlot, failedIdx(randi(numel(failedIdx), 1, min(1, numel(failedIdx))))]; end

% Define matching limits for both plots based on raw grid data
azLims = [min(uAz) - azStep, max(uAz) + azStep];
elLims = [min(uEl_plot) - elStep, max(uEl_plot) + elStep];

for gi = giToPlot
    res = gaussFitResults(gi); % Get the specific result struct
    b = allRFMapping(res.iROI); % Match using the stored original ROI index
    
    % Reconstruct the response grid INLINE
    trialMatrix = b.baselineSubtracted;
    respVec = nan(numel(trialMatrix), 1);
    for posIdx = 1:numel(trialMatrix)
        trials = trialMatrix{posIdx};
        if ~isempty(trials)
            trials = double(trials);
            respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
        end
    end
    respGrid = reshape(respVec, size(AzGrid));
    
    % Generate high-resolution model fit grid matching the exact limits to avoid white space
    [xq, yq] = meshgrid(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100));
    pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
    fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));
    
    figure('Position', [100, 100, 1100, 350]);
    sgtitle(sprintf('ROI %d (Raw R^2 = %.2f | Boot Lower Bound = %.2f | Trusted = %d)', ...
        res.iROI, res.R2, res.bootR2LowerCI, res.isTrusted), 'FontWeight', 'bold');
    
    % Plot 1: Actual Response Data Grid
    subplot(1, 3, 1);
    imagesc(uAz, uEl_plot, respGrid);
    axis image; colormap(gca, 'parula'); colorbar;
    xlim(azLims); ylim(elLims); % Set matched limits
    title('Mean Grid Response (\DeltaF/F)');
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    
    % Plot 2: Fitted 2D Gaussian Model
    subplot(1, 3, 2);
    imagesc(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100), fittedGrid);
    axis image; colormap(gca, 'parula'); colorbar;
    xlim(azLims); ylim(elLims); % Set matched limits
    hold on;
    plot(res.x0, res.y0, 'r+', 'MarkerSize', 12, 'LineWidth', 2); % Fitted Center
    title('Fitted 2D Gaussian Model');
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    
    % Plot 3: Bootstrap R^2 Distribution Histogram
    subplot(1, 3, 3);
    validBoot = res.bootR2(~isnan(res.bootR2));
    histogram(validBoot, linspace(-0.2, 1, 40), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    hold on;
    
    % Plot threshold and calculation lines
    yl = ylim;
    line([res.R2, res.R2], yl, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-'); % Raw R2
    line([res.bootR2LowerCI, res.bootR2LowerCI], yl, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--'); % 5th percentile lower bound
    line([0.2, 0.2], yl, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':'); % Threshold line (e.g. 0.2)
    
    legend('Bootstrap R^2', 'Raw R^2', '5th Pct Lower Bound', 'Threshold (0.2)', 'Location', 'northwest');
    title('1,000 Bootstrap R^2 Distribution');
    xlabel('R^2'); ylabel('Count');
    grid on;
end


%% r2 check 
%% ===================== ELBOW & SENSITIVITY VISUALIZATION =====================
rawR2Vals = [gaussFitResults.R2];
bootCIVals = [gaussFitResults.bootR2LowerCI];
bootCIVals = bootCIVals(~isnan(bootCIVals)); % Drop any NaN bounds for the plot

figure('Position', [150, 150, 1000, 450]);

% Plot 1: Population Histograms with Threshold Lines
subplot(1, 2, 1);
hold on;
h1 = histogram(rawR2Vals, linspace(-0.2, 1, 30), 'FaceColor', [0.2 0.6 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = histogram(bootCIVals, linspace(-0.2, 1, 30), 'FaceColor', [0.8 0.4 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Plot your chosen thresholds
yl = ylim;
line([0.2 0.2], yl, 'Color', [0.1 0.4 0.6], 'LineWidth', 2.5, 'LineStyle', '--');
line([0.1 0.1], yl, 'Color', [0.6 0.2 0.1], 'LineWidth', 2.5, 'LineStyle', ':');

legend('Raw R^2 (Point Estimates)', '95% Boot Lower Bounds', ...
       'Raw Threshold (0.2)', 'Boot CI Threshold (0.1)', 'Location', 'northwest');
title('Population R^2 Distribution');
xlabel('R^2 Value');
ylabel('Number of ROIs');
grid on;

% Plot 2: Cumulative Survival Curves ("The Elbow Plot")
subplot(1, 2, 2);
hold on;

% Sort to create cumulative curves
sortedRaw = sort(rawR2Vals);
yRaw = (numel(sortedRaw):-1:1) / numel(sortedRaw) * 100;

sortedBoot = sort(bootCIVals);
yBoot = (numel(sortedBoot):-1:1) / numel(sortedBoot) * 100;

plot(sortedRaw, yRaw, 'Color', [0.2 0.6 0.8], 'LineWidth', 2.5);
plot(sortedBoot, yBoot, 'Color', [0.8 0.4 0.2], 'LineWidth', 2.5);

% Mark your active selection points
idxRawThresh = find(sortedRaw >= 0.2, 1);
if ~isempty(idxRawThresh)
    plot(0.2, yRaw(idxRawThresh), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerEdgeColor', 'k');
    text(0.22, yRaw(idxRawThresh)+3, sprintf('Raw >= 0.2: Keep %.0f%%', yRaw(idxRawThresh)), 'FontWeight', 'bold');
end

idxBootThresh = find(sortedBoot >= 0.1, 1);
if ~isempty(idxBootThresh)
    plot(0.1, yBoot(idxBootThresh), 's', 'MarkerSize', 10, 'MarkerFaceColor', [0.8 0.4 0.2], 'MarkerEdgeColor', 'k');
    text(0.12, yBoot(idxBootThresh)-4, sprintf('Boot >= 0.1: Keep %.0f%%', yBoot(idxBootThresh)), 'FontWeight', 'bold');
end

title('Cumulative Survival ("Elbow Curve")');
xlabel('Threshold Value');
ylabel('% of ROIs Surviving');
xlim([-0.1, 1]);
ylim([0, 100]);
grid on;

%%
azLims = [min(uAz) - azStep/2, max(uAz) + azStep/2];
elLims = [min(uEl_plot) - elStep/2, max(uEl_plot) + elStep/2];

numFitted = numel(gaussFitResults);
hFig = figure('Position', [50, 100, 1250, 400], 'Name', 'RF Comprehensive Diagnostic Viewer', 'NumberTitle', 'off');

for gi = 1:numFitted
    if ~ishandle(hFig)
        break;
    end
    
    res = gaussFitResults(gi);
    b = allRFMapping(res.iROI);
    
    % Reconstruct raw response grid
    trialMatrix = b.baselineSubtracted;
    respVec = nan(numel(trialMatrix), 1);
    for posIdx = 1:numel(trialMatrix)
        trials = trialMatrix{posIdx};
        if ~isempty(trials)
            trials = double(trials);
            respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
        end
    end
    respGrid = reshape(respVec, size(AzGrid));
    
    % Generate high-resolution model fit grid
    [xq, yq] = meshgrid(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100));
    pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
    fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));
    
    % Prepare boutonData struct dynamically for the trace-plotting function
    boutonData = b;
    boutonData.meanGridResponse = respGrid;
    boutonData.meanTemporalResponse = b.meanTemporalResponse; 
    boutonData.meanBlankResponse = b.meanBlankResponse;
    boutonData.peakAmplitude = max(respGrid(:));
    
    clf(hFig);
    
    sgtitle(sprintf('ROI %d (%d of %d) | Raw R^2 = %.2f | Boot Lower Bound = %.2f | Trusted = %d', ...
        res.iROI, gi, numFitted, res.R2, res.bootR2LowerCI, res.isTrusted), 'FontWeight', 'bold');
    
    % --- Subplot 1: Trace-Overlaid Spatial Heatmap ---
    ax1 = subplot(1, 3, 1);
    
    % Call your custom function on the first axes
    plotRFHeatmapWithTraces(ax1, boutonData, uAz, uEl_plot, timeVector, ...
        'Colormap', 'bone', 'Smooth', false);
    
    % Force matching coordinates and flip vertical axis direction to match standard
    xlim(ax1, azLims); 
    ylim(ax1, elLims);
    set(ax1, 'YDir', 'normal');
    title(ax1, 'Mean Grid Response & Traces');
    
    % --- Subplot 2: Fitted 2D Gaussian Model ---
    ax2 = subplot(1, 3, 2);
    imagesc(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100), fittedGrid);
    axis image; 
    colormap(ax2, 'bone'); 
    colorbar;
    xlim(ax2, azLims); 
    ylim(ax2, elLims);
    set(ax2, 'YDir', 'normal');
    hold on;
    
    % Overlay fitted center and 1-sigma boundary
    plot(ax2, res.x0, res.y0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    t_ellipse = linspace(0, 2*pi, 100);
    plot(ax2, res.x0 + res.sigmaX * cos(t_ellipse), res.y0 + res.sigmaY * sin(t_ellipse), 'r--', 'LineWidth', 1.5);
    
    title(ax2, 'Fitted 2D Gaussian Model');
    xlabel(ax2, 'Azimuth (deg)'); 
    ylabel(ax2, 'Elevation (deg)');
    
    % --- Subplot 3: Bootstrap R^2 Distribution ---
    ax3 = subplot(1, 3, 3);
    validBoot = res.bootR2(~isnan(res.bootR2));
    histogram(ax3, validBoot, linspace(-0.2, 1, 40), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    hold on;
    
    yl = ylim(ax3);
    line([res.R2, res.R2], yl, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-');
    line([res.bootR2LowerCI, res.bootR2LowerCI], yl, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
    line([0.2, 0.2], yl, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':');
    
    legend(ax3, 'Bootstrap R^2', 'Raw R^2', '5th Pct Lower Bound', 'Threshold (0.2)', 'Location', 'northwest');
    title(ax3, '1,000 Bootstrap R^2 Distribution');
    xlabel(ax3, 'R^2'); 
    ylabel(ax3, 'Count');
    grid on;
    
    drawnow;
    
    if gi < numFitted
        waitforbuttonpress;
    end
end

%% count how many categories there are 
%% ===================== COUNT TUNING MORPHOLOGIES =====================
% Define thresholds matching your active gates
rawR2Thresh   = 0.2;
bootCIThresh  = 0.1;

% Initialize category counters
totalFitted   = numel(gaussFitResults);
cntLocalized  = 0;
cntDegenerate = 0;
cntNoise      = 0;

% Lists to store indices for further inspection if needed
idxLocalized  = [];
idxDegenerate = [];

for gi = 1:totalFitted
    res = gaussFitResults(gi);
    
    % Evaluate statistical stability first (Dual-Threshold Gate)
    isRobust = (res.R2 >= rawR2Thresh) && (res.bootR2LowerCI >= bootCIThresh);
    
    if isRobust
        if res.isDegenerateSigma
            cntDegenerate = cntDegenerate + 1;
            idxDegenerate(end+1) = gi; %#ok<SAGROW>
        else
            cntLocalized = cntLocalized + 1;
            idxLocalized(end+1) = gi; %#ok<SAGROW>
        end
    else
        cntNoise = cntNoise + 1;
    end
end

% Calculate percentages relative to the total pool
pctLocalized  = (cntLocalized / totalFitted) * 100;
pctDegenerate = (cntDegenerate / totalFitted) * 100;
pctNoise      = (cntNoise / totalFitted) * 100;

% Display a clean, formatted text table in your Command Window
fprintf('\n=======================================================\n');
fprintf('          POPULATION RECEPTIVE FIELD ANALYSIS          \n');
fprintf('=======================================================\n');
fprintf('Total Fitted ROIs evaluated : %d\n\n', totalFitted);
fprintf('1. Localized RFs (Classic Ovals)     : %d cells (%5.1f%%)\n', cntLocalized, pctLocalized);
fprintf('2. Band/Stripe RFs (Horiz. Stretch)  : %d cells (%5.1f%%)\n', cntDegenerate, pctDegenerate);
fprintf('3. Unstable / Noisy Fits (Discarded) : %d cells (%5.1f%%)\n', cntNoise, pctNoise);
fprintf('-------------------------------------------------------\n');
fprintf('Total Responsive Pool (1 + 2)       : %d cells (%5.1f%%)\n', ...
    (cntLocalized + cntDegenerate), ((cntLocalized + cntDegenerate)/totalFitted)*100);
fprintf('=======================================================\n');


%% plot the horizontal striped ones to see if they are real.. 
rawR2Thresh   = 0.2;
bootCIThresh  = 0.1;

azLims = [min(uAz) - azStep/2, max(uAz) + azStep/2];
elLims = [min(uEl_plot) - elStep/2, max(uEl_plot) + elStep/2];

% Identify the indices of the 54 degenerate/stripe cells
idxDegenerate = [];
for gi = 1:numel(gaussFitResults)
    res = gaussFitResults(gi);
    isRobust = (res.R2 >= rawR2Thresh) && (res.bootR2LowerCI >= bootCIThresh);
    if isRobust && res.isDegenerateSigma
        idxDegenerate(end+1) = gi; %#ok<SAGROW>
    end
end

numStripeCells = numel(idxDegenerate);
fprintf('\n--- Starting viewer for %d Horizontally Stretched/Stripe ROIs ---\n', numStripeCells);
fprintf('Click on the figure window or press ANY KEY to load the next stripe ROI.\n\n');

hFig = figure('Position', [50, 100, 1250, 400], 'Name', 'Stripe/Band RF Diagnostic Viewer', 'NumberTitle', 'off');

for countIdx = 1:numStripeCells
    if ~ishandle(hFig)
        break;
    end
    
    gi = idxDegenerate(countIdx); % Pull out the specific degenerate index
    res = gaussFitResults(gi);
    b = allRFMapping(res.iROI);
    
    % Reconstruct raw response grid
    trialMatrix = b.baselineSubtracted;
    respVec = nan(numel(trialMatrix), 1);
    for posIdx = 1:numel(trialMatrix)
        trials = trialMatrix{posIdx};
        if ~isempty(trials)
            trials = double(trials);
            respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
        end
    end
    respGrid = reshape(respVec, size(AzGrid));
    
    % Generate high-resolution model fit grid
    [xq, yq] = meshgrid(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100));
    pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
    fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));
    
    % Prepare boutonData struct dynamically for the trace-plotting function
    boutonData = b;
    boutonData.meanGridResponse = respGrid;
    boutonData.meanTemporalResponse = b.meanTemporalResponse; 
    boutonData.meanBlankResponse = b.meanBlankResponse;
    boutonData.peakAmplitude = max(respGrid(:));
    
    clf(hFig);
    
    sgtitle(sprintf('STRIPE CELL %d of %d | ROI %d | Raw R^2 = %.2f | Boot Lower Bound = %.2f', ...
        countIdx, numStripeCells, res.iROI, res.R2, res.bootR2LowerCI), 'FontWeight', 'bold');
    
    % --- Subplot 1: Trace-Overlaid Spatial Heatmap ---
    ax1 = subplot(1, 3, 1);
    plotRFHeatmapWithTraces(ax1, boutonData, uAz, uEl_plot, timeVector, ...
        'Colormap', 'bone', 'Smooth', false);
    xlim(ax1, azLims); ylim(ax1, elLims);
    set(ax1, 'YDir', 'normal');
    title(ax1, 'Mean Grid Response & Traces');
    
    % --- Subplot 2: Fitted 2D Gaussian Model ---
    ax2 = subplot(1, 3, 2);
    imagesc(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100), fittedGrid);
    axis image; colormap(ax2, 'bone'); colorbar;
    xlim(ax2, azLims); ylim(ax2, elLims);
    set(ax2, 'YDir', 'normal');
    hold on;
    
    % Overlay fitted center and 1-sigma boundary (will show up very wide/stretched)
    plot(ax2, res.x0, res.y0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    t_ellipse = linspace(0, 2*pi, 100);
    plot(ax2, res.x0 + res.sigmaX * cos(t_ellipse), res.y0 + res.sigmaY * sin(t_ellipse), 'r--', 'LineWidth', 1.5);
    title(ax2, 'Fitted Stretched Gaussian Model');
    xlabel(ax2, 'Azimuth (deg)'); ylabel(ax2, 'Elevation (deg)');
    
    % --- Subplot 3: Bootstrap R^2 Distribution ---
    ax3 = subplot(1, 3, 3);
    validBoot = res.bootR2(~isnan(res.bootR2));
    histogram(ax3, validBoot, linspace(-0.2, 1, 40), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    hold on;
    
    yl = ylim(ax3);
    line([res.R2, res.R2], yl, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-');
    line([res.bootR2LowerCI, res.bootR2LowerCI], yl, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
    line([0.2, 0.2], yl, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':');
    
    legend(ax3, 'Bootstrap R^2', 'Raw R^2', '5th Pct Lower Bound', 'Threshold (0.2)', 'Location', 'northwest');
    title(ax3, '1,000 Bootstrap R^2 Distribution');
    xlabel(ax3, 'R^2'); ylabel(ax3, 'Count');
    grid on;
    
    drawnow;
    
    if countIdx < numStripeCells
        waitforbuttonpress;
    end
end

%% measure population widiths 
rawR2Thresh  = 0.2;
bootCIThresh = 0.1;

pooledSigmaX = [];
pooledSigmaY = [];

for gi = 1:numel(gaussFitResults)
    res = gaussFitResults(gi);
    isRobust = (res.R2 >= rawR2Thresh) && (res.bootR2LowerCI >= bootCIThresh);
    if isRobust
        pooledSigmaX(end+1) = res.sigmaX; 
        pooledSigmaY(end+1) = res.sigmaY; 
    end
end

% Compute statistics for Azimuth Width (Sigma X)
meanX   = mean(pooledSigmaX);
medianX = median(pooledSigmaX);
madX = mad(pooledSigmaX); 
sdX     = std(pooledSigmaX);
iqrX    = iqr(pooledSigmaX);

% Compute statistics for Elevation Width (Sigma Y)
meanY   = mean(pooledSigmaY);
medianY = median(pooledSigmaY);
madY = mad(pooledSigmaY); 

sdY     = std(pooledSigmaY);
iqrY    = iqr(pooledSigmaY);

fprintf('\n=======================================================\n');
fprintf('     RECEPTIVE FIELD WIDTHS (SIGMA) FOR ALL POOLED BOUTONS\n');
fprintf('=======================================================\n');
fprintf('Total Pooled Boutons: %d cells\n\n', numel(pooledSigmaX));
fprintf('--- AZIMUTH WIDTH (\\sigma_x, deg) ---\n');
fprintf('Mean  \\sigma_x : %.2f^\\circ\n', meanX);
fprintf('Median\\sigma_x : %.2f^\\circ\n', medianX);
fprintf('SD    \\sigma_x : %.2f^\\circ\n', sdX);
fprintf('IQR   \\sigma_x : %.2f^\\circ\n\n', iqrX);
fprintf('--- ELEVATION WIDTH (\\sigma_y, deg) ---\n');
fprintf('Mean  \\sigma_y : %.2f^\\circ\n', meanY);
fprintf('Median\\sigma_y : %.2f^\\circ\n', medianY);
fprintf('SD    \\sigma_y : %.2f^\\circ\n', sdY);
fprintf('IQR   \\sigma_y : %.2f^\\circ\n', iqrY);
fprintf('=======================================================\n');

% Calculate the boundaries explicitly
q1_X = prctile(pooledSigmaX, 25);
q3_X = prctile(pooledSigmaX, 75);

q1_Y = prctile(pooledSigmaY, 25);
q3_Y = prctile(pooledSigmaY, 75);

fprintf('\n=======================================================\n');
fprintf('     QUARTILE BOUNDARIES FOR RECEPTIVE FIELD WIDTHS    \n');
fprintf('=======================================================\n');
fprintf('Azimuth (\\sigma_x):\n  [25th Pct (Q1): %.2f^\\circ, 75th Pct (Q3): %.2f^\\circ] (IQR = %.2f^\\circ)\n\n', q1_X, q3_X, q3_X - q1_X);
fprintf('Elevation (\\sigma_y):\n  [25th Pct (Q1): %.2f^\\circ, 75th Pct (Q3): %.2f^\\circ] (IQR = %.2f^\\circ)\n', q1_Y, q3_Y, q3_Y - q1_Y);



%%
% rawR2Thresh = 0.2;
% bootCIThresh = 0.1;
% 
% azLims = [min(uAz) - azStep/2, max(uAz) + azStep/2];
% elLims = [min(uEl_plot) - elStep/2, max(uEl_plot) + elStep/2];
% 
% hFig = figure('Position', [100, 100, 650, 550], 'Name', 'Population RF Coverage Map', 'NumberTitle', 'off');
% ax = axes('Parent', hFig);
% hold(ax, 'on');
% 
% theta = linspace(0, 2*pi, 100);
% 
% for gi = 1:numel(gaussFitResults)
%     res = gaussFitResults(gi);
%     isRobust = (res.R2 >= rawR2Thresh) && (res.bootR2LowerCI >= bootCIThresh);
%     
%     if isRobust
%         if res.isDegenerateSigma
%             if res.sigmaX > res.sigmaY
%                 plot(ax, [azLims(1), azLims(2)], [res.y0, res.y0], 'Color', [0.85, 0.15, 0.15, 0.4], 'LineWidth', 1.5);
%                 fillY = [res.y0 - res.sigmaY, res.y0 - res.sigmaY, res.y0 + res.sigmaY, res.y0 + res.sigmaY];
%                 fillX = [azLims(1), azLims(2), azLims(2), azLims(1)];
%                 fill(ax, fillX, fillY, [0.85, 0.15, 0.15], 'FaceAlpha', 0.04, 'EdgeColor', 'none');
%             else
%                 plot(ax, [res.x0, res.x0], [elLims(1), elLims(2)], 'Color', [0.15, 0.15, 0.85, 0.4], 'LineWidth', 1.5);
%                 fillX = [res.x0 - res.sigmaX, res.x0 - res.sigmaX, res.x0 + res.sigmaX, res.x0 + res.sigmaX];
%                 fillY = [elLims(1), elLims(2), elLims(2), elLims(1)];
%                 fill(ax, fillX, fillY, [0.15, 0.15, 0.85], 'FaceAlpha', 0.04, 'EdgeColor', 'none');
%             end
%         % Replace just the 'else' block inside your main loop with this:
%         else
%             ellX = res.x0 + res.sigmaX * cos(theta);
%             ellY = res.y0 + res.sigmaY * sin(theta);
%             plot(ax, ellX, ellY, 'Color', [0.5, 0.5, 0.5, 0.15], 'LineWidth', 0.8);
%             fill(ax, ellX, ellY, [0.5, 0.5, 0.5], 'FaceAlpha', 0.05, 'EdgeColor', 'none');
%         end
%     end
% end
% 
% xlim(ax, azLims);
% ylim(ax, elLims);
% xlabel(ax, 'Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 11);
% ylabel(ax, 'Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 11);
% title(ax, 'Visual Field Population Coverage Map', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
% 
% dummyGrey = fill(ax, NaN, NaN, [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% dummyRed  = fill(ax, NaN, NaN, [0.85, 0.15, 0.15], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% dummyBlue = fill(ax, NaN, NaN, [0.15, 0.15, 0.85], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% legend([dummyGrey, dummyRed, dummyBlue], {'Localized RFs', 'Horizontal Stripe RFs', 'Vertical Stripe RFs'}, 'Location', 'northeastoutside');
% 
% set(ax, 'Box', 'off', 'TickDir', 'out', 'YDir', 'normal', 'FontName', 'Arial', 'FontSize', 10);
% axis(ax, 'square');
% grid(ax, 'on');
% ax.GridAlpha = 0.1;

%%
responsiveBoutons = [2330,3850, 3189];   %  3189, 
unresponsiveBoutons  = [2637, 69];           % 

exampleSet    = [responsiveBoutons, unresponsiveBoutons];
exampleLabels = {'Responsive', 'Responsive', 'Responsive', ...
                  'Unresponsive', 'Unresponsive'};

nEx  = numel(exampleSet);
colW = 1 / nEx;

figA = figure('Color', 'w', 'Position', [50 50 nEx*550 750], 'Name', 'Panel A: hand-picked examples');

for i = 1:nEx
    b = allRFMapping(exampleSet(i));
    xBase = (i-1)*colW;
    axPanel = axes('Position', [xBase+0.03 0.15 colW-0.05 0.70]);
    plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
    title(axPanel, sprintf('%s\n(Bouton %d)', exampleLabels{i}, exampleSet(i)), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
end


for i = 1:nEx
    b = allRFMapping(exampleSet(i));
    xBase = (i-1)*colW;
    
    axPanel = axes('Position', [xBase+0.03 0.15 colW-0.05 0.70]);
    
    plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
    
    title(axPanel, sprintf('%s\n(Bouton %d)', exampleLabels{i}, exampleSet(i)), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
end

% save
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\eg_rfs';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figA, 'Visible', 'off');
saveFigureFormats(figA, fullfile(outputDir, 'responsive_unresponsive_examplegrid_linesuperimposed'));


%%
trustR2Threshold  = 0.2;
robustCIThreshold = 0.2;   % matched

for gi = 1:numel(gaussFitResults)
    iROI = gaussFitResults(gi).iROI;
    validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));

    if numel(validBoot) >= 800
        lowerBound = prctile(validBoot, lowerPct);
    else
        lowerBound = NaN;
    end

    isRobust = ~isnan(lowerBound) && (lowerBound >= robustCIThreshold);

    r2here = gaussFitResults(gi).R2;
    if isempty(r2here), r2here = -Inf; end

    degenFlag = gaussFitResults(gi).isDegenerateSigma;
    if isempty(degenFlag), degenFlag = false; end
    degenFlag = logical(degenFlag(1)); % force scalar logical

    isTrusted = (r2here >= trustR2Threshold) && isRobust && ~degenFlag;

    gaussFitResults(gi).bootR2LowerCI = lowerBound;
    gaussFitResults(gi).isRobust      = isRobust;
    gaussFitResults(gi).isTrusted     = isTrusted;

    allRFMapping(iROI).gaussFit_bootR2LowerCI = lowerBound;
    allRFMapping(iROI).gaussFit_isRobust      = isRobust;
    allRFMapping(iROI).gaussFit_isTrusted     = isTrusted;
end


%%
fprintf('Number of trusted boutons: %d\n', sum([gaussFitResults.isTrusted]));

trustedIdx = [gaussFitResults.isTrusted];
azVals = [gaussFitResults(trustedIdx).sigmaX];
elVals = [gaussFitResults(trustedIdx).sigmaY];
fprintf('sigAz = %.2f +/- %.2f | sigEl = %.2f +/- %.2f\n', mean(azVals), std(azVals), mean(elVals), std(elVals));
%% Use this for visualising responsive boutons with fits 'Trusted = 1'
azLims = [min(uAz) - azStep/2, max(uAz) + azStep/2];
elLims = [min(uEl_plot) - elStep/2, max(uEl_plot) + elStep/2];

trustedIdxList = find([gaussFitResults.isTrusted]);   % <-- only iterate over trusted fits
numFitted = numel(trustedIdxList);

hFig = figure('Position', [50, 100, 1250, 400], 'Name', 'RF Comprehensive Diagnostic Viewer (Trusted Only)', 'NumberTitle', 'off');

for gi_i = 1:numFitted
    if ~ishandle(hFig)
        break;
    end

    gi = trustedIdxList(gi_i);   % map back to the real index in gaussFitResults
    res = gaussFitResults(gi);
    b = allRFMapping(res.iROI);

    % Reconstruct raw response grid
    trialMatrix = b.baselineSubtracted;
    respVec = nan(numel(trialMatrix), 1);
    for posIdx = 1:numel(trialMatrix)
        trials = trialMatrix{posIdx};
        if ~isempty(trials)
            trials = double(trials);
            respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
        end
    end
    respGrid = reshape(respVec, size(AzGrid));

    % Generate high-resolution model fit grid
    [xq, yq] = meshgrid(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100));
    pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
    fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));

    % Prepare boutonData struct dynamically for the trace-plotting function
    boutonData = b;
    boutonData.meanGridResponse = respGrid;
    boutonData.meanTemporalResponse = b.meanTemporalResponse;
    boutonData.meanBlankResponse = b.meanBlankResponse;
    boutonData.peakAmplitude = max(respGrid(:));

    clf(hFig);

    sgtitle(sprintf('ROI %d (%d of %d TRUSTED) | Raw R^2 = %.2f | Boot Lower Bound = %.2f | Trusted = %d', ...
        res.iROI, gi_i, numFitted, res.R2, res.bootR2LowerCI, res.isTrusted), 'FontWeight', 'bold');

    % Subplot 1: Trace-Overlaid Spatial Heatmap
    ax1 = subplot(1, 3, 1);
    plotRFHeatmapWithTraces(ax1, boutonData, uAz, uEl_plot, timeVector, ...
        'Colormap', 'bone', 'Smooth', false);
    xlim(ax1, azLims);
    ylim(ax1, elLims);
    set(ax1, 'YDir', 'normal');
    title(ax1, 'Mean Grid Response & Traces');

    % Subplot 2: Fitted 2D Gaussian Model
    ax2 = subplot(1, 3, 2);
    imagesc(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100), fittedGrid);
    axis image;
    colormap(ax2, 'bone');
    colorbar;
    xlim(ax2, azLims);
    ylim(ax2, elLims);
    set(ax2, 'YDir', 'normal');
    hold on;
    plot(ax2, res.x0, res.y0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    t_ellipse = linspace(0, 2*pi, 100);
    plot(ax2, res.x0 + res.sigmaX * cos(t_ellipse), res.y0 + res.sigmaY * sin(t_ellipse), 'r--', 'LineWidth', 1.5);
    title(ax2, 'Fitted 2D Gaussian Model');
    xlabel(ax2, 'Azimuth (deg)');
    ylabel(ax2, 'Elevation (deg)');

    % Subplot 3: Bootstrap R^2 Distribution
    ax3 = subplot(1, 3, 3);
    validBoot = res.bootR2(~isnan(res.bootR2));
    histogram(ax3, validBoot, linspace(-0.2, 1, 40), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    hold on;
    yl = ylim(ax3);
    line([res.R2, res.R2], yl, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-');
    line([res.bootR2LowerCI, res.bootR2LowerCI], yl, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
    line([0.2, 0.2], yl, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':');
    legend(ax3, 'Bootstrap R^2', 'Raw R^2', '5th Pct Lower Bound', 'Threshold (0.2)', 'Location', 'northwest');
    title(ax3, '1,000 Bootstrap R^2 Distribution');
    xlabel(ax3, 'R^2');
    ylabel(ax3, 'Count');
    grid on;

    drawnow;

    if gi_i < numFitted
        waitforbuttonpress;
    end
end


%%
%% Rebuild isTrusted: R^2 + bootstrap-CI + beyond-range ONLY (no sigma-degeneracy gate)
r2Vals    = [gaussFitResults.R2];
ciVals    = [gaussFitResults.bootR2LowerCI];
beyondAll = [gaussFitResults.isBeyondRange];

noSigmaGateMask = (r2Vals >= trustR2Threshold) & (ciVals >= robustCIThreshold) & ~beyondAll;

for gi = 1:numel(gaussFitResults)
    gaussFitResults(gi).isTrusted = noSigmaGateMask(gi);
    iROI = gaussFitResults(gi).iROI;
    allRFMapping(iROI).gaussFit_isTrusted = noSigmaGateMask(gi);
end

fprintf('New trusted set (R^2 + bootstrap-CI + beyond-range only): n = %d / %d\n', ...
    sum(noSigmaGateMask), nAboveFloor);
fprintf('sigAz = %.2f +/- %.2f (SD), median = %.2f\n', ...
    mean([gaussFitResults(noSigmaGateMask).sigmaX]), std([gaussFitResults(noSigmaGateMask).sigmaX]), median([gaussFitResults(noSigmaGateMask).sigmaX]));
fprintf('sigEl = %.2f +/- %.2f (SD), median = %.2f\n', ...
    mean([gaussFitResults(noSigmaGateMask).sigmaY]), std([gaussFitResults(noSigmaGateMask).sigmaY]), median([gaussFitResults(noSigmaGateMask).sigmaY]));

%% Shortlist: broadest boutons in the new trusted set, for manual visual review
trustedIdxNew = find(noSigmaGateMask);
combinedSigmaNew = [gaussFitResults(trustedIdxNew).sigmaX] + [gaussFitResults(trustedIdxNew).sigmaY];
[~, sortOrd] = sort(combinedSigmaNew, 'descend');

nToReview = min(15, numel(trustedIdxNew));
reviewList = trustedIdxNew(sortOrd(1:nToReview));

fprintf('\n=== Top %d broadest boutons in trusted set (manual review shortlist) ===\n', nToReview);
fprintf('%-6s %-8s %-8s %-8s %-8s\n', 'iROI', 'sigAz', 'sigEl', 'R2', 'bootCI');
for k = 1:nToReview
    idx = reviewList(k);
    fprintf('%-6d %-8.2f %-8.2f %-8.2f %-8.2f\n', ...
        gaussFitResults(idx).iROI, gaussFitResults(idx).sigmaX, gaussFitResults(idx).sigmaY, ...
        gaussFitResults(idx).R2, gaussFitResults(idx).bootR2LowerCI);
end

fprintf('\nTo review these visually, set:\n');
fprintf('  allIdxList = %s;\n', mat2str(reviewList));
fprintf('and run the diagnostic viewer loop with that list in place of trustedIdxList/allIdxList.\n');


%%
%% Two independent, separately-reportable identifiability checks
% 1. Is the CENTER estimate trustworthy? -> center must fall within the space actually tested
% 2. Is the WIDTH estimate trustworthy?  -> sigma must not exceed the space actually tested
%    (i.e. the curve must have somewhere within your grid to plausibly turn over)

r2Vals    = [gaussFitResults.R2];
ciVals    = [gaussFitResults.bootR2LowerCI];
beyondAll = [gaussFitResults.isBeyondRange];   % CENTER check (already exists)

sigAzAll = [gaussFitResults.sigmaX];
sigElAll = [gaussFitResults.sigmaY];

% WIDTH check: sigma cannot legitimately exceed the actual tested range on that axis
isImplausibleSigmaAz = sigAzAll > range(azRange);
isImplausibleSigmaEl = sigElAll > range(elRange);
isImplausibleSigma   = isImplausibleSigmaAz | isImplausibleSigmaEl;

%% Report each check separately, then combine for isTrusted
fprintf('Center off-screen (beyondRange): %d / %d\n', sum(beyondAll), numel(gaussFitResults));
fprintf('Sigma exceeds tested range (implausible width): %d / %d\n', sum(isImplausibleSigma), numel(gaussFitResults));

newTrustedMask = (r2Vals >= trustR2Threshold) & (ciVals >= robustCIThreshold) & ...
    ~beyondAll & ~isImplausibleSigma;

for gi = 1:numel(gaussFitResults)
    gaussFitResults(gi).isImplausibleSigma = isImplausibleSigma(gi); % store separately, for reporting
    gaussFitResults(gi).isTrusted = newTrustedMask(gi);
    iROI = gaussFitResults(gi).iROI;
    allRFMapping(iROI).gaussFit_isTrusted = newTrustedMask(gi);
end

fprintf('\nTrusted set (R^2 + bootstrap-CI + on-screen center + plausible sigma): n = %d / %d\n', sum(newTrustedMask), nAboveFloor);
fprintf('sigAz = %.2f +/- %.2f (SD), median = %.2f\n', ...
    mean([gaussFitResults(newTrustedMask).sigmaX]), std([gaussFitResults(newTrustedMask).sigmaX]), median([gaussFitResults(newTrustedMask).sigmaX]));
fprintf('sigEl = %.2f +/- %.2f (SD), median = %.2f\n', ...
    mean([gaussFitResults(newTrustedMask).sigmaY]), std([gaussFitResults(newTrustedMask).sigmaY]), median([gaussFitResults(newTrustedMask).sigmaY]));

%%
fprintf('Current azRangeExt = [%.1f, %.1f]\n', azRangeExt(1), azRangeExt(2));
fprintf('Current elRangeExt = [%.1f, %.1f]\n', elRangeExt(1), elRangeExt(2));
fprintf('Current screenAzLimits = [%.1f, %.1f]\n', screenAzLimits(1), screenAzLimits(2));
fprintf('Current screenElLimits = [%.1f, %.1f]\n', screenElLimits(1), screenElLimits(2));

%%
beyondIdx = find([gaussFitResults.isBeyondRange]);
azCenters = [gaussFitResults(beyondIdx).x0];
elCenters = [gaussFitResults(beyondIdx).y0];

fprintf('n beyond-range = %d\n', numel(beyondIdx));
fprintf('Range of x0 among beyond-range boutons: [%.1f, %.1f]\n', min(azCenters), max(azCenters));
fprintf('Range of y0 among beyond-range boutons: [%.1f, %.1f]\n', min(elCenters), max(elCenters));
fprintf('azRangeExt = [%.1f, %.1f], elRangeExt = [%.1f, %.1f]\n', azRangeExt(1), azRangeExt(2), elRangeExt(1), elRangeExt(2));

% How many are beyond on azimuth specifically vs elevation specifically
fprintf('isBeyondRangeAz: %d / %d\n', sum([gaussFitResults.isBeyondRangeAz]), numel(gaussFitResults));
fprintf('isBeyondRangeEl: %d / %d\n', sum([gaussFitResults.isBeyondRangeEl]), numel(gaussFitResults));

% Look at how far outside range these actually are (distance past the nearest edge)
azOvershoot = max(azRangeExt(1) - azCenters, azCenters - azRangeExt(2));
azOvershoot(azOvershoot < 0) = 0; % not actually beyond on this axis
elOvershoot = max(elRangeExt(1) - elCenters, elCenters - elRangeExt(2));
elOvershoot(elOvershoot < 0) = 0;

fprintf('Azimuth overshoot (deg past edge): median=%.1f, max=%.1f\n', median(azOvershoot(azOvershoot>0)), max(azOvershoot));
fprintf('Elevation overshoot (deg past edge): median=%.1f, max=%.1f\n', median(elOvershoot(elOvershoot>0)), max(elOvershoot));