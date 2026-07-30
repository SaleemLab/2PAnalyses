% RFMapping_PlotGaussianFitExamples.m
%
% Visualizes a sample of the 2D Gaussian fits from RFMapping_Gaussian2DFit.m
% against the underlying grid data, so you can SEE what "stable" vs.
% "unstable" and "high R^2" vs. "low R^2" actually look like, rather than
% just reading the summary numbers.
%
% For each selected bouton: plots the raw grid heatmap (like your Image 1
% figures), overlays the fitted Gaussian's 1-sigma ellipse and center, and
% additionally re-runs a SMALL bootstrap (nBootstrapPlot resamples, just
% for these examples) so the scatter of resampled centers is drawn on top
% -- this is the most direct way to SEE what "bootstrap center SD" means:
% a tight cloud of dots = stable; a spread-out cloud = unstable.
%
% UPDATED to match the current field names in RFMapping_Gaussian2DFit.m:
% the old combined 'isStable' field no longer exists on gaussFitResults --
% it's been folded into 'isTrusted' (TRUE only if BOTH axes pass coverage
% + non-degenerate-sigma + not-beyond-range). Categorization below now
% uses 'isTrusted' directly. Also, since the main script now fits ALL
% ROIs (not just responsive ones -- see note in that script), each plot
% title now also shows isResponsive, so you can see at a glance whether
% an example bouton with a passing fit was actually visually responsive
% or just cleared the R^2 floor on noise.
%
% REQUIRES: allRFMapping, uAz, uEl_plot, respIdx, gaussFitResults (from
% RFMapping_Gaussian2DFit.m -- run that first, keep gaussFitResults in
% the workspace).

%% params
nBootstrapPlot = 100;  % smaller than the main run -- just enough to see the spread, for a handful of boutons
nExamplesPerCategory = 3;

if ~exist('gaussFitResults', 'var')
    error('gaussFitResults not found -- run RFMapping_Gaussian2DFit.m first and keep it in the workspace.');
end

azStep  = mean(diff(sort(uAz)));
elStep  = mean(diff(sort(uEl_plot)));
[AzGrid, ElGrid] = meshgrid(uAz, uEl_plot);
xyList = [AzGrid(:), ElGrid(:)];
gaussFit2D = @(p, xy) p(1) * exp( -( (xy(:,1)-p(2)).^2 ./ (2*p(4)^2) + (xy(:,2)-p(3)).^2 ./ (2*p(5)^2) ) );

fitOpts = optimoptions('lsqcurvefit', 'Display', 'off');
azRange = [min(uAz) max(uAz)];
elRange = [min(uEl_plot) max(uEl_plot)];
lb = [0,   azRange(1)-azStep, elRange(1)-elStep, azStep/4,        elStep/4];
ub = [Inf, azRange(2)+azStep, elRange(2)+elStep, 3*range(azRange), 3*range(elRange)];

%% categorize existing fits
% NOTE: 'isTrusted' already combines bootstrap coverage-stability + non-degenerate-sigma +
% not-beyond-range across BOTH axes (see RFMapping_Gaussian2DFit.m), so it's used directly here
% rather than re-deriving it from the old 'isStable'/'isBeyondRange'/'isDegenerateSigma' fields.
allR2vals   = [gaussFitResults.R2];
allTrusted  = [gaussFitResults.isTrusted];

catStableHighR2   = find(allTrusted  & allR2vals >= 0.6);
catStableLowR2    = find(allTrusted  & allR2vals >= 0.3 & allR2vals < 0.5);
catUnstableHighR2 = find(~allTrusted & allR2vals >= 0.6);
catUnstableLowR2  = find(~allTrusted & allR2vals >= 0.3 & allR2vals < 0.5);

categories = {
    'Stable, high R^2',   catStableHighR2;
    'Stable, low R^2',    catStableLowR2;
    'UNSTABLE, high R^2', catUnstableHighR2;
    'UNSTABLE, low R^2',  catUnstableLowR2;
    };

%% pick example boutons from each category
exampleList = struct('idx', {}, 'label', {});
for c = 1:size(categories, 1)
    label   = categories{c, 1};
    poolIdx = categories{c, 2};
    if isempty(poolIdx)
        fprintf('No boutons found for category "%s" -- skipping.\n', label);
        continue;
    end
    nPick = min(nExamplesPerCategory, numel(poolIdx));
    pickedIdx = poolIdx(randperm(numel(poolIdx), nPick)); % random sample within category
    for p = 1:nPick
        exampleList(end+1) = struct('idx', pickedIdx(p), 'label', label); %#ok<SAGROW>
    end
end

nExamples = numel(exampleList);
fprintf('Plotting %d example fits across %d categories...\n', nExamples, size(categories,1));

%% plot
nCols = 3;
nRows = ceil(nExamples / nCols);
figExamples = figure('Color', 'w', 'Position', [50 50 350*nCols 320*nRows], 'Name', 'Example Gaussian fits');

for e = 1:nExamples
    gi   = exampleList(e).idx;
    fitR = gaussFitResults(gi);
    iROI = fitR.iROI;
    b    = allRFMapping(iROI);

    subplot(nRows, nCols, e);
    imagesc(uAz, uEl_plot, double(b.meanGridResponse));
    set(gca, 'YDir', 'normal');
    colormap(gca, gray);
    hold on;

    % 1-sigma ellipse of the fit
    theta = linspace(0, 2*pi, 100);
    ellX  = fitR.x0 + fitR.sigmaX * cos(theta);
    ellY  = fitR.y0 + fitR.sigmaY * sin(theta);
    plot(ellX, ellY, 'r-', 'LineWidth', 1.5);
    plot(fitR.x0, fitR.y0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);

    % small bootstrap just for this bouton, to visualize center spread
    trialMatrix = b.baselineSubtracted;
    if ~isempty(trialMatrix)
        bootCentersPlot = nan(nBootstrapPlot, 2);
        for bIter = 1:nBootstrapPlot
            bootRespVec = nan(numel(xyList(:,1)), 1);
            for posIdx = 1:numel(bootRespVec)
                [rIdx, cIdx] = ind2sub(size(b.meanGridResponse), posIdx);
                trials = trialMatrix{rIdx, cIdx};
                if isempty(trials), continue; end
                trials = double(trials);
                nTrialsHere = size(trials, 1);
                if nTrialsHere == 0, continue; end
                sampIdx = randi(nTrialsHere, nTrialsHere, 1);
                bootRespVec(posIdx) = mean(mean(trials(sampIdx, respIdx), 2, 'omitnan'), 'omitnan');
            end
            if any(isnan(bootRespVec)), continue; end
            try
                pBoot = lsqcurvefit(gaussFit2D, [fitR.A fitR.x0 fitR.y0 fitR.sigmaX fitR.sigmaY], ...
                    xyList, bootRespVec, lb, ub, fitOpts);
                bootCentersPlot(bIter, :) = [pBoot(2), pBoot(3)];
            catch
                continue;
            end
        end
        validRows = ~any(isnan(bootCentersPlot), 2);
        scatter(bootCentersPlot(validRows,1), bootCentersPlot(validRows,2), 10, 'y', 'filled', 'MarkerFaceAlpha', 0.5);
    end

    xlabel('Azimuth (\circ)'); ylabel('Elevation (\circ)');
    % NOTE: 'resp' here shows isResponsive so you can tell apart a real, visually-responsive
    % bouton that happens to fit well from a non-responsive one whose noise cleared the R^2
    % floor by chance -- relevant now that the main script fits ALL ROIs, not just responsive ones.
    title(sprintf('%s\nBouton %d | R^2=%.2f | trusted=%d | resp=%d', ...
        exampleList(e).label, iROI, fitR.R2, fitR.isTrusted, fitR.isResponsive), ...
        'FontSize', 8, 'FontWeight', 'normal');
    axis tight; hold off;
end

sgtitle('Example Gaussian fits: red = fitted 1\sigma ellipse + center, yellow dots = bootstrap-refit centers');