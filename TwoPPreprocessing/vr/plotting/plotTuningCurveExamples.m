function FigH=plotTuningCurveExamples(lapPositionActivity, cvExpVar, roiIdx1, roiIdx2, nFolds, foldToShow)
% Plots train/test tuning curves, raw residuals, and squared residuals
% for two example ROIs, one ROI per ROW (curve | raw residual | squared residual).
% Displays the R^2 for the SPECIFIC fold shown (foldToShow), not the mean across all folds.
if nargin < 5 || isempty(nFolds), nFolds = 5; end
if nargin < 6 || isempty(foldToShow), foldToShow = 1; end
rng(33, 'twister'); % FIXED SEED
exampleROIs = [roiIdx1, roiIdx2];
FigH=figure('Name', 'Tuning curves and residuals', 'Color', 'w', 'Position', [100 100 1000 550]);
for i = 1:2
    roiIdx = exampleROIs(i);
    cellActivity = squeeze(lapPositionActivity(roiIdx, :, :)); % [Laps x Position]
    numTrials = size(cellActivity, 1);
    perms = randperm(numTrials);
    kidx = round(linspace(0, numTrials, nFolds + 1));
    testTrials  = perms(kidx(foldToShow)+1 : kidx(foldToShow+1));
    trainTrials = perms;
    trainTrials(ismember(trainTrials, testTrials)) = [];
    trainCurve = mean(cellActivity(trainTrials, :), 1, 'omitnan');
    testCurve  = mean(cellActivity(testTrials, :), 1, 'omitnan');
    trainMean  = mean(trainCurve, 'omitnan');
    residModel = testCurve - trainCurve;
    residNull  = testCurve - trainMean;
    sqResidModel = residModel.^2;
    sqResidNull  = residNull.^2;

    % R^2 for THIS fold only, computed directly from this fold's own curves
    % (matches the curves/residuals actually plotted, rather than cvExpVar's mean across folds)
    resSS = sum(sqResidModel, 'omitnan');
    totSS = sum(sqResidNull, 'omitnan');
    foldR2 = 1 - resSS/totSS;

% --- Column 1: tuning curves ---
    subplot(2, 3, (i-1)*3 + 1);
    plot(trainCurve, 'b-', 'LineWidth', 1.5); hold on;
    plot(testCurve, 'r-', 'LineWidth', 1.5);
    yline(trainMean, 'k--', 'LineWidth', 1);
    title(sprintf('ROI %d, fold R^2 = %.2f', roiIdx, foldR2), 'FontSize', 9);
    ylabel('Activity');
if i == 2, xlabel('Position bin (cm)'); end
    xticks([40 80 120 160]);
    set(gca, 'TickDir', 'out', 'box', 'off');
if i == 1
        legend({'y_{train}', 'y_{test}', '\bar{y}_{train}'}, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'latex');
        legend boxoff;
end
    defaultAxesProperties(gca, 1)
% --- Column 2: raw residuals ---
    subplot(2, 3, (i-1)*3 + 2);
    plot(residModel, 'Color', [0.2 0.5 0.8], 'LineWidth', 1.2); hold on;
    plot(residNull, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
    yline(0, 'k-', 'LineWidth', 0.5);
    ylabel('Raw residual');
if i == 2, xlabel('Position bin (cm)'); end
    xticks([40 80 120 160]);
    set(gca, 'TickDir', 'out', 'box', 'off');
if i == 1
        legend({'$y_{test} - y_{train}$', '$y_{test} - \bar{y}_{train}$'}, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'latex');
        legend boxoff;
end
    defaultAxesProperties(gca, 1)
% --- Column 3: squared residuals ---
    subplot(2, 3, (i-1)*3 + 3);
    plot(sqResidModel, 'Color', [0.2 0.5 0.8], 'LineWidth', 1.2); hold on;
    plot(sqResidNull, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
    ylabel('Squared residual');
if i == 2, xlabel('Position bin (cm)'); end
    xticks([40 80 120 160]);
    set(gca, 'TickDir', 'out', 'box', 'off');
    ax = gca;
    ax.YAxis.Exponent = 0;
    title(sprintf('resSS=%.3f, totSS=%.3f', resSS, totSS), 'FontSize', 8);
    defaultAxesProperties(gca, 1);
end
end