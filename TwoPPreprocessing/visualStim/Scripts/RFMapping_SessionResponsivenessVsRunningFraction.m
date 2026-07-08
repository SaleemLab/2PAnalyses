% RFMapping_SessionResponsivenessVsRunningFraction.m
%
% Extracts the percentage of responsive boutons per session and plots it
% against the fraction of running trials in that session to test for a
% behavioral bias in your responsive pool yield.

%% Calculate responsiveness per session
uniqueSessAll = unique(sessionLabels, 'stable');
nSessTotal = numel(uniqueSessAll);

sessPctResponsive = zeros(nSessTotal, 1);
sessRunningFraction = zeros(nSessTotal, 1);
sessTotalBoutons = zeros(nSessTotal, 1);

isRespAll = [allRFMapping.isResponsive];

for s = 1:nSessTotal
    thisLabel = uniqueSessAll{s};
    
    % Track responsive vs total boutons for this session
    sessBoutonIdx = strcmp(sessionLabels, thisLabel);
    sessTotalBoutons(s) = sum(sessBoutonIdx);
    sessPctResponsive(s) = (sum(isRespAll(sessBoutonIdx)) / sessTotalBoutons(s)) * 100;
    
    % Match session to find behavioral running fraction from sessionStats
    statIdx = find(strcmp({sessionStats.sessionLabel}, thisLabel), 1);
    if ~isempty(statIdx)
        % Total grid + blank trials classified
        gRun   = sessionStats(statIdx).nGridRun;
        gStat  = sessionStats(statIdx).nGridStat;
        gAmbig = sessionStats(statIdx).nGridAmbig;
        bRun   = sessionStats(statIdx).nBlankRun;
        bStat  = sessionStats(statIdx).nBlankStat;
        bAmbig = sessionStats(statIdx).nBlankAmbig;
        
        totalClassifiedTrials = gRun + gStat + gAmbig + bRun + bStat + bAmbig;
        if totalClassifiedTrials > 0
            sessRunningFraction(s) = ((gRun + bRun) / totalClassifiedTrials) * 100;
        else
            sessRunningFraction(s) = NaN;
        end
    else
        sessRunningFraction(s) = NaN;
    end
end

% Clean out any session that missed behavioral mapping
validSess = ~isnan(sessRunningFraction);
xData = sessRunningFraction(validSess);
yData = sessPctResponsive(validSess);
labelsPreserved = uniqueSessAll(validSess);

%% Plot Session Correlation
figCorr = figure('Color', 'w', 'Position', [100 100 450 450], ...
    'Name', 'Session Yield vs Running Fraction');
hold on;

% Add a linear fit line if there are enough points
if numel(xData) > 2
    ft = polyfit(xData, yData, 1);
    xFit = linspace(0, 100, 100);
    yFit = polyval(ft, xFit);
    plot(xFit, yFit, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit');
    

    [h_x, p_normX] = lillietest(xData);
    [h_y, p_normY] = lillietest(yData);
    fprintf('Normality p-values: X (Running) = %.3f | Y (Yield) = %.3f\n', p_normX, p_normY);

    % Compute correlation coefficient
    [R, pVal] = corrcoef(xData, yData);
    [R_spearman, p_spearman] = corr(xData, yData, 'type', 'Spearman');
    fprintf('Spearman: R = %.2f, p = %.3f\n', R_spearman, p_spearman);
    text(5, 95, sprintf('R = %.2f\np = %.3f', R(1,2), pVal(1,2)), ...
        'FontName', 'Arial', 'FontSize', 10, 'VerticalAlignment', 'top');
end

% Size dots by total number of boutons tracked in that session
scatterSizes = max(20, sessTotalBoutons(validSess) * 0.5); 
hScatter = scatter(xData, yData, scatterSizes, [0.2 0.4 0.6], 'filled', ...
    'MarkerFaceAlpha', 0.8, 'DisplayName', 'Sessions');

if ~isempty(hScatter)
    hScatter.MarkerEdgeColor = [0.1 0.2 0.3];
end

xlabel('% of Total Trials spent Running', 'FontName', 'Arial', 'FontSize', 10);
ylabel('% Responsive Boutons in Session', 'FontName', 'Arial', 'FontSize', 10);
xlim([0 100]); ylim([0 25]);
axis square;
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);
title({'Yield Yield vs. Behavioral Running Bias', ...
       sprintf('(n = %d sessions analyzed)', sum(validSess))}, ...
    'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');

%% Save
% outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\session_behavior_correlation';
% if ~exist(outputDir, 'dir'), mkdir(outputDir); end
% set(figCorr, 'Visible', 'off');
% saveFigureFormats(figCorr, fullfile(outputDir, 'responsiveness_vs_running_fraction_correlation'));