function plotLandmarkIdentityPreference(RSPData, VISpData)
    
    landmarkCentres  = [40, 80, 120, 160];
    gratingLandmarks = [40, 120];
    plaidLandmarks   = [80, 160];
    tolerance = 15;        % Max distance from center to count as a landmark peak (+/- 15cm)
    windowRadius = 15;     % Radius around centers to measure activity (+/- 15cm)
    
    % --- Extract Pure Quantitative Arrays ---
    fprintf('Extracting RSP population profiles...\n');
    [rspSI, rspPeaks, rspSessionID] = extractSplitMetrics(RSPData, ...
        landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius);
        
    fprintf('Extracting VISp population profiles...\n');
    [vispSI, vispPeaks, vispSessionID] = extractSplitMetrics(VISpData, ...
        landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius);
    
    if isempty(rspSI) || isempty(vispSI)
        error('Data extraction returned 0 cells. Check track coordinates or field names.');
    end
    
    p_Curves = ranksum(rspSI, vispSI);

    % median and median absolute deviation (MAD) of each SI distribution
    median_rspSI  = median(rspSI);
    median_vispSI = median(vispSI);
    mad_rspSI     = mad(rspSI, 1);
    mad_vispSI    = mad(vispSI, 1);
    
    % --- Per-session grating/plaid percentages, for the scatter + mean +/- SD bars ---
    rspSession  = computeLandmarkSessionStats(rspPeaks,  rspSessionID,  gratingLandmarks, plaidLandmarks);
    vispSession = computeLandmarkSessionStats(vispPeaks, vispSessionID, gratingLandmarks, plaidLandmarks);
    
    rspMeanG  = mean(rspSession.pctGrating);
    rspMeanP  = mean(rspSession.pctPlaid);
    vispMeanG = mean(vispSession.pctGrating);
    vispMeanP = mean(vispSession.pctPlaid);
    
    rspSDG  = std(rspSession.pctGrating);
    rspSDP  = std(rspSession.pctPlaid);
    vispSDG = std(vispSession.pctGrating);
    vispSDP = std(vispSession.pctPlaid);
    
    % CHANGED: consistent, session-level statistics for all four bar-to-bar
    % comparisons.
    % Within-region comparisons (RSP-G vs RSP-P, VISp-G vs VISp-P) are
    % PAIRED, since the grating and plaid percentages for a given session
    % come from the same session and are not independent (they sum to
    % ~100%), so Wilcoxon signed-rank (signrank) is used, not chi-square
    % on pooled cells, which would incorrectly treat all cells as
    % independent and inflate the effective sample size.
    p_RSP_G_vs_P  = signrank(rspSession.pctGrating,  rspSession.pctPlaid);
    p_VISp_G_vs_P = signrank(vispSession.pctGrating, vispSession.pctPlaid);
    
    % Cross-region comparisons (RSP-G vs VISp-G, RSP-P vs VISp-P) are
    % UNPAIRED, since RSP and VISp sessions are independent (different
    % mice/sessions), so Wilcoxon rank-sum (ranksum) is used.
    p_Grating_RSPvsVISp = ranksum(rspSession.pctGrating, vispSession.pctGrating);
    p_Plaid_RSPvsVISp   = ranksum(rspSession.pctPlaid,   vispSession.pctPlaid);

    fprintf('\n=== Landmark Identity Preference: Session-Level Statistics ===\n');
    fprintf('RSP:  Grating %.1f +/- %.1f%% (SD)  Plaid %.1f +/- %.1f%% (SD)  [n=%d sessions]  signrank p = %.3e\n', ...
        rspMeanG, rspSDG, rspMeanP, rspSDP, length(rspSession.pctGrating), p_RSP_G_vs_P);
    fprintf('VISp: Grating %.1f +/- %.1f%% (SD)  Plaid %.1f +/- %.1f%% (SD)  [n=%d sessions]  signrank p = %.3e\n', ...
        vispMeanG, vispSDG, vispMeanP, vispSDP, length(vispSession.pctGrating), p_VISp_G_vs_P);
    fprintf('Cross-region: Grating RSP vs VISp ranksum p = %.3e | Plaid RSP vs VISp ranksum p = %.3e\n', ...
        p_Grating_RSPvsVISp, p_Plaid_RSPvsVISp);
    fprintf('SI: RSP median = %.3f +/- %.3f MAD, VISp median = %.3f +/- %.3f MAD\n', ...
        median_rspSI, mad_rspSI, median_vispSI, mad_vispSI);
    
    barData = [rspMeanG, vispMeanG; rspMeanP, vispMeanP];
    
    %% Plotting
    fig = figure('Color', 'w', 'Position', [100 100 850 420]);
    
    % LANDMARK IDENTITY PREFERENCE (BAR CHART)
    subplot(1, 2, 1); hold on;
    b = bar(barData, 'EdgeColor', 'none', 'BarWidth', 0.8);
    b(1).FaceColor = 'k';
    b(2).FaceColor = [0.6 0.6 0.6];
    
    jitterWidth = 0.12;
    scatterFaceAlpha = 0.6;
    
    rspX_G   = b(1).XEndPoints(1) + (rand(size(rspSession.pctGrating))  - 0.5) * jitterWidth;
    rspX_P   = b(1).XEndPoints(2) + (rand(size(rspSession.pctPlaid))    - 0.5) * jitterWidth;
    vispX_G  = b(2).XEndPoints(1) + (rand(size(vispSession.pctGrating)) - 0.5) * jitterWidth;
    vispX_P  = b(2).XEndPoints(2) + (rand(size(vispSession.pctPlaid))   - 0.5) * jitterWidth;
    
    scatter(rspX_G,  rspSession.pctGrating,  28, 'k', 'filled', ...
        'MarkerFaceAlpha', scatterFaceAlpha, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    scatter(rspX_P,  rspSession.pctPlaid,    28, 'k', 'filled', ...
        'MarkerFaceAlpha', scatterFaceAlpha, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    scatter(vispX_G, vispSession.pctGrating, 28, [0.6 0.6 0.6], 'filled', ...
        'MarkerFaceAlpha', scatterFaceAlpha, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    scatter(vispX_P, vispSession.pctPlaid,   28, [0.6 0.6 0.6], 'filled', ...
        'MarkerFaceAlpha', scatterFaceAlpha, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    errorbar([b(1).XEndPoints(1), b(1).XEndPoints(2)], ...
             [rspMeanG, rspMeanP], ...
             [rspSDG, rspSDP], ...
             'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.2, 'CapSize', 8, 'HandleVisibility', 'off');
    errorbar([b(2).XEndPoints(1), b(2).XEndPoints(2)], ...
             [vispMeanG, vispMeanP], ...
             [vispSDG, vispSDP], ...
             'LineStyle', 'none', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'CapSize', 8, 'HandleVisibility', 'off');
    
    % CHANGED: significance brackets for the four bar-to-bar comparisons.
    % Heights are stacked upward from the tallest error bar so they don't
    % overlap; adjust stepHeight if your ylim changes.
    topOfBars = max([rspMeanG+rspSDG, rspMeanP+rspSDP, vispMeanG+vispSDG, vispMeanP+vispSDP]);
    stepHeight = 6;
    
    y1 = topOfBars + stepHeight;        % within-region RSP
    y2 = topOfBars + stepHeight;        % within-region VISp
    y3 = topOfBars + stepHeight * 2.2;  % cross-region grating
    y4 = topOfBars + stepHeight * 3.4;  % cross-region plaid
    
    drawSigBracket(b(1).XEndPoints(1), b(1).XEndPoints(2), y1, p_RSP_G_vs_P, 'k');
    drawSigBracket(b(2).XEndPoints(1), b(2).XEndPoints(2), y2, p_VISp_G_vs_P, [0.4 0.4 0.4]);
    drawSigBracket(b(1).XEndPoints(1), b(2).XEndPoints(1), y3, p_Grating_RSPvsVISp, [0.2 0.2 0.2]);
    drawSigBracket(b(1).XEndPoints(2), b(2).XEndPoints(2), y4, p_Plaid_RSPvsVISp, [0.2 0.2 0.2]);
    
    set(gca, 'Box', 'off', 'XTick', 1:2, 'XTickLabel', {'Prefers Grating', 'Prefers Plaid'}, ...
             'FontName', 'Arial', 'FontSize', 11);
    ylabel('% of Landmark Population (mean +/- SD, dots = sessions)', 'FontName', 'Arial', 'FontSize', 12);
    ylim([0, y4 + stepHeight*1.5]);
    title('Landmark Absolute Peak Allocation', 'FontWeight', 'normal');
    legend(b, {sprintf('RSP (n=%d sessions)', length(rspSession.pctGrating)), ...
               sprintf('VISp (n=%d sessions)', length(vispSession.pctGrating))}, ...
        'Location', 'southoutside', 'Box', 'off', 'FontName', 'Arial');
    axis square; box off;
    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
    if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end
    
    % UNIDIRECTIONAL MODULATION STRENGTH (PROPORTION OF ROIs)
    subplot(1, 2, 2); hold on;
    binEdges = 0:0.05:1; 
    
    h2 = histogram(vispSI, 'BinEdges', binEdges, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 2.5);
    
    h1 = histogram(rspSI, 'BinEdges', binEdges, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 2.5);
    
    maxHeight = max([h1.Values, h2.Values]);
    
    plot(median_rspSI, maxHeight * 0.03, 'x', 'Color', 'k', ...
        'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(median_vispSI, maxHeight * 0.03, 'x', 'Color', [0.6 0.6 0.6], ...
        'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    
    text(0.02, maxHeight * 0.80, '\leftarrow Less Selective', 'HorizontalAlignment', 'left', ...
        'FontName', 'Arial', 'Color', [0.5 0.5 0.5], 'FontSize', 9);
    text(0.98, maxHeight * 0.80, 'More Selective \rightarrow', 'HorizontalAlignment', 'right', ...
        'FontName', 'Arial', 'Color', [0.5 0.5 0.5], 'FontSize', 9);
    
    xlim([-0.02, 1.02]); ylim([0, maxHeight * 1.15]);
    xlabel('Feature Selectivity Index (Normalized)', 'FontName', 'Arial', 'FontSize', 11);
    ylabel('Proportion of ROIs', 'FontName', 'Arial', 'FontSize', 11);
    title(sprintf('Identity Tuning Discrimination\nSelectivity Shift p = %.3e\nRSP: %.3f +/- %.3f MAD | VISp: %.3f +/- %.3f MAD', ...
        p_Curves, median_rspSI, mad_rspSI, median_vispSI, mad_vispSI), 'FontWeight', 'normal', 'FontSize', 12);
    
    legend([h1, h2], {sprintf('RSP (n=%d)', length(rspSI)), sprintf('VISp (n=%d)', length(vispSI))}, ...
        'Location', 'southoutside', 'Box', 'off', 'FontName', 'Arial');
        
    axis square; box off;
    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
    if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end
    
    % Save
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\feature_selective_rsp_visp';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    saveFigureFormats(fig, fullfile(outputDir, 'landmark_preference_and_selectivity'));

        p_Curves = ranksum(rspSI, vispSI);

    % median and median absolute deviation (MAD) of each SI distribution
    median_rspSI  = median(rspSI);
    median_vispSI = median(vispSI);
    mad_rspSI     = mad(rspSI, 1);
    mad_vispSI    = mad(vispSI, 1);

    % CHANGED: percentage of cells with SI above 0.25 and above 0.5,
    % matching the same thresholds used for the SMI summary earlier.
    pctRSP_SI025  = (sum(rspSI > 0.25) / length(rspSI)) * 100;
    pctRSP_SI05   = (sum(rspSI > 0.5)  / length(rspSI)) * 100;
    pctVISp_SI025 = (sum(vispSI > 0.25) / length(vispSI)) * 100;
    pctVISp_SI05  = (sum(vispSI > 0.5)  / length(vispSI)) * 100;

        fprintf('SI: RSP median = %.3f +/- %.3f MAD, VISp median = %.3f +/- %.3f MAD\n', ...
        median_rspSI, mad_rspSI, median_vispSI, mad_vispSI);
    fprintf('SI thresholds: RSP %.1f%% (%d/%d) > 0.25, %.1f%% (%d/%d) > 0.5\n', ...
        pctRSP_SI025, sum(rspSI > 0.25), length(rspSI), pctRSP_SI05, sum(rspSI > 0.5), length(rspSI));
    fprintf('SI thresholds: VISp %.1f%% (%d/%d) > 0.25, %.1f%% (%d/%d) > 0.5\n', ...
        pctVISp_SI025, sum(vispSI > 0.25), length(vispSI), pctVISp_SI05, sum(vispSI > 0.5), length(vispSI));
end

%% Draw a significance bracket with stars (or "n.s.") between two x positions
function drawSigBracket(x1, x2, y, pVal, lineColor)
    tickHeight = 0.8;
    plot([x1, x1, x2, x2], [y, y+tickHeight, y+tickHeight, y], ...
        'Color', lineColor, 'LineWidth', 1, 'HandleVisibility', 'off');
    text(mean([x1, x2]), y + tickHeight*1.3, pToStars(pVal), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontName', 'Arial', 'Color', lineColor);
end


%% cross-region grating, cross-region plaid).
function starStr = pToStars(pVal)
    numComparisons = 4;
    alpha001 = 0.001 / numComparisons;
    alpha01  = 0.01  / numComparisons;
    alpha05  = 0.05  / numComparisons;
    
    if pVal <= alpha001
        starStr = '***';
    elseif pVal <= alpha01
        starStr = '**';
    elseif pVal <= alpha05
        starStr = '*';
    else
        starStr = 'n.s.';
    end
end

%% Group pooled peak locations by session and compute per-session grating/plaid percentages
function stats = computeLandmarkSessionStats(peakLocations, sessionID, gratingLandmarks, plaidLandmarks)
    pctGrating = [];
    pctPlaid   = [];
    nCellsArr  = [];
    
    uniqueSessions = unique(sessionID);
    for u = 1:length(uniqueSessions)
        thisSess = uniqueSessions(u);
        sessPeaks = peakLocations(sessionID == thisSess);
        nCells = length(sessPeaks);
        if nCells == 0, continue; end
        
        nG = sum(ismember(sessPeaks, gratingLandmarks));
        nP = sum(ismember(sessPeaks, plaidLandmarks));
        
        pctGrating(end+1) = (nG / nCells) * 100; %#ok<AGROW>
        pctPlaid(end+1)   = (nP / nCells) * 100; %#ok<AGROW>
        nCellsArr(end+1)  = nCells;              %#ok<AGROW>
    end
    
    stats.pctGrating = pctGrating;
    stats.pctPlaid   = pctPlaid;
    stats.nCells     = nCellsArr;
end

function [siVector, peakLocations, sessionID] = extractSplitMetrics(RegionData, ...
    landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius)
    
    siVector = []; peakLocations = []; sessionID = [];
    w_space = gausswin(15); w_space = w_space / sum(w_space);
    
    allowedLandmarkBins = [];
    for c = landmarkCentres
        allowedLandmarkBins = [allowedLandmarkBins, (c - tolerance):(c + tolerance)];
    end
    allowedLandmarkBins = unique(allowedLandmarkBins);
    
    for s = 1:length(RegionData)
        sess = RegionData(s);
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
           ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs), continue; end
        
        lapActivity = sess.ConditionData.Baseline.LapActivity;
        [numROIsTotal, numLaps, numBins] = size(lapActivity);
        
        validSearchBins = allowedLandmarkBins(allowedLandmarkBins >= 1 & allowedLandmarkBins <= numBins);
        
        smoothedActivity = lapActivity;
        for iCell = 1:numROIsTotal
            for iLap = 1:numLaps
                trace = squeeze(lapActivity(iCell, iLap, :));
                if all(isnan(trace)), continue; end
                nanMask = isnan(trace); trace(nanMask) = 0;
                smoothed = filtfilt(w_space, 1, trace); smoothed(nanMask) = NaN;
                smoothedActivity(iCell, iLap, :) = smoothed;
            end
        end
        
        roisToAnalyze = sess.FilteredROIs;
        roiActivity   = smoothedActivity(roisToAnalyze, :, :);
        numROIs       = length(roisToAnalyze);
        if numROIs == 0, continue; end
        
        meanOdd  = squeeze(mean(roiActivity(:, 1:2:end, :), 2, 'omitnan'));
        meanEven = squeeze(mean(roiActivity(:, 2:2:end, :), 2, 'omitnan'));
        
        hasSavedExclusions = isfield(sess, 'SMI') && isfield(sess.SMI, 'ExcludeEdgePeakCells') && ...
                             length(sess.SMI.ExcludeEdgePeakCells) == numROIsTotal;
        
        for i = 1:numROIs
            if hasSavedExclusions
                if sess.SMI.ExcludeEdgePeakCells(roisToAnalyze(i)), continue; end
            end
            
            trainTrace = meanOdd(i, :);
            testTrace  = meanEven(i, :);
            if all(isnan(trainTrace)) || all(isnan(testTrace)), continue; end
            
            minOdd = min(trainTrace, [], 'omitnan');
            maxOdd = max(trainTrace, [], 'omitnan');
            rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end
            
            normTrainOdd  = (trainTrace - minOdd) ./ rangeOdd;
            normTestEven  = (testTrace  - minOdd) ./ rangeOdd;
            
            [~, maxIdxInSearch] = max(normTrainOdd(validSearchBins));
            prefBin = validSearchBins(maxIdxInSearch);
            
            [~, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
            primaryLandmark = landmarkCentres(closestLandmarkIdx);
            
            numBinsThisROI = length(testTrace);
            
            gPeaks = nan(1, length(gratingLandmarks));
            for gi = 1:length(gratingLandmarks)
                thisWindow = (gratingLandmarks(gi) - windowRadius):(gratingLandmarks(gi) + windowRadius);
                thisWindow = thisWindow(thisWindow >= 1 & thisWindow <= numBinsThisROI);
                if isempty(thisWindow), continue; end
                gPeaks(gi) = max(max(0, normTestEven(thisWindow)), [], 'omitnan');
            end
            actG = mean(gPeaks, 'omitnan');
            
            pPeaks = nan(1, length(plaidLandmarks));
            for pi = 1:length(plaidLandmarks)
                thisWindow = (plaidLandmarks(pi) - windowRadius):(plaidLandmarks(pi) + windowRadius);
                thisWindow = thisWindow(thisWindow >= 1 & thisWindow <= numBinsThisROI);
                if isempty(thisWindow), continue; end
                pPeaks(pi) = max(max(0, normTestEven(thisWindow)), [], 'omitnan');
            end
            actP = mean(pPeaks, 'omitnan');
            
            actPeak = max(0, normTestEven(prefBin));
            
            if ismember(primaryLandmark, gratingLandmarks)
                actAlternative = actP;
            else
                actAlternative = actG;
            end
            
            if isnan(actAlternative), continue; end
            
            if (actPeak + actAlternative) > 0
                valSI = (actPeak - actAlternative) / (actPeak + actAlternative);
                if ~isnan(valSI)
                    siVector      = [siVector;      max(0, min(1, valSI))];
                    peakLocations = [peakLocations; primaryLandmark];
                    sessionID     = [sessionID;      s];
                end
            end
        end
    end
end