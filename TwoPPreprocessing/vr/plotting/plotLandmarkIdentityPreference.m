% function plotSMI_LandmarkIdentityPreference_Both(RSPData, VISpData)
%     excludeStart = 30; 
%     excludeEnd = 30;
%     
%     landmarkCentres = [40, 80, 120, 160];
%     gratingLandmarks = [40, 120];
%     plaidLandmarks   = [80, 160];
%     tolerance = 15; 
%     
%     pooledSMI_RSP_Grating = [];
%     pooledSMI_RSP_Plaid   = [];
%     countRSP_Grating = 0;
%     countRSP_Plaid   = 0;
%     
%     fprintf('Processing RSP Landmark Preferences...\n');
%     for s = 1:length(RSPData)
%         sess = RSPData(s);
%         if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
%            ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
%             continue;
%         end
%         
%         lapActivity = sess.ConditionData.Baseline.LapActivity;
%         [~, ~, nBins] = size(lapActivity);
%         validLapRange = (excludeStart + 1) : (nBins - excludeEnd);
%         
%         meanOdd = squeeze(mean(lapActivity(:, 1:2:end, :), 2, 'omitnan'));
%         smiValues = sess.SMI.SMI;
%         
%         for i = 1:length(sess.FilteredROIs)
%             roiIdx = sess.FilteredROIs(i);
%             smiVal = smiValues(roiIdx);
%             if isnan(smiVal); continue; end
%             
%             roiTrace = meanOdd(roiIdx, :);
%             [pks, locs] = findpeaks(roiTrace);
%             validMask = ismember(locs, validLapRange);
%             validPks = pks(validMask); validLocs = locs(validMask);
%             
%             if ~isempty(validLocs)
%                 [~, maxPeakIdx] = max(validPks);
%                 prefBin = validLocs(maxPeakIdx);
%             else
%                 [~, maxRelIdx] = max(roiTrace(validLapRange));
%                 prefBin = validLapRange(maxRelIdx);
%             end
%             
%             [minDist, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
%             
%             if minDist <= tolerance
%                 targetLandmark = landmarkCentres(closestLandmarkIdx);
%                 
%                 if ismember(targetLandmark, gratingLandmarks)
%                     countRSP_Grating = countRSP_Grating + 1;
%                     pooledSMI_RSP_Grating = [pooledSMI_RSP_Grating; smiVal];
%                 elseif ismember(targetLandmark, plaidLandmarks)
%                     countRSP_Plaid = countRSP_Plaid + 1;
%                     pooledSMI_RSP_Plaid = [pooledSMI_RSP_Plaid; smiVal];
%                 end
%             end
%         end
%     end
%     
%     pooledSMI_VISp_Grating = [];
%     pooledSMI_VISp_Plaid   = [];
%     countVISp_Grating = 0;
%     countVISp_Plaid   = 0;
%     
%     fprintf('Processing VISp Landmark Preferences...\n');
%     for s = 1:length(VISpData)
%         sess = VISpData(s);
%         if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
%            ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
%             continue;
%         end
%         
%         lapActivity = sess.ConditionData.Baseline.LapActivity;
%         [~, ~, nBins] = size(lapActivity);
%         validLapRange = (excludeStart + 1) : (nBins - excludeEnd);
%         
%         meanOdd = squeeze(mean(lapActivity(:, 1:2:end, :), 2, 'omitnan'));
%         smiValues = sess.SMI.SMI;
%         
%         for i = 1:length(sess.FilteredROIs)
%             roiIdx = sess.FilteredROIs(i);
%             smiVal = smiValues(roiIdx);
%             if isnan(smiVal); continue; end
%             
%             roiTrace = meanOdd(roiIdx, :);
%             [pks, locs] = findpeaks(roiTrace);
%             validMask = ismember(locs, validLapRange);
%             validPks = pks(validMask); validLocs = locs(validMask);
%             
%             if ~isempty(validLocs)
%                 [~, maxPeakIdx] = max(validPks);
%                 prefBin = validLocs(maxPeakIdx);
%             else
%                 [~, maxRelIdx] = max(roiTrace(validLapRange));
%                 prefBin = validLapRange(maxRelIdx);
%             end
%             
%             [minDist, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
%             
%             if minDist <= tolerance
%                 targetLandmark = landmarkCentres(closestLandmarkIdx);
%                 
%                 if ismember(targetLandmark, gratingLandmarks)
%                     countVISp_Grating = countVISp_Grating + 1;
%                     pooledSMI_VISp_Grating = [pooledSMI_VISp_Grating; smiVal];
%                 elseif ismember(targetLandmark, plaidLandmarks)
%                     countVISp_Plaid = countVISp_Plaid + 1;
%                     pooledSMI_VISp_Plaid = [pooledSMI_VISp_Plaid; smiVal];
%                 end
%             end
%         end
%     end
%     
%     totalRSP = countRSP_Grating + countRSP_Plaid;
%     pctRSP_Grating = (countRSP_Grating / totalRSP) * 100;
%     pctRSP_Plaid   = (countRSP_Plaid / totalRSP) * 100;
%     [pValRSP] = ranksum(pooledSMI_RSP_Grating, pooledSMI_RSP_Plaid);
%     
%     totalVISp = countVISp_Grating + countVISp_Plaid;
%     pctVISp_Grating = (countVISp_Grating / totalVISp) * 100;
%     pctVISp_Plaid   = (countVISp_Plaid / totalVISp) * 100;
%     [pValVISp] = ranksum(pooledSMI_VISp_Grating, pooledSMI_VISp_Plaid);
%     
%     figHandle = figure('Name', 'Landmark Identity Preference: RSP vs VISp', ...
%                        'Color', [1 1 1], 'Position', [150, 150, 600, 550]);
%                    
%     % RSP BAR CHART
%     subplot(2, 2, 1); hold on;
%     grid off;
%     b1 = bar([pctRSP_Grating, pctRSP_Plaid], 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.5);
%     b1.CData(1,:) = [0.15 0.15 0.15]; 
%     b1.CData(2,:) = [0.60 0.60 0.60]; 
%     
%     set(gca, 'Box', 'off', 'Visible', 'off'); 
%     xlim([0.5, 2.5]); ylim([0, 100]); axis square;
%     
%     plot([0.4, 0.4], [0, 100], 'k', 'LineWidth', 0.5, 'Clipping', 'off'); 
%     
%     for yTick = 0:20:100
%         plot([0.33, 0.4], [yTick, yTick], 'k', 'LineWidth', 0.5, 'Clipping', 'off');
%         text(0.23, yTick, num2str(yTick), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
%     end
%     
%     text(1, -8, 'Grating', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
%     text(2, -8, 'Plaid', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
%     
%     text(-0.05, 50, '% of selective ROIs', 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
%     title(sprintf('RSP Preferred Feature\n(n = %d)', totalRSP));
%     if exist('offsetAxes(gca);', 'file') == 2; offsetAxes(gca); end
%     
%     % RSP CDF CHART
%     subplot(2, 2, 2); hold on;
%     grid off;
%     [fRG, xRG] = ecdf(pooledSMI_RSP_Grating);
%     [fRP, xRP] = ecdf(pooledSMI_RSP_Plaid);
%     plot(xRG, fRG, 'LineWidth', 1.5, 'Color', [0.15 0.15 0.15], 'DisplayName', sprintf('Grating (n=%d)', length(pooledSMI_RSP_Grating)));
%     plot(xRP, fRP, 'LineWidth', 1.5, 'Color', [0.60 0.60 0.60], 'DisplayName', sprintf('Plaid (n=%d)', length(pooledSMI_RSP_Plaid)));
%     xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
%     yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
%     xlabel('Spatial modulation index'); ylabel('Cumulative probability');
%     title(sprintf('RSP Modulation \n(p = %.4e)', pValRSP));
%     xlim([-1.1, 1.1]); ylim([0, 1.02]); 
%     box off;
%     axis square;
%     defaultAxesProperties(gca);
%     offsetAxes(gca);
%     
%     % VISp SOMAS BAR CHART
%     subplot(2, 2, 3); hold on;
%     grid off;
%     b2 = bar([pctVISp_Grating, pctVISp_Plaid], 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.5);
%     b2.CData(1,:) = [0.15 0.15 0.15]; 
%     b2.CData(2,:) = [0.60 0.60 0.60]; 
%     
%     set(gca, 'Box', 'off', 'Visible', 'off'); 
%     xlim([0.5, 2.5]); ylim([0, 100]); axis square;
%     
%     plot([0.4, 0.4], [0, 100], 'k', 'LineWidth', 0.5, 'Clipping', 'off'); 
%     
%     for yTick = 0:20:100
%         plot([0.33, 0.4], [yTick, yTick], 'k', 'LineWidth', 0.5, 'Clipping', 'off');
%         text(0.23, yTick, num2str(yTick), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
%     end
%     
%     text(1, -8, 'Grating', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
%     text(2, -8, 'Plaid', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
%     
%     text(-0.05, 50, '% of selective ROIs', 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
%     title(sprintf('VISp Preferred Feature\n(n = %d)', totalVISp));
%     if exist('offsetAxes(gca);', 'file') == 2; offsetAxes(gca); end
%     
%     % VISp SOMAS CDF CHART
%     subplot(2, 2, 4); hold on;
%     grid off;
%     [fVG, xVG] = ecdf(pooledSMI_VISp_Grating);
%     [fVP, xVP] = ecdf(pooledSMI_VISp_Plaid);
%     plot(xVG, fVG, 'LineWidth', 1.5, 'Color', [0.15 0.15 0.15], 'DisplayName', sprintf('Grating (n=%d)', length(pooledSMI_VISp_Grating)));
%     plot(xVP, fVP, 'LineWidth', 1.5, 'Color', [0.60 0.60 0.60], 'DisplayName', sprintf('Plaid (n=%d)', length(pooledSMI_VISp_Plaid)));
%     xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
%     yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
%     xlabel('Spatial modulation index'); ylabel('Cumulative probability');
%     title(sprintf('VISp Modulation \n(p = %.4e)', pValVISp));
%     xlim([-1.1, 1.1]); ylim([0, 1.02]);
%     box off;
%     axis square;
%     defaultAxesProperties(gca);
%     offsetAxes(gca);
%     
%     outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1\RSPVsVISp\';
%     if ~exist(outputDir, 'dir')
%         mkdir(outputDir);
%     end
%     
%     baseFileName = 'landmark_identity_preference_both_regions';
%     fullSavePath = fullfile(outputDir, baseFileName);
%     
%     saveFigureFormats(figHandle, fullSavePath);
%     fprintf('Landmark identity preference analysis completed and plots saved for both regions.\n');
% end


function plotLandmarkIdentityPreference(RSPData, VISpData)
    
    landmarkCentres  = [40, 80, 120, 160];
    gratingLandmarks = [40, 120];
    plaidLandmarks   = [80, 160];
    tolerance = 15;        % Max distance from center to count as a landmark peak (+/- 15cm)
    windowRadius = 15;     % Radius around centers to measure activity (+/- 15cm)
    
    % --- Extract Pure Quantitative Arrays ---
    fprintf('Extracting RSP population profiles...\n');
    [rspSI, rspPeaks] = extractSplitMetrics(RSPData, ...
        landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius);
        
    fprintf('Extracting VISp population profiles...\n');
    [vispSI, vispPeaks] = extractSplitMetrics(VISpData, ...
        landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius);
    
    if isempty(rspSI) || isempty(vispSI)
        error('Data extraction returned 0 cells. Check track coordinates or field names.');
    end
    
    rspG  = sum(ismember(rspPeaks,  gratingLandmarks)) / length(rspPeaks)  * 100;
    rspP  = sum(ismember(rspPeaks,  plaidLandmarks))   / length(rspPeaks)  * 100;
    vispG = sum(ismember(vispPeaks, gratingLandmarks)) / length(vispPeaks) * 100;
    vispP = sum(ismember(vispPeaks, plaidLandmarks))   / length(vispPeaks) * 100;
    
    barData = [rspG, vispG; rspP, vispP];
    p_Bars   = ranksum(rspPeaks, vispPeaks);
    p_Curves = ranksum(rspSI, vispSI);

    % Within-region: does RSP significantly prefer grating over plaid?
    rspGCount  = sum(ismember(rspPeaks, gratingLandmarks));
    rspPCount  = sum(ismember(rspPeaks, plaidLandmarks));
    p_RSP_G_vs_P = chi2GoodnessOfFit5050(rspGCount, rspPCount);

    % Within-region: does VISp significantly prefer plaid over grating?
    vispGCount = sum(ismember(vispPeaks, gratingLandmarks));
    vispPCount = sum(ismember(vispPeaks, plaidLandmarks));
    p_VISp_G_vs_P = chi2GoodnessOfFit5050(vispGCount, vispPCount);

    fprintf('\n=== Within-Region Grating vs Plaid Preference ===\n');
    fprintf('RSP:  Grating n=%d (%.1f%%)  Plaid n=%d (%.1f%%)  p = %.3e\n', ...
        rspGCount, rspG, rspPCount, rspP, p_RSP_G_vs_P);
    fprintf('VISp: Grating n=%d (%.1f%%)  Plaid n=%d (%.1f%%)  p = %.3e\n', ...
        vispGCount, vispG, vispPCount, vispP, p_VISp_G_vs_P);
    
    %% Plotting
    fig = figure('Color', 'w', 'Position', [100 100 850 420]);
    
    % LANDMARK IDENTITY PREFERENCE (BAR CHART)
    subplot(1, 2, 1); hold on;
    b = bar(barData, 'EdgeColor', 'none', 'BarWidth', 0.8);
    b(1).FaceColor = 'k';
    b(2).FaceColor = [0.6 0.6 0.6];
    
    set(gca, 'Box', 'off', 'XTick', 1:2, 'XTickLabel', {'Prefers Grating', 'Prefers Plaid'}, ...
             'FontName', 'Arial', 'FontSize', 11);
    ylabel('% of Landmark Population', 'FontName', 'Arial', 'FontSize', 12);
    title(sprintf('Landmark Absolute Peak Allocation\nPopulation Shift p = %.3e\nRSP G vs P p = %.3e | VISp G vs P p = %.3e', ...
        p_Bars, p_RSP_G_vs_P, p_VISp_G_vs_P), 'FontWeight', 'normal');
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
    
    text(0.02, maxHeight * 0.80, '\leftarrow Less Selective', 'HorizontalAlignment', 'left', ...
        'FontName', 'Arial', 'Color', [0.5 0.5 0.5], 'FontSize', 9);
    text(0.98, maxHeight * 0.80, 'More Selective \rightarrow', 'HorizontalAlignment', 'right', ...
        'FontName', 'Arial', 'Color', [0.5 0.5 0.5], 'FontSize', 9);
    
    xlim([-0.02, 1.02]); ylim([0, maxHeight * 1.15]);
    xlabel('Feature Selectivity Index (Normalized)', 'FontName', 'Arial', 'FontSize', 11);
    ylabel('Proportion of ROIs', 'FontName', 'Arial', 'FontSize', 11);
    title(sprintf('Identity Tuning Discrimination\nSelectivity Shift p = %.3e', p_Curves), ...
        'FontWeight', 'normal', 'FontSize', 12);
    
    legend([h1, h2], {sprintf('RSP (n=%d)', length(rspSI)), sprintf('VISp (n=%d)', length(vispSI))}, ...
        'Location', 'southoutside', 'Box', 'off', 'FontName', 'Arial');
        
    axis square; box off;
    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
    if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end
    
    % Save
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4\feature_selective_rsp_visp';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    saveFigureFormats(fig, fullfile(outputDir, 'landmark_preference_and_selectivity'));
end

%% Chi-Square Goodness of Fit (50/50 split)
function p = chi2GoodnessOfFit5050(count1, count2)
    observed = [count1, count2];
    total = sum(observed);
    if total == 0, p = 1; return; end
    expected = [total/2, total/2];
    chi2Val = sum((observed - expected).^2 ./ expected);
    df = 1;
    p = 1 - chi2cdf(chi2Val, df);
end

%%
function [siVector, peakLocations] = extractSplitMetrics(RegionData, ...
    landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius)
    siVector = []; peakLocations = [];
    w_space = gausswin(15); w_space = w_space / sum(w_space);
    
    gratingBins = []; for g = gratingLandmarks, gratingBins = [gratingBins, (g - windowRadius):(g + windowRadius)]; end
    plaidBins   = []; for p = plaidLandmarks,   plaidBins   = [plaidBins,   (p - windowRadius):(p + windowRadius)]; end
    
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
            
            actG    = mean(max(0, normTestEven(gratingBins)), 'omitnan');
            actP    = mean(max(0, normTestEven(plaidBins)),   'omitnan');
            actPeak = max(0, normTestEven(prefBin));
            
            if ismember(primaryLandmark, gratingLandmarks)
                actAlternative = actP;
            else
                actAlternative = actG;
            end
            
            if (actPeak + actAlternative) > 0
                valSI = (actPeak - actAlternative) / (actPeak + actAlternative);
                if ~isnan(valSI)
                    siVector      = [siVector;      max(0, min(1, valSI))];
                    peakLocations = [peakLocations; primaryLandmark];
                end
            end
        end
    end
end

