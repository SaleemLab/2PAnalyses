function plotSMI_LandmarkIdentityPreference_Both(RSPData, VISpData)
    excludeStart = 30; 
    excludeEnd = 30;
    
    landmarkCentres = [40, 80, 120, 160];
    gratingLandmarks = [40, 120];
    plaidLandmarks   = [80, 160];
    tolerance = 15; 
    
    pooledSMI_RSP_Grating = [];
    pooledSMI_RSP_Plaid   = [];
    countRSP_Grating = 0;
    countRSP_Plaid   = 0;
    
    fprintf('Processing RSP Landmark Preferences...\n');
    for s = 1:length(RSPData)
        sess = RSPData(s);
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
           ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end
        
        lapActivity = sess.ConditionData.Baseline.LapActivity;
        [~, ~, nBins] = size(lapActivity);
        validLapRange = (excludeStart + 1) : (nBins - excludeEnd);
        
        meanOdd = squeeze(mean(lapActivity(:, 1:2:end, :), 2, 'omitnan'));
        smiValues = sess.SMI.SMI;
        
        for i = 1:length(sess.FilteredROIs)
            roiIdx = sess.FilteredROIs(i);
            smiVal = smiValues(roiIdx);
            if isnan(smiVal); continue; end
            
            roiTrace = meanOdd(roiIdx, :);
            [pks, locs] = findpeaks(roiTrace);
            validMask = ismember(locs, validLapRange);
            validPks = pks(validMask); validLocs = locs(validMask);
            
            if ~isempty(validLocs)
                [~, maxPeakIdx] = max(validPks);
                prefBin = validLocs(maxPeakIdx);
            else
                [~, maxRelIdx] = max(roiTrace(validLapRange));
                prefBin = validLapRange(maxRelIdx);
            end
            
            [minDist, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
            
            if minDist <= tolerance
                targetLandmark = landmarkCentres(closestLandmarkIdx);
                
                if ismember(targetLandmark, gratingLandmarks)
                    countRSP_Grating = countRSP_Grating + 1;
                    pooledSMI_RSP_Grating = [pooledSMI_RSP_Grating; smiVal];
                elseif ismember(targetLandmark, plaidLandmarks)
                    countRSP_Plaid = countRSP_Plaid + 1;
                    pooledSMI_RSP_Plaid = [pooledSMI_RSP_Plaid; smiVal];
                end
            end
        end
    end
    
    pooledSMI_VISp_Grating = [];
    pooledSMI_VISp_Plaid   = [];
    countVISp_Grating = 0;
    countVISp_Plaid   = 0;
    
    fprintf('Processing VISp Landmark Preferences...\n');
    for s = 1:length(VISpData)
        sess = VISpData(s);
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
           ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end
        
        lapActivity = sess.ConditionData.Baseline.LapActivity;
        [~, ~, nBins] = size(lapActivity);
        validLapRange = (excludeStart + 1) : (nBins - excludeEnd);
        
        meanOdd = squeeze(mean(lapActivity(:, 1:2:end, :), 2, 'omitnan'));
        smiValues = sess.SMI.SMI;
        
        for i = 1:length(sess.FilteredROIs)
            roiIdx = sess.FilteredROIs(i);
            smiVal = smiValues(roiIdx);
            if isnan(smiVal); continue; end
            
            roiTrace = meanOdd(roiIdx, :);
            [pks, locs] = findpeaks(roiTrace);
            validMask = ismember(locs, validLapRange);
            validPks = pks(validMask); validLocs = locs(validMask);
            
            if ~isempty(validLocs)
                [~, maxPeakIdx] = max(validPks);
                prefBin = validLocs(maxPeakIdx);
            else
                [~, maxRelIdx] = max(roiTrace(validLapRange));
                prefBin = validLapRange(maxRelIdx);
            end
            
            [minDist, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
            
            if minDist <= tolerance
                targetLandmark = landmarkCentres(closestLandmarkIdx);
                
                if ismember(targetLandmark, gratingLandmarks)
                    countVISp_Grating = countVISp_Grating + 1;
                    pooledSMI_VISp_Grating = [pooledSMI_VISp_Grating; smiVal];
                elseif ismember(targetLandmark, plaidLandmarks)
                    countVISp_Plaid = countVISp_Plaid + 1;
                    pooledSMI_VISp_Plaid = [pooledSMI_VISp_Plaid; smiVal];
                end
            end
        end
    end
    
    totalRSP = countRSP_Grating + countRSP_Plaid;
    pctRSP_Grating = (countRSP_Grating / totalRSP) * 100;
    pctRSP_Plaid   = (countRSP_Plaid / totalRSP) * 100;
    [pValRSP] = ranksum(pooledSMI_RSP_Grating, pooledSMI_RSP_Plaid);
    
    totalVISp = countVISp_Grating + countVISp_Plaid;
    pctVISp_Grating = (countVISp_Grating / totalVISp) * 100;
    pctVISp_Plaid   = (countVISp_Plaid / totalVISp) * 100;
    [pValVISp] = ranksum(pooledSMI_VISp_Grating, pooledSMI_VISp_Plaid);
    
    figHandle = figure('Name', 'Landmark Identity Preference: RSP vs VISp', ...
                       'Color', [1 1 1], 'Position', [150, 150, 600, 550]);
                   
    % RSP BAR CHART
    subplot(2, 2, 1); hold on;
    grid off;
    b1 = bar([pctRSP_Grating, pctRSP_Plaid], 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.5);
    b1.CData(1,:) = [0.15 0.15 0.15]; 
    b1.CData(2,:) = [0.60 0.60 0.60]; 
    
    set(gca, 'Box', 'off', 'Visible', 'off'); 
    xlim([0.5, 2.5]); ylim([0, 100]); axis square;
    
    plot([0.4, 0.4], [0, 100], 'k', 'LineWidth', 0.5, 'Clipping', 'off'); 
    
    for yTick = 0:20:100
        plot([0.33, 0.4], [yTick, yTick], 'k', 'LineWidth', 0.5, 'Clipping', 'off');
        text(0.23, yTick, num2str(yTick), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
    
    text(1, -8, 'Grating', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(2, -8, 'Plaid', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    
    text(-0.05, 50, '% of selective ROIs', 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    title(sprintf('RSP Preferred Feature\n(n = %d)', totalRSP));
    if exist('offsetAxes(gca);', 'file') == 2; offsetAxes(gca); end
    
    % RSP CDF CHART
    subplot(2, 2, 2); hold on;
    grid off;
    [fRG, xRG] = ecdf(pooledSMI_RSP_Grating);
    [fRP, xRP] = ecdf(pooledSMI_RSP_Plaid);
    plot(xRG, fRG, 'LineWidth', 1.5, 'Color', [0.15 0.15 0.15], 'DisplayName', sprintf('Grating (n=%d)', length(pooledSMI_RSP_Grating)));
    plot(xRP, fRP, 'LineWidth', 1.5, 'Color', [0.60 0.60 0.60], 'DisplayName', sprintf('Plaid (n=%d)', length(pooledSMI_RSP_Plaid)));
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlabel('Spatial modulation index'); ylabel('Cumulative probability');
    title(sprintf('RSP Modulation \n(p = %.4e)', pValRSP));
    xlim([-1.1, 1.1]); ylim([0, 1.02]); 
    box off;
    axis square;
    defaultAxesProperties(gca);
    offsetAxes(gca);
    
    % VISp SOMAS BAR CHART
    subplot(2, 2, 3); hold on;
    grid off;
    b2 = bar([pctVISp_Grating, pctVISp_Plaid], 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.5);
    b2.CData(1,:) = [0.15 0.15 0.15]; 
    b2.CData(2,:) = [0.60 0.60 0.60]; 
    
    set(gca, 'Box', 'off', 'Visible', 'off'); 
    xlim([0.5, 2.5]); ylim([0, 100]); axis square;
    
    plot([0.4, 0.4], [0, 100], 'k', 'LineWidth', 0.5, 'Clipping', 'off'); 
    
    for yTick = 0:20:100
        plot([0.33, 0.4], [yTick, yTick], 'k', 'LineWidth', 0.5, 'Clipping', 'off');
        text(0.23, yTick, num2str(yTick), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end
    
    text(1, -8, 'Grating', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(2, -8, 'Plaid', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    
    text(-0.05, 50, '% of selective ROIs', 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    title(sprintf('VISp Preferred Feature\n(n = %d)', totalVISp));
    if exist('offsetAxes(gca);', 'file') == 2; offsetAxes(gca); end
    
    % VISp SOMAS CDF CHART
    subplot(2, 2, 4); hold on;
    grid off;
    [fVG, xVG] = ecdf(pooledSMI_VISp_Grating);
    [fVP, xVP] = ecdf(pooledSMI_VISp_Plaid);
    plot(xVG, fVG, 'LineWidth', 1.5, 'Color', [0.15 0.15 0.15], 'DisplayName', sprintf('Grating (n=%d)', length(pooledSMI_VISp_Grating)));
    plot(xVP, fVP, 'LineWidth', 1.5, 'Color', [0.60 0.60 0.60], 'DisplayName', sprintf('Plaid (n=%d)', length(pooledSMI_VISp_Plaid)));
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlabel('Spatial modulation index'); ylabel('Cumulative probability');
    title(sprintf('VISp Modulation \n(p = %.4e)', pValVISp));
    xlim([-1.1, 1.1]); ylim([0, 1.02]);
    box off;
    axis square;
    defaultAxesProperties(gca);
    offsetAxes(gca);
    
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1\RSPVsVISp\';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    baseFileName = 'landmark_identity_preference_both_regions';
    fullSavePath = fullfile(outputDir, baseFileName);
    
    saveFigureFormats(figHandle, fullSavePath);
    fprintf('Landmark identity preference analysis completed and plots saved for both regions.\n');
end