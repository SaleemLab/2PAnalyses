function plotPeaksCorridorSplit(RSPData, VISpData)
    % 
    trackStart = 30;
    trackEnd = 170;
    midPoint = 100; % Exact middle of the 30-170cm valid corridor
    
    fprintf('Extracting RSP peaks from SMI_Metrics fields...\n');
    rspPeaks = extractSavedStructPeaks(RSPData, trackStart, trackEnd);
    
    fprintf('Extracting VISp peaks from SMI_Metrics fields...\n');
    vispPeaks = extractSavedStructPeaks(VISpData, trackStart, trackEnd);
    
    % Count Cells in first half vs. second half 
    rspFirstHalf  = sum(rspPeaks >= trackStart & rspPeaks <= midPoint);
    rspSecondHalf = sum(rspPeaks > midPoint & rspPeaks <= trackEnd);
    
    vispFirstHalf  = sum(vispPeaks >= trackStart & vispPeaks <= midPoint);
    vispSecondHalf = sum(vispPeaks > midPoint & vispPeaks <= trackEnd);
    
    % Convert to percentages 
    rspProportions  = [rspFirstHalf, rspSecondHalf] / length(rspPeaks) * 100;
    vispProportions = [vispFirstHalf, vispSecondHalf] / length(vispPeaks) * 100;
    
    %  Rows = Halves, Columns = Regions
    barData = [rspProportions(1), vispProportions(1); ...
               rspProportions(2), vispProportions(2)];
           
    %  cross-regional comparison (Is RSP different from VISp overall?)
    observed = [rspFirstHalf, vispFirstHalf; rspSecondHalf, vispSecondHalf];
    p_SpatialShift = chi2TestIndependence(observed);
    
    % within-region comparison (Does each region prefer one half over the other?)
    p_RSP_1st_vs_2nd = chi2GoodnessOfFit5050(rspFirstHalf, rspSecondHalf);
    p_VISp_1st_vs_2nd = chi2GoodnessOfFit5050(vispFirstHalf, vispSecondHalf);
    
    %% plotting
    fig = figure('Color', 'w', 'Position', [150 150 550 480]);
    hold on;
    
    b = bar(barData, 'EdgeColor', 'none', 'BarWidth', 0.8);
    b(1).FaceColor = 'k';           % RSP = Black
    b(2).FaceColor = [0.6 0.6 0.6]; % VISp = Gray
    
    set(gca, 'Box', 'off', 'XTick', 1:2, ...
             'XTickLabel', {sprintf('First Half (%d-%d cm)', trackStart, midPoint), ...
                            sprintf('Second Half (%d-%d cm)', midPoint, trackEnd)}, ...
             'FontName', 'Arial', 'FontSize', 11);
    ylabel('% of ROIs', 'FontName', 'Arial', 'FontSize', 12);
    
    title(sprintf('SMI RespP bin \nRegional Diff p = %.2e\nWithin-RSP p = %.2e | Within-VISp p = %.2e', ...
        p_SpatialShift, p_RSP_1st_vs_2nd, p_VISp_1st_vs_2nd), ...
        'FontWeight', 'normal', 'FontSize', 11, 'FontName', 'Arial');
    
    legend(b, {sprintf('RSP (n=%d)', length(rspPeaks)), sprintf('VISp (n=%d)', length(vispPeaks))}, ...
        'Location', 'southoutside', 'Box', 'off', 'FontName', 'Arial');
        
    axis square; box off;
    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
    if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end
    hold off;

    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4\peak_lap-halves_rsp_visp';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    saveFigureFormats(fig, fullfile(outputDir, 'peakbin_halves_within_lap_rsp_visp'));

    %% Per-session breakdown — check if bias is driven by one session/animal
    fprintf('\n=== Per-Session Peak Half Analysis ===\n');

    for region = 1:2
        if region == 1
            Data = RSPData;
            regionName = 'RSP';
        else
            Data = VISpData;
            regionName = 'VISp';
        end

        fprintf('\n--- %s ---\n', regionName);
        fprintf('%-20s %-12s %-12s %-12s %-12s\n', 'Session', 'N Cells', 'First Half', 'Second Half', '% First');

        for s = 1:length(Data)
            sess = Data(s);
            if ~isfield(sess, 'SMI') || ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
                continue;
            end

            RespP        = sess.SMI.RpBin;
            excludeFlags = sess.SMI.ExcludeEdgePeakCells;
            roisToAnalyze = sess.FilteredROIs;

            sessionPeaks = [];
            for i = 1:length(roisToAnalyze)
                roiIdx = roisToAnalyze(i);
                if excludeFlags(roiIdx), continue; end
                cellPeakBin = RespP(roiIdx);
                if cellPeakBin >= trackStart && cellPeakBin <= trackEnd
                    sessionPeaks = [sessionPeaks; cellPeakBin];
                end
            end

            if isempty(sessionPeaks), continue; end

            nCells  = length(sessionPeaks);
            nFirst  = sum(sessionPeaks <= midPoint);
            nSecond = sum(sessionPeaks > midPoint);
            pctFirst = (nFirst / nCells) * 100;

            % Chi-square test for this session alone vs 50/50
            p_sess = chi2GoodnessOfFit5050(nFirst, nSecond);

            if isfield(sess, 'sessionName')
                sessLabel = sess.sessionName;
            elseif isfield(sess, 'animal')
                sessLabel = sprintf('%s_s%d', sess.animal, s);
            else
                sessLabel = sprintf('Session_%d', s);
            end

            fprintf('%-20s %-12d %-12d %-12d %-12.1f  (p = %.3f)\n', ...
                sessLabel, nCells, nFirst, nSecond, pctFirst, p_sess);
        end
    end

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

%% Use the RespP bin already computed for SMI analyses 
function [validPeaks] = extractSavedStructPeaks(RegionData, trackStart, trackEnd)
    validPeaks = [];
    for s = 1:length(RegionData)
        sess = RegionData(s);
        if ~isfield(sess, 'SMI') || ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end
        RespP = sess.SMI.RpBin;
        excludeFlags = sess.SMI.ExcludeEdgePeakCells;
        roisToAnalyze = sess.FilteredROIs;
        for i = 1:length(roisToAnalyze)
            roiIdx = roisToAnalyze(i);
            if excludeFlags(roiIdx)
                continue;
            end
            cellPeakBin = RespP(roiIdx);
            if cellPeakBin >= trackStart && cellPeakBin <= trackEnd
                validPeaks = [validPeaks; cellPeakBin];
            end
        end
    end
end

%% Chi-Square Test of Independence
function p = chi2TestIndependence(observed)
    rowTotals = sum(observed, 2);
    colTotals = sum(observed, 1);
    grandTotal = sum(observed, 'all');
    if grandTotal == 0, p = 1; return; end
    expected = (rowTotals * colTotals) / grandTotal;
    expected(expected == 0) = eps; 
    chi2Val = sum((observed - expected).^2 ./ expected, 'all');
    df = (size(observed, 1) - 1) * (size(observed, 2) - 1);
    p = 1 - chi2cdf(chi2Val, df);
end