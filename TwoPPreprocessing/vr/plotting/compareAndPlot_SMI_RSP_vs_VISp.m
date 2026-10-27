function compareAndPlot_SMI_RSP_vs_VISp(RSPData, VISpData, outputDir)
    
    pooledSMI_RSP = [];
    for s = 1:length(RSPData)
        sess = RSPData(s);
        if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
            rawSMI = sess.SMI.SMI;
            cleanSMI = rawSMI(sess.FilteredROIs);
            cleanSMI = cleanSMI(~isnan(cleanSMI));
            pooledSMI_RSP = [pooledSMI_RSP; cleanSMI(:)];
        end
    end
    
    pooledSMI_VISp = [];
    for s = 1:length(VISpData)
        sess = VISpData(s);
        if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
            rawSMI = sess.SMI.SMI;
            cleanSMI = rawSMI(sess.FilteredROIs);
            cleanSMI = cleanSMI(~isnan(cleanSMI));
            pooledSMI_VISp = [pooledSMI_VISp; cleanSMI(:)];
        end
    end
    
    if isempty(pooledSMI_RSP) || isempty(pooledSMI_VISp)
        error('one or both brain regions do not contain any valid filtered data.');
    end
    
    % Calculate Medians and Median Absolute Deviations
    median_RSP  = median(pooledSMI_RSP);
    median_VISp = median(pooledSMI_VISp);
    mad_RSP     = mad(pooledSMI_RSP, 1);   % flag = 1 -> median absolute deviation
    mad_VISp    = mad(pooledSMI_VISp, 1);
    
    % ADDED: 25th/50th/75th percentiles for each region
    pct_RSP  = prctile(pooledSMI_RSP,  [25, 50, 75]);
    pct_VISp = prctile(pooledSMI_VISp, [25, 50, 75]);
    
    [pVal, hStat] = ranksum(pooledSMI_RSP, pooledSMI_VISp);
    fprintf('RSP Boutons pooled ROIs:  %d (median SMI: %.3f +/- %.3f MAD)\n', length(pooledSMI_RSP), median_RSP, mad_RSP);
    fprintf('VISp Somas pooled ROIs:   %d (median SMI: %.3f +/- %.3f MAD)\n', length(pooledSMI_VISp), median_VISp, mad_VISp);
    fprintf('p-value (rank-sum test):  %.4e\n', pVal);
    if hStat
        disp('result: significantly different distributions.');
    else
        disp('result: no significant difference found.');
    end
    
    % plotting 
    figHandle = figure('Name', 'Spatial Modulation Index: RSP vs VISp Comparison', ...
                       'Color', [1 1 1], 'Position', [100 100 700 700]);
    
    rspColor  = 'k'; 
    vispColor = [0.6 0.6 0.6]; 
    
    %% --- Subplot 1: ECDF ---
    subplot(1, 2, 1); hold on;
    
    [fRSP, xRSP] = ecdf(pooledSMI_RSP);
    [fVISp, xVISp] = ecdf(pooledSMI_VISp);
    
    plot(xRSP, fRSP, 'LineWidth',2, 'Color', rspColor, ...
         'DisplayName', sprintf('rsp boutons (n=%d)', length(pooledSMI_RSP)));
    plot(xVISp, fVISp, 'LineWidth', 2, 'Color', vispColor, ...
         'DisplayName', sprintf('visp somas (n=%d)', length(pooledSMI_VISp)));
     
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    
    % ADDED: mark 25th/50th/75th percentiles directly on the ECDF
    plot(pct_RSP,  [0.25, 0.5, 0.75], 'x', 'Color', rspColor,  'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(pct_VISp, [0.25, 0.5, 0.75], 'x', 'Color', vispColor, 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    
    xlabel('Spatial modulation index');
    ylabel('Cumul. Probability');
    title(sprintf('Distribution of SMI \n(rank-sum p = %.4e)', pVal));
    
    xlim([-1.1, 1.1]);
    ylim([0, 1.02]);
    legend('Location', 'best', 'Interpreter', 'none', 'Box','off');
    
    yticks([0, 1]); 
    
    set(gca, 'Box', 'off'); 
    defaultAxesProperties(gca)
    offsetAxes(gca)
    axis square; 
    
    %% --- Subplot 2: Back-to-back filled histogram (RSP up, VISp flipped down) ---
    subplot(1, 2, 2); hold on;
    
    binEdges = linspace(-1.1, 1.1, 55);
    binCentres = binEdges(1:end-1) + diff(binEdges)/2;
    
    [counts_RSP, ~]  = histcounts(pooledSMI_RSP, binEdges, 'Normalization', 'probability');
    [counts_VISp, ~] = histcounts(pooledSMI_VISp, binEdges, 'Normalization', 'probability');
    
    combined_counts = counts_RSP + counts_VISp;
    first_idx = find(combined_counts > 0, 1, 'first');
    last_idx  = find(combined_counts > 0, 1, 'last');
    
    cropped_centres = binCentres(first_idx : last_idx);
    cropped_RSP     = counts_RSP(first_idx : last_idx);
    cropped_VISp    = -counts_VISp(first_idx : last_idx);   % negated so it plots below the axis
    binWidth        = diff(binEdges(1:2));
    
    % Filled bars: RSP above zero, VISp flipped below zero
    h1 = bar(cropped_centres, cropped_RSP, 1, 'FaceColor', rspColor, 'EdgeColor', 'none', ...
             'FaceAlpha', 0.85, 'DisplayName', 'rsp boutons');
    h2 = bar(cropped_centres, cropped_VISp, 1, 'FaceColor', vispColor, 'EdgeColor', 'none', ...
             'FaceAlpha', 0.85, 'DisplayName', 'visp somas');
    
    yline(0, 'Color', [0 0 0], 'LineWidth', 1, 'HandleVisibility', 'off');
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    
    xlabel('Spatial modulation index');
    ylabel('Probability (VISp flipped below zero)');
    
    xlim([-1, 1]);
    
    max_density = max([max(counts_RSP), max(counts_VISp)]);
    y_max = ceil(max_density * 5) / 5; 
    y_max = 0.15;
    ylim([-y_max, y_max]); 
    yticks([-y_max, 0, y_max]);
    yticklabels({num2str(y_max), '0', num2str(y_max)});   % show positive labels on both sides
    
    % ADDED: crosses marking 25th/50th/75th percentiles on each histogram.
    % RSP crosses sit above zero (its own baseline); VISp crosses sit
    % below zero, at the NEGATIVE percentile value, matching how its bars
    % are drawn negated.
    plot(pct_RSP,  repmat(y_max * 0.05, 1, 3), 'x', 'Color', rspColor, ...
        'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(pct_VISp, repmat(-y_max * 0.05, 1, 3), 'x', 'Color', vispColor, ...
        'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    
    text_x = cropped_centres(end) * 0.45; 
    text(text_x, y_max * 0.85, sprintf('RSP Med: %.3f +/- %.3f', median_RSP, mad_RSP), ...
         'Color', 'k', 'FontWeight', 'bold', 'FontSize', 9);
    text(text_x, -y_max * 0.85, sprintf('VISp Med: %.3f +/- %.3f', median_VISp, mad_VISp), ...
         'Color', [0.4 0.4 0.4], 'FontWeight', 'bold', 'FontSize', 9);
    
    legend([h1, h2], 'Location', 'best', 'Box','off');
    
    set(gca, 'Box', 'off');
    defaultAxesProperties(gca)
    offsetAxes(gca)
    axis square;
    
    
    %% --- Save Block ---
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    baseFileName = 'rsp_vs_vispdFFNeu_smi_comparison';
    fullSavePath = fullfile(outputDir, baseFileName);
    
    saveFigureFormats(figHandle, fullSavePath);

    printSMISummaryStats(RSPData, VISpData)
end

function printSMISummaryStats(RSPData, VISpData)
% printSMISummaryStats: Prints median SMI, median absolute deviation (MAD),
% a one-sample Wilcoxon signed-rank test against zero, a cross-region
% Wilcoxon rank-sum test, and the percentage (and count) of cells with
% SMI > 0.25 and SMI > 0.5, for both RSP boutons and VISp somata.

    [rspSMI, rspN]   = poolSMIValues(RSPData);
    [vispSMI, vispN] = poolSMIValues(VISpData);

    if isempty(rspSMI) || isempty(vispSMI)
        error('One or both regions returned no valid SMI values.');
    end

    % RSP stats
    medianRSP = median(rspSMI);
    madRSP    = mad(rspSMI, 1); % flag = 1 -> median absolute deviation
    nAboveRSP_025 = sum(rspSMI > 0.25);
    pctAboveRSP_025 = (nAboveRSP_025 / length(rspSMI)) * 100;
    nAboveRSP_05 = sum(rspSMI > 0.5);
    pctAboveRSP_05 = (nAboveRSP_05 / length(rspSMI)) * 100;
    pRSP_vs_zero = signrank(rspSMI, 0, 'tail', 'right');

    % VISp stats
    medianVISp = median(vispSMI);
    madVISp    = mad(vispSMI, 1);
    nAboveVISp_025 = sum(vispSMI > 0.25);
    pctAboveVISp_025 = (nAboveVISp_025 / length(vispSMI)) * 100;
    nAboveVISp_05 = sum(vispSMI > 0.5);
    pctAboveVISp_05 = (nAboveVISp_05 / length(vispSMI)) * 100;
    pVISp_vs_zero = signrank(vispSMI, 0, 'tail', 'right');

    % Cross-region comparison
    pRankSum = ranksum(rspSMI, vispSMI, 'tail', 'left');

    fprintf('\n=== SMI Summary Statistics ===\n');
    fprintf('RSP boutons:  n = %d cells\n', rspN);
    fprintf('  Median SMI = %.3f +/- %.3f MAD\n', medianRSP, madRSP);
    fprintf('  Signed-rank test vs 0 (right-tailed): p = %.3e\n', pRSP_vs_zero);
    fprintf('  %d/%d cells (%.1f%%) had SMI > 0.25\n', nAboveRSP_025, length(rspSMI), pctAboveRSP_025);
    fprintf('  %d/%d cells (%.1f%%) had SMI > 0.5\n', nAboveRSP_05, length(rspSMI), pctAboveRSP_05);

    fprintf('\nVISp somas:  n = %d cells\n', vispN);
    fprintf('  Median SMI = %.3f +/- %.3f MAD\n', medianVISp, madVISp);
    fprintf('  Signed-rank test vs 0 (right-tailed): p = %.3e\n', pVISp_vs_zero);
    fprintf('  %d/%d cells (%.1f%%) had SMI > 0.25\n', nAboveVISp_025, length(vispSMI), pctAboveVISp_025);
    fprintf('  %d/%d cells (%.1f%%) had SMI > 0.5\n', nAboveVISp_05, length(vispSMI), pctAboveVISp_05);

    fprintf('\nCross-region comparison:\n');
    fprintf('  Wilcoxon rank-sum test (RSP vs VISp SMI, left-tailed): p = %.3e\n', pRankSum);
    fprintf('\n');
end

%% Pool SMI values (with NaNs removed) across all sessions for a region
function [pooledSMI, nCells] = poolSMIValues(RegionData)
    pooledSMI = [];
    for s = 1:length(RegionData)
        sess = RegionData(s);
        if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
            rawSMI = sess.SMI.SMI;
            cleanSMI = rawSMI(sess.FilteredROIs);
            cleanSMI = cleanSMI(~isnan(cleanSMI));
            pooledSMI = [pooledSMI; cleanSMI(:)]; %#ok<AGROW>
        end
    end
    nCells = length(pooledSMI);
end