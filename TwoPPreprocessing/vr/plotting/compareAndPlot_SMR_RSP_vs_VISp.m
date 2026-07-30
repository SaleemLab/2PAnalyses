function compareAndPlot_SMR_RSP_vs_VISp(RSPData, VISpData, outputDir)
    
    % Same landmark geometry as computeSpatialModulationRatio, used to
    % derive landmark-locked status directly from GlobalPeakBin (which is
    % already saved), rather than requiring a pre-saved IsLandmarkLocked
    % field.
    landmarkCentres = [40, 80, 120, 160];
    tolerance = 10;
    allowedLandmarkBins = [];
    for c = landmarkCentres
        allowedLandmarkBins = [allowedLandmarkBins, (c - tolerance):(c + tolerance)]; %#ok<AGROW>
    end
    allowedLandmarkBins = unique(allowedLandmarkBins);

    pooledSMR_RSP = [];
    for s = 1:length(RSPData)
        sess = RSPData(s);
        if isfield(sess, 'SMR') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
            rawSMR = sess.SMR.SMR;
            totalROIs = numel(rawSMR);
            if islogical(sess.FilteredROIs)
                filteredMask = sess.FilteredROIs(:);
            else
                filteredMask = false(totalROIs, 1);
                filteredMask(sess.FilteredROIs) = true;
            end
            keepIdx = filteredMask;
            if isfield(sess.SMR, 'GlobalPeakBin')
                isLandmarkLocked = ismember(sess.SMR.GlobalPeakBin(:), allowedLandmarkBins);
                keepIdx = keepIdx & isLandmarkLocked;
            end
            cleanSMR = rawSMR(keepIdx);
            cleanSMR = cleanSMR(~isnan(cleanSMR));
            pooledSMR_RSP = [pooledSMR_RSP; cleanSMR(:)];
        end
    end
    
    pooledSMR_VISp = [];
    for s = 1:length(VISpData)
        sess = VISpData(s);
        if isfield(sess, 'SMR') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
            rawSMR = sess.SMR.SMR;
            totalROIs = numel(rawSMR);
            if islogical(sess.FilteredROIs)
                filteredMask = sess.FilteredROIs(:);
            else
                filteredMask = false(totalROIs, 1);
                filteredMask(sess.FilteredROIs) = true;
            end
            keepIdx = filteredMask;
            if isfield(sess.SMR, 'GlobalPeakBin')
                isLandmarkLocked = ismember(sess.SMR.GlobalPeakBin(:), allowedLandmarkBins);
                keepIdx = keepIdx & isLandmarkLocked;
            end
            cleanSMR = rawSMR(keepIdx);
            cleanSMR = cleanSMR(~isnan(cleanSMR));
            pooledSMR_VISp = [pooledSMR_VISp; cleanSMR(:)];
        end
    end
    
    if isempty(pooledSMR_RSP) || isempty(pooledSMR_VISp)
        error('one or both brain regions do not contain any valid filtered data.');
    end
    
    % Calculate Medians and Mean +/- SD (Saleem et al., 2018 report mean +/- SD, e.g. "0.62 +/- 0.26")
    median_RSP  = median(pooledSMR_RSP);
    median_VISp = median(pooledSMR_VISp);
    mean_RSP  = mean(pooledSMR_RSP);
    sd_RSP    = std(pooledSMR_RSP);
    mean_VISp = mean(pooledSMR_VISp);
    sd_VISp   = std(pooledSMR_VISp);
    
    [pVal, hStat] = ranksum(pooledSMR_RSP, pooledSMR_VISp);
    fprintf('RSP Boutons pooled ROIs:  %d (median SMR: %.3f, mean +/- SD: %.3f +/- %.3f)\n', length(pooledSMR_RSP), median_RSP, mean_RSP, sd_RSP);
    fprintf('VISp Somas pooled ROIs:   %d (median SMR: %.3f, mean +/- SD: %.3f +/- %.3f)\n', length(pooledSMR_VISp), median_VISp, mean_VISp, sd_VISp);
    fprintf('p-value (rank-sum test):  %.4e\n', pVal);
    if hStat
        disp('result: significantly different distributions.');
    else
        disp('result: no significant difference found.');
    end
    
    % plotting 
    figHandle = figure('Name', 'Spatial Modulation Ratio: RSP vs VISp Comparison', ...
                       'Color', [1 1 1], 'Position', [100 100 700 700]);
    
    rspColor  = 'k'; 
    vispColor = [0.6 0.6 0.6]; 
    
    %% --- Subplot 1: ECDF ---
    subplot(1, 2, 1); hold on;
    
    [fRSP, xRSP] = ecdf(pooledSMR_RSP);
    [fVISp, xVISp] = ecdf(pooledSMR_VISp);
    
    plot(xRSP, fRSP, 'LineWidth',2, 'Color', rspColor, ...
         'DisplayName', sprintf('rsp boutons (n=%d)', length(pooledSMR_RSP)));
    plot(xVISp, fVISp, 'LineWidth', 2, 'Color', vispColor, ...
         'DisplayName', sprintf('visp somas (n=%d)', length(pooledSMR_VISp)));
     
    xline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    
    xlabel('Spatial modulation ratio');
    ylabel('Cumul. Probability');
    title(sprintf('Distribution of SMR \n(rank-sum p = %.4e)', pVal));
    
    % NOTE: underlying data (pooledSMR_*) stays linear -- only the AXIS
    % DISPLAY is log-scaled. Stats (mean/median/SD/ranksum above) are all
    % computed on the untransformed linear values, matching Saleem et al.
    % A log-scaled axis cannot show values <= 0, so those points (if any)
    % simply won't render on this subplot; they are NOT removed from the
    % underlying data or from the stats.
    set(gca, 'XScale', 'log');
    % Show log2-style ticks (powers of 2) on this log-scaled axis. The
    % axis itself is still a standard log (base-10) scale under the
    % hood -- log spacing looks identical regardless of base, since the
    % on-screen distance between x and 2x is fixed either way. This just
    % relabels/repositions the ticks at powers of 2 rather than 1/10/100.
    dataMinRSP = min([pooledSMR_RSP(pooledSMR_RSP>0); pooledSMR_VISp(pooledSMR_VISp>0)]);
    dataMaxRSP = max([pooledSMR_RSP; pooledSMR_VISp]);
    lowPow  = floor(log2(dataMinRSP));
    highPow = ceil(log2(dataMaxRSP));
    pow2Ticks = 2.^(lowPow:highPow);
    xticks(pow2Ticks);
    xticklabels(arrayfun(@(v) sprintf('%g', v), pow2Ticks, 'UniformOutput', false));
    ylim([0, 1.02]);
    legend('Location', 'best', 'Interpreter', 'none', 'Box','off');
    
    yticks([0, 1]); 
    
    set(gca, 'Box', 'off'); 
    defaultAxesProperties(gca)
    offsetAxes(gca)
    axis square; 
    
    %% --- Subplot 2: Histogram ---
    subplot(1, 2, 2); hold on;
    
    % NOTE: log-spaced bins (via logspace) so bin widths look even once
    % the x-axis is displayed on a log scale below. Values <= 0 can't be
    % placed on a log axis and are excluded from this histogram only --
    % they remain in pooledSMR_RSP/VISp for all stats above.
    allSMR = [pooledSMR_RSP; pooledSMR_VISp];
    posSMR = allSMR(allSMR > 0);
    nNonPositive = sum(allSMR <= 0);
    if nNonPositive > 0
        fprintf('Note: %d/%d pooled ROIs have SMR <= 0 and are omitted from the log-axis histogram only (not from stats).\n', nNonPositive, length(allSMR));
    end
    binEdges = 2.^linspace(log2(min(posSMR)), log2(max(posSMR)), 100);
    
    [counts_RSP, ~]  = histcounts(pooledSMR_RSP, binEdges, 'Normalization', 'pdf');
    [counts_VISp, ~] = histcounts(pooledSMR_VISp, binEdges, 'Normalization', 'pdf');
    
    combined_counts = counts_RSP + counts_VISp;
    first_idx = find(combined_counts > 0, 1, 'first');
    last_idx  = find(combined_counts > 0, 1, 'last');
    
    cropped_edges  = binEdges(first_idx : last_idx+1);
    cropped_RSP    = [counts_RSP(first_idx : last_idx), 0];   
    cropped_VISp   = [counts_VISp(first_idx : last_idx), 0];  
    
    h1 = stairs(cropped_edges, cropped_RSP, 'Color', [0 0 0], 'LineWidth', 1.2, ...
                'DisplayName', 'rsp boutons');
    h2 = stairs(cropped_edges, cropped_VISp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
                'DisplayName', 'visp somas');
          
    xline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlabel('Spatial modulation ratio');
    ylabel('Probability');
    set(gca, 'XScale', 'log');
    lowPowHist  = floor(log2(min(posSMR)));
    highPowHist = ceil(log2(max(posSMR)));
    pow2TicksHist = 2.^(lowPowHist:highPowHist);
    xticks(pow2TicksHist);
    xticklabels(arrayfun(@(v) sprintf('%g', v), pow2TicksHist, 'UniformOutput', false));
    
    max_density = max([max(counts_RSP), max(counts_VISp)]);
    y_max = ceil(max_density * 5) / 5; 
    ylim([0, y_max]); 
    yticks([0, y_max/2, y_max]); 
    
    % --- PRINT MEDIANS ON THE FIGURE ---
    % Places text in the upper-right sector of the plot coordinates
    text_x = cropped_edges(end) * 0.45; 
    text(text_x, y_max * 0.85, sprintf('RSP Med: %.3f', median_RSP), ...
         'Color', 'k', 'FontWeight', 'bold', 'FontSize', 9);
    text(text_x, y_max * 0.75, sprintf('VISp Med: %.3f', median_VISp), ...
         'Color', [0.4 0.4 0.4], 'FontWeight', 'bold', 'FontSize', 9);
    % 
    
    legend('Location', 'best', 'Box','off');
    
    set(gca, 'Box', 'off');
    defaultAxesProperties(gca)
    offsetAxes(gca)
    axis square;
    
    
    %% --- Save Block ---
   
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    baseFileName = 'rsp_vs_vispdFFNeu_smr_comparison';
    fullSavePath = fullfile(outputDir, baseFileName);
    
    saveFigureFormats(figHandle, fullSavePath);
end
% function compareAndPlot_SMR_RSP_vs_VISp(RSPData, VISpData, outputDir)
%     
%     pooledSMR_RSP = [];
%     for s = 1:length(RSPData)
%         sess = RSPData(s);
%         if isfield(sess, 'SMR') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
%             rawSMR = sess.SMR.SMR;
%             keepIdx = sess.FilteredROIs;
%             if isfield(sess.SMR, 'IsLandmarkLocked')
%                 keepIdx = keepIdx & sess.SMR.IsLandmarkLocked(:);
%             end
%             cleanSMR = rawSMR(keepIdx);
%             cleanSMR = cleanSMR(~isnan(cleanSMR));
%             pooledSMR_RSP = [pooledSMR_RSP; cleanSMR(:)];
%         end
%     end
%     
%     pooledSMR_VISp = [];
%     for s = 1:length(VISpData)
%         sess = VISpData(s);
%         if isfield(sess, 'SMR') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
%             rawSMR = sess.SMR.SMR;
%             keepIdx = sess.FilteredROIs;
%             if isfield(sess.SMR, 'IsLandmarkLocked')
%                 keepIdx = keepIdx & sess.SMR.IsLandmarkLocked(:);
%             end
%             cleanSMR = rawSMR(keepIdx);
%             cleanSMR = cleanSMR(~isnan(cleanSMR));
%             pooledSMR_VISp = [pooledSMR_VISp; cleanSMR(:)];
%         end
%     end
%     
%     if isempty(pooledSMR_RSP) || isempty(pooledSMR_VISp)
%         error('one or both brain regions do not contain any valid filtered data.');
%     end
%     
%     % Calculate Medians and Mean +/- SD (Saleem et al., 2018 report mean +/- SD, e.g. "0.62 +/- 0.26")
%     median_RSP  = median(pooledSMR_RSP);
%     median_VISp = median(pooledSMR_VISp);
%     mean_RSP  = mean(pooledSMR_RSP);
%     sd_RSP    = std(pooledSMR_RSP);
%     mean_VISp = mean(pooledSMR_VISp);
%     sd_VISp   = std(pooledSMR_VISp);
%     
%     [pVal, hStat] = ranksum(pooledSMR_RSP, pooledSMR_VISp);
%     fprintf('RSP Boutons pooled ROIs:  %d (median SMR: %.3f, mean +/- SD: %.3f +/- %.3f)\n', length(pooledSMR_RSP), median_RSP, mean_RSP, sd_RSP);
%     fprintf('VISp Somas pooled ROIs:   %d (median SMR: %.3f, mean +/- SD: %.3f +/- %.3f)\n', length(pooledSMR_VISp), median_VISp, mean_VISp, sd_VISp);
%     fprintf('p-value (rank-sum test):  %.4e\n', pVal);
%     if hStat
%         disp('result: significantly different distributions.');
%     else
%         disp('result: no significant difference found.');
%     end
%     
%     % plotting 
%     figHandle = figure('Name', 'Spatial Modulation Ratio: RSP vs VISp Comparison', ...
%                        'Color', [1 1 1], 'Position', [100 100 700 700]);
%     
%     rspColor  = 'k'; 
%     vispColor = [0.6 0.6 0.6]; 
%     
%     %% --- Subplot 1: ECDF ---
%     subplot(1, 2, 1); hold on;
%     
%     [fRSP, xRSP] = ecdf(pooledSMR_RSP);
%     [fVISp, xVISp] = ecdf(pooledSMR_VISp);
%     
%     plot(xRSP, fRSP, 'LineWidth',2, 'Color', rspColor, ...
%          'DisplayName', sprintf('rsp boutons (n=%d)', length(pooledSMR_RSP)));
%     plot(xVISp, fVISp, 'LineWidth', 2, 'Color', vispColor, ...
%          'DisplayName', sprintf('visp somas (n=%d)', length(pooledSMR_VISp)));
%      
%     xline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
%     yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
%     
%     xlabel('Spatial modulation ratio');
%     ylabel('Cumul. Probability');
%     title(sprintf('Distribution of SMR \n(rank-sum p = %.4e)', pVal));
%     
%     % NOTE: underlying data (pooledSMR_*) stays linear -- only the AXIS
%     % DISPLAY is log-scaled. Stats (mean/median/SD/ranksum above) are all
%     % computed on the untransformed linear values, matching Saleem et al.
%     % A log-scaled axis cannot show values <= 0, so those points (if any)
%     % simply won't render on this subplot; they are NOT removed from the
%     % underlying data or from the stats.
%     set(gca, 'XScale', 'log');
%     % Show log2-style ticks (powers of 2) on this log-scaled axis. The
%     % axis itself is still a standard log (base-10) scale under the
%     % hood -- log spacing looks identical regardless of base, since the
%     % on-screen distance between x and 2x is fixed either way. This just
%     % relabels/repositions the ticks at powers of 2 rather than 1/10/100.
%     dataMinRSP = min([pooledSMR_RSP(pooledSMR_RSP>0); pooledSMR_VISp(pooledSMR_VISp>0)]);
%     dataMaxRSP = max([pooledSMR_RSP; pooledSMR_VISp]);
%     lowPow  = floor(log2(dataMinRSP));
%     highPow = ceil(log2(dataMaxRSP));
%     pow2Ticks = 2.^(lowPow:highPow);
%     xticks(pow2Ticks);
%     xticklabels(arrayfun(@(v) sprintf('%g', v), pow2Ticks, 'UniformOutput', false));
%     ylim([0, 1.02]);
%     legend('Location', 'best', 'Interpreter', 'none', 'Box','off');
%     
%     yticks([0, 1]); 
%     
%     set(gca, 'Box', 'off'); 
%     defaultAxesProperties(gca)
%     offsetAxes(gca)
%     axis square; 
%     
%     %% --- Subplot 2: Histogram ---
%     subplot(1, 2, 2); hold on;
%     
%     % NOTE: log-spaced bins (via logspace) so bin widths look even once
%     % the x-axis is displayed on a log scale below. Values <= 0 can't be
%     % placed on a log axis and are excluded from this histogram only --
%     % they remain in pooledSMR_RSP/VISp for all stats above.
%     allSMR = [pooledSMR_RSP; pooledSMR_VISp];
%     posSMR = allSMR(allSMR > 0);
%     nNonPositive = sum(allSMR <= 0);
%     if nNonPositive > 0
%         fprintf('Note: %d/%d pooled ROIs have SMR <= 0 and are omitted from the log-axis histogram only (not from stats).\n', nNonPositive, length(allSMR));
%     end
%     binEdges = 2.^linspace(log2(min(posSMR)), log2(max(posSMR)), 100);
%     
%     [counts_RSP, ~]  = histcounts(pooledSMR_RSP, binEdges, 'Normalization', 'pdf');
%     [counts_VISp, ~] = histcounts(pooledSMR_VISp, binEdges, 'Normalization', 'pdf');
%     
%     combined_counts = counts_RSP + counts_VISp;
%     first_idx = find(combined_counts > 0, 1, 'first');
%     last_idx  = find(combined_counts > 0, 1, 'last');
%     
%     cropped_edges  = binEdges(first_idx : last_idx+1);
%     cropped_RSP    = [counts_RSP(first_idx : last_idx), 0];   
%     cropped_VISp   = [counts_VISp(first_idx : last_idx), 0];  
%     
%     h1 = stairs(cropped_edges, cropped_RSP, 'Color', [0 0 0], 'LineWidth', 1.2, ...
%                 'DisplayName', 'rsp boutons');
%     h2 = stairs(cropped_edges, cropped_VISp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
%                 'DisplayName', 'visp somas');
%           
%     xline(1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
%     xlabel('Spatial modulation ratio');
%     ylabel('Probability');
%     set(gca, 'XScale', 'log');
%     lowPowHist  = floor(log2(min(posSMR)));
%     highPowHist = ceil(log2(max(posSMR)));
%     pow2TicksHist = 2.^(lowPowHist:highPowHist);
%     xticks(pow2TicksHist);
%     xticklabels(arrayfun(@(v) sprintf('%g', v), pow2TicksHist, 'UniformOutput', false));
%     
%     max_density = max([max(counts_RSP), max(counts_VISp)]);
%     y_max = ceil(max_density * 5) / 5; 
%     ylim([0, y_max]); 
%     yticks([0, y_max/2, y_max]); 
%     
%     % --- PRINT MEDIANS ON THE FIGURE ---
%     % Places text in the upper-right sector of the plot coordinates
%     text_x = cropped_edges(end) * 0.45; 
%     text(text_x, y_max * 0.85, sprintf('RSP Med: %.3f', median_RSP), ...
%          'Color', 'k', 'FontWeight', 'bold', 'FontSize', 9);
%     text(text_x, y_max * 0.75, sprintf('VISp Med: %.3f', median_VISp), ...
%          'Color', [0.4 0.4 0.4], 'FontWeight', 'bold', 'FontSize', 9);
%     % 
%     
%     legend('Location', 'best', 'Box','off');
%     
%     set(gca, 'Box', 'off');
%     defaultAxesProperties(gca)
%     offsetAxes(gca)
%     axis square;
%     
%     
%     %% --- Save Block ---
%    
%     if ~exist(outputDir, 'dir')
%         mkdir(outputDir);
%     end
%     
%     baseFileName = 'rsp_vs_vispdFFNeu_smr_comparison';
%     fullSavePath = fullfile(outputDir, baseFileName);
%     
%     saveFigureFormats(figHandle, fullSavePath);
% end

