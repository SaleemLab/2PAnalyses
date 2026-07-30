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
    
    % Calculate Medians
    median_RSP  = median(pooledSMI_RSP);
    median_VISp = median(pooledSMI_VISp);
    
    [pVal, hStat] = ranksum(pooledSMI_RSP, pooledSMI_VISp);
    fprintf('RSP Boutons pooled ROIs:  %d (median SMI: %.3f)\n', length(pooledSMI_RSP), median_RSP);
    fprintf('VISp Somas pooled ROIs:   %d (median SMI: %.3f)\n', length(pooledSMI_VISp), median_VISp);
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
    
    %% --- Subplot 2: Histogram ---
    subplot(1, 2, 2); hold on;
    
    binEdges = linspace(-1.1, 1.1, 55);
    
    [counts_RSP, ~]  = histcounts(pooledSMI_RSP, binEdges, 'Normalization', 'pdf');
    [counts_VISp, ~] = histcounts(pooledSMI_VISp, binEdges, 'Normalization', 'pdf');
    
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
          
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlabel('Spatial modulation index');
    ylabel('Probability');
    
    xlim([-1, 1]);
    
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
    
    baseFileName = 'rsp_vs_vispdFFNeu_smi_comparison';
    fullSavePath = fullfile(outputDir, baseFileName);
    
    saveFigureFormats(figHandle, fullSavePath);
end