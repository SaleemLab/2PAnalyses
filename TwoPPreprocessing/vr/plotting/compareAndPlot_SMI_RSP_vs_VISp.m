function compareAndPlot_SMI_RSP_vs_VISp(RSPData, VISpData)

    
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
    
    % 
    [pVal, hStat] = ranksum(pooledSMI_RSP, pooledSMI_VISp);
    fprintf('RSP Boutons pooled ROIs:  %d (median SMI: %.3f)\n', length(pooledSMI_RSP), median(pooledSMI_RSP));
    fprintf('VISp Somas pooled ROIs:   %d (median SMI: %.3f)\n', length(pooledSMI_VISp), median(pooledSMI_VISp));
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

    set(gca, 'Box', 'off'); 
    defaultAxesProperties(gca)
    offsetAxes()
    axis square; 

    % 
    subplot(1, 2, 2); hold on;
    
    binEdges = linspace(-1.1, 1.1, 55);
    

          
    histogram(pooledSMI_VISp, 'BinEdges', binEdges, 'Normalization', 'pdf', ...
        'FaceColor', vispColor, 'EdgeColor', 'w', 'FaceAlpha', 0.3, ...
        'DisplayName', 'VISp somas');
    histogram(pooledSMI_RSP, 'BinEdges', binEdges, 'Normalization', 'pdf', ...
        'FaceColor', rspColor, 'EdgeColor', 'w', 'FaceAlpha', 0.4, ...
        'DisplayName', 'RSP boutons');
%     histogram(pooledSMI_RSP, 'BinEdges', binEdges, 'Normalization', 'pdf', ...
%               'DisplayStyle', 'stairs', 'EdgeColor', [0 0 0], 'LineWidth', 1, ...
%               'DisplayName', 'rsp boutons');
%           
%     % VISp Somas - Dashed Dark Gray Line
%     histogram(pooledSMI_VISp, 'BinEdges', binEdges, 'Normalization', 'pdf', ...
%               'DisplayStyle', 'stairs', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1, ...
%               'LineStyle', '--', 'DisplayName', 'visp somas');
          
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xlabel('Spatial modulation index');
    ylabel('Probability');

    
    
    xlim([-1.1, 1.1]);
   
    legend('Location', 'best', 'Box','off');
    
    % 
    set(gca, 'Box', 'off');
    defaultAxesProperties(gca)
    offsetAxes(gca)
    axis square;
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1\RSPVsVISp\';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    baseFileName = 'rsp_vs_visp_smi_comparison';
    fullSavePath = fullfile(outputDir, baseFileName);
    
    saveFigureFormats(figHandle, fullSavePath);

end