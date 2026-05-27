function plotFilteredSMIDistributions(sessionMetrics)
% plotFilteredSMIDistributions Pools filtered SMI values across all sessions
% and visualizes the aggregate dataset using a probability density histogram
% and a unified cumulative distribution function (CDF).

    numSessions = length(sessionMetrics);
    pooledSMI = []; % Vector to accumulate all valid ROIs across sessions
    
    % Extract and pool valid ROIs using the pre-computed FilteredROIs index
    for s = 1:numSessions
        sess = sessionMetrics(s);
        
        % Verify SMI data and FilteredROIs exist for this session
        if ~isfield(sess, 'SMI') || isempty(sess.SMI) || ~isfield(sess.SMI, 'SMI')
            warning('No valid SMI metadata found for session %d. Skipping.', s);
            continue;
        end
        if ~isfield(sess, 'FilteredROIs')
            warning('Session %d is missing the ''FilteredROIs'' index field. Skipping.', s);
            continue;
        end
        
        rawSMI = sess.SMI.SMI;
        activeIndices = sess.FilteredROIs;
        
        if isempty(activeIndices)
            continue;
        end
        
        % Filter to  indices and remove any NaNs
        cleanSMI = rawSMI(activeIndices);
        cleanSMI = cleanSMI(~isnan(cleanSMI));
        
        % append to pool 
        pooledSMI = [pooledSMI; cleanSMI(:)]; 
    end
    
    if isempty(pooledSMI)
        error('No valid filtered ROI data was found across any sessions to plot.');
    end
    
%     plotMask = (pooledSMI > -1) & (pooledSMI < 1);
%     plotSMI = pooledSMI(plotMask);

    totalPooledROIs = length(pooledSMI);
    fprintf('Successfully pooled %d filtered ROIs across %d sessions.\n', ...
        totalPooledROIs, numSessions);
    
    % Generate Figures
    figure('Name', 'Pooled Spatial Modulation Index Metrics (Filtered)', ...
           'Color', [1 1 1], 'Position', [150 150 1100 500]);
    
    % Histogram 
    subplot(1, 2, 1);
    % Setting 'Normalization' to 'pdf' generates the probability density profile
    histogram(pooledSMI, 60, 'Normalization', 'pdf', ...
              'FaceColor', [0.200, 0.600, 0.800], 'EdgeColor', 'w', 'FaceAlpha', 0.8);
    
    xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    xlabel('smi value');
    ylabel('probability density');
    title(sprintf('pooled smi density distribution\n(n = %d rois)', totalPooledROIs));
    xlim([-1.1, 1.1]);
    grid on;
    
    %
    subplot(1, 2, 2); hold on;
    
    % Compute cumulative distribution parameters for the aggregate pool
    [f, x] = ecdf(pooledSMI);
    
    plot(x, f, 'LineWidth', 3, 'Color', [0.000, 0.400, 1.000], ...
         'DisplayName', sprintf('pooled population (n=%d)', totalPooledROIs));
    
    xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(0.5, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    xlabel('smi value');
    ylabel('cumulative proportion');
    title('pooled cumulative distribution function');
    xlim([-1.1, 1.1]);
    ylim([0, 1.02]);

    legend('Location', 'best', 'Interpreter', 'none');
    
    drawnow;
end