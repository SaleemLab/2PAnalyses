function filterAndPlotSMIsByEVs(sessionMetrics)
    evValueRanges  = 0:0.1:1;
    numSessions = length(sessionMetrics);
    SMIs = []; % Vector to accumulate all ROIs across sessions
    medianEVs = []; % Vector to accumulate all valid ROIs across sessions
    
    % Extract and pool valid ROIs using the pre-computed FilteredROIs index
    for s = 1:numSessions
        sess = sessionMetrics(s);
        
        % SMI data 
        if ~isfield(sess, 'SMI') || isempty(sess.SMI) || ~isfield(sess.SMI, 'SMI')
            warning('No valid SMI data found for session %d. Skipping.', s);
            continue;
        end
        
        % EV data 
        if ~isfield(sess, 'cvExpVar') || isempty(sess.cvExpVar) || ~isfield(sess.cvExpVar, 'cvExpVar')
            warning('No valid EV data found for session %d. Skipping.', s);
            continue;
        end
      
        % Load smis and evs for all rois recorded
        sessionsSMIValues = sess.SMI.SMI;
        sessionEVValues = sess.cvExpVar.medianExpVar;
        
        % append to pool 
        SMIs = [SMIs; sessionsSMIValues]; 
        medianEVs = [medianEVs; sessionEVValues]; 
    end
    
fprintf('SMIs length: %d, medianEVs length: %d\n', length(SMIs), length(medianEVs));
% filter the medianEVs based on the evValueRangesand plt as superimposed
% cumultiave probability hisrograms as sudo code below:
% Generate Figures
    figure('Name', 'Pooled Spatial Modulation Index Metrics (Filtered)', ...
           'Color', [1 1 1], 'Position', [150 150 500 500]);
    
    hold on;
    
    % 
    numColors = length(evValueRanges) - 1;
    colors = colorcube(numColors + 2); % Add padding to avoid pure white/black if needed
    colors = colors(1:numColors, :);   % Take just the ones you need
    
    % Loop through each bin defined by evValueRanges
    for i = 1:(length(evValueRanges) - 1)
        lowBound = evValueRanges(i);
        highBound = evValueRanges(i+1);
        
        % Filter SMIs falling within the current EV range
        idx = (medianEVs >= lowBound) & (medianEVs < highBound);
        filteredSMI = SMIs(idx);
        nROIs = length(filteredSMI);
        
        % Only plot if there are elements in the bin
        if nROIs > 0
            [f, x] = ecdf(filteredSMI);
            plot(x, f, 'LineWidth', 2, 'Color', colors(i, :), ...
                 'DisplayName', sprintf('EV [%.1f, %.1f) (n=%d)', lowBound, highBound, nROIs));
        end
    end
    
    xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(0.5, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, 'HandleVisibility', 'off');
    
    xlabel('smi value');
    ylabel('cumulative proportion');
    title('pooled cumulative distribution function');
    xlim([-1.1, 1.1]);
    ylim([0, 1.02]);
    legend('Location', 'best', 'Interpreter', 'none');
    
    drawnow;
    
end