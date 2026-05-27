function plotSpatialModulationSummary(allData, sessionIdx)
    % plotSpatialModulationSummary: Visualizes population tiling and SMI distribution
    % for the Baseline condition of a specific session.
    
    if sessionIdx > length(allData)
        error('Session index out of bounds.');
    end
    
    sess = allData(sessionIdx);
    if ~isfield(sess.ConditionData, 'Baseline') || ~isfield(sess.ConditionData.Baseline, 'SMI')
        error('Baseline SMI data not found. Run computeSMIForBaseline first.');
    end
    
    % Extract Data
    data = sess.ConditionData.Baseline;
    activity = data.Activity; % [Cells x Laps x Position]
    smi = data.SMI;
    rho = sess.RawStabilityRho; % Original stability scores
    
    % Calculate Mean Tuning for Heatmap
    meanTuning = squeeze(mean(activity, 2, 'omitnan'));
    normTuning = meanTuning ./ max(meanTuning, [], 2);
    
    % Order by peak position for the heatmap
    [~, peakPos] = max(normTuning, [], 2);
    [~, sortIdx] = sort(peakPos);
    sortedTuning = normTuning(sortIdx, :);
    
    figure('Color', 'w', 'Position', [100 100 1200 400]);
    
    % --- Panel 1: Population Heatmap (Tiling) ---
    subplot(1, 3, 1);
    imagesc(sortedTuning);
    colormap(jet);
    hold on;
    % Draw vertical lines at 40cm repeats [cite: 96]
    for x = 40:40:160
        line([x x], [1 size(sortedTuning,1)], 'Color', 'w', 'LineStyle', '--');
    end
    xlabel('Position (cm)');
    ylabel('Stable Neurons (Sorted)');
    title(sprintf('V1 Tiling (Day %d)', sess.Day));
    
    % --- Panel 2: SMI Cumulative Distribution ---
    subplot(1, 3, 2);
    validSMI = smi(~isnan(smi));
    [f, x] = ecdf(validSMI);
    plot(x, f, 'k', 'LineWidth', 2);
    hold on;
    medSMI = median(validSMI);
    line([medSMI medSMI], [0 0.5], 'Color', 'r', 'LineStyle', ':');
    line([-1 medSMI], [0.5 0.5], 'Color', 'r', 'LineStyle', ':');
    xlim([-1 1]); ylim([0 1]); grid on;
    xlabel('Spatial Modulation Index (SMI)');
    ylabel('Cumulative Prob.');
    title(['SMI Distribution (Med: ', num2str(medSMI, '%.2f'), ')']);
    
    % --- Panel 3: Stability vs SMI Correlation ---
    subplot(1, 3, 3);
    % We need to filter rho to match the filtered stable neurons used for SMI
    stableRho = rho(sess.RawStabilityRho >= sess.StabilityThreshold);
    scatter(stableRho, smi, 20, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('Stability (\rho)');
    ylabel('SMI');
    title('Stability vs Spatial Modulation');
    grid on;
end