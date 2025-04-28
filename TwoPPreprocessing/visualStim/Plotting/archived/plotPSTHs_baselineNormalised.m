function plotPSTHs_baselineNormalised(psthData, preStimDuration, neuronIdx)
    % Function to plot baseline-normalized PSTHs for a single neuron 
    % with different stimuli superimposed on the same plot.
    %
    % Inputs:
    %   - psthData: Structure array containing PSTH results for different stimuli
    %   - preStimDuration: Time (in seconds) before stimulus onset for baseline correction
    %   - neuronIdx: Index of the neuron to plot
    %
    
    numStimuli = length(psthData);
    figure;
    hold on;

    colors = lines(numStimuli); % Generate distinct colors for different stimuli
    
    allMeans = []; % Store all means to adjust y-axis
    
    for thisStim = 1:numStimuli
        timeVector = psthData(thisStim).timeVector;
        
        % Find indices corresponding to the pre-stimulus period
        preStimIndices = timeVector >= -preStimDuration & timeVector < 0;
        
        % Extract responses
        alignedResponses = psthData(thisStim).alignedResponses;
        numNeurons = size(alignedResponses, 1);
        
        % Check if neuron index is valid
        if neuronIdx > numNeurons || neuronIdx < 1
            warning('Neuron index %d is out of range for stimulus %d. Skipping...', neuronIdx, psthData(thisStim).stimValue);
            continue;
        end
        
        % Compute baseline: Mean response across the pre-stimulus period
        baselineMean = nanmean(alignedResponses(neuronIdx, preStimIndices, :), [2, 3]); % Across trials
        
        % Normalize response by subtracting baseline
        alignedResponses(neuronIdx, :, :) = alignedResponses(neuronIdx, :, :) - baselineMean;
        
        % Compute mean and SEM for this neuron
        meanResponse = squeeze(nanmean(alignedResponses(neuronIdx, :, :), 3)); % Average across trials
        semResponse = squeeze(nanstd(alignedResponses(neuronIdx, :, :), 0, 3)) ./ sqrt(size(alignedResponses, 3));
        
        % Store mean responses for dynamic Y-limits
        allMeans = [allMeans; meanResponse];

        % Shaded error bars
        x = [timeVector, fliplr(timeVector)];
        y = [meanResponse - semResponse, fliplr(meanResponse + semResponse)];
        fill(x, y, colors(thisStim, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        
        % Plot mean response line
        plot(timeVector, meanResponse, 'Color', colors(thisStim, :), 'LineWidth', 4, ...
             'DisplayName', sprintf('%d', psthData(thisStim).stimValue));
    end
    
    % Compute Y-axis limits dynamically
    yMin = min(allMeans(:));
    yMax = max(allMeans(:));
    yBuffer = 0.1 * (yMax - yMin); % Add 10% buffer

    % Grey shaded stimulus period (e.g., 0-1 sec)
    grey_x = [0 1 1 0];
    grey_y = [yMin - yBuffer, yMin - yBuffer, yMax + yBuffer, yMax + yBuffer];
    fill(grey_x, grey_y, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Labels and title
    xlabel('Time (s)');
    ylabel('Neural Response (Baseline Normalized)');
    title(sprintf('Superimposed PSTHs for Neuron %d', neuronIdx));
    legend('show');
    hold off;
end
