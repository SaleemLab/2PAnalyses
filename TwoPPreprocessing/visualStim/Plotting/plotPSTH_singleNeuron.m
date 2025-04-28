function plotPSTH_singleNeuron(psthData, preStimDuration, neuronIdx)
    % Function to plot baseline-normalized PSTH for a selected neuron
    %
    % Inputs:
    %   - psthData: Structure array containing PSTH results per stimulus type
    %   - preStimDuration: Time (in seconds) before stimulus onset for baseline correction
    %   - neuronIdx: Index of the neuron to plot
    %
    
    numStimuli = length(psthData);
    
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
        
        % Create figure
        figure;
        hold on;
        
        % Define the grey shaded stimulus period (e.g., 0-1 sec)
        grey_x = [0 1 1 0];
        grey_y = [min(meanResponse) min(meanResponse) max(meanResponse+100) max(meanResponse+100)];
        fill(grey_x, grey_y, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        
        % Shaded error bars
        x = [timeVector, fliplr(timeVector)];
        y = [meanResponse - semResponse, fliplr(meanResponse + semResponse)];
        fill(x, y, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        
        % Plot mean response line
        plot(timeVector, meanResponse, 'k', 'LineWidth', 1.5);
        
        % Labels and title
        title(sprintf('Neuron %d - Stimulus %d', neuronIdx, psthData(thisStim).stimValue));
        xlabel('Time (s)');
        ylabel('Neural Response (Baseline Normalized)');
        grid on;
        hold off;
    end
end
