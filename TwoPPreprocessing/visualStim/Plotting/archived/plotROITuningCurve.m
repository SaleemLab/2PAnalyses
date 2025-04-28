function plotROITuningCurve(psthData, neuronIndex)
    % plotNeuronTuningCurve - Plots a neuron's tuning curve across all stimulus types.
    %
    % Inputs:
    %   - psthData: Structure containing PSTH responses per stimulus condition.
    %   - neuronIndex: Index of the neuron to analyze.
    %
    % The function extracts the mean response of the given neuron for each stimulus
    % condition, plots its tuning curve, and includes the time vector for reference.

    % Number of stimulus conditions
    numStimuli = length(psthData);
    
    % Initialize arrays for tuning curve
    stimValues = zeros(numStimuli, 1);  % Store stimulus values
    neuronResponses = zeros(numStimuli, 1);  % Store mean responses
    neuronSEM = zeros(numStimuli, 1);  % Store SEM
    timeVector = psthData(1).timeVector;  % Extract time vector from any condition (assumed same across)

    % Extract mean response for each stimulus type
    for i = 1:numStimuli
        stimValues(i) = psthData(i).stimValue;  % Get stimulus value
        meanResponse = psthData(i).meanResponse(neuronIndex, :);  % Extract mean response for this neuron
        neuronResponses(i) = nanmean(meanResponse);  % Compute mean across time
        stdResponse = nanstd(meanResponse);  % Compute standard deviation
        neuronSEM(i) = stdResponse / sqrt(sum(~isnan(meanResponse)));  % Compute SEM
    end

    % Sort by stimulus values (in case they are unordered)
    [stimValues, sortIdx] = sort(stimValues);
    neuronResponses = neuronResponses(sortIdx);
    neuronSEM = neuronSEM(sortIdx);

    % Create figure
    figure;
    
    % Plot Tuning Curve (Mean Response vs Stimulus Type)
    subplot(2, 1, 1);  % First subplot for tuning curve
    hold on;
    errorbar(stimValues, neuronResponses, neuronSEM, '-o', 'MarkerFaceColor', 'b', ...
        'MarkerEdgeColor', 'k', 'Color', 'b', 'LineWidth', 2);
    xlabel('Stimulus Type');
    ylabel('Mean Neural Response');
    title(sprintf('Tuning Curve for Neuron %d', neuronIndex));
    legend({'Mean ± SEM'});
    grid on;
    hold off;
    
    % Plot Time Course of Mean Response for Each Stimulus
    subplot(2, 1, 2);  % Second subplot for time series response
    hold on;
    cmap = lines(numStimuli); % Generate different colors for each stimulus
    for i = 1:numStimuli
        plot(timeVector, psthData(i).meanResponse(neuronIndex, :), 'Color', cmap(i, :), 'LineWidth', 2);
    end
    xlabel('Time (s)');
    ylabel('Response');
    title(sprintf('Time Course of Neuron %d Across Stimuli', neuronIndex));
    legend(arrayfun(@(x) sprintf('%d', x), stimValues, 'UniformOutput', false));
    grid on;
    hold off;
end
