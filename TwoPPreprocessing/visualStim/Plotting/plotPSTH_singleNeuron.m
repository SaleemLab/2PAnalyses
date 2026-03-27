function plotPSTH_singleNeuron(psthData, preStimDuration, neuronIdx)
    % plotPSTH_singleNeuron plots baseline-normalized PSTH for a selected neuron.
    % 
    % Inputs:
    %   psthData: Structure with fields 'timeVector', 'alignedResponses', 'stimValue'
    %   preStimDuration: Scalar (e.g., 0.5) representing seconds before 0 to use as baseline
    %   neuronIdx: Integer index of the neuron
    
    numStimuli = length(psthData);
    
    for thisStim = 1:numStimuli
        % 1. Setup Data
        timeVector = psthData(thisStim).timeVector;
        % Force timeVector to be a row vector to prevent 'fill' errors
        tRow = timeVector(:)'; 
        
        % Extract responses [Neurons x Time x Trials]
        resp = psthData(thisStim).alignedResponses;
        
        % Safety check for neuron index
        if neuronIdx > size(resp, 1) || neuronIdx < 1
            warning('Neuron %d is out of range for Stimulus %d.', neuronIdx, thisStim);
            continue;
        end
        
        % 2. Baseline Correction
        % Note: We look for time between -preStimDuration and 0
        preStimIndices = tRow >= -abs(preStimDuration) & tRow < 0;
        
        if ~any(preStimIndices)
            warning('No time points found between -%.2f and 0. Check your timeVector!', preStimDuration);
            baselineMean = 0; % Fallback
        else
            % Average across Time (dim 2) and Trials (dim 3)
            baselineMean = nanmean(nanmean(resp(neuronIdx, preStimIndices, :), 2), 3);
        end
        
        % 3. Calculate Mean and SEM (Across Trials)
        % Get data for this neuron, subtract baseline
        neuronData = squeeze(resp(neuronIdx, :, :)) - baselineMean;
        
        meanResponse = nanmean(neuronData, 2)'; % Force to row
        semResponse = (nanstd(neuronData, 0, 2) ./ sqrt(size(neuronData, 2)))'; % Force to row
        
        % 4. Plotting
        figure('Color', 'w', 'Name', sprintf('Neuron %d', neuronIdx));
        hold on;
        
        % Draw Stimulus Shade (Assume stimulus is from 0 to 1 second)
        % We use 'ylim' later to make this box full-height
        yl = [min(meanResponse - semResponse), max(meanResponse + semResponse)];
        fill([0 1 1 0], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        
        % Draw Shaded SEM
        xPath = [tRow, fliplr(tRow)];
        yPath = [(meanResponse - semResponse), fliplr(meanResponse + semResponse)];
        fill(xPath, yPath, [0.7 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        
        % Plot Mean Line
        plot(tRow, meanResponse, 'b', 'LineWidth', 2);
        
        % Formatting
        title(sprintf('Neuron %d | Stimulus: %s', neuronIdx, num2str(psthData(thisStim).stimValue)));
        xlabel('Time (s)');
        ylabel('Baseline Subtracted Activity');
        grid on;
        set(gca, 'TickDir', 'out');
    end
end