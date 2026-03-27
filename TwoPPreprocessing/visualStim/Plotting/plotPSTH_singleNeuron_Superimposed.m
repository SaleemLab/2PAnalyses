function plotPSTH_singleNeuron_Superimposed(psthData, preStimDuration, neuronIdx)
    % Plots all stimuli for a single neuron on one set of axes.
    
    numStimuli = length(psthData);
    
    % Create one figure for the neuron
    figure('Color', 'w', 'Name', ['Neuron ' num2str(neuronIdx)]);
    hold on;
    
    % Generate a set of colors (one for each stimulus)
    colors = lines(numStimuli); 
    
    % Track min/max for the grey stimulus box
    allMin = 0; allMax = 0;

    for thisStim = 1:numStimuli
        tRow = psthData(thisStim).timeVector(:)';
        resp = psthData(thisStim).alignedResponses;
        
        if neuronIdx > size(resp, 1), continue; end
        
        % Baseline corrected
        preStimIndices = tRow >= -abs(preStimDuration) & tRow < 0;
        baselineMean = 0;
        if any(preStimIndices)
            baselineMean = nanmean(nanmean(resp(neuronIdx, preStimIndices, :), 2), 3);
        end
        
        % 2. Calculate Stats
        neuronData = squeeze(resp(neuronIdx, :, :)) - baselineMean;
        meanResponse = nanmean(neuronData, 2)';
        semResponse = (nanstd(neuronData, 0, 2) ./ sqrt(size(neuronData, 2)))';
        
        % Update Y-limits for the stimulus box later
        allMin = min([allMin, meanResponse - semResponse]);
        allMax = max([allMax, meanResponse + semResponse]);

        % 3. Plot Shaded SEM
        xPath = [tRow, fliplr(tRow)];
        yPath = [(meanResponse - semResponse), fliplr(meanResponse + semResponse)];
        fill(xPath, yPath, colors(thisStim, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
        
        % 4. Plot Mean Line (with Label for Legend)
        plot(tRow, meanResponse, 'Color', colors(thisStim, :), 'LineWidth', 2, ...
            'DisplayName', ['Stim: ' num2str(psthData(thisStim).stimValue)]);
    end
    
    % 5. Final Formatting
    % Add grey stimulus box from 0 to 1s
    fill([0 2 2 0], [allMin allMin allMax allMax], [0.9 0.9 0.9], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
    uistack(gca, 'child', 1); % Move grey box to background
    
    title(['Neuron ' num2str(neuronIdx) ' - Comparison Across Stimuli']);
    xlabel('Time (s)');
    ylabel('Baseline Subtracted Activity');
    legend('Location', 'best');
    set(gca, 'TickDir', 'out');
end