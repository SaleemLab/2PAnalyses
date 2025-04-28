function plotPSTHs(psthData, outputFilename)
    % plotPSTHs - Plot PSTH data from getTrialResponsePSTHs for all neurons and save as PDF.
    %
    % Inputs:
    %   psthData (struct) - Output from getTrialResponsePSTHs
    %   outputFilename (char) - Name of the output PDF file
    %
    % Example Usage:
    %   plotPSTHs(psthData, 'psth_plots.pdf')
    %
    
    numNeurons = size(psthData(1).meanResponse, 1);
    numStimuli = length(psthData);
    
    pdfFilename = outputFilename;
    
    % Create figure for plotting
    hFig = figure('Visible', 'off');
    
    for neuronIdx = 1:numNeurons
        clf(hFig); % Clear figure for each neuron
        hold on;
        
        % Iterate through all stimuli and plot PSTHs
        for stimIdx = 1:numStimuli
            stimValue = psthData(stimIdx).stimValue;
            psthMean = psthData(stimIdx).meanResponse(neuronIdx, :);
            psthSEM = psthData(stimIdx).semResponse(neuronIdx, :);
            timeAxis = psthData(stimIdx).timeVector;
            
            % Define shaded error bars
            x = [timeAxis, fliplr(timeAxis)];
            y = [psthMean - psthSEM, fliplr(psthMean + psthSEM)];
            fill(x, y, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            
            % Plot mean response line
            plot(timeAxis, psthMean, 'LineWidth', 1.5, 'DisplayName', sprintf('Stim %d', stimValue));
        end
        
        title(sprintf('PSTH for Neuron %d', neuronIdx));
        xlabel('Time (s)');
        ylabel('Fluorescence Response');
        legend show;
        grid on;
        hold off;
        
        % Save each figure to PDF
        if neuronIdx == 1
            exportgraphics(hFig, pdfFilename, 'Append', false);
        else
            exportgraphics(hFig, pdfFilename, 'Append', true);
        end
    end
    
    close(hFig); % Close figure after saving
    fprintf('PSTHs saved to %s\n', pdfFilename);
end
