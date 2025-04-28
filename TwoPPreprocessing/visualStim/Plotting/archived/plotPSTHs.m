function plotPSTHs(psthData, preStimTime, postStimTime)
    % Function to plot the overall PSTHs for different stimulus conditions
    % with baseline normalization and a specified time window.
    %
    % Inputs:
    %   - psthData: Structure array containing PSTH results per stimulus type
    %   - preStimTime: Time (in seconds) before stimulus onset to include
    %   - postStimTime: Time (in seconds) after stimulus onset to include

    numStimuli = length(psthData);
    figure('Position', [100, 100, 1000, 500]);
    hold on;
    
    colors = lines(numStimuli); % Generate distinct colors for each stimulus
    
    timeVector = psthData(1).timeVector;
    timeIndices = (timeVector >= -preStimTime) & (timeVector <= postStimTime);
    
    for thisStim = 1:numStimuli
        % Extract relevant data
        meanResponse = nanmean(psthData(thisStim).meanResponse, 1);
        semResponse = nanmean(psthData(thisStim).semResponse, 1);
        
        % Baseline normalization
        preStimIndices = timeVector >= -preStimTime & timeVector < 0;
        baselineMean = nanmean(meanResponse(preStimIndices));
        meanResponse = meanResponse - baselineMean;
        
        % Restrict to specified time range
        meanResponse = meanResponse(timeIndices);
        semResponse = semResponse(timeIndices);
        timePlot = timeVector(timeIndices);
        
%         Plot mean response with shaded error bars
        fill([timePlot, fliplr(timePlot)], ...
             [meanResponse - semResponse, fliplr(meanResponse + semResponse)], ...
             colors(thisStim, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        plot(timePlot, meanResponse, 'Color', colors(thisStim, :), 'LineWidth', 2.5, ...
             'DisplayName', sprintf('%d', psthData(thisStim).stimValue));
    end
    
    % Add gray shaded stimulus period (0 to postStimTime)
    yLimits = ylim;
    hFill = fill([0 2 2 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
             [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    set(get(get(hFill, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); % Exclude from legend

    
    hold off;
    xlabel('Time (s)');
    ylabel('Mean raw F');
    title('Overall PSTHs for Direction Tuning');
     legend('show', 'Location', 'northeastoutside');
%     grid on;
end