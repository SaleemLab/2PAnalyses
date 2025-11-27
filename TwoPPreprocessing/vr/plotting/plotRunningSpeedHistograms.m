function plotRunningSpeedHistograms(sessionSpeed, lapSpeed, binSpeed)
    
    % number of laps
    nLaps = length(lapSpeed);

    % mean speed for each individual lap
    meanSpeedPerLap = cellfun(@mean, lapSpeed); % applied function to each cell 
    
    figure('Position', [100, 100, 1500, 400]);

    % --- Plot 1: Full Session Speed Distribution ---
    subplot(1, 4, 1);
    histogram(sessionSpeed, 50, 'Normalization', 'probability', 'FaceColor', [0 0.447 0.741]);
    title('Full Session Speed Distribution');
    xlabel('Running Speed (cm/s)');
    ylabel('Probability Density');
    grid on;

    % Combined Lap Speed Distribution ---
    subplot(1, 4, 2);
    if ~isempty(lapSpeed)
        allLapSpeeds = cell2mat(lapSpeed);
        histogram(allLapSpeeds, 50, 'Normalization', 'probability', 'FaceColor', [0.85 0.325 0.098]);
        title( 'Combined Lap Speed Distribution');
        xlabel('Running Speed (cm/s)');
        ylabel('Probability Density');
        grid on;
    else
        text(0.5, 0.5, 'Lap data unavailable', 'HorizontalAlignment', 'center');
    end

    %  Mean Speed Across Laps
    subplot(1, 4, 3);
    if nLaps > 0
        % Plot the mean speed for each lap as a scatter or line plot
        plot(1:nLaps, meanSpeedPerLap, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5);
        
        % Add the overall average as a reference line
        hold on;
        yline(mean(meanSpeedPerLap), 'r--', 'Overall Mean');
        hold off;
        
        title('Mean Running Speed Per Lap');
        xlabel('Lap Number');
        ylabel('Mean Speed (cm/s)');
        xlim([0.5, nLaps + 0.5]);
        grid on;
    else
        text(0.5, 0.5, 'Lap data unavailable', 'HorizontalAlignment', 'center');
    end

    %  Mean Speed Profile Across Position Bins
    subplot(1, 4, 4);
    if ~isempty(binSpeed)
        % Mean speed for each bin, averaged across all laps
        meanSpeedPerBin = mean(binSpeed, 1, 'omitnan');
        
        % Plot the speed profile across the track
        plot(meanSpeedPerBin, 'LineWidth', 2, 'Color', [0.494 0.184 0.556]);
        hold on;
        % Add error bars (SEM)
        semSpeedPerBin = std(binSpeed, 0, 1, 'omitnan') / sqrt(size(binSpeed, 1));
        errorbar(1:size(binSpeed, 2), meanSpeedPerBin, semSpeedPerBin, 'LineStyle', 'none', 'Color', [0.5 0.5 0.5]);
        
        title('🛣️ Mean Speed Profile Across Position Bins');
        xlabel('Position Bin (1 cm bins)');
        ylabel('Mean Running Speed (cm/s)');
        grid on;
    else
        text(0.5, 0.5, 'Binned data unavailable', 'HorizontalAlignment', 'center');
    end
end