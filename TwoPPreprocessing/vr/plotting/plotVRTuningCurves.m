function plotVRTuningCurves(lapPositionActivity, roiToPlot)
    % plotTuningCurves Visualizes neural activity as a function of position.
    %
    %   plotTuningCurves(lap_position_activity)
    %   plotTuningCurves(lap_position_activity, neuron_to_plot)
    %
    %   This function takes a 3D array of neural activity (neurons x laps x
    %   position bins) and generates two plots:
    %   1. A sorted heatmap of the entire neural population's average activity.
    %   2. A line plot of a single neuron's tuning curve if 'neuron_to_plot'
    %      is specified.
    %
    % INPUTS:
    %   lap_position_activity - A 3D array where dimensions are
    %                           (number of neurons, number of laps, number of position bins).
    %   neuron_to_plot        - (Optional) The index of a single neuron you
    %                           wish to visualize in a separate line plot.
    
    %
    % October 2025
    
    %% Input Validation and Data Preparation ---
    if nargin < 1
        error('A 3D activity array (neurons x laps x bins) must be provided.');
    end
    if ndims(lapPositionActivity) ~= 3
        error('Input must be a 3D array of size (neurons x laps x position bins).');
    end

    % Calculate the mean activity across all laps to get the average tuning curve
    % This results in a 2D array: (neurons x position bins)
    mean_activity = squeeze(mean(lapPositionActivity, 2));

    % Get dimensions for plotting
    [num_neurons, num_bins] = size(mean_activity);
    position_axis = 1:num_bins;

    %%Plot Population Heatmap ---
    
    % Create a new figure for the heatmap
    figure('Name', 'Population Activity Heatmap', 'NumberTitle', 'off');

    % Normalize each neuron's activity from 0 to 1 for better visualization
    min_vals = min(mean_activity, [], 2);
    max_vals = max(mean_activity, [], 2);
    range_vals = max_vals - min_vals;
    range_vals(range_vals == 0) = 1; % Avoid division by zero for non-active neurons
    normalized_activity = (mean_activity - min_vals) ./ range_vals;

    % Find the position of peak activity for each neuron
    [~, peak_locations] = max(mean_activity, [], 2);

    % Sort the neurons based on their peak activity location
    [~, sorted_indices] = sort(peak_locations);
    sorted_normalized_activity = normalized_activity(sorted_indices, :);

    % Display the sorted, normalized data as an image
    imagesc(position_axis, 1:num_neurons, sorted_normalized_activity);
    
    % Add labels and a color bar
    colormap('hot'); % 'jet' or 'hot' are good choices
    colorbar;
    xlabel('Position (cm)');
    ylabel('Sorted Neuron ID');
    title(sprintf('Normalized Tuning Curves for %d Neurons', num_neurons));
    
    %%  Plot Single Neuron Tuning Curve 

    % Check if the user specified a neuron to plot
    if nargin > 1 && ~isempty(roiToPlot)
        if roiToPlot > 0 && roiToPlot <= num_neurons
            
            % Create a new figure for the single neuron plot
            figure('Name', ['Tuning Curve for Neuron ' num2str(roiToPlot)], 'NumberTitle', 'off');
            
            % Extract the activity for the specified neuron
            single_neuron_activity = mean_activity(roiToPlot, :);
            
            % Create the line plot
            plot(position_axis, single_neuron_activity, 'b-', 'LineWidth', 2);
            
            % Add labels, title, and formatting
            xlabel('Position (cm)');
            ylabel('Average Activity (\DeltaF/F)');
            title(sprintf('Tuning Curve for Neuron %d', roiToPlot));
            grid on;
            box off;
            xlim([position_axis(1) position_axis(end)]); % Ensure x-axis is tight
            
        else
            warning('Invalid neuron index provided. Must be between 1 and %d.', num_neurons);
        end
    end

end