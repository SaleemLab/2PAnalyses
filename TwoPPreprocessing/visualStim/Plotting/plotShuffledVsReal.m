function plotShuffledVsReal(response, sigData, roiIdx, fname, wheelSpeed, sIdx, binIndices, numShiftsToPlot)
    % Arguments:
    % response: the output struct from your main function
    % sigData: the raw signal matrix [ROIs x Frames] for the current field (e.g., dFF)
    % roiIdx: the index of the ROI you want to inspect
    % fname: string name of the field (e.g., 'spks')
    % wheelSpeed, sIdx, binIndices: indices used during the main function

    if nargin < 8, numShiftsToPlot = 20; end
    
    % Setup X-axis (0 for stationary, then bin centers)
    edges = response.speedBins;
    binCenters = edges(1:end-1) + diff(edges)/2;
    xVals = [0, binCenters];
    
    figure('Color', 'w', 'Position', [100 100 1000 400]);
    t = tiledlayout(1, 2, 'TileSpacing', 'Loose');
    title(t, sprintf('ROI %d: %s Tuning vs. Shuffle', roiIdx, fname));

    %% Left Plot: Tuning Curves
    ax1 = nexttile; hold on;
    
    % 1. Plot a subset of shuffled curves in light grey
    rng(1); % Match your main function's seed
    minShift = 600;
    maxShift = length(wheelSpeed) - minShift;
    
    for i = 1:numShiftsToPlot
        thisShift = randi([minShift, maxShift]);
        sigS = circshift(sigData(roiIdx, :), thisShift);
        
        % Calculate shuffled means
        shuffStat = mean(sigS(sIdx), 'omitnan');
        shuffMove = zeros(1, length(binIndices));
        for b = 1:length(binIndices)
            shuffMove(b) = mean(sigS(binIndices{b}), 'omitnan');
        end
        
        plot(ax1, xVals, [shuffStat, shuffMove], 'Color', [0.8 0.8 0.8, 0.5], 'LineWidth', 0.5);
    end
    
    % 2. Plot Real Tuning Curve in Bold Black/Red
    realCurve = [response.(fname).statMean(roiIdx), response.(fname).moveMean(roiIdx,:)];
    realErr   = [response.(fname).statSEM(roiIdx), response.(fname).moveSEM(roiIdx,:)];
    
    errorbar(ax1, xVals, realCurve, realErr, 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'Real');
    
    xlabel(ax1, 'Speed (cm/s)'); ylabel(ax1, 'Response');
    title(ax1, 'Tuning Shape');
    grid on;

    %% Right Plot: Variance Distribution (The "Why" of the p-Value)
    ax2 = nexttile; hold on;
    
    % We need the full distribution of shuffled variances for this ROI
    % (Note: If you didn't save shuffVars in 'response', you'd recalculate here)
    % For illustration, let's assume we are visualizing the pVal result:
    
    text(0.5, 0.5, {['p-Value (Full): ' num2str(response.(fname).pValFull(roiIdx))], ...
                   ['Significant: ' num2str(response.(fname).isSignificant_999(roiIdx))]}, ...
                   'HorizontalAlignment', 'center', 'FontSize', 12);
    title(ax2, 'Significance Statistics');
    axis off;
end