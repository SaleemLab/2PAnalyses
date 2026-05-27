function plotRoiActivity_AcrossTrials(response, roiIdx, varargin)
    defaultSignal = 'dFFNeuropilCorrected';
    expectedSignals = {'dFF', 'dFFNeuropilCorrected', 'spks'};
    
    p = inputParser;
    addRequired(p, 'response', @isstruct);
    addRequired(p, 'roiIdx', @isnumeric);
    addParameter(p, 'SignalType', defaultSignal, @(x) any(validatestring(x, expectedSignals)));
    parse(p, response, roiIdx, varargin{:});
    
    signalType = p.Results.SignalType;
    fullData = response.lapPositionActivity.(signalType);
    
    roiData = squeeze(fullData(roiIdx, :, :));
    if isvector(roiData)
        roiData = roiData(:)'; 
    end
    
    meanActivity = nanmean(roiData, 1);
    
    figure('Color', 'w', 'Position', [100, 100, 800, 450]);
    hold on;
    
    numBins = size(roiData, 2);
    xBins = 1:numBins;
    
    hTrials = plot(xBins, roiData', 'Color', [0.5, 0.5, 0.5, 0.25], 'LineWidth', 0.75);
    hMean = plot(xBins, meanActivity, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2.5);
    
    title(sprintf('ROI %d (%s)', roiIdx, signalType), 'Interpreter', 'none');
    xlabel('Position (cm)');
    ylabel('\DeltaF/F'); xline([40 80 120 160], 'k--'); 
    xlim([1, numBins]);
    
    legend([hTrials(1), hMean], {'Individual Laps', 'Mean Activity'}, 'Location', 'best');
    
   
    ax = gca;
    ax.Box = 'off';
    hold off;
end