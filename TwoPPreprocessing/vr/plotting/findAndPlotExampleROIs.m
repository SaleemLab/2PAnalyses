function findAndPlotExampleROIs(RegionData)
    % --- Configuration Options ---
    excludeStart = 30; 
    excludeEnd = 30;
    landmarkCentres  = [40, 80, 120, 160];
    gratingLandmarks = [40, 120];
    plaidLandmarks   = [80, 160];
    tolerance = 15;
    windowRadius = 15;
    
    gratingBins = []; for g = gratingLandmarks, gratingBins = [gratingBins, (g - windowRadius):(g + windowRadius)]; end
    plaidBins   = []; for p = plaidLandmarks,   plaidBins   = [plaidBins,   (p - windowRadius):(p + windowRadius)]; end

    % Storage master matrix to log all valid cells
    % Columns: [sessionIdx, roiIdx, indexValue, prefLandmark]
    cellLog = [];  
    
    %% --- Scan Engine ---
    for s = 1:length(RegionData)
        sess = RegionData(s);
        if ~isfield(sess, 'ConditionData') || ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs), continue; end
        
        condNames = fieldnames(sess.ConditionData);
        activeCond = condNames{1}; if ismember('Baseline', condNames), activeCond = 'Baseline'; end
        
        lapActivity = sess.ConditionData.(activeCond).LapActivity;
        nBins = size(lapActivity, 3);
        meanOdd = squeeze(mean(lapActivity(:, 1:2:end, :), 2, 'omitnan'));
        
        for i = 1:length(sess.FilteredROIs)
            roiIdx = sess.FilteredROIs(i);
            rawTrace = meanOdd(roiIdx, :);
            if all(isnan(rawTrace)), continue; end
            
            % Min-Max Normalization
            minVal = min(rawTrace((excludeStart+1):(nBins-excludeEnd)));
            maxVal = max(rawTrace((excludeStart+1):(nBins-excludeEnd)));
            rangeVal = maxVal - minVal; if rangeVal == 0, rangeVal = eps; end 
            roiTrace = (rawTrace - minVal) / rangeVal;
            
            % Locate Peak
            [~, maxRelIdx] = max(roiTrace((excludeStart+1):(nBins-excludeEnd)));
            prefBin = maxRelIdx + excludeStart;
            [minDist, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
            
            if minDist <= tolerance
                primaryLandmark = landmarkCentres(closestLandmarkIdx);
                
                actG = mean(roiTrace(gratingBins), 'omitnan');
                actP = mean(roiTrace(plaidBins), 'omitnan');
                actPeak = roiTrace(prefBin);
                
                if ismember(primaryLandmark, gratingLandmarks)
                    actAlternative = actP;
                else
                    actAlternative = actG;
                end
                
                if (actPeak + actAlternative) > 0
                    valSI = (actPeak - actAlternative) / (actPeak + actAlternative);
                    if ~isnan(valSI)
                        cellLog = [cellLog; s, roiIdx, valSI, primaryLandmark];
                    end
                end
            end
        end
    end
    
    %% --- Dynamic Selection Engine ---
    if isempty(cellLog)
        error('Could not find any valid landmark-responsive cells in this dataset.');
    end
    
    % Sort all logged cells from least selective to most selective
    [~, sortIdx] = sort(cellLog(:, 3));
    sortedCellLog = cellLog(sortIdx, :);
    
    % Dynamically grab the absolute best examples available in the dataset
    exLow  = sortedCellLog(1, :);            % The absolute lowest selectivity cell
    exHigh = sortedCellLog(end, :);          % The absolute highest selectivity cell
    
    %% --- Plotting Engine ---
    fig = figure('Color', 'w', 'Position', [150 150 900 380]);
    xBins = 1:nBins;
    
    % --- Plot Panel 1: Highly Selective ROI ---
    subplot(1, 2, 1); hold on;
    sessIdx1 = exHigh(1); roiIdx1 = exHigh(2); indexVal1 = exHigh(3);
    trace1 = squeeze(mean(RegionData(sessIdx1).ConditionData.(activeCond).LapActivity(roiIdx1, 1:2:end, :), 2, 'omitnan'));
    trace1 = (trace1 - min(trace1)) / (max(trace1) - min(trace1) + eps);
    
    % Layering shading boxes behind the neural data lines
    fillLandmarkZones(landmarkCentres, gratingLandmarks, plaidLandmarks, 1);
    plot(xBins, trace1, 'Color', 'k', 'LineWidth', 2.5);
    
    xlim([excludeStart, nBins-excludeEnd]); ylim([0 1.1]);
    xlabel('Track Position (cm)', 'FontName', 'Arial', 'FontSize', 10); 
    ylabel('Normalized Activity', 'FontName', 'Arial', 'FontSize', 10);
    title(sprintf('Example: Highly Selective ROI\nFeature Selectivity Index = %.3f', indexVal1), 'FontWeight', 'normal', 'FontSize', 11);
    axis square; box off;
    
    % --- Plot Panel 2: Less Selective ROI ---
    subplot(1, 2, 2); hold on;
    sessIdx2 = exLow(1); roiIdx2 = exLow(2); indexVal2 = exLow(3);
    trace2 = squeeze(mean(RegionData(sessIdx2).ConditionData.(activeCond).LapActivity(roiIdx2, 1:2:end, :), 2, 'omitnan'));
    trace2 = (trace2 - min(trace2)) / (max(trace2) - min(trace2) + eps);
    
    fillLandmarkZones(landmarkCentres, gratingLandmarks, plaidLandmarks, 1);
    plot(xBins, trace2, 'Color', 'k', 'LineWidth', 2.5);
    
    xlim([excludeStart, nBins-excludeEnd]); ylim([0 1.1]);
    xlabel('Track Position (cm)', 'FontName', 'Arial', 'FontSize', 10); 
    ylabel('Normalized Activity', 'FontName', 'Arial', 'FontSize', 10);
    title(sprintf('Example: Less Selective ROI\nFeature Selectivity Index = %.3f', indexVal2), 'FontWeight', 'normal', 'FontSize', 11);
    axis square; box off;
    
    % Print details to the terminal window
    fprintf('\n================ EXAMPLE ROI IDENTIFIED =================\n');
    fprintf('HIGH SELECTIVITY EXAMPLE: Session #%d, ROI #%d (Index: %.4f)\n', sessIdx1, roiIdx1, indexVal1);
    fprintf('LOW SELECTIVITY EXAMPLE:  Session #%d, ROI #%d (Index: %.4f)\n', sessIdx2, roiIdx2, indexVal2);
    fprintf('=========================================================\n\n');
end

%% --- Helper function to shade your track landmark locations ---
function fillLandmarkZones(landmarkCentres, gratingLandmarks, plaidLandmarks, yMax)
    for l = landmarkCentres
        xStart = l - 15;
        xEnd = l + 15;
        xPatch = [xStart, xEnd, xEnd, xStart];
        yPatch = [0, 0, yMax*1.05, yMax*1.05];
        
        if ismember(l, gratingLandmarks)
            % Light translucent gray for Grating regions
            fill(xPatch, yPatch, [0.85 0.85 0.85], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); 
        else
            % Darker translucent gray for Plaid regions
            fill(xPatch, yPatch, [0.65 0.65 0.65], 'FaceAlpha', 0.25, 'EdgeColor', 'none'); 
        end
    end
end