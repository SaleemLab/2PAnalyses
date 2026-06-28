pairs=struct; 
pairs.M26005 = ['20260305', '20260306', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M26004 = ['20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs
pairs.M25131 = ['20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M25126 = ['20260311', '20260312', '20260313']; % unique fovs 
 

VISpSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

% edit this function to include spikes
VISpDataSpks = getTuningDataByCondition(VISpSessions, 'signalToUse', 'spks');

% these critera were included to plot all rois 
VISpDataSpks = appendFilteredROIs(VISpDataSpks, 'UseExpVar_SigNullDist', true, 'cvExpvarThreshold', 0.1, 'ExpVarSigThreshold', 0.01);
% filter used to plot the 3 conditions summary for aman
VISpDataSpks = appendFilteredROIs(VISpDataSpks, 'UseExpVar_SigNullDist', true, 'cvExpvarThreshold', 0.1, 'ExpVarSigThreshold', 0.01, 'UseHalves', true, 'RhoHalvesThreshold', 0.8, 'FilterSomasByRF', true);


% plotPooledPopulation_AcrossConditions(VISpDataSpks, 'V1', ...
%     'TypeToPlot', 'Somas', ...
%     'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\conditionsPooled\conditions_pooledAcrossMice');

% re-ran the crtiera two times to generate the two different figures 
plotThreeConditions(VISpDataSpks, 'V1', ...
    'TypeToPlot', 'Somas', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\conditionsPooled\conditions_pooledAcrossMice_3Conditions_withoutHalvesAndRF');

% this includes the difference plots
plotThreeConditions_DifferenceIncluded(VISpDataSpks, 'V1', ...
    'TypeToPlot', 'Somas', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\conditionsPooled\conditionsAndDiff_pooledAcrossMice_3Conditions_withoutHalvesAndRF')
%%
% find background rois with the smaller pool 
plotBackgroundROIs(VISpDataSpks, 'V1', ...
    'TypeToPlot', 'Somas', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\conditionsPooled\backroundROIs_3Conditions_withoutHalvesAndRF_new');

% run this with all rois 
plotBackgroundROIs_DifferenceOnly(VISpDataSpks, 'V1', ...
    'TypeToPlot', 'Somas', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\conditionsPooled\backroundROIs_Difference_5Conditions_withoutHalvesAndRF_new');


%% 
plotConditions_DifferenceIncluded(VISpDataSpks, 'V1', ...
    'TypeToPlot', 'Somas', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\conditionsPooled\conditionsAndDiff_pooledAcrossMice_5Conditions_withoutHalvesAndRF');

%% reformat julien's space order maps figure for session 20260318 m26004
fig = openfig("Z:\ibn-vision\USERS\Sonali\JulienModelFigs_M26004_20260318\SpaceOrderedMaps_M26004_20260318.fig");
ax = findobj(fig, 'Type', 'axes');

axSpace     = [];
axLandmarks = [];
axBG        = [];
meanKernelAxes = [];

for thisAxes = 1:length(ax)
    thisAx = ax(thisAxes);
    
    axTitle = get(thisAx.Title, 'String');
    axLabel = get(thisAx.YLabel, 'String');
    
    if contains(axTitle, 'Space') || contains(axLabel, 'Space')
        axSpace = thisAx;
    elseif contains(axTitle, 'Landmarks') || contains(axLabel, 'Landmarks')
        axLandmarks = thisAx;
    elseif contains(axTitle, 'BG') || contains(axLabel, 'BG')
        axBG = thisAx;
    else
        meanKernelAxes = [meanKernelAxes; thisAx];
    end
end

if isempty(axSpace) || isempty(axLandmarks) || isempty(axBG)
    axSpace     = ax(13); 
    axLandmarks = ax(12); 
    axBG        = ax(11); 
end

if isempty(meanKernelAxes)
    meanKernelAxes = ax(1:10);
end

saveDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\ForAS\';
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
targetFont = 'Arial';
landmarks = [40, 80, 120, 160];
tickLabels = {'40', '80', '120', '160'};

%% --- FIGURE 1: VERTICALLY ELONGATED HEATMAPS ---
newFig = figure('Color', 'w', 'Position', [30 30 1420 680]);
copiedHeatmaps = copyobj([axSpace, axLandmarks, axBG], newFig);
axSpaceNew     = copiedHeatmaps(1);
axLandmarksNew = copiedHeatmaps(2);
axBGNew        = copiedHeatmaps(3);

oldCBs = findobj(newFig, 'Type', 'colorbar');
if ~isempty(oldCBs), delete(oldCBs); end

targetHeatmaps = {axSpaceNew, axLandmarksNew, axBGNew};
displayTitles = {'Space', 'Landmarks', 'BG'};
widths      = [0.27, 0.17, 0.17];  
leftOffsets = [0.10, 0.46, 0.72];  
plotBottom  = 0.12;
plotHeight  = 0.75;

for aIdx = 1:3
    axTarget = targetHeatmaps{aIdx};
    set(axTarget, 'Position', [leftOffsets(aIdx), plotBottom, widths(aIdx), plotHeight]);
    hold(axTarget, 'on');
    
    imgObj = findobj(axTarget, 'Type', 'image');
    if ~isempty(imgObj)
        numN = size(imgObj.CData, 1);
    else
        numN = 186;
    end
    
    set(axTarget, 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out', ...
                  'XColor', 'k', 'YColor', 'k');
              
    set(axTarget, 'XTick', landmarks, 'XTickLabel', tickLabels);
    set(axTarget, 'YTick', [1, numN], 'YTickLabel', {'1', num2str(numN)});
    
    for lVal = landmarks
        plot(axTarget, [lVal, lVal], [numN * 0.96, numN], 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        plot(axTarget, [lVal, lVal], [1, numN * 0.04], 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    end
    
    xlabel(axTarget, 'Position (cm)', 'FontName', targetFont, 'FontSize', 11, 'FontWeight', 'bold');
    title(axTarget, displayTitles{aIdx}, 'FontName', targetFont, 'FontSize', 12, 'FontWeight', 'bold');
    
    if aIdx == 1
        yLabelText = {['\bfSorted Neurons '], ['\rm\it(n = ', num2str(numN), ')']};
        ylabel(axTarget, yLabelText, 'FontName', targetFont, 'FontSize', 11);
    else
        set(axTarget, 'YTickLabel', {'', ''});
    end
    
    set(axTarget, 'FontName', targetFont, 'FontSize', 10, 'LineWidth', 0.8);
    defaultAxesProperties(axTarget, false); 
    set(axTarget, 'XMinorTick', 'off', 'YMinorTick', 'off', 'TickDir', 'out', 'Box', 'off');
    
    cb = colorbar(axTarget, 'eastoutside');
    set(cb, 'Units', 'normalized', 'Position', [leftOffsets(aIdx) + widths(aIdx) + 0.012, plotBottom, 0.012, plotHeight]);
    cb.Ticks = [0 0.5 1]; 
    cb.TickLabels = {'0', '0.5', '1'};
    cb.TickDirection = 'out'; 
    cb.Box = 'off';
    cb.FontName = targetFont;
    cb.FontSize = 10;
    cb.Label.String = 'log (gain)';
    cb.Label.FontName = targetFont;
    cb.Label.FontSize = 12;
    cb.Label.FontWeight = 'bold';
end

exportgraphics(newFig, fullfile(saveDir, 'SpacLandBGeOrderedMaps.pdf'), 'ContentType', 'vector');


% %% --- FIGURE 2: TWO-ROW BALANCED MEAN TUNING CURVES ---
% curvesFig = figure('Color', 'w', 'Position', [30 30 1100 460]);
% copiedCurves = copyobj(meanKernelAxes, curvesFig);
% numCurves = length(copiedCurves);
% 
% oldCBsCurves = findobj(curvesFig, 'Type', 'colorbar');
% if ~isempty(oldCBsCurves), delete(oldCBsCurves); end
% 
% colsPerRow = ceil(numCurves / 2); 
% startPosCurve = 0.09;
% endPosCurve   = 0.96;
% curveWidth    = (endPosCurve - startPosCurve) / colsPerRow - 0.04; 
% gapCurve      = (endPosCurve - startPosCurve - (colsPerRow * curveWidth)) / (colsPerRow - 1);
% 
% curveHeight = 0.31; 
% rowBottoms  = [0.55, 0.14]; 
% 
% for cIdx = 1:numCurves
%     rowIdx = ceil(cIdx / colsPerRow); 
%     colIdx = mod(cIdx - 1, colsPerRow) + 1;
%     
%     curveX = startPosCurve + (colIdx - 1) * (curveWidth + gapCurve);
%     curveY = rowBottoms(rowIdx);
%     
%     axCurve = copiedCurves(cIdx);
%     set(axCurve, 'Position', [curveX, curveY, curveWidth, curveHeight]);
%     hold(axCurve, 'on');
%     
%     title(axCurve, '');
%     xlabel(axCurve, '');
%     ylabel(axCurve, '');
%     legend(axCurve, 'off');
%     delete(findobj(axCurve, 'Type', 'text')); 
%     
%     set(axCurve, 'XLim', [0 200], 'YLim', [-30 200]);
%     set(axCurve, 'XTick', landmarks, 'XTickLabel', tickLabels);
%     set(axCurve, 'YTick', [-30, 0, 100, 200], 'YTickLabel', {'-30', '0', '100', '200'});
%     
%     set(axCurve, 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out', ...
%                  'XColor', 'k', 'YColor', 'k');
%     
%     if colIdx ~= 1
%         set(axCurve, 'YTickLabel', {});
%     end
%     
%     set(axCurve, 'FontName', targetFont, 'FontSize', 10, 'LineWidth', 0.8, 'XMinorTick', 'off', 'YMinorTick', 'off');
%     set(axCurve, 'TickDir', 'out', 'Box', 'off');
%     
%     axCurve.XAxis.TickLabelRotation = 90;
%     axCurve.YAxis.TickLabelRotation = 90;
%     
%     defaultAxesProperties(axCurve, false); 
%     if exist('offsetAxes', 'file') == 2
%         offsetAxes(axCurve); 
%     end
% end
% 
% exportgraphics(curvesFig, fullfile(saveDir, 'SpaceOrderedMaps_MeanCurves.pdf'), 'ContentType', 'vector');
% 
% close(fig);

%% extract figure examples from the other figure 
fig = openfig("Z:\ibn-vision\USERS\Sonali\JulienModelFigs_M26004_20260318\CellExamples_GLM_M26004_20260318.fig");
ax = findobj(fig, 'Type', 'axes');

saveDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\ForAS\';
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

targetFont = 'Arial';
landmarks = [40, 80, 120, 160];
tickLabels = {'40', '80', '120', '160'};

% Loop through every axis found in this massive overview figure
for i = 1:length(ax)
    thisAx = ax(i);
    axTitle = get(thisAx.Title, 'String');
    
    % Strictly target only the Spatial Position columns
    if contains(axTitle, 'Spatial position kernel')
        hold(thisAx, 'on');
        
        % Force upper bound to 200 while letting the minimum adjust dynamically to the data range
        currentYLim = get(thisAx, 'YLim');
        set(thisAx, 'YLim', [currentYLim(1), 200]);
        
        % FIX: Forces uniform X-axis limits, ticks, and labels for alignment consistency
        set(thisAx, 'XLim', [0 200]);
        set(thisAx, 'XTick', landmarks, 'XTickLabel', tickLabels);
        
        % Ensure clean, vector-friendly properties before export
        set(thisAx, 'Box', 'off', 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');
        set(thisAx, 'XMinorTick', 'off', 'YMinorTick', 'off');
        
        % Keep the tick label values standing straight up at 90 degrees
        thisAx.XAxis.TickLabelRotation = 90;
        thisAx.YAxis.TickLabelRotation = 90;
        
        % Run your lab formatting tools on these specific axes
        defaultAxesProperties(thisAx, false); 
        if exist('offsetAxes', 'file') == 2
            offsetAxes(thisAx); 
        end
        
        % Re-enforce tick directions and box overrides that formatting tools can break
        set(thisAx, 'TickDir', 'out', 'Box', 'off');
    end
end

% Export the entire layout safely as a clean vector PDF for Affinity
exportgraphics(fig, fullfile(saveDir, 'SpatialPositionKernels_Overview.pdf'), 'ContentType', 'vector');