function [sessionFileInfo, response] = computeSpatialModulationRatio(sessionFileInfo, VRStimName, applySmoothing, excludeFlaggedLaps, doPlot)
if nargin < 2, error('Must provide the response struct.'); end
if nargin < 3, applySmoothing = true; end
if nargin < 4, excludeFlaggedLaps = true; end
if nargin < 5, doPlot = true; end

stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName ''%s'' not found.', VRStimName); end

response = load(sessionFileInfo.stimFiles(stimIdx).Response, 'lapPositionActivity', 'trialIndicesByCondition', 'stimName', 'flaggedLaps');
partnerLandmarkPosition = 80;
excludeStart = 30; excludeEnd = 30;  
landmarkCentres = [40, 80, 120, 160];
tolerance = 10; % +/- 10cm search window around landmarks

signalNames = fieldnames(response.lapPositionActivity);
disp('Computing Landmark-Restricted SMR for this session...\n');

SMR = struct(); SMR_Metrics = struct();

for iSignal = 1:length(signalNames)
    currentSignalName = signalNames{iSignal};
    lapActivity = response.lapPositionActivity.(currentSignalName);
    baseTrials = response.trialIndicesByCondition.Baseline;
    
    if excludeFlaggedLaps && isfield(response, 'flaggedLaps') && ~isempty(response.flaggedLaps)
        baseTrials = setdiff(baseTrials, response.flaggedLaps);
    end
    
    baseLapActivity = lapActivity(:, baseTrials, :);
    [numROIs, numLaps, numBins] = size(baseLapActivity);
    
    % landmark search window
    allowedLandmarkBins = [];
    for c = landmarkCentres
        allowedLandmarkBins = [allowedLandmarkBins, (c - tolerance):(c + tolerance)];
    end
    allowedLandmarkBins = unique(allowedLandmarkBins);
    validInnerRange = (excludeStart + 1) : (numBins - excludeEnd);
    allowedSearchBins = intersect(allowedLandmarkBins, validInnerRange);
    
    if applySmoothing
        w = gausswin(15); w = w / sum(w);
        for iCell = 1:numROIs
            for iLap = 1:numLaps
                trace = squeeze(baseLapActivity(iCell, iLap, :));
                if all(isnan(trace)), continue; end
                nanMask = isnan(trace); trace(nanMask) = 0;
                smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
                baseLapActivity(iCell, iLap, :) = smoothed;
            end
        end
    end
    
    % Train on Odd (Odds), Test on Even (Evens)
    meanOdd  = squeeze(mean(baseLapActivity(:, 1:2:end, :), 2, 'omitnan'));
    meanEven = squeeze(mean(baseLapActivity(:, 2:2:end, :), 2, 'omitnan'));
    
    smrValues = NaN(numROIs, 1); rpValues = NaN(numROIs, 1); rnValues = NaN(numROIs, 1);
    peakBins = NaN(numROIs, 1); rnBins = NaN(numROIs, 1);
    globalPeakBins = NaN(numROIs, 1); excludeEdgePeakCells = false(numROIs, 1); 
    allPks = cell(numROIs, 1); allLocs = cell(numROIs, 1);
    
    for thsROI = 1:numROIs
        trainTrace = meanOdd(thsROI, :);
        testTrace  = meanEven(thsROI, :);
        if all(isnan(trainTrace)) || all(isnan(testTrace)), continue; end

        % CROSS-NORMALIZATION (@changed aman)
        minOdd = min(trainTrace, [], 'omitnan'); 
        maxOdd = max(trainTrace, [], 'omitnan');
        rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end
        
        normTrain = (trainTrace - minOdd) ./ rangeOdd; 
        normTest  = (testTrace - minOdd) ./ rangeOdd;  
        
        % Find Global Peak across the entire track for the edge-check flag
        [globalMaxPk, globalPrefBin] = max(normTrain);
        globalPeakBins(thsROI) = globalPrefBin;
        allPks{thsROI} = globalMaxPk; allLocs{thsROI} = globalPrefBin;
        
        if (globalPrefBin <= excludeStart) || (globalPrefBin > (numBins - excludeEnd))
            excludeEdgePeakCells(thsROI) = true;
        end
        
        % RESTRICTED PEAK SEARCH: Look ONLY inside the landmark search windows
        [~, maxIdxInSearch] = max(normTrain(allowedSearchBins));
        prefBin = allowedSearchBins(maxIdxInSearch); 
        peakBins(thsROI) = prefBin;
        
        % Extract peak magnitude from independent Test Set (Evens)
        Rp = normTest(prefBin);
        rpValues(thsROI) = Rp;
        
        % Calculate spatial partner location (+/- 80 cm)
        if prefBin <= 100
            targetPartner = prefBin + partnerLandmarkPosition;
        else
            targetPartner = prefBin - partnerLandmarkPosition;
        end
        
        if targetPartner >= 1 && targetPartner <= numBins
            matchingBin = targetPartner;
            rnBins(thsROI) = matchingBin;
            %raw (possibly negative) normalised value used directly.
            Rn = normTest(matchingBin);
            rnValues(thsROI) = Rn;
        else
            continue;
        end
        
        % Spatial Modulation Ratio: ratio = peak_nonpreferred /
        % peak_preferred; Turns balues to nans instead of inf 
        % Only guard against exact division by zero (Rp == 0); no flooring applied.
        if Rp ~= 0
            smrValues(thsROI) = Rn / Rp;
        end
    end
    
    SMR.(currentSignalName) = smrValues;
    SMR_Metrics.(currentSignalName).SMR = smrValues;
    SMR_Metrics.(currentSignalName).Rp = rpValues;
    SMR_Metrics.(currentSignalName).Rn = rnValues;
    SMR_Metrics.(currentSignalName).RpBin = peakBins;
    SMR_Metrics.(currentSignalName).RnBin = rnBins;
    SMR_Metrics.(currentSignalName).PeaksAcrossFullTrack = allPks;
    SMR_Metrics.(currentSignalName).PeakBinAcrossFullTrack = allLocs;
    SMR_Metrics.(currentSignalName).GlobalPeakBin = globalPeakBins;
    SMR_Metrics.(currentSignalName).ExcludeEdgePeakCells = excludeEdgePeakCells;
end

disp('Landmark-Restricted SMR calculation complete.');
disp(['appending variable SMR_Metrics to file: ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, 'SMR_Metrics', '-append');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Done.');

if doPlot
    figure('Name', 'spatial modulation ratio summary', 'Color', [1 1 1]);
    numSignals = length(signalNames);
    
    subplot(1, 2, 1); hold on;
    plotGroupData = []; plotGroupLabels = {};
    for iSignal = 1:numSignals
        sigName = signalNames{iSignal};
        vals = SMR.(sigName)(~isnan(SMR.(sigName)));
        if ~isempty(vals)
            vals = SMR.(sigName)(~SMR_Metrics.(sigName).ExcludeEdgePeakCells & ~isnan(SMR.(sigName)));
            jitter = (rand(size(vals)) - 0.5) * 0.15 + iSignal;
            scatter(jitter, vals, 25, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.4);
            plotGroupData = [plotGroupData; vals];
            labels = repmat({sigName}, length(vals), 1);
            plotGroupLabels = [plotGroupLabels; labels];
        end
    end
    if ~isempty(plotGroupData), boxplot(plotGroupData, plotGroupLabels, 'Widths', 0.4, 'Colors', 'k'); end
    ylabel('spatial modulation ratio (smr)'); title('smr distribution per signal type'); grid on;
    
    subplot(1, 2, 2); hold on;
    for iSignal = 1:numSignals
        sigName = signalNames{iSignal};
        vals = SMR.(sigName)(~isnan(SMR.(sigName)));
        if ~isempty(vals)
            [f, x] = ecdf(vals);
            plot(x, f, 'LineWidth', 2.5, 'DisplayName', sigName);
        end
    end
    xlim([-2, 5]);
    xlabel('smr value'); ylabel('cumulative proportion'); title('cumulative distribution function');
    legend('Location', 'best', 'Interpreter', 'none'); grid on;
    drawnow;
end
end