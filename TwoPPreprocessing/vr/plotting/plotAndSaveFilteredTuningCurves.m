function plotAndSaveFilteredTuningCurves(sessionMetrics, applySmoothing)
    if nargin < 2
        applySmoothing = true;
    end
    
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1\@Aman\EVFilter_SMILandBoundary\VISp';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    numSessions = length(sessionMetrics);
    for s = 1:numSessions
        sess = sessionMetrics(s);
        
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline')
            warning('session %d has no baseline dataset profile. skipping.', s);
            continue;
        end
        if ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            fprintf('session %d (%s) contains no valid filtered rois to save.\n', s, sess.Session);
            continue;
        end
        
        lapActivityFull = sess.ConditionData.Baseline.LapActivity;
        [~, nLaps, nBins] = size(lapActivityFull);
        validLapRange = 31 : (nBins - 30); 
        partnerLandmarkPosition = 80;
        landmarkCentres = [40, 80, 120, 160];
        peakWindow = 15;
        
        pdfName = sprintf('session_%d_%s_day%d_all_filtered_profiles.pdf', s, sess.MouseID, sess.Day);
        pdfPath = fullfile(outputDir, pdfName);
                      
        if exist(pdfPath, 'file')
            delete(pdfPath);
        end
        
        targetROIs = sess.FilteredROIs;
        
        % Calculate median cross-validated explained variance across folds
        % medianEV = median(sess.cvExpVar, 2, 'omitnan');
        
        % Safely extract p-values matching the ROIs
        if isfield(sess, 'cvExpVar')
            pValues = sess.cvExpVar.pValues;
        else
            pValues = NaN(size(targetROIs)); % Fallback if missing
        end

        if isfield(sess, 'SMI')
            smi = sess.SMI.SMI; 

        else
            smi = nan(length(targetROIs));
        end 
        
        % include all filtered rois and also add the median ev to the title
        % and also the pval  
        validROIs = targetROIs; 
        
        baseLapActivity = lapActivityFull(:, :, :);
        if applySmoothing
            w = gausswin(15); w = w / sum(w);
            for iCell = 1:size(baseLapActivity, 1)
                for iLap = 1:nLaps
                    trace = squeeze(baseLapActivity(iCell, iLap, :));
                    if all(isnan(trace)), continue; end
                    nanMask = isnan(trace);
                    trace(nanMask) = 0;
                    smoothed = filtfilt(w, 1, trace);
                    smoothed(nanMask) = NaN;
                    baseLapActivity(iCell, iLap, :) = smoothed;
                end
            end
        end
        
        meanOdd = squeeze(mean(baseLapActivity(:, 1:2:end, :), 2, 'omitnan'));
        meanEven = squeeze(mean(baseLapActivity(:, 2:2:end, :), 2, 'omitnan'));
        
        fprintf('generating complete for %d filtered rois in session %d (%s)...\n', ...
                length(validROIs), s, sess.Session);
            
        for iROI = 1:length(validROIs)
            targetROI = validROIs(iROI);
            
            % 
            fig = figure('Visible', 'off', 'Position', [100 100 1200 450]);
            
            %% Lines plot odd and even 
            subplot(1, 2, 1);
            hold on;
            
            yLims = selectYLim(meanOdd(targetROI,:), meanEven(targetROI,:));
            
            % Plot original landmark shading
            for iLand = 1:length(landmarkCentres)
                lowB = max(1, landmarkCentres(iLand) - peakWindow);
                highB = min(nBins, landmarkCentres(iLand) + peakWindow);
                shadeLow = max(lowB, validLapRange(1));
                shadeHigh = min(highB, validLapRange(end));
                if shadeLow <= shadeHigh
                    fill([shadeLow shadeHigh shadeHigh shadeLow], [yLims(1) yLims(1) yLims(2) yLims(2)], ...
                        [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
                end
            end
            
            % 
            plot(1:nBins, meanOdd(targetROI, :), 'Color', [0.824, 0.016, 0.176], 'LineWidth', 2, 'DisplayName', 'mean odd laps');
            plot(1:nBins, meanEven(targetROI, :), 'Color', [0.000, 0.400, 1.000], 'LineWidth', 2, 'DisplayName', 'mean even laps');
            
            % Calculate preference parameters
            landmarkBins = [];
            for iLand = 1:length(landmarkCentres)
                lowBound = max(1, landmarkCentres(iLand) - peakWindow);
                highBound = min(nBins, landmarkCentres(iLand) + peakWindow);
                landmarkBins = [landmarkBins, lowBound:highBound];
            end
            validPeakBins = intersect(unique(landmarkBins), validLapRange);
            
            [~, maxRelIdx] = max(meanOdd(targetROI, validPeakBins));
            prefBin = validPeakBins(maxRelIdx);
            
            if prefBin <= 100
                partBin = prefBin + partnerLandmarkPosition; 
            else
                partBin = prefBin - partnerLandmarkPosition; 
            end
            
            % Superimpose scatter targets evaluated against Even steps
            plot(prefBin, meanEven(targetROI, prefBin), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Rp (pref bin)');
            if partBin >= 1 && partBin <= nBins
                plot(partBin, meanEven(targetROI, partBin), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Rn (partner bin)');
            end
            
            % Title updated to display ROI, Median EV, and Null Dist p-value
            title(sprintf('ROI %d Tuning Curve (median EV: %.2f | p: %.3f and SMI: %.2f )', ...
                  targetROI, sess.cvExpVar.medianExpVar(targetROI), pValues(targetROI)), smi(targetROI), 'Interpreter', 'none');
            xlabel('position bins');
            ylabel('activity');
            ylim(yLims);
            grid on;
            legend('Location', 'best');
            
            %% 
            subplot(1, 2, 2);
            
            % Extract matrix layout of [Laps x Bins] for the active ROI
            roiAllLapsData = squeeze(baseLapActivity(targetROI, :, :));
            
            % Generate the structural heatmap image
            imagesc(1:nBins, 1:nLaps, roiAllLapsData);
            colormap(gca, 'parula'); 
            colorbar;
            
            title('All Laps Activity Heatmap');
            xlabel('spatial position bins');
            ylabel('laps');
            
            % Export comprehensive graphic layout
            exportgraphics(fig, pdfPath, 'Append', true);
            close(fig);
        end
        fprintf('session %d full documentation complete -> output saved to: %s\n', s, pdfPath);
    end
end

function yLimits = selectYLim(traceA, traceB)
    mx = max([traceA, traceB]);
    mn = min([traceA, traceB]);
    rangeVal = mx - mn;
    if rangeVal <= 1e-9, rangeVal = 1; end
    yLimits = [mn - 0.1*rangeVal, mx + 0.1*rangeVal];
end