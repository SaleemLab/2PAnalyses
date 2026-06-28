function response = classifySpeedTuningFromCurves_80_20Split(sessionFileInfo, response, useField)
% Parameterizes and classifies neurons into 4 main tuning classes using
% a Multi-Template Competitive Gaussian Fit framework.
%
% SALEEM et al 2013
%   1. Segments the continuous timeline into 20-second contiguous blocks.
%   2. Performs a random 80/20 block split (Train/Test).
%   3. Dynamically calculates adaptive quantile edges strictly using 
%      the 80% training data pool, matching your exact original bin count.
%   4. Fits Edd's symmetric templates to the Training Curves.
%   5. Evaluates model quality (R2) against the unseen 20% Testing Curves.
    if nargin < 3, useField = 'dFFNeuropilCorrected'; end
    
    stimFileName = sprintf('%s_%s_Response_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
    fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
    
    if ~exist(fileFullPath, 'file')
        error('Internal path recovery failed. File does not exist: %s', fileFullPath);
    end
    
    % Load raw time-series data required for the block slice
    stimIdx = find(strcmp(response.stimName, {sessionFileInfo.stimFiles.name}), 1);
    data2P = load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData);
    dataPeriph = load(sessionFileInfo.stimFiles(stimIdx).processedPeripheralData, 'peripheralData');
    
    sourceStructName = response.tuningCurve.source;
    
    if contains(useField, 'spks')
        sigData = data2P.spks;
    else 
        signals = data2P.(sourceStructName);
        sigData = signals.(useField);
    end 
    
    fs = 60; 
    tickToCm = 3.1415 * 20 / 1024; 
    wheelSpeed = [0; diff(dataPeriph.peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(dataPeriph.peripheralData.Wheel.sampleTimes)];
    wheelSpeed(abs(wheelSpeed) > 150) = NaN; 
    totalFrames = length(wheelSpeed);

    % 80/20 CONTIGUOUS TIME-BLOCK SPLIT (DIMENSION ORIENTED)
    blockSizeFrames = 20 * fs; 
    numBlocks = floor(totalFrames / blockSizeFrames);
    
    rng(42); 
    shuffledBlockIdx = randperm(numBlocks);
    numTrainBlocks = round(0.80 * numBlocks);
    
    trainBlocks = shuffledBlockIdx(1:numTrainBlocks);
    testBlocks  = shuffledBlockIdx(numTrainBlocks+1:end);
    
    trainMask = false(size(wheelSpeed)); 
    for b = trainBlocks
        startF = ((b-1) * blockSizeFrames) + 1;
        endF   = b * blockSizeFrames;
        trainMask(startF:endF) = true;
    end
    testMask = ~trainMask;

    % Only use the running bins 
    movingTrainMask = (wheelSpeed > 1) & ~isnan(wheelSpeed) & trainMask; 
    trainSpeeds = wheelSpeed(movingTrainMask);
    
    % FIX: Read the exact baseline bin count your pipeline originally used
    % This stops 'edgesQuantileTrain' from shifting lengths and breaking your plots!
    originalBinCount = length(response.tuningCurve.speedBins) - 1;
    
    if ~isempty(trainSpeeds) && originalBinCount > 0
        edgesQuantileTrain = quantile(trainSpeeds, linspace(0, 1, originalBinCount + 1));
        edgesQuantileTrain = unique(edgesQuantileTrain);
        
        % If unique dropped a bin boundary due to identical values, restore length 
        if length(edgesQuantileTrain) ~= (originalBinCount + 1)
            edgesQuantileTrain = linspace(min(trainSpeeds), max(trainSpeeds), originalBinCount + 1);
        end
    else
        edgesQuantileTrain = linspace(1, 50, originalBinCount + 1); 
    end
    
    % Standard fixed rigid boundaries for cross-session tracking alignment
    edgesFixed = [1.5, 3.5, 7.0, 14.0, 25.0, 45.0];
    
    % Reassign matching speed axes
    response.tuningCurve.speedBins = edgesQuantileTrain;
    response.tuningCurveFixedBins.speedBins = edgesFixed;
    frameworks = {'tuningCurve', 'tuningCurveFixedBins'};
    edgeDefinitions = {edgesQuantileTrain, edgesFixed};
    
    % STEP 3: CROSS-VALIDATED PARAMETRIC COMPETING FITS
    for f = 1:2
        targetStruct = frameworks{f};
        edges = edgeDefinitions{f};
        movingCenters = (edges(1:end-1) + diff(edges)/2)'; 
        numBins = length(movingCenters);
        numROIs = size(sigData, 1);
        
        cls = struct();
        cls.tuningType        = cell(numROIs, 1);
        cls.tuningCode        = nan(numROIs, 1);
        cls.R2                = nan(numROIs, 1); 
        cls.fitParams         = nan(numROIs, 4);
        cls.preferredSpeed    = nan(numROIs, 1);
        cls.usedMovingOnly    = true;
        cls.classifierVersion = "80_20_cm_s";
        
        gaussFun = @(params, xdata) params(1) + params(2) .* exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
        optsFit = optimset('Display', 'off');
        
        minSpeed = min(movingCenters); maxSpeed = max(movingCenters); speedRange = maxSpeed - minSpeed; customOffset = 0.5; 
        
        % --- GEOMETRIC AXIS CONSTRAINTS ---
        minAllowableCenter = minSpeed + (0.15 * speedRange);
        maxAllowableCenter = maxSpeed - (0.15 * speedRange);
        
        % Enforce lower and upper bounds on width separately to prevent overrides
        minAllowableSigma  = 1.0;
        maxAllowableSigma  = 0.25 * speedRange;
        
        fprintf('Executing Saleem 80/20 CV fits on %d ROIs [%s Framework]...\n', numROIs, targetStruct);
        
        for r = 1:numROIs
            yTrain = nan(numBins, 1);
            yTest  = nan(numBins, 1);
            
            isMoving = (wheelSpeed > 1) & ~isnan(wheelSpeed);
            
            for b = 1:numBins
                % Restrict bin strictly to active moving frames
                movbinIdx = wheelSpeed >= edges(b) & wheelSpeed < edges(b+1) & isMoving;
                
                yTrain(b) = mean(sigData(r, movbinIdx & trainMask), 'omitnan');
                yTest(b)  = mean(sigData(r, movbinIdx & testMask), 'omitnan');
            end

            validTrain = ~isnan(yTrain);
            validTest  = ~isnan(yTest);
            
            if sum(validTrain) < 4 || sum(validTest) < 4
                cls.tuningType{r} = 'untuned'; cls.tuningCode(r) = 0; cls.R2(r) = 0; continue;
            end
            
            xtrainvals = movingCenters(validTrain); ytrainvals = yTrain(validTrain);
            xtestvals  = movingCenters(validTest);  ytestvals  = yTest(validTest);
            
            [maxVal, maxidx] = max(ytrainvals); [minVal, minidx] = min(ytrainvals);
            rangeVal = maxVal - minVal; if rangeVal <= 0, rangeVal = 1e-5; end
            
            % Compute Pearson correlation coefficient for selection gating
            speedCorrelation = corr(xtrainvals, ytrainvals, 'rows', 'complete', 'type', 'Pearson');
            param_out = nan(6, 4); resnorm_train = inf(6, 1);
            
            % Setup matrices natively on the physical scale
            lb_low1  = [-200,   0.1, minSpeed - customOffset*10, minAllowableSigma];  ub_low1  = [200,  200, minSpeed + customOffset,    speedRange*5]; x0_low1  = [0,  rangeVal, minSpeed - customOffset*2, speedRange/2];
            lb_low2  = [-200, -200.0, maxSpeed - customOffset,   minAllowableSigma];  ub_low2  = [200, -0.1, maxSpeed + customOffset*10, speedRange*5]; x0_low2  = [0, -rangeVal, maxSpeed + customOffset*2, speedRange/2];
            lb_high1 = [-200,   0.1, maxSpeed - customOffset,   minAllowableSigma];  ub_high1 = [200,  200, maxSpeed + customOffset*10, speedRange*5]; x0_high1 = [0,  rangeVal, maxSpeed + customOffset*2, speedRange/2];
            lb_high2 = [-200, -200.0, minSpeed - customOffset*10, minAllowableSigma]; ub_high2 = [200, -0.1, minSpeed + customOffset,    speedRange*5]; x0_high2 = [0, -rangeVal, minSpeed - customOffset*2, speedRange/2];
            lb_band  = [-200,   0.1, minAllowableCenter,         minAllowableSigma];  ub_band  = [200,  200, maxAllowableCenter,         maxAllowableSigma]; x0_band  = [0,  rangeVal, xtrainvals(maxidx), speedRange/4];
            lb_inv   = [-200, -200.0, minAllowableCenter,         minAllowableSigma];  ub_inv   = [200, -0.1, maxAllowableCenter,         maxAllowableSigma]; x0_inv   = [0, -rangeVal, xtrainvals(minidx), speedRange/4];
            x0_band(3) = min(max(x0_band(3), lb_band(3)), ub_band(3));
            x0_inv(3)  = min(max(x0_inv(3),  lb_inv(3)),  ub_inv(3));
            
            try [param_out(1,:), resnorm_train(1)] = lsqcurvefit(gaussFun, x0_low1,  xtrainvals, ytrainvals, lb_low1,  ub_low1,  optsFit); catch, end
            try [param_out(2,:), resnorm_train(2)] = lsqcurvefit(gaussFun, x0_low2,  xtrainvals, ytrainvals, lb_low2,  ub_low2,  optsFit); catch, end
            try [param_out(3,:), resnorm_train(3)] = lsqcurvefit(gaussFun, x0_high1, xtrainvals, ytrainvals, lb_high1, ub_high1, optsFit); catch, end
            try [param_out(4,:), resnorm_train(4)] = lsqcurvefit(gaussFun, x0_high2, xtrainvals, ytrainvals, lb_high2, ub_high2, optsFit); catch, end
            try [param_out(5,:), resnorm_train(5)] = lsqcurvefit(gaussFun, x0_band,  xtrainvals, ytrainvals, lb_band,  ub_band,  optsFit); catch, end
            try [param_out(6,:), resnorm_train(6)] = lsqcurvefit(gaussFun, x0_inv,   xtrainvals, ytrainvals, lb_inv,   ub_inv,   optsFit); catch, end
            
            [~, fitIdx] = sort(resnorm_train, 'ascend');
            chosenIdx = nan; bestParams = nan(1,4); typeStr = 'untuned'; charCode = 0;
            idx_idx = 1;
            
            % --- SELECTION PROCESSOR ---
            while idx_idx <= 6
                bestIdx = fitIdx(idx_idx); candidateParams = param_out(bestIdx, :);
                if any(~isfinite(candidateParams)), idx_idx = idx_idx + 1; continue; end
                
                if bestIdx == 1 || bestIdx == 2
                    % Low-pass: correlates negatively AND maximum must be at first bin
                    firstBinIsMax = (maxidx == 1);
                    if ~isnan(speedCorrelation) && speedCorrelation < -0.35 && firstBinIsMax
                        chosenIdx = bestIdx; bestParams = candidateParams; typeStr = 'lowpass'; charCode = 1; break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                elseif bestIdx == 3 || bestIdx == 4
                    % High-pass: must correlate positively AND last bin near maximum
                    lastBinIsMax = (ytrainvals(end) >= 0.80 * max(ytrainvals));
                    if ~isnan(speedCorrelation) && speedCorrelation > 0.50 && lastBinIsMax
                        chosenIdx = bestIdx; bestParams = candidateParams; typeStr = 'highpass'; charCode = 2; break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                elseif bestIdx == 5
                    % Band-pass: strict geometric + raw data criteria
                    fitY = gaussFun(candidateParams, linspace(minSpeed, maxSpeed, 200));
                    % Gaussian peak must sit in the interior 70% of the speed range
                    centerInside      = candidateParams(3) > minAllowableCenter && candidateParams(3) < maxAllowableCenter;
                    % Width of the Gaussian can't be so broad it's essentially flat
                    sigmaAcceptable   = candidateParams(4) <= maxAllowableSigma;
                    % The actual data maximum must not be at the first or last bin.
                    rawPeakInterior   = maxidx > 1 && maxidx < numel(ytrainvals);
                    % The fit must drop by at least 15% of the range on both sides. Prevents a monotonic curve being called band-pass.
                    dropsBothSidesFit = (max(fitY) - fitY(1) > 0.15 * rangeVal) && (max(fitY) - fitY(end) > 0.15 * rangeVal);
                    % Same check but on raw data — 10% drop on both sides.
                    dropsBothSidesRaw = (maxVal - ytrainvals(1) > 0.10 * rangeVal) && (maxVal - ytrainvals(end) > 0.10 * rangeVal);
                    % Peak must stand at least 15% above the higher of the two endpoints.
                    fitProminence     = max(fitY) - max(fitY(1), fitY(end));
                    if centerInside && sigmaAcceptable && rawPeakInterior && dropsBothSidesFit && dropsBothSidesRaw && (fitProminence > 0.15 * rangeVal)
                        chosenIdx = bestIdx; bestParams = candidateParams; typeStr = 'bandpass'; charCode = 3; break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                elseif bestIdx == 6
                    % Trough-inverted: symmetric strict criteria to band-pass
                    fitY = gaussFun(candidateParams, linspace(minSpeed, maxSpeed, 200));
                    centerInside      = candidateParams(3) > minAllowableCenter && candidateParams(3) < maxAllowableCenter;
                    sigmaAcceptable   = candidateParams(4) <= maxAllowableSigma;
                    rawTroughInterior = minidx > 1 && minidx < numel(ytrainvals);
                    dropsBothSidesFit = (fitY(1) - min(fitY) > 0.15 * rangeVal) && (fitY(end) - min(fitY) > 0.15 * rangeVal);
                    dropsBothSidesRaw = (ytrainvals(1) - minVal > 0.10 * rangeVal) && (ytrainvals(end) - minVal > 0.10 * rangeVal);
                    fitProminence     = min(fitY(1), fitY(end)) - min(fitY);
                    if centerInside && sigmaAcceptable && rawTroughInterior && dropsBothSidesFit && dropsBothSidesRaw && (fitProminence > 0.15 * rangeVal)
                        chosenIdx = bestIdx; bestParams = candidateParams; typeStr = 'trough_inverted'; charCode = 4; break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                end
            end
            
            if isnan(chosenIdx)
                typeStr = 'untuned'; charCode = 0; bestParams = [nan nan nan nan]; cvR2 = 0;
            else
                yModelPred = gaussFun(bestParams, xtestvals);
                SS_res = sum((ytestvals - yModelPred).^2);
                SS_tot = sum((ytestvals - mean(ytestvals)).^2);
                if SS_tot < 1e-9, SS_tot = 1; end
                cvR2 = 1 - (SS_res / SS_tot);
                if cvR2 < -2, cvR2 = -2; end 
            end
            
            if charCode == 1
                prefSpeedVal = minSpeed;
            elseif charCode == 2
                prefSpeedVal = maxSpeed;
            elseif charCode == 3 || charCode == 4
                pVal = bestParams(3);
                if pVal < minSpeed, pVal = minSpeed; end
                if pVal > maxSpeed, pVal = maxSpeed; end
                prefSpeedVal = pVal;
            else
                prefSpeedVal = nan;
            end
            
            cls.tuningType{r}     = typeStr;
            cls.tuningCode(r)     = charCode;
            cls.R2(r)             = cvR2; 
            cls.fitParams(r, :)   = bestParams;
            cls.preferredSpeed(r) = prefSpeedVal;
        end
        
        response.(targetStruct).(useField).classification_80_20Split = cls;
        allTypes = cls.tuningType;
        fprintf('\n--- Multi-Template Competitive Speed Summary [%s] ---\n', targetStruct);
        fprintf('  Total Functional ROIs evaluated: %d\n', numROIs);
        fprintf('  ---------------------------------------------------\n');
        fprintf('  Template Group 1 [Low-Pass / Decay] : %d\n', sum(strcmp(allTypes, 'lowpass')));
        fprintf('  Template Group 2 [High-Pass / Rise] : %d\n', sum(strcmp(allTypes, 'highpass')));
        fprintf('  Template Group 3 [Band-Pass / Peak] : %d\n', sum(strcmp(allTypes, 'bandpass')));
        fprintf('  Template Group 4 [Trough / Inverted]: %d\n', sum(strcmp(allTypes, 'trough_inverted')));
        fprintf('  Untuned Failsafes (Insufficient Data): %d\n', sum(strcmp(allTypes, 'untuned')));
        fprintf('-----------------------------------------------------\n\n');
    end
    save(fileFullPath, 'response', '-append');
    fprintf('Successfully finished Saleem-style cross-validated execution loop.\n');
end
