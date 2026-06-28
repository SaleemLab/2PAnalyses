function response = classifySpeedTuningFromCurves_8020_split(sessionFileInfo, response, useField)
% CLASSIFYSPEEDTUNINGFROMCURVES_V3_8020
%
% Edd-style Multi-Template Competitive Gaussian Fit with 80/20 CV.
%
% 
% ─────────────────────────────────────────────────────────────────────────
%   Fits on 80% of frames (training), evaluates CV-R² on held-out 20%.
%   Low-pass:  maxidx==1 (data max must be at first fitted bin)
%                + firstBinWasFit (bin 1 must not be NaN)
%   High-pass: last bin >= 80% of max (data must peak at end)
%                + lastBinWasFit (last bin must not be NaN)
%   Band-pass / trough: prominence gate only (rangeVal/3)
%   candidateR2 > 0: fit must beat a flat line
%   Tight mu bounds + sigma cap on monotonic templates
%   CV-R² saved to cls.R2 — threshold in your population script
%   Result saved to .classification_8020
% ─────────────────────────────────────────────────────────────────────────

    if nargin < 3 || isempty(useField), useField = 'dFFNeuropilCorrected'; end

    stimFileName = sprintf('%s_%s_Response_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
    fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

    if ~exist(fileFullPath, 'file')
        error('Response file not found: %s', fileFullPath);
    end

    % ── Load raw signals and wheel ───────────────────────────────────────
    stimIdx    = find(strcmp(response.stimName, {sessionFileInfo.stimFiles.name}), 1);
    data2P     = load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData);
    dataPeriph = load(sessionFileInfo.stimFiles(stimIdx).processedPeripheralData, 'peripheralData');

    sourceStructName = response.tuningCurve.source;

    if contains(useField, 'spks')
        sigData = data2P.spks;
    else
        signals = data2P.(sourceStructName);
        sigData = signals.(useField);
    end

    fs         = 60;
    tickToCm   = 3.1415 * 20 / 1024;
    wheelSpeed = [0; diff(dataPeriph.peripheralData.Wheel.Value * tickToCm)] ./ ...
                 [1; diff(dataPeriph.peripheralData.Wheel.sampleTimes)];
    wheelSpeed(abs(wheelSpeed) > 150) = NaN;
    totalFrames = length(wheelSpeed);

    % ── 80/20 contiguous block split ────────────────────────────────────
    blockSizeFrames = 20 * fs;
    numBlocks       = floor(totalFrames / blockSizeFrames);

    rng(42);
    shuffledBlockIdx = randperm(numBlocks);
    numTrainBlocks   = round(0.80 * numBlocks);
    trainBlocks      = shuffledBlockIdx(1:numTrainBlocks);

    trainMask = false(size(wheelSpeed));
    for b = trainBlocks
        startF = (b-1) * blockSizeFrames + 1;
        endF   = b * blockSizeFrames;
        trainMask(startF:endF) = true;
    end
    testMask = ~trainMask;

    isMoving = (wheelSpeed > 1) & ~isnan(wheelSpeed);

    frameworks      = {'tuningCurve', 'tuningCurveFixedBins'};
    edgeDefinitions = { response.tuningCurve.speedBins, ...
                        response.tuningCurveFixedBins.speedBins };

    for f = 1:2
        targetStruct = frameworks{f};
        if ~isfield(response, targetStruct) || ...
           ~isfield(response.(targetStruct), useField), continue; end

        edges         = edgeDefinitions{f};
        movingCenters = (edges(1:end-1) + diff(edges)/2)';
        numBins       = length(movingCenters);
        numROIs       = size(sigData, 1);

        cls = struct();
        cls.tuningType        = cell(numROIs, 1);
        cls.tuningCode        = nan(numROIs, 1);
        cls.R2                = nan(numROIs, 1);
        cls.fitParams         = nan(numROIs, 4);
        cls.preferredSpeed    = nan(numROIs, 1);
        cls.usedMovingOnly    = true;
        cls.classifierVersion = "edd_style_v3_8020_final2";

        gaussFun = @(params, xdata) params(1) + params(2) .* ...
            exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
        optsFit = optimset('Display', 'off');

        minSpeed   = min(movingCenters);
        maxSpeed   = max(movingCenters);
        speedRange = maxSpeed - minSpeed;

        customOffset               = mean(diff(movingCenters));
        minAllowableCenter         = minSpeed  + (0.15 * speedRange);
        maxAllowableCenter         = maxSpeed  - (0.15 * speedRange);
        minAllowableSigma          = 1.0;
        maxAllowableSigma          = 0.25 * speedRange;
        maxAllowableSigmaMonotonic = speedRange;

        xDense = linspace(minSpeed, maxSpeed, 200);

        fprintf('Running Edd-style v3 80/20 final fits on %d ROIs [%s]...\n', numROIs, targetStruct);

        for r = 1:numROIs

            % ── Build tuning curves from train / test frames ─────────────
            yTrain = nan(numBins, 1);
            yTest  = nan(numBins, 1);

            for b = 1:numBins
                binIdx    = wheelSpeed >= edges(b) & wheelSpeed < edges(b+1) & isMoving;
                yTrain(b) = mean(sigData(r, binIdx & trainMask), 'omitnan');
                yTest(b)  = mean(sigData(r, binIdx & testMask),  'omitnan');
            end

            firstBinWasFit = ~isnan(yTrain(1));
            lastBinWasFit  = ~isnan(yTrain(end));

            validTrain = ~isnan(yTrain);
            validTest  = ~isnan(yTest);

            if sum(validTrain) < 4 || sum(validTest) < 2
                cls.tuningType{r} = 'untuned'; cls.tuningCode(r) = 0;
                cls.R2(r) = 0; continue;
            end

            xtrainvals = movingCenters(validTrain);
            ytrainvals = yTrain(validTrain);
            xtestvals  = movingCenters(validTest);
            ytestvals  = yTest(validTest);

            [maxVal, maxidx] = max(ytrainvals);
            [minVal, minidx] = min(ytrainvals);
            rangeVal = maxVal - minVal;
            if rangeVal <= 0, rangeVal = 1e-5; end

            % ── 6 competing templates ────────────────────────────────────
            param_out = nan(6, 4);
            resnorm   = inf(6, 1);

            lb_low1  = [-200,   0.1, minSpeed - customOffset, minAllowableSigma         ];
            ub_low1  = [ 200,   200, minSpeed + customOffset, maxAllowableSigmaMonotonic ];
            x0_low1  = [   0, rangeVal, minSpeed,             speedRange/2              ];

            lb_low2  = [-200,  -200, maxSpeed - customOffset, minAllowableSigma         ];
            ub_low2  = [ 200,  -0.1, maxSpeed + customOffset, maxAllowableSigmaMonotonic];
            x0_low2  = [   0,-rangeVal, maxSpeed,             speedRange/2              ];

            lb_high1 = [-200,   0.1, maxSpeed - customOffset, minAllowableSigma         ];
            ub_high1 = [ 200,   200, maxSpeed + customOffset, maxAllowableSigmaMonotonic];
            x0_high1 = [   0, rangeVal, maxSpeed,             speedRange/2              ];

            lb_high2 = [-200,  -200, minSpeed - customOffset, minAllowableSigma         ];
            ub_high2 = [ 200,  -0.1, minSpeed + customOffset, maxAllowableSigmaMonotonic];
            x0_high2 = [   0,-rangeVal, minSpeed,             speedRange/2              ];

            lb_band  = [-200,   0.1, minAllowableCenter, minAllowableSigma];
            ub_band  = [ 200,   200, maxAllowableCenter, maxAllowableSigma ];
            x0_band  = [   0, rangeVal, xtrainvals(maxidx), speedRange/4  ];
            x0_band(3) = min(max(x0_band(3), lb_band(3)), ub_band(3));

            lb_inv   = [-200,  -200, minAllowableCenter, minAllowableSigma];
            ub_inv   = [ 200,  -0.1, maxAllowableCenter, maxAllowableSigma];
            x0_inv   = [   0,-rangeVal, xtrainvals(minidx), speedRange/4  ];
            x0_inv(3) = min(max(x0_inv(3), lb_inv(3)), ub_inv(3));

            try [param_out(1,:),resnorm(1)] = lsqcurvefit(gaussFun,x0_low1, xtrainvals,ytrainvals,lb_low1, ub_low1, optsFit); catch, end
            try [param_out(2,:),resnorm(2)] = lsqcurvefit(gaussFun,x0_low2, xtrainvals,ytrainvals,lb_low2, ub_low2, optsFit); catch, end
            try [param_out(3,:),resnorm(3)] = lsqcurvefit(gaussFun,x0_high1,xtrainvals,ytrainvals,lb_high1,ub_high1,optsFit); catch, end
            try [param_out(4,:),resnorm(4)] = lsqcurvefit(gaussFun,x0_high2,xtrainvals,ytrainvals,lb_high2,ub_high2,optsFit); catch, end
            try [param_out(5,:),resnorm(5)] = lsqcurvefit(gaussFun,x0_band, xtrainvals,ytrainvals,lb_band, ub_band, optsFit); catch, end
            try [param_out(6,:),resnorm(6)] = lsqcurvefit(gaussFun,x0_inv,  xtrainvals,ytrainvals,lb_inv,  ub_inv,  optsFit); catch, end

            [~, fitIdx] = sort(resnorm, 'ascend');

            chosenIdx  = nan;
            bestParams = nan(1,4);
            typeStr    = 'untuned';
            charCode   = 0;
            cvR2       = 0;
            idx_idx    = 1;

            while idx_idx <= 6
                bestIdx         = fitIdx(idx_idx);
                candidateParams = param_out(bestIdx, :);
                if any(~isfinite(candidateParams)), idx_idx = idx_idx + 1; continue; end

                SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
                if SS_tot < 1e-9, SS_tot = 1; end
                candidateR2 = 1 - (resnorm(bestIdx) / SS_tot);

                % Must beat a flat line
                if candidateR2 <= 0, idx_idx = idx_idx + 1; continue; end

                % ── LOW-PASS ────────────────────────────────────────────
                if bestIdx == 1 || bestIdx == 2
                    % Check for secondary peak in bins 2:end
                    [~, ~, ~, secProm] = findpeaks(ytrainvals(3:end));
                    noSecondaryPeak = isempty(secProm) || max(secProm) < rangeVal / 3;
                    if firstBinWasFit && (maxidx == 1) && noSecondaryPeak
                        chosenIdx = bestIdx; bestParams = candidateParams;
                        typeStr = 'lowpass'; charCode = 1;
                        break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end

                % ── HIGH-PASS ───────────────────────────────────────────
                elseif bestIdx == 3 || bestIdx == 4
                    if lastBinWasFit && (ytrainvals(end) >= 0.80 * maxVal)
                        chosenIdx = bestIdx; bestParams = candidateParams;
                        typeStr = 'highpass'; charCode = 2;
                        break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end

                % ── BAND-PASS ───────────────────────────────────────────
                elseif bestIdx == 5
                    fitY = gaussFun(candidateParams, xDense);
                    [~, ~, ~, prom] = findpeaks(fitY);
                    if ~isempty(prom) && prom(1) >= rangeVal / 3
                        chosenIdx = bestIdx; bestParams = candidateParams;
                        typeStr = 'bandpass'; charCode = 3;
                        break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end

                % ── TROUGH-INVERTED ─────────────────────────────────────
                elseif bestIdx == 6
                    fitY = gaussFun(candidateParams, xDense);
                    [~, ~, ~, prom] = findpeaks(-fitY);
                    if ~isempty(prom) && prom(1) >= rangeVal / 3
                        chosenIdx = bestIdx; bestParams = candidateParams;
                        typeStr = 'trough_inverted'; charCode = 4;
                        break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                end
            end % while

            % ── CV-R² on held-out 20% ───────────────────────────────────
            if ~isnan(chosenIdx)
                yPred  = gaussFun(bestParams, xtestvals);
                SS_res = sum((ytestvals - yPred).^2);
                SS_tot = sum((ytestvals - mean(ytestvals)).^2);
                if SS_tot < 1e-9, SS_tot = 1; end
                cvR2 = 1 - (SS_res / SS_tot);
                cvR2 = max(cvR2, -2);
            end

            if charCode == 1
                prefSpeedVal = minSpeed;
            elseif charCode == 2
                prefSpeedVal = maxSpeed;
            elseif charCode == 3 || charCode == 4
                pVal = bestParams(3);
                pVal = max(min(pVal, maxSpeed), minSpeed);
                prefSpeedVal = pVal;
            else
                prefSpeedVal = nan;
                bestParams   = [nan nan nan nan];
                cvR2         = 0;
                typeStr      = 'untuned';
                charCode     = 0;
            end

            cls.tuningType{r}     = typeStr;
            cls.tuningCode(r)     = charCode;
            cls.R2(r)             = cvR2;
            cls.fitParams(r, :)   = bestParams;
            cls.preferredSpeed(r) = prefSpeedVal;

        end % ROI loop

        response.(targetStruct).(useField).classification_8020 = cls;

        allTypes = cls.tuningType;
        nLP  = sum(strcmp(allTypes, 'lowpass'));
        nHP  = sum(strcmp(allTypes, 'highpass'));
        nBP  = sum(strcmp(allTypes, 'bandpass'));
        nTR  = sum(strcmp(allTypes, 'trough_inverted'));
        nUN  = sum(strcmp(allTypes, 'untuned'));
        nCls = nLP + nHP + nBP + nTR;

        fprintf('\n--- Edd-style v3 80/20 final Summary [%s] ---\n', targetStruct);
        fprintf('  Total ROIs evaluated : %d\n', numROIs);
        fprintf('  ─────────────────────────────────\n');
        fprintf('  Low-Pass   (decay)   : %d  (%.1f%%)\n', nLP, 100*nLP/max(nCls,1));
        fprintf('  High-Pass  (rise)    : %d  (%.1f%%)\n', nHP, 100*nHP/max(nCls,1));
        fprintf('  Band-Pass  (peak)    : %d  (%.1f%%)\n', nBP, 100*nBP/max(nCls,1));
        fprintf('  Trough     (inverted): %d  (%.1f%%)\n', nTR, 100*nTR/max(nCls,1));
        fprintf('  Untuned              : %d\n', nUN);
        fprintf('─────────────────────────────────────\n\n');

    end % framework loop

    save(fileFullPath, 'response', '-append');
    fprintf('v3 80/20 final classifier finished.\n');
end
