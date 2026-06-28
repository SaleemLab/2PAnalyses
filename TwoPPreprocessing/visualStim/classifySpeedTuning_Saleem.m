function response = classifySpeedTuning_Saleem(sessionFileInfo, response, useField)
% CLASSIFYSPEEDTUNING_SALEEM
%
% Saleem et al. 2013 - style speed tuning classification.
%
% DESCRIPTIVE FUNCTION:
%   Asymmetric (split) Gaussian — different sigma on each side of peak:
%   y(s) = ymax * exp( -(s - smax)^2 / sigma(s) )
%   where sigma(s) = sigma_minus if s < smax, sigma_plus if s >= smax
%   Free parameters: ymax, smax, sigma_minus, sigma_plus
%   NOTE: no baseline parameter — function decays to 0 at extremes.
%
% THREE COMPETING CURVES:
%   1. High-pass:  smax constrained >= 30 cm/s  (monotonically increasing)
%   2. Low-pass:   smax constrained <= 1  cm/s  (monotonically decreasing)
%   3. Band-pass:  smax unconstrained
%
% MODEL SELECTION (80/20 split, seed 42):
%   Fit on 80% of frames, evaluate variance explained on held-out 20%.
%   Band-pass wins only if:
%     (a) its CV variance explained > both monotonic curves, AND
%     (b) 2 < smax < 25 cm/s
%   Otherwise the better-fitting monotonic curve wins.
%
% CLASSIFICATION:
%   highpass  — high-pass curve wins
%   lowpass   — low-pass curve wins
%   bandpass  — band-pass curve wins (with peak location gate)
%   untuned   — fewer than 4 valid bins

    if nargin < 3 || isempty(useField), useField = 'dFFNeuropilCorrected'; end

    stimFileName = sprintf('%s_%s_Response_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
    fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

    % ── Load raw signals and wheel for 80/20 split ──────────────────────
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

    fs        = 60;
    tickToCm  = 3.1415 * 20 / 1024;
    wheelSpeed = [0; diff(dataPeriph.peripheralData.Wheel.Value * tickToCm)] ./ ...
                 [1; diff(dataPeriph.peripheralData.Wheel.sampleTimes)];
    wheelSpeed(abs(wheelSpeed) > 150) = NaN;
    totalFrames = length(wheelSpeed);

    % ── 80/20 contiguous block split (seed 42) ──────────────────────────
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

    % ── Asymmetric Gaussian (split Gaussian) ────────────────────────────
    % y(s) = ymax * exp( -(s-smax)^2 / sigma(s) )
    % sigma(s) = sigma_minus if s < smax, sigma_plus if s >= smax
    % params = [ymax, smax, sigma_minus, sigma_plus]
    splitGauss = @(params, xdata) params(1) .* exp( ...
        -( (xdata - params(2)).^2 ) ./ ...
        ( (xdata < params(2)) * params(3) + (xdata >= params(2)) * params(4) ) );

    optsFit = optimset('Display', 'off');

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

        minSpeed   = min(movingCenters);
        maxSpeed   = max(movingCenters);
        speedRange = maxSpeed - minSpeed;

        cls = struct();
        cls.tuningType        = cell(numROIs, 1);
        cls.tuningCode        = nan(numROIs, 1);
        cls.R2                = nan(numROIs, 1);   % CV variance explained on 20%
        cls.fitParams         = nan(numROIs, 4);   % [ymax, smax, sigma_minus, sigma_plus]
        cls.preferredSpeed    = nan(numROIs, 1);   % smax in cm/s
        cls.usedMovingOnly    = true;
        cls.classifierVersion = "saleem2013_splitGauss_8020";

        fprintf('Running Saleem 2013 split-Gaussian classifier on %d ROIs [%s]...\n', numROIs, targetStruct);

        for r = 1:numROIs

            % ── Build train/test tuning curves ───────────────────────────
            yTrain = nan(numBins, 1);
            yTest  = nan(numBins, 1);

            for b = 1:numBins
                binIdx    = wheelSpeed >= edges(b) & wheelSpeed < edges(b+1) & isMoving;
                yTrain(b) = mean(sigData(r, binIdx & trainMask), 'omitnan');
                yTest(b)  = mean(sigData(r, binIdx & testMask),  'omitnan');
            end

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

            [maxVal, ~] = max(ytrainvals);
            [~, minidx] = min(ytrainvals);
            rangeVal    = maxVal - min(ytrainvals);
            if rangeVal <= 0, rangeVal = 1e-5; end

            % ── Fit 3 curves on training data ────────────────────────────
            % params = [ymax, smax, sigma_minus, sigma_plus]
            % Bounds chosen to allow flexible but physiologically sensible fits
            % sigma in units of cm/s^2 (note: no sqrt in exponent per Saleem eq)

            sigMin = 1;
            sigMax = speedRange^2 * 4;   % generous upper bound

            % 1. High-pass: smax >= 30 cm/s
            lb_hp  = [0,    30,         sigMin, sigMin];
            ub_hp  = [maxVal*3, maxSpeed*3,  sigMax, sigMax];
            x0_hp  = [maxVal, maxSpeed,   speedRange, speedRange];
            try
                [p_hp, rn_hp] = lsqcurvefit(splitGauss, x0_hp, xtrainvals, ytrainvals, lb_hp, ub_hp, optsFit);
            catch
                p_hp = nan(1,4); rn_hp = inf;
            end

            % 2. Low-pass: smax <= 1 cm/s
            lb_lp  = [0,    -speedRange,  sigMin, sigMin];
            ub_lp  = [maxVal*3, 1,         sigMax, sigMax];
            x0_lp  = [maxVal, minSpeed,   speedRange, speedRange];
            try
                [p_lp, rn_lp] = lsqcurvefit(splitGauss, x0_lp, xtrainvals, ytrainvals, lb_lp, ub_lp, optsFit);
            catch
                p_lp = nan(1,4); rn_lp = inf;
            end

            % 3. Band-pass: smax unconstrained
            lb_bp  = [0,    minSpeed,    sigMin, sigMin];
            ub_bp  = [maxVal*3, maxSpeed, sigMax, sigMax];
            x0_bp  = [maxVal, xtrainvals(find(ytrainvals == maxVal, 1)), speedRange/2, speedRange/2];
            x0_bp(2) = min(max(x0_bp(2), lb_bp(2)), ub_bp(2));
            try
                [p_bp, rn_bp] = lsqcurvefit(splitGauss, x0_bp, xtrainvals, ytrainvals, lb_bp, ub_bp, optsFit);
            catch
                p_bp = nan(1,4); rn_bp = inf;
            end

            % ── Evaluate CV variance explained on held-out 20% ──────────
            SS_tot = sum((ytestvals - mean(ytestvals)).^2);
            if SS_tot < 1e-9, SS_tot = 1; end

            cvR2_hp = -inf; cvR2_lp = -inf; cvR2_bp = -inf;
            if all(isfinite(p_hp))
                SS_res = sum((ytestvals - splitGauss(p_hp, xtestvals)).^2);
                cvR2_hp = max(1 - SS_res/SS_tot, -2);
            end
            if all(isfinite(p_lp))
                SS_res = sum((ytestvals - splitGauss(p_lp, xtestvals)).^2);
                cvR2_lp = max(1 - SS_res/SS_tot, -2);
            end
            if all(isfinite(p_bp))
                SS_res = sum((ytestvals - splitGauss(p_bp, xtestvals)).^2);
                cvR2_bp = max(1 - SS_res/SS_tot, -2);
            end

            % ── Classification (Saleem 2013 rules) ──────────────────────
            % Band-pass wins ONLY if:
            %   (a) cvR2_bp > cvR2_hp AND cvR2_bp > cvR2_lp
            %   (b) 2 < smax < 25 cm/s
            % Otherwise the better monotonic curve wins.

            bandpassPeakOK = all(isfinite(p_bp)) && p_bp(2) > 2 && p_bp(2) < 25;
            bandpassWins   = bandpassPeakOK && cvR2_bp > cvR2_hp && cvR2_bp > cvR2_lp;

            if bandpassWins
                typeStr  = 'bandpass';  charCode = 3;
                bestP    = p_bp;        cvR2     = cvR2_bp;
                prefSpeedVal = p_bp(2);
            elseif cvR2_hp >= cvR2_lp
                typeStr  = 'highpass';  charCode = 2;
                bestP    = p_hp;        cvR2     = cvR2_hp;
                prefSpeedVal = maxSpeed;
            else
                typeStr  = 'lowpass';   charCode = 1;
                bestP    = p_lp;        cvR2     = cvR2_lp;
                prefSpeedVal = minSpeed;
            end

            % If no fit was valid, mark untuned
            if ~any(isfinite([cvR2_hp, cvR2_lp, cvR2_bp]))
                typeStr = 'untuned'; charCode = 0;
                bestP   = [nan nan nan nan]; cvR2 = 0;
                prefSpeedVal = nan;
            end

            cls.tuningType{r}     = typeStr;
            cls.tuningCode(r)     = charCode;
            cls.R2(r)             = cvR2;
            cls.fitParams(r, :)   = bestP;
            cls.preferredSpeed(r) = prefSpeedVal;

        end % ROI loop

        response.(targetStruct).(useField).classification_saleem = cls;

        allTypes = cls.tuningType;
        nLP  = sum(strcmp(allTypes, 'lowpass'));
        nHP  = sum(strcmp(allTypes, 'highpass'));
        nBP  = sum(strcmp(allTypes, 'bandpass'));
        nUN  = sum(strcmp(allTypes, 'untuned'));
        nCls = nLP + nHP + nBP;

        fprintf('\n--- Saleem 2013 Split-Gaussian Summary [%s] ---\n', targetStruct);
        fprintf('  Total ROIs evaluated : %d\n', numROIs);
        fprintf('  ─────────────────────────────────\n');
        fprintf('  Low-Pass   (decay)   : %d  (%.1f%%)\n', nLP, 100*nLP/max(nCls,1));
        fprintf('  High-Pass  (rise)    : %d  (%.1f%%)\n', nHP, 100*nHP/max(nCls,1));
        fprintf('  Band-Pass  (peak)    : %d  (%.1f%%)\n', nBP, 100*nBP/max(nCls,1));
        fprintf('  Untuned              : %d\n', nUN);
        fprintf('─────────────────────────────────────\n\n');

    end % framework loop

    save(fileFullPath, 'response', '-append');
    fprintf('Saleem 2013 classifier finished.\n');
end
