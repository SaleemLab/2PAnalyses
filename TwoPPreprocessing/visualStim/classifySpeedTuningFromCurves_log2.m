function response = classifySpeedTuningFromCurves_log2(sessionFileInfo, response, useField)
% CLASSIFYSPEEDTUNINGFROMCURVES_LOG2
%
% Multi-Template Competitive Gaussian Fit in log2 speed space.
%
% KEY FEATURES:
% ─────────────────────────────────────────────────────────────────────────
%   • Speeds log2-transformed before fitting — symmetric templates
%   • Low-pass:  maxidx==1 + firstBinWasFit + no secondary peak (bins 3:end)
%   • High-pass: lastBinWasFit + last bin >= 80% of max
%   • Band-pass: prominence gate on fit (rangeVal/3) AND raw data (rangeVal/2)
%   • Trough:    same as band-pass but inverted
%   • candidateR2 > 0: fit must beat flat line
%   • preferredSpeed saved in linear cm/s (2^mu)
%   • fitParams saved in log2 space
% ─────────────────────────────────────────────────────────────────────────

    if nargin < 3 || isempty(useField), useField = 'dFFNeuropilCorrected'; end

    stimFileName = sprintf('%s_%s_Response_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
    fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

    frameworks = {'tuningCurve', 'tuningCurveFixedBins'};

    for f = 1:numel(frameworks)
        targetStruct = frameworks{f};
        if ~isfield(response, targetStruct) || ...
           ~isfield(response.(targetStruct), useField), continue; end

        edges         = response.(targetStruct).speedBins;
        movingCenters = (edges(1:end-1) + diff(edges)/2)';
        y_session     = response.(targetStruct).(useField).moveMean;
        numROIs       = size(y_session, 1);

        % ── Log2-transform the speed axis ────────────────────────────────
        movingCenters_log2 = log2(max(movingCenters, 0.01));

        cls = struct();
        cls.tuningType        = cell(numROIs, 1);
        cls.tuningCode        = nan(numROIs, 1);
        cls.R2                = nan(numROIs, 1);
        cls.fitParams         = nan(numROIs, 4);   % in log2 space
        cls.preferredSpeed    = nan(numROIs, 1);   % in linear cm/s
        cls.usedMovingOnly    = true;
        cls.classifierVersion = "log2_v3_fullgates";

        gaussFun = @(params, xdata) params(1) + params(2) .* ...
            exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
        optsFit = optimset('Display', 'off');

        minSpeed   = min(movingCenters_log2);
        maxSpeed   = max(movingCenters_log2);
        speedRange = maxSpeed - minSpeed;

        customOffset               = mean(diff(movingCenters_log2));
        minAllowableCenter         = minSpeed  + (0.15 * speedRange);
        maxAllowableCenter         = maxSpeed  - (0.15 * speedRange);
        minAllowableSigma          = 0.1;
        maxAllowableSigma          = 0.25 * speedRange;
        maxAllowableSigmaMonotonic = speedRange;

        xDense = linspace(minSpeed, maxSpeed, 200);

        fprintf('Running log2 v3 classifier on %d ROIs [%s]...\n', numROIs, targetStruct);

        for r = 1:numROIs

            yraw  = y_session(r, :)';
            xall  = movingCenters_log2;
            valid = ~isnan(yraw);

            % Track whether boundary bins were actually fitted
            firstBinWasFit = valid(1);
            lastBinWasFit  = valid(end);

            if sum(valid) < 4
                cls.tuningType{r} = 'untuned'; cls.tuningCode(r) = 0;
                cls.R2(r) = 0; continue;
            end

            xvals = xall(valid);
            yvals = yraw(valid);

            [maxVal, maxidx] = max(yvals);
            [minVal, minidx] = min(yvals);
            rangeVal = maxVal - minVal;
            if rangeVal <= 0, rangeVal = 1e-5; end

            % ── 6 competing templates in log2 space ─────────────────────
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
            x0_band  = [   0, rangeVal, xvals(maxidx),  speedRange/4      ];
            x0_band(3) = min(max(x0_band(3), lb_band(3)), ub_band(3));

            lb_inv   = [-200,  -200, minAllowableCenter, minAllowableSigma];
            ub_inv   = [ 200,  -0.1, maxAllowableCenter, maxAllowableSigma];
            x0_inv   = [   0,-rangeVal, xvals(minidx),  speedRange/4      ];
            x0_inv(3) = min(max(x0_inv(3), lb_inv(3)), ub_inv(3));

            try [param_out(1,:),resnorm(1)] = lsqcurvefit(gaussFun,x0_low1, xvals,yvals,lb_low1, ub_low1, optsFit); catch, end
            try [param_out(2,:),resnorm(2)] = lsqcurvefit(gaussFun,x0_low2, xvals,yvals,lb_low2, ub_low2, optsFit); catch, end
            try [param_out(3,:),resnorm(3)] = lsqcurvefit(gaussFun,x0_high1,xvals,yvals,lb_high1,ub_high1,optsFit); catch, end
            try [param_out(4,:),resnorm(4)] = lsqcurvefit(gaussFun,x0_high2,xvals,yvals,lb_high2,ub_high2,optsFit); catch, end
            try [param_out(5,:),resnorm(5)] = lsqcurvefit(gaussFun,x0_band, xvals,yvals,lb_band, ub_band, optsFit); catch, end
            try [param_out(6,:),resnorm(6)] = lsqcurvefit(gaussFun,x0_inv,  xvals,yvals,lb_inv,  ub_inv,  optsFit); catch, end

            [~, fitIdx] = sort(resnorm, 'ascend');

            bestParams = nan(1,4);
            typeStr    = 'untuned';
            charCode   = 0;
            fullR2     = 0;
            idx_idx    = 1;

            while idx_idx <= 6
                bestIdx         = fitIdx(idx_idx);
                candidateParams = param_out(bestIdx, :);
                if any(~isfinite(candidateParams)), idx_idx = idx_idx + 1; continue; end

                SS_tot = sum((yvals - mean(yvals)).^2);
                if SS_tot < 1e-9, SS_tot = 1; end
                candidateR2 = 1 - (resnorm(bestIdx) / SS_tot);

                % Must beat a flat line
                if candidateR2 <= 0, idx_idx = idx_idx + 1; continue; end

                % ── LOW-PASS ────────────────────────────────────────────
                if bestIdx == 1 || bestIdx == 2
                    [~, ~, ~, secProm] = findpeaks(yvals(3:end));
                    noSecondaryPeak = isempty(secProm) || max(secProm) < rangeVal / 3;
                    if firstBinWasFit && (maxidx == 1) && noSecondaryPeak
                        bestParams = candidateParams;
                        typeStr = 'lowpass'; charCode = 1; fullR2 = candidateR2;
                        break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end

                % ── HIGH-PASS ───────────────────────────────────────────
                elseif bestIdx == 3 || bestIdx == 4
                    if lastBinWasFit && (yvals(end) >= 0.80 * maxVal)
                        bestParams = candidateParams;
                        typeStr = 'highpass'; charCode = 2; fullR2 = candidateR2;
                        break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end

          
                    % ── BAND-PASS ───────────────────────────────────────────
                elseif bestIdx == 5
                    fitY = gaussFun(candidateParams, xDense);
                    [~, ~, ~, prom]    = findpeaks(fitY);
                    [~, ~, ~, rawProm] = findpeaks(yvals);
                    fitOK         = ~isempty(prom)    && prom(1)      >= rangeVal / 3;
                    rawOK         = ~isempty(rawProm) && max(rawProm) >= rangeVal / 2;
                    sigmaOK       = candidateParams(4) >= 2 * mean(diff(xvals));
                    sufficientAmp = rangeVal > 0.03;
                    if fitOK && rawOK && sigmaOK && sufficientAmp
                        bestParams = candidateParams;
                        typeStr = 'bandpass'; charCode = 3; fullR2 = candidateR2;
                        break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end

                    % ── TROUGH-INVERTED ─────────────────────────────────────
                elseif bestIdx == 6
                    fitY = gaussFun(candidateParams, xDense);
                    [~, ~, ~, prom]    = findpeaks(-fitY);
                    [~, ~, ~, rawProm] = findpeaks(-yvals);
                    fitOK         = ~isempty(prom)    && prom(1)      >= rangeVal / 3;
                    rawOK         = ~isempty(rawProm) && max(rawProm) >= rangeVal / 2;
                    sigmaOK       = candidateParams(4) >= 2 * mean(diff(xvals));
                    sufficientAmp = rangeVal > 0.03;
                    if fitOK && rawOK && sigmaOK && sufficientAmp
                        bestParams = candidateParams;
                        typeStr = 'trough_inverted'; charCode = 4; fullR2 = candidateR2;
                        break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                end
            end % while

            % ── Preferred speed — convert back to linear cm/s ────────────
            if charCode == 1
                prefSpeedVal = 2^minSpeed;
            elseif charCode == 2
                prefSpeedVal = 2^maxSpeed;
            elseif charCode == 3 || charCode == 4
                pVal_log2 = bestParams(3);
                pVal_log2 = max(min(pVal_log2, maxSpeed), minSpeed);
                prefSpeedVal = 2^pVal_log2;
            else
                prefSpeedVal = nan;
                bestParams   = [nan nan nan nan];
                fullR2       = 0;
                typeStr      = 'untuned';
                charCode     = 0;
            end

            cls.tuningType{r}     = typeStr;
            cls.tuningCode(r)     = charCode;
            cls.R2(r)             = fullR2;
            cls.fitParams(r, :)   = bestParams;
            cls.preferredSpeed(r) = prefSpeedVal;

        end % ROI loop

        response.(targetStruct).(useField).classification_log2 = cls;

        allTypes = cls.tuningType;
        nLP  = sum(strcmp(allTypes, 'lowpass'));
        nHP  = sum(strcmp(allTypes, 'highpass'));
        nBP  = sum(strcmp(allTypes, 'bandpass'));
        nTR  = sum(strcmp(allTypes, 'trough_inverted'));
        nUN  = sum(strcmp(allTypes, 'untuned'));
        nCls = nLP + nHP + nBP + nTR;

        fprintf('\n--- Log2 v3 Speed Tuning Summary [%s] ---\n', targetStruct);
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
    fprintf('Log2 v3 classifier finished.\n');
end
