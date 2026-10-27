function response = classifySpeedTuningFromCurves(sessionFileInfo, response, useField, doSave)
% CLASSIFYSPEEDTUNINGFROMCURVES  (gaussCompete_linearSpace_v5)  [FINAL]
%
% v5 (Sep 2026): Gaussian templates fitted to the mean response in each running
% speed bin, with speed in cm/s (bin centres), using the acceptance rule of
% Horrocks et al. (2024):
%   - templates are ranked by residual error
%   - the best-fitting low-pass or high-pass template is accepted directly
%   - band-pass and trough templates are accepted only if the peak/trough of the
%     fitted curve has a prominence of at least a third of the data range
% No additional shape checks and no R2 threshold.
% Selection of speed-tuned boutons (pValMoving) is applied in the population script.
%
% Labels:
%   'lowpass', 'highpass', 'bandpass', 'trough_inverted'
%   'untuned'    - fewer than 4 valid speed bins (cannot be fitted)
%   'fit_failed' - all six template fits failed (check the printed error)
% Previous version: gaussCompete_linearSpace_v4

    if nargin < 3 || isempty(useField), useField = 'dFFNeuropilCorrected'; end
    if nargin < 4 || isempty(doSave), doSave = true; end

    stimFileName = sprintf('%s_%s_Response_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
    fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

    frameworks = {'tuningCurve', 'tuningCurveFixedBins'};   % remove the second to skip fixed bins

    gaussFun = @(params, xdata) params(1) + params(2) .* ...
        exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
    optsFit = optimset('Display', 'off');

    for f = 1:numel(frameworks)
        targetStruct = frameworks{f};
        if ~isfield(response, targetStruct) || ~isfield(response.(targetStruct), useField)
            continue;
        end

        edges         = response.(targetStruct).speedBins;
        movingCenters = edges(1:end-1) + diff(edges) / 2;
        movingCenters = movingCenters(:);
        y_session     = response.(targetStruct).(useField).moveMean;
        numROIs       = size(y_session, 1);
        nBinsTotal    = numel(movingCenters);

        cls = struct();
        cls.tuningType        = repmat({'untuned'}, numROIs, 1);
        cls.tuningCode        = zeros(numROIs, 1);   % 0 untuned, 1 low, 2 high, 3 band, 4 trough, -1 fit failed
        cls.R2                = zeros(numROIs, 1);
        cls.fitParams         = nan(numROIs, 4);
        cls.preferredSpeed    = nan(numROIs, 1);
        cls.usedMovingOnly    = true;
        cls.classifierVersion = "gaussCompete_linearSpace_v5";

        nFitErrors = 0;
        firstError = '';

        fprintf('Running v5 fits on %d ROIs [%s], %d speed bins...\n', numROIs, targetStruct, nBinsTotal);
        if nBinsTotal < 4
            fprintf('  Fewer than 4 speed bins: all ROIs left untuned.\n');
        end

        for r = 1:numROIs
            yraw  = y_session(r, :)';
            valid = ~isnan(yraw);
            if sum(valid) < 4
                continue;   % stays 'untuned'
            end

            xvals = movingCenters(valid);
            yvals = yraw(valid);

            minSpeed   = min(xvals);
            maxSpeed   = max(xvals);
            speedRange = maxSpeed - minSpeed;
            binStep    = mean(diff(xvals));

            minAllowableCenter         = minSpeed + 0.15 * speedRange;
            maxAllowableCenter         = maxSpeed - 0.15 * speedRange;
            minAllowableSigma          = 1.0;
            maxAllowableSigma          = 0.25 * speedRange;
            maxAllowableSigmaMonotonic = speedRange;
            xDense = linspace(minSpeed, maxSpeed, 200);

            [maxVal, maxidx] = max(yvals);
            [minVal, minidx] = min(yvals);
            rangeVal = maxVal - minVal;
            if rangeVal <= 0, rangeVal = 1e-5; end

            %        [base,      amp,               mu,                  sigma]
            lb = { [-200,       0.1, minSpeed - binStep, minAllowableSigma], ...   % 1 low-pass (+)
                   [-200,      -200, maxSpeed - binStep, minAllowableSigma], ...   % 2 low-pass (-)
                   [-200,       0.1, maxSpeed - binStep, minAllowableSigma], ...   % 3 high-pass (+)
                   [-200,      -200, minSpeed - binStep, minAllowableSigma], ...   % 4 high-pass (-)
                   [-200,       0.1, minAllowableCenter, minAllowableSigma], ...   % 5 band-pass
                   [-200,      -200, minAllowableCenter, minAllowableSigma] };     % 6 trough
            ub = { [ 200,       200, minSpeed + binStep, maxAllowableSigmaMonotonic], ...
                   [ 200,      -0.1, maxSpeed + binStep, maxAllowableSigmaMonotonic], ...
                   [ 200,       200, maxSpeed + binStep, maxAllowableSigmaMonotonic], ...
                   [ 200,      -0.1, minSpeed + binStep, maxAllowableSigmaMonotonic], ...
                   [ 200,       200, maxAllowableCenter, maxAllowableSigma], ...
                   [ 200,      -0.1, maxAllowableCenter, maxAllowableSigma] };
            x0 = { [0,  rangeVal, minSpeed,      speedRange / 2], ...
                   [0, -rangeVal, maxSpeed,      speedRange / 2], ...
                   [0,  rangeVal, maxSpeed,      speedRange / 2], ...
                   [0, -rangeVal, minSpeed,      speedRange / 2], ...
                   [0,  rangeVal, xvals(maxidx), speedRange / 4], ...
                   [0, -rangeVal, xvals(minidx), speedRange / 4] };

            param_out = nan(6, 4);
            resnorm   = inf(6, 1);
            for t = 1:6
                x0t = min(max(x0{t}, lb{t}), ub{t});   % keep starting point within bounds
                try
                    [param_out(t, :), resnorm(t)] = lsqcurvefit(gaussFun, x0t, xvals, yvals, lb{t}, ub{t}, optsFit);
                catch ME
                    nFitErrors = nFitErrors + 1;
                    if isempty(firstError)
                        firstError = sprintf('ROI %d, template %d: %s', r, t, ME.message);
                    end
                end
            end

            if all(isinf(resnorm))
                cls.tuningType{r} = 'fit_failed';
                cls.tuningCode(r) = -1;
                continue;
            end

            [~, fitIdx] = sort(resnorm, 'ascend');
            SS_tot = sum((yvals - mean(yvals)).^2);
            if SS_tot < 1e-9, SS_tot = 1; end

            typeStr = 'untuned'; charCode = 0; bestParams = nan(1, 4); fullR2 = 0;

            for idx_idx = 1:6
                bestIdx         = fitIdx(idx_idx);
                candidateParams = param_out(bestIdx, :);
                if any(~isfinite(candidateParams)), continue; end
                candidateR2 = 1 - resnorm(bestIdx) / SS_tot;

                if bestIdx == 1 || bestIdx == 2            % low-pass: accepted if best-fitting
                    typeStr = 'lowpass'; charCode = 1;
                elseif bestIdx == 3 || bestIdx == 4        % high-pass: accepted if best-fitting
                    typeStr = 'highpass'; charCode = 2;
                elseif bestIdx == 5                        % band-pass: prominence check
                    [~, ~, ~, prom] = findpeaks(gaussFun(candidateParams, xDense));
                    if isempty(prom) || prom(1) < rangeVal / 3, continue; end
                    typeStr = 'bandpass'; charCode = 3;
                elseif bestIdx == 6                        % trough: prominence check
                    [~, ~, ~, prom] = findpeaks(-gaussFun(candidateParams, xDense));
                    if isempty(prom) || prom(1) < rangeVal / 3, continue; end
                    typeStr = 'trough_inverted'; charCode = 4;
                end
                bestParams = candidateParams;
                fullR2     = candidateR2;
                break;
            end

            if charCode == 3 || charCode == 4
                prefSpeedVal = max(min(bestParams(3), maxSpeed), minSpeed);
            elseif charCode == 1
                prefSpeedVal = minSpeed;
            elseif charCode == 2
                prefSpeedVal = maxSpeed;
            else
                prefSpeedVal = NaN;
            end

            cls.tuningType{r}     = typeStr;
            cls.tuningCode(r)     = charCode;
            cls.R2(r)             = fullR2;
            cls.fitParams(r, :)   = bestParams;
            cls.preferredSpeed(r) = prefSpeedVal;
        end % ROI loop

        response.(targetStruct).(useField).classification = cls;

        % --- Summary ---
        allTypes = cls.tuningType;
        nLP = sum(strcmp(allTypes, 'lowpass'));
        nHP = sum(strcmp(allTypes, 'highpass'));
        nBP = sum(strcmp(allTypes, 'bandpass'));
        nTR = sum(strcmp(allTypes, 'trough_inverted'));
        nUN = sum(strcmp(allTypes, 'untuned'));
        nFF = sum(strcmp(allTypes, 'fit_failed'));
        nCls = nLP + nHP + nBP + nTR;

        fprintf('\n--- v5 Summary [%s] ---\n', targetStruct);
        fprintf('  Total ROIs evaluated : %d (%d speed bins)\n', numROIs, nBinsTotal);
        fprintf('  Low-Pass   (decay)   : %d  (%.1f%%)\n', nLP, 100 * nLP / max(nCls, 1));
        fprintf('  High-Pass  (rise)    : %d  (%.1f%%)\n', nHP, 100 * nHP / max(nCls, 1));
        fprintf('  Band-Pass  (peak)    : %d  (%.1f%%)\n', nBP, 100 * nBP / max(nCls, 1));
        fprintf('  Trough     (inverted): %d  (%.1f%%)\n', nTR, 100 * nTR / max(nCls, 1));
        fprintf('  Untuned (<4 bins)    : %d\n', nUN);
        fprintf('  Fit failed           : %d\n', nFF);
        if nFitErrors > 0
            fprintf('  WARNING: %d individual template fits errored. First error:\n    %s\n', nFitErrors, firstError);
        end
        fprintf('-----------------------------\n\n');
    end % framework loop

%     save(fileFullPath, 'response', '-append');
%     fprintf('v5 classifier finished.\n');
    if doSave
       save(fileFullPath, 'response', '-append');
       fprintf('v5 classifier finished.\n');
   end
end