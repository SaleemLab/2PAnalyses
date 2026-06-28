% function response = classifySpeedTuningFromCurves(sessionFileInfo, response, useField)
% % CLASSIFYSPEEDTUNINGFROMCURVES
% % Classifies speed tuning phenotypes by fitting Gaussian templates in
% % log-speed space, following Saleem et al. (2013).
% % Preferred speeds are back-transformed and stored in cm/s.
% 
% if nargin < 3 || isempty(useField), useField = 'dFFNeuropilCorrected'; end
% 
% stimFileName = sprintf('%s_%s_Response_%s.mat', ...
%     sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
% fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
% frameworks = {'tuningCurve', 'tuningCurveFixedBins'};
% 
% for f = 1:numel(frameworks)
%     targetStruct = frameworks{f};
%     if ~isfield(response, targetStruct) || ~isfield(response.(targetStruct), useField), continue; end
% 
%     edges = response.(targetStruct).speedBins;
%     movingCenters_linear = (edges(1:end-1) + diff(edges)/2)';
% 
%     % LOG-TRANSFORM SPEED AXIS ---
%     movingCenters = log(movingCenters_linear);
% 
%     y_session = response.(targetStruct).(useField).moveMean;
%     numROIs   = size(y_session, 1);
% 
%     % Preallocate outputs
%     cls = struct();
%     cls.tuningType        = cell(numROIs, 1);
%     cls.tuningCode        = nan(numROIs, 1);
%     cls.R2                = nan(numROIs, 1);
%     cls.fitParams         = nan(numROIs, 4);
%     cls.preferredSpeed    = nan(numROIs, 1);
%     cls.usedMovingOnly    = true;
%     cls.classifierVersion = "log_space_gaussian_v1";
% 
%     gaussFun = @(params, xdata) params(1) + params(2) .* ...
%         exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
%     optsFit = optimset('Display', 'off');
% 
%     % --- BOUNDS IN LOG-SPACE ---
%     minSpeed   = min(movingCenters);
%     maxSpeed   = max(movingCenters);
%     speedRange = maxSpeed - minSpeed;
%     customOffset = 0.5 * (speedRange / range(movingCenters_linear)) ;
% 
%     minAllowableCenter = minSpeed + (0.15 * speedRange);
%     maxAllowableCenter = maxSpeed - (0.15 * speedRange);
%     minAllowableSigma  = 0.05;                % log-space sigma lower bound
%     maxAllowableSigma  = 0.5 * speedRange;    % log-space sigma upper bound
% 
%     fprintf('Running Log-Space Gaussian Fits on %d ROIs [%s]...\n', numROIs, targetStruct);
% 
%     for r = 1:numROIs
%         ytrainvals = y_session(r, :)';
%         xtrainvals = movingCenters;
%         validTrain = ~isnan(ytrainvals);
% 
%         if sum(validTrain) < 4
%             cls.tuningType{r} = 'untuned'; cls.tuningCode(r) = 0;
%             cls.R2(r) = 0; continue;
%         end
% 
%         xtrainvals = xtrainvals(validTrain);
%         ytrainvals = ytrainvals(validTrain);
% 
%         [maxVal, maxidx] = max(ytrainvals);
%         [minVal, minidx] = min(ytrainvals);
%         rangeVal = maxVal - minVal;
%         if rangeVal <= 0, rangeVal = 1e-5; end
% 
%         % Pearson correlation in log-space
%         speedCorrelation = corr(xtrainvals, ytrainvals, ...
%             'rows', 'complete', 'type', 'Pearson');
% 
%         param_out = nan(6, 4);
%         resnorm   = inf(6, 1);
%         co        = customOffset;
% 
%         % --- BOUNDS (all in log-space) ---
%         lb_low1  = [-200,   0.1, minSpeed - co*10,  minAllowableSigma];
%         ub_low1  = [ 200,   200, minSpeed + co,      speedRange*5     ];
%         x0_low1  = [   0, rangeVal, minSpeed - co*2, speedRange/2     ];
% 
%         lb_low2  = [-200,  -200, maxSpeed - co,      minAllowableSigma];
%         ub_low2  = [ 200,  -0.1, maxSpeed + co*10,   speedRange*5     ];
%         x0_low2  = [   0, -rangeVal, maxSpeed + co*2, speedRange/2    ];
% 
%         lb_high1 = [-200,   0.1, maxSpeed - co,      minAllowableSigma];
%         ub_high1 = [ 200,   200, maxSpeed + co*10,   speedRange*5     ];
%         x0_high1 = [   0, rangeVal, maxSpeed + co*2, speedRange/2     ];
% 
%         lb_high2 = [-200,  -200, minSpeed - co*10,   minAllowableSigma];
%         ub_high2 = [ 200,  -0.1, minSpeed + co,      speedRange*5     ];
%         x0_high2 = [   0, -rangeVal, minSpeed - co*2, speedRange/2    ];
% 
%         lb_band  = [-200,   0.1, minAllowableCenter, minAllowableSigma];
%         ub_band  = [ 200,   200, maxAllowableCenter,  maxAllowableSigma];
%         x0_band  = [   0, rangeVal, xtrainvals(maxidx), speedRange/4  ];
% 
%         lb_inv   = [-200,  -200, minAllowableCenter, minAllowableSigma];
%         ub_inv   = [ 200,  -0.1, maxAllowableCenter,  maxAllowableSigma];
%         x0_inv   = [   0, -rangeVal, xtrainvals(minidx), speedRange/4 ];
% 
%         x0_band(3) = min(max(x0_band(3), lb_band(3)), ub_band(3));
%         x0_inv(3)  = min(max(x0_inv(3),  lb_inv(3)),  ub_inv(3));
% 
%         try [param_out(1,:), resnorm(1)] = lsqcurvefit(gaussFun, x0_low1,  xtrainvals, ytrainvals, lb_low1,  ub_low1,  optsFit); catch, end
%         try [param_out(2,:), resnorm(2)] = lsqcurvefit(gaussFun, x0_low2,  xtrainvals, ytrainvals, lb_low2,  ub_low2,  optsFit); catch, end
%         try [param_out(3,:), resnorm(3)] = lsqcurvefit(gaussFun, x0_high1, xtrainvals, ytrainvals, lb_high1, ub_high1, optsFit); catch, end
%         try [param_out(4,:), resnorm(4)] = lsqcurvefit(gaussFun, x0_high2, xtrainvals, ytrainvals, lb_high2, ub_high2, optsFit); catch, end
%         try [param_out(5,:), resnorm(5)] = lsqcurvefit(gaussFun, x0_band,  xtrainvals, ytrainvals, lb_band,  ub_band,  optsFit); catch, end
%         try [param_out(6,:), resnorm(6)] = lsqcurvefit(gaussFun, x0_inv,   xtrainvals, ytrainvals, lb_inv,   ub_inv,   optsFit); catch, end
% 
%         [~, fitIdx] = sort(resnorm, 'ascend');
%         chosenIdx = nan; bestParams = nan(1,4);
%         typeStr = 'untuned'; charCode = 0;
% 
%         idx_idx = 1;
%         while idx_idx <= 6
%             bestIdx = fitIdx(idx_idx);
%             candidateParams = param_out(bestIdx, :);
%             if any(~isfinite(candidateParams)), idx_idx = idx_idx + 1; continue; end
% 
%             if bestIdx == 1 || bestIdx == 2
%                 if ~isnan(speedCorrelation) && speedCorrelation < -0.35
%                     chosenIdx = bestIdx; bestParams = candidateParams;
%                     typeStr = 'lowpass'; charCode = 1; break;
%                 else
%                     idx_idx = idx_idx + 1; continue;
%                 end
% 
%             elseif bestIdx == 3 || bestIdx == 4
%                 % High-Pass: Must correlate positively AND last bin must be near maximum
%                 lastBinIsMax = (ytrainvals(end) >= 0.80 * max(ytrainvals));
%                 if ~isnan(speedCorrelation) && speedCorrelation > 0.50 && lastBinIsMax
%                     chosenIdx = bestIdx; bestParams = candidateParams;
%                     typeStr = 'highpass'; charCode = 2; break;
%                 else
%                     idx_idx = idx_idx + 1; continue;
%                 end
% 
%             elseif bestIdx == 5
%                 evalOfFun = feval(gaussFun, candidateParams, ...
%                     linspace(minSpeed, maxSpeed, 100));
%                 [~, ~, ~, prom] = findpeaks(evalOfFun);
%                 if isempty(prom) || prom < (rangeVal / 3)
%                     idx_idx = idx_idx + 1; continue;
%                 end
%                 chosenIdx = bestIdx; bestParams = candidateParams;
%                 typeStr = 'bandpass'; charCode = 3; break;
% 
%             elseif bestIdx == 6
%                 evalOfFun = -1 .* feval(gaussFun, candidateParams, ...
%                     linspace(minSpeed, maxSpeed, 100));
%                 [~, ~, ~, prom] = findpeaks(evalOfFun);
%                 if isempty(prom) || prom < (rangeVal / 3)
%                     idx_idx = idx_idx + 1; continue;
%                 end
%                 chosenIdx = bestIdx; bestParams = candidateParams;
%                 typeStr = 'trough_inverted'; charCode = 4; break;
%             end
%         end
% 
%         % --- R2 ---
%         if isnan(chosenIdx)
%             typeStr = 'untuned'; charCode = 0;
%             bestParams = [nan nan nan nan]; fullR2 = 0;
%         else
%             SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
%             if SS_tot < 1e-9, SS_tot = 1; end
%             fullR2 = 1 - (resnorm(chosenIdx) / SS_tot);
%         end
% 
%         % --- PREFERRED SPEED: back-transform to cm/s ---
%         if charCode == 1
%             prefSpeedVal = exp(minSpeed);
%         elseif charCode == 2
%             prefSpeedVal = exp(maxSpeed);
%         elseif charCode == 3 || charCode == 4
%             pVal = bestParams(3);
%             pVal = min(max(pVal, minSpeed), maxSpeed);
%             prefSpeedVal = exp(pVal);         % back to cm/s
%         else
%             prefSpeedVal = nan;
%         end
% 
%         cls.tuningType{r}     = typeStr;
%         cls.tuningCode(r)     = charCode;
%         cls.R2(r)             = fullR2;
%         cls.fitParams(r, :)   = bestParams;   % stored in log-space
%         cls.preferredSpeed(r) = prefSpeedVal; % stored in cm/s
%     end
% 
%     response.(targetStruct).(useField).classification = cls;
% end
% 
% save(fileFullPath, 'response', '-append');
% end
function response = classifySpeedTuningFromCurves(sessionFileInfo, response, useField)
% CLASSIFYSPEEDTUNINGFROMCURVES
% Classifies speed tuning phenotypes by fitting Gaussian templates to physical speeds.
% - High-pass: stricter gate requiring last bin near maximum and r > 0.50
% - Band-pass: strict geometric + raw data criteria replacing simple findpeaks prominence
% - Trough-inverted: symmetric strict criteria to band-pass
    if nargin < 3 || isempty(useField), useField = 'dFFNeuropilCorrected'; end
    stimFileName = sprintf('%s_%s_Response_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
    fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
    frameworks = {'tuningCurve', 'tuningCurveFixedBins'};
    for f = 1:numel(frameworks)
        targetStruct = frameworks{f};
        if ~isfield(response, targetStruct) || ~isfield(response.(targetStruct), useField), continue; end
        edges         = response.(targetStruct).speedBins;
        movingCenters = (edges(1:end-1) + diff(edges)/2)';
        y_session = response.(targetStruct).(useField).moveMean;
        numROIs   = size(y_session, 1);
        cls = struct();
        cls.tuningType        = cell(numROIs, 1);
        cls.tuningCode        = nan(numROIs, 1);
        cls.R2                = nan(numROIs, 1);
        cls.fitParams         = nan(numROIs, 4);
        cls.preferredSpeed    = nan(numROIs, 1);
        cls.usedMovingOnly    = true;
        cls.classifierVersion = "strict_bandpass_highpass_v2";
        gaussFun = @(params, xdata) params(1) + params(2) .* ...
            exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
        optsFit = optimset('Display', 'off');
        minSpeed   = min(movingCenters);
        maxSpeed   = max(movingCenters);
        speedRange = maxSpeed - minSpeed;
        customOffset = 0.5;
        minAllowableCenter = minSpeed + (0.15 * speedRange);
        maxAllowableCenter = maxSpeed - (0.15 * speedRange);
        minAllowableSigma  = 1.0;
        maxAllowableSigma  = 0.25 * speedRange;
        fprintf('Running Strict-Gated Physical Fits on %d ROIs [%s]...\n', numROIs, targetStruct);
        for r = 1:numROIs
            ytrainvals = y_session(r, :)';
            xtrainvals = movingCenters;
            validTrain = ~isnan(ytrainvals);
            if sum(validTrain) < 4
                cls.tuningType{r} = 'untuned'; cls.tuningCode(r) = 0;
                cls.R2(r) = 0; continue;
            end
            xtrainvals = xtrainvals(validTrain);
            ytrainvals = ytrainvals(validTrain);
            [maxVal, maxidx] = max(ytrainvals);
            [minVal, minidx] = min(ytrainvals);
            rangeVal = maxVal - minVal;
            if rangeVal <= 0, rangeVal = 1e-5; end
            speedCorrelation = corr(xtrainvals, ytrainvals, ...
                'rows', 'complete', 'type', 'Pearson');
            param_out = nan(6, 4);
            resnorm   = inf(6, 1);
            lb_low1  = [-200,   0.1, minSpeed - customOffset*10, minAllowableSigma];  ub_low1  = [200,  200, minSpeed + customOffset,    speedRange*5]; x0_low1  = [0,  rangeVal, minSpeed - customOffset*2, speedRange/2];
            lb_low2  = [-200, -200.0, maxSpeed - customOffset,   minAllowableSigma];  ub_low2  = [200, -0.1, maxSpeed + customOffset*10, speedRange*5]; x0_low2  = [0, -rangeVal, maxSpeed + customOffset*2, speedRange/2];
            lb_high1 = [-200,   0.1, maxSpeed - customOffset,   minAllowableSigma];  ub_high1 = [200,  200, maxSpeed + customOffset*10, speedRange*5]; x0_high1 = [0,  rangeVal, maxSpeed + customOffset*2, speedRange/2];
            lb_high2 = [-200, -200.0, minSpeed - customOffset*10, minAllowableSigma]; ub_high2 = [200, -0.1, minSpeed + customOffset,    speedRange*5]; x0_high2 = [0, -rangeVal, minSpeed - customOffset*2, speedRange/2];
            lb_band  = [-200,   0.1, minAllowableCenter,         minAllowableSigma];  ub_band  = [200,  200, maxAllowableCenter,         maxAllowableSigma]; x0_band  = [0,  rangeVal, xtrainvals(maxidx), speedRange/4];
            lb_inv   = [-200, -200.0, minAllowableCenter,         minAllowableSigma];  ub_inv   = [200, -0.1, maxAllowableCenter,         maxAllowableSigma]; x0_inv   = [0, -rangeVal, xtrainvals(minidx), speedRange/4];
            x0_band(3) = min(max(x0_band(3), lb_band(3)), ub_band(3));
            x0_inv(3)  = min(max(x0_inv(3),  lb_inv(3)),  ub_inv(3));
            try [param_out(1,:), resnorm(1)] = lsqcurvefit(gaussFun, x0_low1,  xtrainvals, ytrainvals, lb_low1,  ub_low1,  optsFit); catch, end
            try [param_out(2,:), resnorm(2)] = lsqcurvefit(gaussFun, x0_low2,  xtrainvals, ytrainvals, lb_low2,  ub_low2,  optsFit); catch, end
            try [param_out(3,:), resnorm(3)] = lsqcurvefit(gaussFun, x0_high1, xtrainvals, ytrainvals, lb_high1, ub_high1, optsFit); catch, end
            try [param_out(4,:), resnorm(4)] = lsqcurvefit(gaussFun, x0_high2, xtrainvals, ytrainvals, lb_high2, ub_high2, optsFit); catch, end
            try [param_out(5,:), resnorm(5)] = lsqcurvefit(gaussFun, x0_band,  xtrainvals, ytrainvals, lb_band,  ub_band,  optsFit); catch, end
            try [param_out(6,:), resnorm(6)] = lsqcurvefit(gaussFun, x0_inv,   xtrainvals, ytrainvals, lb_inv,   ub_inv,   optsFit); catch, end
            [~, fitIdx] = sort(resnorm, 'ascend');
            chosenIdx = nan; bestParams = nan(1,4);
            typeStr = 'untuned'; charCode = 0;
            idx_idx = 1;
            while idx_idx <= 6
                bestIdx         = fitIdx(idx_idx);
                candidateParams = param_out(bestIdx, :);
                if any(~isfinite(candidateParams)), idx_idx = idx_idx + 1; continue; end
                if bestIdx == 1 || bestIdx == 2
                    % Low-pass: correlates negatively AND maximum must be at first bin
                    firstBinIsMax = (maxidx == 1);
                    if ~isnan(speedCorrelation) && speedCorrelation < -0.35 && firstBinIsMax
                        chosenIdx = bestIdx; bestParams = candidateParams;
                        typeStr = 'lowpass'; charCode = 1; break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                elseif bestIdx == 3 || bestIdx == 4
                    % High-pass: must correlate positively AND last bin near maximum
                    lastBinIsMax = (ytrainvals(end) >= 0.80 * max(ytrainvals));
                    if ~isnan(speedCorrelation) && speedCorrelation > 0.50 && lastBinIsMax
                        chosenIdx = bestIdx; bestParams = candidateParams;
                        typeStr = 'highpass'; charCode = 2; break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                elseif bestIdx == 5
                    % Band-pass: strict geometric + raw data criteria
                    fitY = gaussFun(candidateParams, linspace(minSpeed, maxSpeed, 200));
                    % Gaussian peak must sit in the interior 70% of the speed range
                    centerInside      = candidateParams(3) > minAllowableCenter && ...
                                        candidateParams(3) < maxAllowableCenter;
                    % Width of the Gaussian can't be so broad it's essentially flat
                    sigmaAcceptable   = candidateParams(4) <= maxAllowableSigma;
                    % The actual data maximum must not be at the first or last bin.
                    rawPeakInterior   = maxidx > 1 && maxidx < numel(ytrainvals);
                    % The fit must drop by at least 15% of the range on both sides. 
                    % Prevents a monotonic curve being called band-pass.
                    dropsBothSidesFit = (max(fitY) - fitY(1)   > 0.15 * rangeVal) && ...
                                        (max(fitY) - fitY(end) > 0.15 * rangeVal);
                    % Same check but on raw data — 10% drop on both sides.
                    dropsBothSidesRaw = (maxVal - ytrainvals(1)   > 0.10 * rangeVal) && ...
                                        (maxVal - ytrainvals(end) > 0.10 * rangeVal);
                    % Peak must stand at least 15% above the higher of the two endpoints.
                    fitProminence     = max(fitY) - max(fitY(1), fitY(end));
                    if centerInside && sigmaAcceptable && rawPeakInterior && ...
                       dropsBothSidesFit && dropsBothSidesRaw && (fitProminence > 0.15 * rangeVal)
                        chosenIdx = bestIdx; bestParams = candidateParams;
                        typeStr = 'bandpass'; charCode = 3; break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                elseif bestIdx == 6
                    % Trough-inverted: symmetric strict criteria to band-pass
                    fitY = gaussFun(candidateParams, linspace(minSpeed, maxSpeed, 200));
                    
                    centerInside      = candidateParams(3) > minAllowableCenter && ...
                                        candidateParams(3) < maxAllowableCenter;
                    sigmaAcceptable   = candidateParams(4) <= maxAllowableSigma;
                    rawTroughInterior = minidx > 1 && minidx < numel(ytrainvals);
                    dropsBothSidesFit = (fitY(1)   - min(fitY) > 0.15 * rangeVal) && ...
                                        (fitY(end) - min(fitY) > 0.15 * rangeVal);
                    dropsBothSidesRaw = (ytrainvals(1)   - minVal > 0.10 * rangeVal) && ...
                                        (ytrainvals(end) - minVal > 0.10 * rangeVal);
                    fitProminence     = min(fitY(1), fitY(end)) - min(fitY);
                    if centerInside && sigmaAcceptable && rawTroughInterior && ...
                       dropsBothSidesFit && dropsBothSidesRaw && (fitProminence > 0.15 * rangeVal)
                        chosenIdx = bestIdx; bestParams = candidateParams;
                        typeStr = 'trough_inverted'; charCode = 4; break;
                    else
                        idx_idx = idx_idx + 1; continue;
                    end
                end
            end
            if isnan(chosenIdx)
                typeStr = 'untuned'; charCode = 0;
                bestParams = [nan nan nan nan]; fullR2 = 0;
            else
                SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
                if SS_tot < 1e-9, SS_tot = 1; end
                fullR2 = 1 - (resnorm(chosenIdx) / SS_tot);
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
            cls.R2(r)             = fullR2;
            cls.fitParams(r, :)   = bestParams;
            cls.preferredSpeed(r) = prefSpeedVal;
        end
        response.(targetStruct).(useField).classification = cls;
        
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
    fprintf('Successfully finished Saleem-style execution loop.\n');
end


% % function response = classifySpeedTuningFromCurves(sessionFileInfo, response, useField)
% % % CLASSIFYSPEEDTUNINGFROMCURVES
% % % Classifies speed tuning curves using moving-speed data bins only.
% % % 
% % % Phenotype Classifications:
% % %   0 = Untuned
% % %   1 = Low-Pass / Decay
% % %   2 = High-Pass / Rise
% % %   3 = Band-Pass / Selective Peak
% % %   4 = Trough / Inverted Peak
% % 
% %     if nargin < 3 || isempty(useField)
% %         useField = 'dFFNeuropilCorrected';
% %     end
% % 
% %     stimFileName = sprintf('%s_%s_Response_%s.mat', ...
% %         sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
% %     fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
% %     frameworks = {'tuningCurve', 'tuningCurveFixedBins'};
% % 
% %     for f = 1:numel(frameworks)
% %         targetStruct = frameworks{f};
% %         if ~isfield(response, targetStruct) || ~isfield(response.(targetStruct), useField)
% %             continue;
% %         end
% % 
% %         edges = response.(targetStruct).speedBins;
% %         movingCenters = (edges(1:end-1) + diff(edges)/2)';
% % 
% %         y_session = response.(targetStruct).(useField).moveMean;
% %         numROIs = size(y_session, 1);
% % 
% %         % Preallocate primary classification outputs
% %         cls = struct();
% %         cls.tuningType          = cell(numROIs, 1);
% %         cls.tuningCode          = nan(numROIs, 1);
% %         cls.R2                  = nan(numROIs, 1);
% %         cls.fitParams           = nan(numROIs, 4);
% %         cls.preferredSpeed      = nan(numROIs, 1);
% %         cls.usedMovingOnly      = true;
% %         cls.classifierVersion   = "running_only_streamlined_v1";
% % 
% %         gaussFun = @(params, xdata) ...
% %             params(1) + params(2) .* exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
% %         optsFit = optimset('Display', 'off');
% % 
% %         minSpeed   = min(movingCenters);
% %         maxSpeed   = max(movingCenters);
% %         speedRange = maxSpeed - minSpeed;
% %         customOffset = 0.5;
% % 
% %         % --- HARDCODED CRITERIA BOUNDARIES ---
% %         % Bins must peak away from the raw edges, and widths must stay bounded
% %         minAllowableCenter = minSpeed + (0.15 * speedRange);
% %         maxAllowableCenter = maxSpeed - (0.15 * speedRange);
% %         maxAllowableSigma  = 0.25 * speedRange;
% % 
% %         fprintf('Executing Shape-Constrained Template Fits on %d ROIs [%s]...\n', numROIs, targetStruct);
% % 
% %         for r = 1:numROIs
% %             ytrainvals = y_session(r, :)';
% %             xtrainvals = movingCenters;
% %             validTrain = ~isnan(ytrainvals);
% % 
% %             % Gate 1: Enforce minimum raw data points
% %             if sum(validTrain) < 4
% %                 cls.tuningType{r} = 'untuned';
% %                 cls.tuningCode(r) = 0;
% %                 cls.R2(r) = 0;
% %                 continue;
% %             end
% % 
% %             xtrainvals = xtrainvals(validTrain);
% %             ytrainvals = ytrainvals(validTrain);
% % 
% %             [maxVal, maxidx] = max(ytrainvals);
% %             [minVal, minidx] = min(ytrainvals);
% %             rangeVal = maxVal - minVal;
% %             if rangeVal <= 0, rangeVal = 1e-5; end
% % 
% %             % Extract linear trajectory metrics
% %             speedCorrelation = corr(xtrainvals, ytrainvals, 'rows', 'complete', 'type', 'Pearson');
% %             netRise = ytrainvals(end) - ytrainvals(1);
% %             netDrop = ytrainvals(1) - ytrainvals(end);
% % 
% %             % Pre-build template bound configurations
% %             param_out = nan(6, 4);
% %             resnorm   = inf(6, 1);
% % 
% %             lb_low1  = [-200,   0.1, minSpeed - customOffset*10, 0.7];     ub_low1  = [ 200, 200.0, minSpeed + customOffset,    speedRange*5]; x0_low1  = [0, rangeVal, minSpeed - customOffset*2, speedRange/2];
% %             lb_low2  = [-200, -200.0, maxSpeed - customOffset,   0.7];     ub_low2  = [ 200,  -0.1, maxSpeed + customOffset*10, speedRange*5]; x0_low2  = [0, -rangeVal, maxSpeed + customOffset*2, speedRange/2];
% %             lb_high1 = [-200,   0.1, maxSpeed - customOffset,   0.7];     ub_high1 = [ 200, 200.0, maxSpeed + customOffset*10, speedRange*5]; x0_high1 = [0, rangeVal, maxSpeed + customOffset*2, speedRange/2];
% %             lb_high2 = [-200, -200.0, minSpeed - customOffset*10, 0.7];    ub_high2 = [ 200,  -0.1, minSpeed + customOffset,     speedRange*5]; x0_high2 = [0, -rangeVal, minSpeed - customOffset*2, speedRange/2];
% %             lb_band  = [-200,   0.1, minAllowableCenter, 0.7];             ub_band  = [ 200, 200.0, maxAllowableCenter, maxAllowableSigma];    x0_band  = [0, rangeVal, xtrainvals(maxidx), max(speedRange/6, 0.7)];
% %             lb_inv   = [-200, -200.0, minAllowableCenter, 0.7];             ub_inv   = [ 200,   -0.1, maxAllowableCenter, maxAllowableSigma];    x0_inv   = [0, -rangeVal, xtrainvals(minidx), max(speedRange/6, 0.7)];
% % 
% %             x0_band(3) = min(max(x0_band(3), lb_band(3)), ub_band(3));
% %             x0_inv(3)  = min(max(x0_inv(3),  lb_inv(3)),  ub_inv(3));
% % 
% %             % Compute mathematical curve optimizations
% %             try [param_out(1,:), resnorm(1)] = lsqcurvefit(gaussFun, x0_low1,  xtrainvals, ytrainvals, lb_low1,  ub_low1,  optsFit); catch, end
% %             try [param_out(2,:), resnorm(2)] = lsqcurvefit(gaussFun, x0_low2,  xtrainvals, ytrainvals, lb_low2,  ub_low2,  optsFit); catch, end
% %             try [param_out(3,:), resnorm(3)] = lsqcurvefit(gaussFun, x0_high1, xtrainvals, ytrainvals, lb_high1, ub_high1, optsFit); catch, end
% %             try [param_out(4,:), resnorm(4)] = lsqcurvefit(gaussFun, x0_high2, xtrainvals, ytrainvals, lb_high2, ub_high2, optsFit); catch, end
% %             try [param_out(5,:), resnorm(5)] = lsqcurvefit(gaussFun, x0_band,  xtrainvals, ytrainvals, lb_band,  ub_band,  optsFit); catch, end
% %             try [param_out(6,:), resnorm(6)] = lsqcurvefit(gaussFun, x0_inv,   xtrainvals, ytrainvals, lb_inv,   ub_inv,   optsFit); catch, end
% % 
% %             [~, fitIdx] = sort(resnorm, 'ascend');
% % 
% %             chosenIdx  = nan;
% %             bestParams = nan(1,4);
% %             typeStr    = 'untuned';
% %             charCode   = 0;
% % 
% %             % --- COMPETITIVE TEMPLATE EVALUATION LOOP ---
% %             idx_idx = 1;
% %             while idx_idx <= 6
% %                 bestIdx = fitIdx(idx_idx);
% %                 candidateParams = param_out(bestIdx, :);
% %                 if any(~isfinite(candidateParams)), idx_idx = idx_idx + 1; continue; end
% % 
% %                 if bestIdx == 1 || bestIdx == 2
% %                     % CRITERIA: Low-Pass Validation (Net drop > 12%, Pearson correlation < -0.35, no interior peak)
% %                     isMonotonicDecay = (netDrop > 0.12 * rangeVal) && (~isnan(speedCorrelation) && speedCorrelation < -0.35);
% %                     hasStrongInteriorPeak = maxidx > 1 && maxidx < numel(ytrainvals) && ((maxVal - max(ytrainvals(1), ytrainvals(end))) > 0.15 * rangeVal);
% % 
% %                     if isMonotonicDecay && ~hasStrongInteriorPeak
% %                         chosenIdx = bestIdx; bestParams = candidateParams; typeStr = 'lowpass'; charCode = 1; break;
% %                     else
% %                         idx_idx = idx_idx + 1; continue;
% %                     end
% % 
% %                 elseif bestIdx == 3 || bestIdx == 4
% %                     % CRITERIA: High-Pass Validation (Net rise > 12%, Pearson correlation > 0.35, no interior peak)
% %                     isMonotonicRise = (netRise > 0.12 * rangeVal) && (~isnan(speedCorrelation) && speedCorrelation > 0.35);
% %                     hasStrongInteriorPeak = maxidx > 1 && maxidx < numel(ytrainvals) && ((maxVal - max(ytrainvals(1), ytrainvals(end))) > 0.15 * rangeVal);
% % 
% %                     if isMonotonicRise && ~hasStrongInteriorPeak
% %                         chosenIdx = bestIdx; bestParams = candidateParams; typeStr = 'highpass'; charCode = 2; break;
% %                     else
% %                         idx_idx = idx_idx + 1; continue;
% %                     end
% % 
% %                 elseif bestIdx == 5
% %                     % CRITERIA: Band-Pass Validation
% %                     fitY = gaussFun(candidateParams, linspace(minSpeed, maxSpeed, 200));
% % 
% %                     centerInside    = candidateParams(3) > minAllowableCenter && candidateParams(3) < maxAllowableCenter;
% %                     sigmaAcceptable = candidateParams(4) <= maxAllowableSigma;
% %                     rawPeakInterior = maxidx > 1 && maxidx < numel(ytrainvals);
% % 
% %                     dropsBothSidesFit = (max(fitY) - fitY(1) > 0.15 * rangeVal) && (max(fitY) - fitY(end) > 0.15 * rangeVal);
% %                     dropsBothSidesRaw = (maxVal - ytrainvals(1) > 0.10 * rangeVal) && (maxVal - ytrainvals(end) > 0.10 * rangeVal);
% %                     fitProminence     = max(fitY) - max(fitY(1), fitY(end));
% % 
% %                     if centerInside && sigmaAcceptable && rawPeakInterior && ...
% %                        dropsBothSidesFit && dropsBothSidesRaw && (fitProminence > 0.15 * rangeVal)
% %                         chosenIdx = bestIdx; bestParams = candidateParams; typeStr = 'bandpass'; charCode = 3; break;
% %                     else
% %                         idx_idx = idx_idx + 1; continue;
% %                     end
% % 
% %                 elseif bestIdx == 6
% %                     % CRITERIA: Trough / Inverted Validation
% %                     fitY = gaussFun(candidateParams, linspace(minSpeed, maxSpeed, 200));
% % 
% %                     centerInside      = candidateParams(3) > minAllowableCenter && candidateParams(3) < maxAllowableCenter;
% %                     sigmaAcceptable   = candidateParams(4) <= maxAllowableSigma;
% %                     rawTroughInterior = minidx > 1 && minidx < numel(ytrainvals);
% % 
% %                     risesBothSidesFit = (fitY(1) - min(fitY) > 0.15 * rangeVal) && (fitY(end) - min(fitY) > 0.15 * rangeVal);
% %                     risesBothSidesRaw = (ytrainvals(1) - minVal > 0.10 * rangeVal) && (ytrainvals(end) - minVal > 0.10 * rangeVal);
% %                     fitProminence     = min(fitY(1), fitY(end)) - min(fitY);
% % 
% %                     if centerInside && sigmaAcceptable && rawTroughInterior && ...
% %                        risesBothSidesFit && risesBothSidesRaw && (fitProminence > 0.15 * rangeVal)
% %                         chosenIdx = bestIdx; bestParams = candidateParams; typeStr = 'trough_inverted'; charCode = 4; break;
% %                     else
% %                         idx_idx = idx_idx + 1; continue;
% %                     end
% %                 end
% %                 idx_idx = idx_idx + 1;
% %             end
% % 
% %             % Compute final R2 scores
% %             if isnan(chosenIdx)
% %                 typeStr = 'untuned'; charCode = 0; bestParams = [nan nan nan nan]; fullR2 = 0;
% %             else
% %                 SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
% %                 if SS_tot < 1e-9, SS_tot = 1; end
% %                 fullR2 = 1 - (resnorm(chosenIdx) / SS_tot);
% %             end
% % 
% %             % Target preferred speed mapping coordinates
% %             if charCode == 3 || charCode == 4
% %                 prefSpeedVal = bestParams(3);
% %             elseif charCode == 1 || charCode == 2
% %                 prefSpeedVal = xtrainvals(maxidx);
% %             else
% %                 prefSpeedVal = nan;
% %             end
% % 
% %             cls.tuningType{r}         = typeStr;
% %             cls.tuningCode(r)         = charCode;
% %             cls.R2(r)                 = fullR2;
% %             cls.fitParams(r, :)       = bestParams;
% %             cls.preferredSpeed(r)     = prefSpeedVal;
% %         end
% % 
% %         response.(targetStruct).(useField).classification = cls;
% % 
% %         % Output simple performance metrics to execution terminal
% %         allTypes = cls.tuningType;
% %         fprintf('\n--- Streamlined Running-only Summary [%s] ---\n', targetStruct);
% %         fprintf('  Total Functional ROIs evaluated: %d\n', numROIs);
% %         fprintf('  Group 1 [Low-Pass / Decay]  : %d\n', sum(strcmp(allTypes, 'lowpass')));
% %         fprintf('  Group 2 [High-Pass / Rise]  : %d\n', sum(strcmp(allTypes, 'highpass')));
% %         fprintf('  Group 3 [Band-Pass / Peak]  : %d\n', sum(strcmp(allTypes, 'bandpass')));
% %         fprintf('  Group 4 [Trough / Inverted] : %d\n', sum(strcmp(allTypes, 'trough_inverted')));
% %         fprintf('  Untuned Failsafes           : %d\n', sum(strcmp(allTypes, 'untuned')));
% %         fprintf('----------------------------------------------\n\n');
% %     end
% % 
% %     save(fileFullPath, 'response', '-append');
% % end
% % 
% % 



% function response = classifySpeedTuningFromCurves(sessionFileInfo, response, useField, opts)
% % CLASSIFYSPEEDTUNINGFROMCURVES_RUNNINGONLY : this was clauds version 
% % Classifies speed tuning shape using moving-speed bins only.
% %
% % Shapes:
% %   0 = untuned
% %   1 = lowpass
% %   2 = highpass
% %   3 = bandpass
% %   4 = trough_inverted
% %
% % This function:
% %   - fits competitive Gaussian templates on moveMean only
% %   - does NOT use stationary baseline to decide bandpass/highpass/lowpass
% %   - stores descriptive fit R2 and diagnostics
% %   - uses shape-validation rules to avoid obvious false labels
% %   - does not hard-apply a descriptive R2 threshold unless opts.minFitR2 is set
% %
% % Usage:
% %   response = classifySpeedTuningFromCurves_runningOnly(sessionFileInfo, response)
% %   response = classifySpeedTuningFromCurves_runningOnly(sessionFileInfo, response, 'dFFNeuropilCorrected')
% %   response = classifySpeedTuningFromCurves_runningOnly(sessionFileInfo, response, 'dFFNeuropilCorrected', opts)
% 
%     if nargin < 3 || isempty(useField)
%         useField = 'dFFNeuropilCorrected';
%     end
% 
%     if nargin < 4 || isempty(opts)
%         opts = struct();
%     end
% 
%     % ----------------------------
%     % Defaults
%     % ----------------------------
%     if ~isfield(opts, 'minValidBins'),            opts.minValidBins = 4; end
%     if ~isfield(opts, 'centerMarginFrac'),        opts.centerMarginFrac = 0.15; end
%     if ~isfield(opts, 'maxBandSigmaFrac'),        opts.maxBandSigmaFrac = 0.25; end
%     if ~isfield(opts, 'fitEdgeDropFrac'),         opts.fitEdgeDropFrac = 0.15; end
%     if ~isfield(opts, 'rawEdgeDropFrac'),         opts.rawEdgeDropFrac = 0.10; end
%     if ~isfield(opts, 'prominenceFrac'),          opts.prominenceFrac = 0.15; end
%     if ~isfield(opts, 'monoNetFrac'),             opts.monoNetFrac = 0.12; end
%     if ~isfield(opts, 'monoCorrThresh'),          opts.monoCorrThresh = 0.35; end
%     if ~isfield(opts, 'minDynamicRangeAbs'),      opts.minDynamicRangeAbs = 0; end
%     if ~isfield(opts, 'minDynamicRangeFracMean'), opts.minDynamicRangeFracMean = 0; end
%     if ~isfield(opts, 'minFitR2'),                opts.minFitR2 = []; end
%     if ~isfield(opts, 'saveToDisk'),              opts.saveToDisk = true; end
%     if ~isfield(opts, 'verbose'),                 opts.verbose = true; end
% 
%     stimFileName = sprintf('%s_%s_Response_%s.mat', ...
%         sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
% 
%     fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
% 
%     frameworks = {'tuningCurve', 'tuningCurveFixedBins'};
% 
%     for f = 1:numel(frameworks)
%         targetStruct = frameworks{f};
% 
%         if ~isfield(response, targetStruct)
%             continue;
%         end
%         if ~isfield(response.(targetStruct), useField)
%             continue;
%         end
% 
%         edges = response.(targetStruct).speedBins;
%         movingCenters = (edges(1:end-1) + diff(edges)/2)';
%         numBins = numel(movingCenters);
% 
%         y_session = response.(targetStruct).(useField).moveMean;
%         y_sem     = response.(targetStruct).(useField).moveSEM;
%         y_count   = response.(targetStruct).(useField).moveCount;
% 
%         numROIs = size(y_session, 1);
% 
%         % Preallocate
%         cls = struct();
%         cls.tuningType          = cell(numROIs, 1);
%         cls.tuningCode          = nan(numROIs, 1);
%         cls.R2                  = nan(numROIs, 1);
%         cls.fitParams           = nan(numROIs, 4);
%         cls.preferredSpeed      = nan(numROIs, 1);
%         cls.dynamicRange        = nan(numROIs, 1);
%         cls.meanLevel           = nan(numROIs, 1);
%         cls.rangeOverMean       = nan(numROIs, 1);
%         cls.speedCorr           = nan(numROIs, 1);
%         cls.netRise             = nan(numROIs, 1);
%         cls.netDrop             = nan(numROIs, 1);
%         cls.fitProminence       = nan(numROIs, 1);
%         cls.leftEdgeDeltaFit    = nan(numROIs, 1);
%         cls.rightEdgeDeltaFit   = nan(numROIs, 1);
%         cls.leftEdgeDeltaRaw    = nan(numROIs, 1);
%         cls.rightEdgeDeltaRaw   = nan(numROIs, 1);
%         cls.centerInside        = false(numROIs, 1);
%         cls.sigmaAcceptable     = false(numROIs, 1);
%         cls.rawPeakInterior     = false(numROIs, 1);
%         cls.rawTroughInterior   = false(numROIs, 1);
%         cls.passesShapeCheck    = false(numROIs, 1);
%         cls.passesDynamicRange  = false(numROIs, 1);
%         cls.passesFitR2         = false(numROIs, 1);
%         cls.passesMonotonicity  = false(numROIs, 1);
%         cls.bestTemplateIdx     = nan(numROIs, 1);
%         cls.candidateOrder      = nan(numROIs, 6);
%         cls.usedMovingOnly      = true;
%         cls.classifierVersion   = "running_only_v1";
% 
%         gaussFun = @(params, xdata) ...
%             params(1) + params(2) .* exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
% 
%         optsFit = optimset('Display', 'off');
% 
%         minSpeed   = min(movingCenters);
%         maxSpeed   = max(movingCenters);
%         speedRange = maxSpeed - minSpeed;
%         customOffset = 0.5;
% 
%         minAllowableCenter = minSpeed + (opts.centerMarginFrac * speedRange);
%         maxAllowableCenter = maxSpeed - (opts.centerMarginFrac * speedRange);
%         maxAllowableSigma  = opts.maxBandSigmaFrac * speedRange;
% 
%         if opts.verbose
%             fprintf('Running-only template fits on %d ROIs [%s | %s]...\n', ...
%                 numROIs, targetStruct, useField);
%         end
% 
%         for r = 1:numROIs
%             ytrainvals = y_session(r, :)';
%             xtrainvals = movingCenters;
% 
%             validTrain = ~isnan(ytrainvals);
%             if sum(validTrain) < opts.minValidBins
%                 cls.tuningType{r} = 'untuned';
%                 cls.tuningCode(r) = 0;
%                 cls.R2(r) = 0;
%                 continue;
%             end
% 
%             xtrainvals = xtrainvals(validTrain);
%             ytrainvals = ytrainvals(validTrain);
% 
%             [maxVal, maxidx] = max(ytrainvals);
%             [minVal, minidx] = min(ytrainvals);
% 
%             rangeVal = maxVal - minVal;
%             if rangeVal <= 0
%                 rangeVal = 1e-5;
%             end
% 
%             meanLevel = mean(ytrainvals, 'omitnan');
%             if abs(meanLevel) < 1e-9
%                 rangeOverMean = inf;
%             else
%                 rangeOverMean = rangeVal / abs(meanLevel);
%             end
% 
%             cls.dynamicRange(r)  = rangeVal;
%             cls.meanLevel(r)     = meanLevel;
%             cls.rangeOverMean(r) = rangeOverMean;
% 
%             passesDynamicRange = (rangeVal >= opts.minDynamicRangeAbs) && ...
%                                  (rangeOverMean >= opts.minDynamicRangeFracMean);
%             cls.passesDynamicRange(r) = passesDynamicRange;
% 
%             if numel(unique(xtrainvals)) > 1 && numel(unique(ytrainvals)) > 1
%                 rho = corr(xtrainvals, ytrainvals, 'rows', 'complete', 'type', 'Pearson');
%             else
%                 rho = NaN;
%             end
%             cls.speedCorr(r) = rho;
% 
%             netRise = ytrainvals(end) - ytrainvals(1);
%             netDrop = ytrainvals(1) - ytrainvals(end);
%             cls.netRise(r) = netRise;
%             cls.netDrop(r) = netDrop;
% 
%             param_out = nan(6, 4);
%             resnorm   = inf(6, 1);
% 
%             lb_low1  = [-200,   0.1, minSpeed - customOffset*10, 0.7];
%             ub_low1  = [ 200, 200.0, minSpeed + customOffset,    speedRange*5];
%             x0_low1  = [0, rangeVal, minSpeed - customOffset*2, speedRange/2];
% 
%             lb_low2  = [-200, -200.0, maxSpeed - customOffset,   0.7];
%             ub_low2  = [ 200,  -0.1, maxSpeed + customOffset*10, speedRange*5];
%             x0_low2  = [0, -rangeVal, maxSpeed + customOffset*2, speedRange/2];
% 
%             lb_high1 = [-200,   0.1, maxSpeed - customOffset,   0.7];
%             ub_high1 = [ 200, 200.0, maxSpeed + customOffset*10, speedRange*5];
%             x0_high1 = [0, rangeVal, maxSpeed + customOffset*2, speedRange/2];
% 
%             lb_high2 = [-200, -200.0, minSpeed - customOffset*10, 0.7];
%             ub_high2 = [ 200,  -0.1, minSpeed + customOffset,     speedRange*5];
%             x0_high2 = [0, -rangeVal, minSpeed - customOffset*2, speedRange/2];
% 
%             lb_band  = [-200,   0.1, minAllowableCenter, 0.7];
%             ub_band  = [ 200, 200.0, maxAllowableCenter, maxAllowableSigma];
%             x0_band  = [0, rangeVal, xtrainvals(maxidx), max(speedRange/6, 0.7)];
% 
%             lb_inv   = [-200, -200.0, minAllowableCenter, 0.7];
%             ub_inv   = [ 200,   -0.1, maxAllowableCenter, maxAllowableSigma];
%             x0_inv   = [0, -rangeVal, xtrainvals(minidx), max(speedRange/6, 0.7)];
% 
%             x0_band(3) = min(max(x0_band(3), lb_band(3)), ub_band(3));
%             x0_inv(3)  = min(max(x0_inv(3),  lb_inv(3)),  ub_inv(3));
% 
%             try
%                 [param_out(1,:), resnorm(1)] = lsqcurvefit(gaussFun, x0_low1,  xtrainvals, ytrainvals, lb_low1,  ub_low1,  optsFit);
%             catch
%             end
%             try
%                 [param_out(2,:), resnorm(2)] = lsqcurvefit(gaussFun, x0_low2,  xtrainvals, ytrainvals, lb_low2,  ub_low2,  optsFit);
%             catch
%             end
%             try
%                 [param_out(3,:), resnorm(3)] = lsqcurvefit(gaussFun, x0_high1, xtrainvals, ytrainvals, lb_high1, ub_high1, optsFit);
%             catch
%             end
%             try
%                 [param_out(4,:), resnorm(4)] = lsqcurvefit(gaussFun, x0_high2, xtrainvals, ytrainvals, lb_high2, ub_high2, optsFit);
%             catch
%             end
%             try
%                 [param_out(5,:), resnorm(5)] = lsqcurvefit(gaussFun, x0_band,  xtrainvals, ytrainvals, lb_band,  ub_band,  optsFit);
%             catch
%             end
%             try
%                 [param_out(6,:), resnorm(6)] = lsqcurvefit(gaussFun, x0_inv,   xtrainvals, ytrainvals, lb_inv,   ub_inv,   optsFit);
%             catch
%             end
% 
%             [~, fitIdx] = sort(resnorm, 'ascend');
%             cls.candidateOrder(r, :) = fitIdx(:)';
% 
%             chosenIdx  = nan;
%             bestParams = nan(1,4);
%             typeStr    = 'untuned';
%             charCode   = 0;
%             passesShapeCheck = false;
%             passesMonotonicity = false;
% 
%             idx_idx = 1;
%             while idx_idx <= 6
%                 bestIdx = fitIdx(idx_idx);
%                 candidateParams = param_out(bestIdx, :);
% 
%                 if any(~isfinite(candidateParams))
%                     idx_idx = idx_idx + 1;
%                     continue;
%                 end
% 
%                 SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
%                 if SS_tot < 1e-9
%                     SS_tot = 1;
%                 end
%                 currentR2 = 1 - (resnorm(bestIdx) / SS_tot);
% 
%                 if bestIdx == 1 || bestIdx == 2
%                     monoOK = (netDrop > opts.monoNetFrac * rangeVal) && ...
%                              (~isnan(rho) && rho < -opts.monoCorrThresh);
% 
%                     % reject if obvious interior peak stronger than monotonic explanation
%                     [~, rawPeakIdx] = max(ytrainvals);
%                     hasStrongInteriorPeak = rawPeakIdx > 1 && rawPeakIdx < numel(ytrainvals) && ...
%                                             ((maxVal - max(ytrainvals(1), ytrainvals(end))) > opts.prominenceFrac * rangeVal);
% 
%                     if monoOK && ~hasStrongInteriorPeak
%                         chosenIdx = bestIdx;
%                         bestParams = candidateParams;
%                         typeStr = 'lowpass';
%                         charCode = 1;
%                         passesShapeCheck = true;
%                         passesMonotonicity = true;
%                         break;
%                     else
%                         idx_idx = idx_idx + 1;
%                         continue;
%                     end
% 
%                 elseif bestIdx == 3 || bestIdx == 4
%                     monoOK = (netRise > opts.monoNetFrac * rangeVal) && ...
%                              (~isnan(rho) && rho > opts.monoCorrThresh);
% 
%                     [~, rawPeakIdx] = max(ytrainvals);
%                     hasStrongInteriorPeak = rawPeakIdx > 1 && rawPeakIdx < numel(ytrainvals) && ...
%                                             ((maxVal - max(ytrainvals(1), ytrainvals(end))) > opts.prominenceFrac * rangeVal);
% 
%                     if monoOK && ~hasStrongInteriorPeak
%                         chosenIdx = bestIdx;
%                         bestParams = candidateParams;
%                         typeStr = 'highpass';
%                         charCode = 2;
%                         passesShapeCheck = true;
%                         passesMonotonicity = true;
%                         break;
%                     else
%                         idx_idx = idx_idx + 1;
%                         continue;
%                     end
% 
%                 elseif bestIdx == 5
%                     peakCenter = candidateParams(3);
%                     sigmaVal   = candidateParams(4);
% 
%                     fineX = linspace(minSpeed, maxSpeed, 200);
%                     fitY  = gaussFun(candidateParams, fineX);
% 
%                     peakValFit   = max(fitY);
%                     leftEdgeFit  = fitY(1);
%                     rightEdgeFit = fitY(end);
% 
%                     rawPeakVal  = max(ytrainvals);
%                     rawLeftVal  = ytrainvals(1);
%                     rawRightVal = ytrainvals(end);
%                     [~, rawPeakIdx] = max(ytrainvals);
% 
%                     centerInside = peakCenter > minAllowableCenter && peakCenter < maxAllowableCenter;
%                     sigmaAcceptable = sigmaVal <= maxAllowableSigma;
% 
%                     leftDropFit  = peakValFit - leftEdgeFit;
%                     rightDropFit = peakValFit - rightEdgeFit;
%                     leftDropRaw  = rawPeakVal - rawLeftVal;
%                     rightDropRaw = rawPeakVal - rawRightVal;
%                     fitProminence = peakValFit - max(leftEdgeFit, rightEdgeFit);
% 
%                     rawPeakInterior = rawPeakIdx > 1 && rawPeakIdx < numel(ytrainvals);
% 
%                     cls.leftEdgeDeltaFit(r)  = leftDropFit;
%                     cls.rightEdgeDeltaFit(r) = rightDropFit;
%                     cls.leftEdgeDeltaRaw(r)  = leftDropRaw;
%                     cls.rightEdgeDeltaRaw(r) = rightDropRaw;
%                     cls.fitProminence(r)     = fitProminence;
%                     cls.centerInside(r)      = centerInside;
%                     cls.sigmaAcceptable(r)   = sigmaAcceptable;
%                     cls.rawPeakInterior(r)   = rawPeakInterior;
% 
%                     dropsBothSidesFit = (leftDropFit  > opts.fitEdgeDropFrac * rangeVal) && ...
%                                         (rightDropFit > opts.fitEdgeDropFrac * rangeVal);
% 
%                     dropsBothSidesRaw = (leftDropRaw  > opts.rawEdgeDropFrac * rangeVal) && ...
%                                         (rightDropRaw > opts.rawEdgeDropFrac * rangeVal);
% 
%                     if centerInside && sigmaAcceptable && rawPeakInterior && ...
%                        dropsBothSidesFit && dropsBothSidesRaw && ...
%                        (fitProminence > opts.prominenceFrac * rangeVal)
% 
%                         chosenIdx = bestIdx;
%                         bestParams = candidateParams;
%                         typeStr = 'bandpass';
%                         charCode = 3;
%                         passesShapeCheck = true;
%                         break;
%                     else
%                         idx_idx = idx_idx + 1;
%                         continue;
%                     end
% 
%                 elseif bestIdx == 6
%                     troughCenter = candidateParams(3);
%                     sigmaVal     = candidateParams(4);
% 
%                     fineX = linspace(minSpeed, maxSpeed, 200);
%                     fitY  = gaussFun(candidateParams, fineX);
% 
%                     troughValFit  = min(fitY);
%                     leftEdgeFit   = fitY(1);
%                     rightEdgeFit  = fitY(end);
% 
%                     rawTroughVal = min(ytrainvals);
%                     rawLeftVal   = ytrainvals(1);
%                     rawRightVal  = ytrainvals(end);
%                     [~, rawTroughIdx] = min(ytrainvals);
% 
%                     centerInside = troughCenter > minAllowableCenter && troughCenter < maxAllowableCenter;
%                     sigmaAcceptable = sigmaVal <= maxAllowableSigma;
% 
%                     leftRiseFit  = leftEdgeFit - troughValFit;
%                     rightRiseFit = rightEdgeFit - troughValFit;
%                     leftRiseRaw  = rawLeftVal - rawTroughVal;
%                     rightRiseRaw = rawRightVal - rawTroughVal;
%                     fitProminence = min(leftEdgeFit, rightEdgeFit) - troughValFit;
% 
%                     rawTroughInterior = rawTroughIdx > 1 && rawTroughIdx < numel(ytrainvals);
% 
%                     cls.leftEdgeDeltaFit(r)  = leftRiseFit;
%                     cls.rightEdgeDeltaFit(r) = rightRiseFit;
%                     cls.leftEdgeDeltaRaw(r)  = leftRiseRaw;
%                     cls.rightEdgeDeltaRaw(r) = rightRiseRaw;
%                     cls.fitProminence(r)     = fitProminence;
%                     cls.centerInside(r)      = centerInside;
%                     cls.sigmaAcceptable(r)   = sigmaAcceptable;
%                     cls.rawTroughInterior(r) = rawTroughInterior;
% 
%                     risesBothSidesFit = (leftRiseFit  > opts.fitEdgeDropFrac * rangeVal) && ...
%                                         (rightRiseFit > opts.fitEdgeDropFrac * rangeVal);
% 
%                     risesBothSidesRaw = (leftRiseRaw  > opts.rawEdgeDropFrac * rangeVal) && ...
%                                         (rightRiseRaw > opts.rawEdgeDropFrac * rangeVal);
% 
%                     if centerInside && sigmaAcceptable && rawTroughInterior && ...
%                        risesBothSidesFit && risesBothSidesRaw && ...
%                        (fitProminence > opts.prominenceFrac * rangeVal)
% 
%                         chosenIdx = bestIdx;
%                         bestParams = candidateParams;
%                         typeStr = 'trough_inverted';
%                         charCode = 4;
%                         passesShapeCheck = true;
%                         break;
%                     else
%                         idx_idx = idx_idx + 1;
%                         continue;
%                     end
%                 end
% 
%                 idx_idx = idx_idx + 1;
%             end
% 
%             if isnan(chosenIdx)
%                 typeStr = 'untuned';
%                 charCode = 0;
%                 bestParams = [nan nan nan nan];
%                 fullR2 = 0;
%             else
%                 SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
%                 if SS_tot < 1e-9
%                     SS_tot = 1;
%                 end
%                 fullR2 = 1 - (resnorm(chosenIdx) / SS_tot);
%             end
% 
%             passesFitR2 = true;
%             if ~isempty(opts.minFitR2)
%                 passesFitR2 = fullR2 >= opts.minFitR2;
%             end
% 
%             if ~passesDynamicRange || ~passesFitR2
%                 typeStr = 'untuned';
%                 charCode = 0;
%             end
% 
%             if charCode == 3 || charCode == 4
%                 prefSpeedVal = bestParams(3);
%             elseif charCode == 1 || charCode == 2
%                 [~, peakLocationIdx] = max(ytrainvals);
%                 prefSpeedVal = xtrainvals(peakLocationIdx);
%             else
%                 prefSpeedVal = nan;
%             end
% 
%             cls.tuningType{r}         = typeStr;
%             cls.tuningCode(r)         = charCode;
%             cls.R2(r)                 = fullR2;
%             cls.fitParams(r, :)       = bestParams;
%             cls.preferredSpeed(r)     = prefSpeedVal;
%             cls.passesShapeCheck(r)   = passesShapeCheck;
%             cls.passesFitR2(r)        = passesFitR2;
%             cls.passesMonotonicity(r) = passesMonotonicity;
%             if ~isnan(chosenIdx)
%                 cls.bestTemplateIdx(r) = chosenIdx;
%             end
%         end
% 
%         response.(targetStruct).(useField).classification = cls;
% 
%         if opts.verbose
%             allTypes = cls.tuningType;
%             fprintf('\n--- Running-only Speed Summary [%s | %s] ---\n', targetStruct, useField);
%             fprintf('  Total ROIs evaluated: %d\n', numROIs);
%             fprintf('  Low-Pass          : %d\n', sum(strcmp(allTypes, 'lowpass')));
%             fprintf('  High-Pass         : %d\n', sum(strcmp(allTypes, 'highpass')));
%             fprintf('  Band-Pass         : %d\n', sum(strcmp(allTypes, 'bandpass')));
%             fprintf('  Trough-Inverted   : %d\n', sum(strcmp(allTypes, 'trough_inverted')));
%             fprintf('  Untuned           : %d\n', sum(strcmp(allTypes, 'untuned')));
%             fprintf('----------------------------------------------\n\n');
%         end
%     end
% 
%     if opts.saveToDisk
%         save(fileFullPath, 'response', '-append');
%         if opts.verbose
%             fprintf('Running-only classifications saved to disk.\n');
%         end
%     end
% end


% function response = classifySpeedTuningFromCurves(sessionFileInfo, response, useField)
% % Parameterizes and classifies neurons into 4 main tuning classes using
% % a Multi-Template Competitive Gaussian Fit framework on 100% of session data.
% %   1. Reads your exact master speed bins and stable session mean curves.
% %   2. Fits 4 symmetric Gaussian templates natively on the physical speed axis (cm/s).
% %   3. Evaluates standard full-session R2.
% %   4. Downstream noise filtering is handled by your 1,000-permutation shuffles.
% 
%     if nargin < 3, useField = 'dFFNeuropilCorrected'; end
% 
%     % 
%     stimFileName = sprintf('%s_%s_Response_%s.mat', ...
%         sessionFileInfo.animal_name, sessionFileInfo.session_name, response.stimName);
%     fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
% 
%     if ~exist(fileFullPath, 'file')
%         error('Internal path recovery failed. File does not exist: %s', fileFullPath);
%     end
% 
%     frameworks = {'tuningCurve', 'tuningCurveFixedBins'};
% 
%     for f = 1:2
%         targetStruct = frameworks{f};
%         if ~isfield(response, targetStruct) || ~isfield(response.(targetStruct), useField), continue; end
% 
%         % Recover master session-wide velocity axis boundaries
%         edges = response.(targetStruct).speedBins;
%         movingCenters = (edges(1:end-1) + diff(edges)/2)'; 
%         numBins = length(movingCenters);
% 
%         % Pull clean, stable full-session mean response curves
%         y_session = response.(targetStruct).(useField).moveMean;
%         numROIs = size(y_session, 1);
% 
%         % Pre-allocate classification subfields
%         response.(targetStruct).(useField).classification.tuningType     = cell(numROIs, 1);
%         response.(targetStruct).(useField).classification.tuningCode     = nan(numROIs, 1);
%         response.(targetStruct).(useField).classification.R2             = nan(numROIs, 1); 
%         response.(targetStruct).(useField).classification.fitParams      = nan(numROIs, 4);
%         response.(targetStruct).(useField).classification.preferredSpeed = nan(numROIs, 1);
% 
%         gaussFun = @(params, xdata) params(1) + params(2) .* exp(-(((xdata - params(3)).^2) / (2 * (params(4).^2))));
%         opts = optimset('Display', 'off');
% 
%         minSpeed = min(movingCenters); maxSpeed = max(movingCenters); speedRange = maxSpeed - minSpeed; customOffset = 0.5; 
% 
%         fprintf('Executing Full-Session Template Fits on %d ROIs [%s Framework]...\n', numROIs, targetStruct);
% 
%         for r = 1:numROIs
%             % Extract the full-session data vector for this specific ROI
%             ytrainvals = y_session(r, :)';
%             xtrainvals = movingCenters;
% 
%             validTrain = ~isnan(ytrainvals);
%             if sum(validTrain) < 4
%                 response.(targetStruct).(useField).classification.tuningType{r} = 'untuned';
%                 response.(targetStruct).(useField).classification.tuningCode(r) = 0;
%                 response.(targetStruct).(useField).classification.R2(r)         = 0;
%                 continue;
%             end
% 
%             xtrainvals = xtrainvals(validTrain);
%             ytrainvals = ytrainvals(validTrain);
% 
%             [maxVal, maxidx] = max(ytrainvals); [minVal, minidx] = min(ytrainvals);
%             rangeVal = maxVal - minVal; if rangeVal <= 0, rangeVal = 1e-5; end
% 
%             param_out = zeros(6, 4); resnorm = zeros(6, 1);
% 
%             % Boundary boxes constraint matrices (Symmetric fits bounded in true physical cm/s)
%             lb_low1 = [-200, 0.1, minSpeed - customOffset*10, 0.7]; ub_low1 = [200, 200, minSpeed + customOffset, speedRange*5]; x0_low1 = [0, rangeVal, minSpeed - customOffset*2, speedRange/2];
%             lb_low2 = [-200, -200, maxSpeed - customOffset, 0.7]; ub_low2 = [200, -0.1, maxSpeed + customOffset*10, speedRange*5]; x0_low2 = [0, -rangeVal, maxSpeed + customOffset*2, speedRange/2];
%             lb_high1 = [-200, 0.1, maxSpeed - customOffset, 0.7]; ub_high1 = [200, 200, maxSpeed + customOffset*10, speedRange*5]; x0_high1 = [0, rangeVal, maxSpeed + customOffset*2, speedRange/2];
%             lb_high2 = [-200, -200, minSpeed - customOffset*10, 0.7]; ub_high2 = [200, -0.1, minSpeed + customOffset, speedRange*5]; x0_high2 = [0, -rangeVal, minSpeed - customOffset*2, speedRange/2];
%             lb_band = [-200, 0.1, minSpeed + customOffset, 0.7]; ub_band = [200, 200, maxSpeed - customOffset, speedRange]; x0_band = [0, rangeVal, xtrainvals(maxidx), speedRange/4];
%             lb_inv = [-200, -200, minSpeed + customOffset, 0.7]; ub_inv = [200, -0.1, maxSpeed - customOffset, speedRange]; x0_inv = [0, -rangeVal, xtrainvals(minidx), speedRange/4];
% 
%             % Run non-linear curve fitting natively on the physical axis coordinates
%             [param_out(1,:), resnorm(1)] = lsqcurvefit(gaussFun, x0_low1,  xtrainvals, ytrainvals, lb_low1,  ub_low1,  opts);
%             [param_out(2,:), resnorm(2)] = lsqcurvefit(gaussFun, x0_low2,  xtrainvals, ytrainvals, lb_low2,  ub_low2,  opts);
%             [param_out(3,:), resnorm(3)] = lsqcurvefit(gaussFun, x0_high1, xtrainvals, ytrainvals, lb_high1, ub_high1, opts);
%             [param_out(4,:), resnorm(4)] = lsqcurvefit(gaussFun, x0_high2, xtrainvals, ytrainvals, lb_high2, ub_high2, opts);
%             [param_out(5,:), resnorm(5)] = lsqcurvefit(gaussFun, x0_band,  xtrainvals, ytrainvals, lb_band,  ub_band,  opts);
%             [param_out(6,:), resnorm(6)] = lsqcurvefit(gaussFun, x0_inv,   xtrainvals, ytrainvals, lb_inv,   ub_inv,   opts);
% 
% 
% 
%             [~, fitIdx] = sort(resnorm);
%             idx_idx = 1; done = false;
% 
%             while ~done && idx_idx <= 6
%                 bestIdx = fitIdx(idx_idx); bestParams = param_out(bestIdx, :);
%                 if bestIdx == 1 || bestIdx == 2, char = 1; typeStr = 'lowpass';
%                 elseif bestIdx == 3 || bestIdx == 4, char = 2; typeStr = 'highpass';
%                 elseif bestIdx == 5
%                     char = 3; typeStr = 'bandpass';
%                     fineX = linspace(minSpeed, maxSpeed, 100); evalOfFun = feval(gaussFun, bestParams, fineX); [~, ~, ~, prom] = findpeaks(evalOfFun);
%                     if isempty(prom) || prom < (rangeVal / 3), idx_idx = idx_idx + 1; continue; end
%                 elseif bestIdx == 6
%                     char = 4; typeStr = 'trough_inverted';
%                     fineX = linspace(minSpeed, maxSpeed, 100); evalOfFun = -1 .* feval(gaussFun, bestParams, fineX); [~, ~, ~, prom] = findpeaks(evalOfFun);
%                     if isempty(prom) || prom < (rangeVal / 3), idx_idx = idx_idx + 1; continue; end
%                 end
%                 done = true;
%             end
% 
%             % Standard R2 Evaluation over 100% data points
%             SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
%             if SS_tot < 1e-9, SS_tot = 1; end
%             fullR2 = 1 - (resnorm(bestIdx) / SS_tot);
% 
%             if char == 4, [~, peakLocationIdx] = min(ytrainvals);
%             else, [~, peakLocationIdx] = max(ytrainvals); end
%             prefSpeedVal = xtrainvals(peakLocationIdx);
% 
%             %
%             response.(targetStruct).(useField).classification.tuningType{r}      = typeStr;
%             response.(targetStruct).(useField).classification.tuningCode(r)      = char;
%             response.(targetStruct).(useField).classification.R2(r)              = fullR2;
%             response.(targetStruct).(useField).classification.fitParams(r, :)    = bestParams;
%             response.(targetStruct).(useField).classification.preferredSpeed(r)  = prefSpeedVal;
%         end
% 
%         %
%         allTypes = response.(targetStruct).(useField).classification.tuningType;
%         fprintf('\n--- Multi-Template Speed Summary [%s] ---\n', targetStruct);
%         fprintf('  Total Functional ROIs evaluated: %d\n', numROIs);
%         fprintf('  ---------------------------------------------------\n');
%         fprintf('  Template Group 1 [Low-Pass / Decay] : %d\n', sum(strcmp(allTypes, 'lowpass')));
%         fprintf('  Template Group 2 [High-Pass / Rise] : %d\n', sum(strcmp(allTypes, 'highpass')));
%         fprintf('  Template Group 3 [Band-Pass / Peak] : %d\n', sum(strcmp(allTypes, 'bandpass')));
%         fprintf('  Template Group 4 [Trough / Inverted]: %d\n', sum(strcmp(allTypes, 'trough_inverted')));
%         fprintf('  Untuned Failsafes (Insufficient Data): %d\n', sum(strcmp(allTypes, 'untuned')));
%         fprintf('-----------------------------------------------------\n\n');
%     end
% 
%     save(fileFullPath, 'response', '-append');
%     fprintf('Full-session classifications successfully saved to disk.\n');
% end