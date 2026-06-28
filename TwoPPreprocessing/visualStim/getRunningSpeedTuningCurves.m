function [response, sessionFileInfo] = getRunningSpeedTuningCurves(sessionFileInfo, stimName, useZScoredSignals, shuffle)
% Computes 1D speed tuning curves using BOTH dynamic quantiles (Saleem et al style)
% and a standardized fixed grid for cross-stim type comparisons (gray-darkness).
%
% SHUFFLE SIGNIFICANCE CRITERION:
%   Matches Saleem et al. (2013): Real tuning curve variance must exceed 99.9% 
%   of circular-shifted null variances (p <= 0.001) across 1,000 permutations.

    if nargin < 3 || isempty(useZScoredSignals), useZScoredSignals = false; end
    if nargin < 4 || isempty(shuffle), shuffle = true; end 
    
    stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
    if isempty(stimIdx), error('Stimulus %s not found.', stimName); end
    
    data2P = load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData);
    dataPeriph = load(sessionFileInfo.stimFiles(stimIdx).processedPeripheralData, 'peripheralData');
    
    sourceStructName = 'processedSignals';
    if useZScoredSignals, sourceStructName = 'zScoredProcessedSignals'; end
    signals = data2P.(sourceStructName);
    if isfield(data2P, 'spks') && ~isempty(data2P.spks)
        signals.spks = data2P.spks; 
    end 
    
    fs = 60; 
    tickToCm = 3.1415 * 20 / 1024; 
    wheelSpeed = [0; diff(dataPeriph.peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(dataPeriph.peripheralData.Wheel.sampleTimes)];
    wheelSpeed(abs(wheelSpeed) > 150) = NaN; 
    
    mIdx = find(wheelSpeed > 1 & ~isnan(wheelSpeed)); 
    sIdx = find(wheelSpeed <= 1 & ~isnan(wheelSpeed)); 
    
    % 
    allFields = fieldnames(signals);
    firstFieldName = allFields{1};
    numROIs = size(signals.(firstFieldName), 1);
    
    %  Define quantile boundaries (Saleem et al. 2013: 7% frames/bin) ---
    ptsPerBin = floor(0.07 * length(wheelSpeed)); 
    nB = floor(length(mIdx) / ptsPerBin); 
    edgesQuantile = quantile(wheelSpeed(mIdx), linspace(0, 1, max(1, nB) + 1));
    edgesQuantile = unique(edgesQuantile);
    
    % ddefine fixed boundaries for gray darkness comparison (For gray-darkness comparisons) ---
    edgesFixed = [1.5, 3.5, 7.0, 14.0, 25.0, 45.0];
    
    
    response.tuningCurve.source = sourceStructName;
    response.tuningCurve.occupancy.stationary = length(sIdx) / fs;
    response.tuningCurve.speedBins = edgesQuantile;
    
    response.tuningCurveFixedBins.source = sourceStructName;
    response.tuningCurveFixedBins.occupancy.stationary = length(sIdx) / fs;
    response.tuningCurveFixedBins.speedBins = edgesFixed;
    
    % Configure Time-Shift Permutation Offsets 
    numShifts = 1000;  
    minShift = 600; % 10-second 
    maxShift = length(wheelSpeed) - minShift; 
    rng(1); 
    shiftValues = randi([minShift, maxShift], [1, numShifts]);
    
    for thisField = 1:length(allFields)
        fname = allFields{thisField};
        sigData = signals.(fname);
        if size(sigData, 2) ~= length(wheelSpeed), continue; end
        
        statMean = mean(sigData(:, sIdx), 2, 'omitnan');
        statSEM  = std(sigData(:, sIdx), 0, 2, 'omitnan') / sqrt(length(sIdx));
        statCount = repmat(length(sIdx), numROIs, 1);
        runFrameworks = {'tuningCurve', 'tuningCurveFixedBins'};
        edgeDefinitions = {edgesQuantile, edgesFixed};
        
        % Separately compute the complete metric per framework
        for f = 1:2
            targetStruct = runFrameworks{f};
            edges = edgeDefinitions{f};
            numBins = length(edges) - 1;
            
            mMean  = zeros(numROIs, numBins);
            mSEM   = zeros(numROIs, numBins);
            mCount = zeros(numROIs, numBins);
            
            for b = 1:numBins
                idx = find(wheelSpeed >= edges(b) & wheelSpeed < edges(b+1));
                
                mMean(:, b)  = mean(sigData(:, idx), 2, 'omitnan');
                mSEM(:, b)   = std(sigData(:, idx), 0, 2, 'omitnan') / sqrt(length(idx));
                mCount(:, b) = length(idx);
                if thisField == 1
                    response.(targetStruct).occupancy.moving(b) = length(idx) / fs;
                end
            end
            
            moveMax = max(mMean, [], 2); moveMin = min(mMean, [], 2);
            modRatio_move = 2 * (moveMax - moveMin) ./ (moveMax + moveMin);
            
            fullCurve = [statMean, mMean];
            fullMax = max(fullCurve, [], 2); fullMin = min(fullCurve, [], 2);
            modRatio_full = 2 * (fullMax - fullMin) ./ (fullMax + fullMin);
            
            % Measure true curve variance 
            realVarFull = var(fullCurve, 0, 2, 'omitnan');
            realVarMovingOnly = var(mMean, 0, 2, 'omitnan');
            
            pValFull = ones(numROIs, 1);
            pValMoving = ones(numROIs, 1);
            
            %  shuffle
            if shuffle
                fprintf('Shuffling %s [%s Framework]...\n', fname, targetStruct);
                shuffVarsFull = zeros(numROIs, numShifts);
                shuffVarsMoving = zeros(numROIs, numShifts);
                
                % entire session shuffled together for full curves 
                for iperm = 1:numShifts
                    sigS_full = circshift(sigData, shiftValues(iperm), 2);
                    
                    sM_shuff = mean(sigS_full(:, sIdx), 2, 'omitnan');
                    mM_shuff_full = zeros(numROIs, numBins);
                    for b = 1:numBins
                        idxB = wheelSpeed >= edges(b) & wheelSpeed < edges(b+1);
                        mM_shuff_full(:, b) = mean(sigS_full(:, idxB), 2, 'omitnan');
                    end
                    shuffVarsFull(:, iperm) = var([sM_shuff, mM_shuff_full], 0, 2, 'omitnan');
                end
                
                % moving only shuffled together 
                sigDataMovingOnly = sigData(:, mIdx);
                wheelSpeedMovingOnly = wheelSpeed(mIdx);
                
                % Set bounds safe for the size of the moving-only matrix
                maxShiftMove = length(mIdx) - minShift;
                if maxShiftMove > minShift
                    rng(1); shiftValuesMove = randi([minShift, maxShiftMove], [1, numShifts]);
                else
                    shiftValuesMove = shiftValues; % Fallback if running data is short
                end
                
                for iperm = 1:numShifts
                    % Shift ONLY within the active running timeseries block
                    sigS_move = circshift(sigDataMovingOnly, shiftValuesMove(iperm), 2);
                    
                    mM_shuff_move = zeros(numROIs, numBins);
                    for b = 1:numBins
                        % Align the shifted running traces with the running speed vector
                        idxB = wheelSpeedMovingOnly >= edges(b) & wheelSpeedMovingOnly < edges(b+1);
                        mM_shuff_move(:, b) = mean(sigS_move(:, idxB), 2, 'omitnan');
                    end
                    shuffVarsMoving(:, iperm) = var(mM_shuff_move, 0, 2, 'omitnan');
                end
                
                % --- COMPUTE SIGNIFICANCE ---
                pValFull = sum(shuffVarsFull >= realVarFull, 2) / numShifts;
                pValFull(realVarFull == 0) = 1.0; 
                
                pValMoving = sum(shuffVarsMoving >= realVarMovingOnly, 2) / numShifts;
                pValMoving(realVarMovingOnly == 0) = 1.0; 
            end
            
            % Map data 
            response.(targetStruct).(fname).statMean = statMean;
            response.(targetStruct).(fname).statSEM  = statSEM;
            response.(targetStruct).(fname).statCount = statCount;
            response.(targetStruct).(fname).moveMean = mMean;
            response.(targetStruct).(fname).moveSEM  = mSEM;
            response.(targetStruct).(fname).moveCount = mCount;
            
            response.(targetStruct).(fname).modulationRatio_move = modRatio_move;
            response.(targetStruct).(fname).percentModulated_move = modRatio_move * 100;
            response.(targetStruct).(fname).modulationRatio_full = modRatio_full;
            response.(targetStruct).(fname).percentModulated_full = modRatio_full * 100;
            
            response.(targetStruct).(fname).tuningVarFull = realVarFull;
            response.(targetStruct).(fname).tuningVarMoving = realVarMovingOnly;
            response.(targetStruct).(fname).pValFull = pValFull;
            response.(targetStruct).(fname).pValMoving = pValMoving;
            response.(targetStruct).(fname).isSignificant_999 = pValFull <= 0.001;
            response.(targetStruct).(fname).isSignificantMoving_999 = pValMoving <= 0.001;
            response.stimName = stimName; 
        end
    end
    
   % save response structure
    stimFileName = sprintf('%s_%s_Response_%s.mat', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, stimName);
    savePath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
    
    % update file path in sessionsfileinfo
    sessionFileInfo.stimFiles(stimIdx).Response = savePath; 
    
    % save sfi
    save(savePath, 'response', '-v7.3');
    save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end
