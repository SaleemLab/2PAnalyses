function [response] = getRunningSpeedTuningCurves_FixedBins(sessionFileInfo, stimName, useZScoredSignals, shuffle, numBins)
% From Saleem et al, 2013 Significance: Real Variance > 99.9% of Shuffled
% Variance. Fig 2 Darkness 
% Can you this for plotting: plotRunningTuningCurves(respGray, respDark, pdfPath)
if nargin < 3 || isempty(useZScoredSignals), useZScoredSignals = false; end
if nargin < 4 || isempty(shuffle), shuffle = true; end 
if nargin < 5 || isempty(numBins), numBins = 10; end 

%% Load Data
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

%% Speed & Binning
fs = 60; 
tickToCm = 3.1415 * 20 / 1024; 
wheelSpeed = [0; diff(dataPeriph.peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(dataPeriph.peripheralData.Wheel.sampleTimes)];
wheelSpeed(abs(wheelSpeed) > 150) = NaN; 

mIdx = find(wheelSpeed > 1 & ~isnan(wheelSpeed)); % running
sIdx = find(wheelSpeed <= 1 & ~isnan(wheelSpeed)); % stationary

% Create speed bin edges using quantiles for exactly numBins
edges = quantile(wheelSpeed(mIdx), linspace(0, 1, numBins + 1)); 
edges = unique(edges); 
actualNumBins = length(edges)-1;
speedAxis = [0, edges(1:end-1) + diff(edges)/2];

% Save Binning and Occupancy Information
response.tuningCurve.speedBins = edges;
response.tuningCurve.source = sourceStructName;
response.tuningCurve.occupancy.stationary = length(sIdx) / fs;
samplesPerMoveBin = length(mIdx) / actualNumBins; 
response.tuningCurve.occupancy.moving = repmat(samplesPerMoveBin / fs, 1, actualNumBins);

%% Shuffle Parameters
numShifts = 1000;  
minShift = 600; 
maxShift = length(wheelSpeed) - minShift; 
rng(1); 
shiftValues = randi([minShift, maxShift], [1, numShifts]);

%% Processing
allFields = fieldnames(signals);
for thisField = 1:length(allFields)
    fname = allFields{thisField};
    sigData = signals.(fname);
    if size(sigData, 2) ~= length(wheelSpeed), continue; end
    
    numROIs = size(sigData, 1);
    mMean = zeros(numROIs, actualNumBins);
    mSEM  = zeros(numROIs, actualNumBins);
    mCount = zeros(numROIs, actualNumBins);
    
    % Stationary Stats (SEM)
    statMean = mean(sigData(:, sIdx), 2, 'omitnan');
    statSEM  = std(sigData(:, sIdx), 0, 2, 'omitnan') / sqrt(length(sIdx));
    statCount = repmat(length(sIdx), numROIs, 1);
    
    for b = 1:actualNumBins
        idx = find(wheelSpeed >= edges(b) & wheelSpeed < edges(b+1));
        mMean(:, b) = mean(sigData(:, idx), 2, 'omitnan');
        mSEM(:, b)  = std(sigData(:, idx), 0, 2, 'omitnan') / sqrt(length(idx));
        mCount(:, b) = length(idx);
    end
    
    % Move Only modulation
    moveMax = max(mMean, [], 2);
    moveMin = min(mMean, [], 2);
    % Saleem et al 2013
    modRatio_move = 2 * (moveMax - moveMin) ./ (moveMax + moveMin);
    
    % Full (Move + Stat) modulation
    fullCurve = [statMean, mMean];
    fullMax = max(fullCurve, [], 2);
    fullMin = min(fullCurve, [], 2);
    modRatio_full = 2 * (fullMax - fullMin) ./ (fullMax + fullMin);
    
    %% sig + variance
    realVarFull = var(fullCurve, 0, 2, 'omitnan');
    realVarMovingOnly = var(mMean, 0, 2, 'omitnan');
    [~, peakIdx] = max(fullCurve, [], 2);
    peakSpeed = speedAxis(peakIdx);
    
    pValFull = ones(numROIs, 1);
    pValMoving = ones(numROIs, 1);
    
    if shuffle
        fprintf('Shuffling %s...\n', fname);
        shuffVarsFull = zeros(numROIs, numShifts);
        shuffVarsMoving = zeros(numROIs, numShifts);
        for iperm = 1:numShifts
            sigS = circshift(sigData, shiftValues(iperm), 2);
            sM_shuff = mean(sigS(:, sIdx), 2, 'omitnan');
            mM_shuff = zeros(numROIs, actualNumBins);
            for b = 1:actualNumBins
                idxB = wheelSpeed >= edges(b) & wheelSpeed < edges(b+1);
                mM_shuff(:, b) = mean(sigS(:, idxB), 2, 'omitnan');
            end
            shuffVarsFull(:, iperm) = var([sM_shuff, mM_shuff], 0, 2, 'omitnan');
            shuffVarsMoving(:, iperm) = var(mM_shuff, 0, 2, 'omitnan');
        end
        pValFull = sum(shuffVarsFull >= realVarFull, 2) / numShifts;
        pValMoving = sum(shuffVarsMoving >= realVarMovingOnly, 2) / numShifts;
    end
    
    % Store
    response.tuningCurve.(fname).statMean = statMean;
    response.tuningCurve.(fname).statSEM  = statSEM;
    response.tuningCurve.(fname).statCount = statCount;
    response.tuningCurve.(fname).moveMean = mMean;
    response.tuningCurve.(fname).moveSEM  = mSEM;
    response.tuningCurve.(fname).moveCount = mCount;
    
    response.tuningCurve.(fname).modulationRatio_move = modRatio_move;
    response.tuningCurve.(fname).percentModulated_move = modRatio_move * 100;
    response.tuningCurve.(fname).modulationRatio_full = modRatio_full;
    response.tuningCurve.(fname).percentModulated_full = modRatio_full * 100;
    
    response.tuningCurve.(fname).peakSpeed = peakSpeed;
    response.tuningCurve.(fname).tuningVarFull = realVarFull;
    response.tuningCurve.(fname).tuningVarMoving = realVarMovingOnly;
    response.tuningCurve.(fname).pValFull = pValFull;
    response.tuningCurve.(fname).pValMoving = pValMoving;
    response.tuningCurve.(fname).isSignificant_999 = pValFull <= 0.001;
    response.tuningCurve.(fname).isSignificantMoving_999 = pValMoving <= 0.001;
end

stimFileName = sprintf('%s_%s_Response_%s_FixedBins.mat', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, stimName);
savePath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
save(savePath, 'response', '-v7.3');
end