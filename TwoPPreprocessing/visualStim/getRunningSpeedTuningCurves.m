function [response] = getRunningSpeedTuningCurves(sessionFileInfo, stimName, useZScoredSignals, shuffle)
% From Saleem et al, 2014 Significance: Real Variance > 99.9% of Shuffled
% Variance. Fig 2 Darkness 
% Can you this for plotting: plotRunningTuningCurves(respGray, respDark, pdfPath)
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

%% 
fs = 60; 
tickToCm = 3.1415 * 20 / 1024; 
wheelSpeed = [0; diff(dataPeriph.peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(dataPeriph.peripheralData.Wheel.sampleTimes)];
wheelSpeed(abs(wheelSpeed) > 150) = NaN; 
mIdx = find(wheelSpeed > 1 & ~isnan(wheelSpeed)); % running
sIdx = find(wheelSpeed <= 1 & ~isnan(wheelSpeed)); % stationary
ptsPerBin = floor(0.07 * length(wheelSpeed)); % Set points per bin to 7% of total session length
nB = floor(length(mIdx) / ptsPerBin); % Calculate number of moving bins based on available moving data 
edges = quantile(wheelSpeed(mIdx), linspace(0, 1, max(1, nB) + 1)); % create speed bin edges using quantiles (equally distributed data points per bin)
edges = unique(edges); 
numBins = length(edges)-1;
speedAxis = [0, edges(1:end-1) + diff(edges)/2];

response.tuningCurve.speedBins = edges;
response.tuningCurve.source = sourceStructName;
response.tuningCurve.occupancy.stationary = length(sIdx) / fs;
response.tuningCurve.occupancy.moving = zeros(1, numBins);

%% Shuffle 
numShifts = 1000;  
minShift = 600; % 10s
maxShift = length(wheelSpeed) - minShift; 
rng(1); 
shiftValues = randi([minShift, maxShift], [1, numShifts]);

allFields = fieldnames(signals);
for thisField = 1:length(allFields)
    fname = allFields{thisField};
    sigData = signals.(fname);
    if size(sigData, 2) ~= length(wheelSpeed), continue; end
    
    numROIs = size(sigData, 1);
    mMean = zeros(numROIs, numBins);
    mSEM  = zeros(numROIs, numBins);
    mCount = zeros(numROIs, numBins);
    
    % Stationary Stats (Back to SEM)
    statMean = mean(sigData(:, sIdx), 2, 'omitnan');
    statSEM  = std(sigData(:, sIdx), 0, 2, 'omitnan') / sqrt(length(sIdx));
    statCount = repmat(length(sIdx), numROIs, 1);
    
    for b = 1:numBins
        idx = find(wheelSpeed >= edges(b) & wheelSpeed < edges(b+1));
        mMean(:, b) = mean(sigData(:, idx), 2, 'omitnan');
        mSEM(:, b)  = std(sigData(:, idx), 0, 2, 'omitnan') / sqrt(length(idx));
        mCount(:, b) = length(idx);
        if thisField == 1, response.tuningCurve.occupancy.moving(b) = length(idx) / fs; end
    end
    
    % Move Only modulation
    moveMax = max(mMean, [], 2);
    moveMin = min(mMean, [], 2);
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
            mM_shuff = zeros(numROIs, numBins);
            for b = 1:numBins
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

stimFileName = sprintf('%s_%s_Response_%s.mat', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, stimName);
savePath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
save(savePath, 'response', '-v7.3');
end

%%
% could use for plotting: 
% plotRunningTuningCurves(respGray, respDark, pdfPath)
%%
% To check variance during debugging: not saving the shiffle distribution
% if needed 
% debugROI = 4; 
% debugField = fname; % or 'dFF'
% 
% %
% realStat = mean(sigData(debugROI, sIdx), 'omitnan');
% realMove = zeros(1, numBins);
% for b = 1:numBins
%     realMove(b) = mean(sigData(debugROI, binIndices{b}), 'omitnan');
% end
% realCurve = [realStat, realMove];
% realVar = var(realCurve, 0, 2, 'omitnan');
% 
% % Run a quick mini-shuffle (50 perms) just for this ROI to see the "cloud"
% numDebugShifts = 50;
% debugShuffVars = zeros(1, numDebugShifts);
% figure('Color', 'w', 'Name', ['Debug ROI ' num2str(debugROI)]);
% subplot(1,2,1); hold on;
% 
% for i = 1:numDebugShifts
%     shift = randi([600, length(wheelSpeed)-600]);
%     sigS = circshift(sigData(debugROI, :), shift);
% 
%     sM = mean(sigS(sIdx), 'omitnan');
%     mM = zeros(1, numBins);
%     for b = 1:numBins
%         mM(b) = mean(sigS(binIndices{b}), 'omitnan');
%     end
%     shuffCurve = [sM, mM];
%     debugShuffVars(i) = var(shuffCurve, 0, 2, 'omitnan');
% 
%     % Plot the "Null" grey lines
%     plot([0, 1:numBins], shuffCurve, 'Color', [0.8 0.8 0.8]);
% end
% 
% % 3. Plot the Real curve on top
% plot([0, 1:numBins], realCurve, 'r', 'LineWidth', 2);
% title('Tuning: Real (Red) vs Shuffled (Grey)');
% xticks([0, 1:numBins]); xticklabels(['Stat', string(1:numBins)]);
% grid on;
% 
% % 4. Plot the Variance Distribution
% subplot(1,2,2);
% histogram(debugShuffVars, 'FaceColor', [.5 .5 .5]); hold on;
% xline(realVar, 'r', 'LineWidth', 2, 'Label', 'Real Var');
% title('Variance Distribution');
% xlabel('Variance'); ylabel('Count');