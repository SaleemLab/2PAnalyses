function [response] = Copy_of_binAllNeuralSignalsByRunningSpeeds(sessionFileInfo, stimName, useZScoredSignals, shuffle,plotflag)
% Saleem et al, 2014 Significance: Real Variance > 99.9% of Shuffled Variance.
% added 99% (p < 0.01) threshold.

if nargin < 3 || isempty(useZScoredSignals), useZScoredSignals = false; end
if nargin < 4, plotflag = true; end

%% Load Data
stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
if isempty(stimIdx), error('Stimulus %s not found.', stimName); end

data2P = load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData);
dataPeriph = load(sessionFileInfo.stimFiles(stimIdx).processedPeripheralData, 'peripheralData');

sourceStructName = 'processedSignals';
if useZScoredSignals, sourceStructName = 'zScoredProcessedSignals'; end
signals = data2P.(sourceStructName);

%% Speed & Binning (7% Rule)
fs = 60; 
tickToCm = 3.1415 * 20 / 1024; 
wheelSpeed = [0; diff(dataPeriph.peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(dataPeriph.peripheralData.Wheel.sampleTimes)];
wheelSpeed(abs(wheelSpeed) > 150) = NaN; 

% Define speed indices
mIdx = wheelSpeed > 1 & ~isnan(wheelSpeed); % running
sIdx = wheelSpeed <= 1 & ~isnan(wheelSpeed); % stationary

 % each speed bin contains at least 7% of the condition
ptsPerBin = floor(0.07 * length(wheelSpeed));
nB = floor(sum(mIdx) / ptsPerBin);

% Quantile binning makes the count of frames in each bin equal.
edges = quantile(wheelSpeed(mIdx), linspace(0, 1, max(1, nB) + 1));
edges = unique(edges); 

%% Shuffle setp

maxShift = round(length(data2P.TwoPFrameTime)/2); % maximum number of elements by which the signal can be circularly shifted 
minShift = 600; % 10s ? 
numShifts = 1000;  % number of times the randomization loop will run
% Generates a vector of numShifts random integers.
rng(1) % the same shift value will apply to all rois. 
shiftVal = randi([minShift maxShift] ,[1 numShifts]); 
%% Currently processing dff/dffneu/spks
allFields = fieldnames(signals);
response.speedBins = edges;
response.source = sourceStructName;

% Initialize Occupancy
response.occupancy.stationary = sum(sIdx) / fs;
response.occupancy.moving = zeros(1, length(edges)-1);

% For circular shifts:

for f = 1:length(allFields)
    fname = allFields{f};
    sigData = signals.(fname);
    if size(sigData, 2) ~= length(wheelSpeed), continue; end
    
    numROIs = size(sigData, 1);
    numBins = length(edges)-1;
    
    % means/SEMs
    statMean = mean(sigData(:, sIdx), 2, 'omitnan');
    statSEM  = std(sigData(:, sIdx), 0, 2, 'omitnan') / sqrt(sum(sIdx));
    mMean = zeros(numROIs, numBins);
    mSEM  = zeros(numROIs, numBins);
    
    binIndices = cell(numBins, 1);
    for b = 1:numBins
        binIndices{b} = find(wheelSpeed >= edges(b) & wheelSpeed < edges(b+1));
        mMean(:, b) = mean(sigData(:, binIndices{b}), 2, 'omitnan');
        mSEM(:, b)  = std(sigData(:, binIndices{b}), 0, 2, 'omitnan') / sqrt(length(binIndices{b}));
        
        % Fill occupancy once
        if f == 1, response.occupancy.moving(b) = length(binIndices{b}) / fs; end
    end
    
    %Saleem et al. 2014
    fprintf('Calculating 99.9%% Significance for %s...\n', fname);
    realVar = var([statMean, mMean], 0, 2, 'omitnan');
    shuffVars = zeros(numROIs, nPerms);
    
    for thisShift = 1:numShifts

        sigS = circshift(sigData, shiftVal(thisShift), 2);
        
        sM = mean(sigS(:, sIdx), 2, 'omitnan');
        mM = zeros(numROIs, numBins);
        for b = 1:numBins
            mM(:, b) = mean(sigS(:, binIndices{b}), 2, 'omitnan');
        end
        shuffVars(:, thisShift) = var([sM, mM], 0, 2, 'omitnan');
    end
    
    % Calculate p-values
    pVals = sum(shuffVars >= realVar, 2) / nPerms;
    
    % Store everything
    response.(fname).statMean = statMean;
    response.(fname).statSEM  = statSEM;
    response.(fname).moveMean = mMean;
    response.(fname).moveSEM  = mSEM;
    response.(fname).pVal     = pVals;
    response.(fname).isSignificant_999 = pVals <= 0.001; % Saleem et al.
    response.(fname).isSignificant_99  = pVals <= 0.01;  % 99% Threshold
end

%% 5. Save & Plot
stimFileName = sprintf('%s_%s_Response_%s.mat', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, stimName);
savePath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
save(savePath, 'response', '-v7.3');
if plotflag
    roiIdx = 1; 
    fieldsInResp = fieldnames(response);
    plotCandidates = {'spks', 'spikes', 'dFF', 'dFFNeuropilCorrected'};
    fieldsToPlot = intersect(fieldsInResp, plotCandidates, 'stable');
    numPlots = length(fieldsToPlot);
    
    if numPlots > 0
        fig = figure('Color', 'w', 'Name', stimName, 'Position', [100 100 800 900]);
        % No extra tile for occupancy anymore, just the signal plots
        tlo = tiledlayout(numPlots, 1, 'TileSpacing', 'compact', 'Padding', 'loose'); 
        
        % Get occupancy values once
        occ_vals = [response.occupancy.stationary, response.occupancy.moving(:)'];
        x_vals = [0.5, response.speedBins(1:end-1)];

        for thisShift = 1:numPlots
            fn = fieldsToPlot{thisShift};
            ax = nexttile(tlo); hold(ax, 'on');
            
            y_mean = [response.(fn).statMean(roiIdx); response.(fn).moveMean(roiIdx,:)'];
            y_err  = [response.(fn).statSEM(roiIdx);  response.(fn).moveSEM(roiIdx,:)'] * 2;
            
            % 1. Plot the Tuning Curve
            errorbar(ax, x_vals, y_mean, y_err, 'o-k', 'MarkerFaceColor', 'k', 'LineWidth', 1.2);
            
            % 2. Add Occupancy Labels above each point
            for i = 1:length(x_vals)
                % Place text slightly above the error bar
                text_y = y_mean(i) + y_err(i); 
                text(ax, x_vals(i), text_y, sprintf(' %.1fs', occ_vals(i)), ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 7, 'Color', [0.4 0.4 0.4]);
            end
            
            % Formatting
            set(ax, 'XScale', 'log'); grid(ax, 'on');
            xticks(ax, [0.5, 1, 5, 10, 20, 40, 60, 100]);
            xticklabels(ax, {'0','1','5','10','20','40','60','100'});
            xlim(ax, [0.4, 110]);
            
            % Add some headroom for the text labels
            curr_ylim = ylim(ax);
            ylim(ax, [curr_ylim(1), curr_ylim(2) * 1.2]);
            
            ylabel(ax, [fn ' (\pm 2 s.e.)'], 'Interpreter', 'none');
            if thisShift == 1, title(ax, sprintf('ROI %d Speed Tuning with Occupancy (s)', roiIdx)); end
        end
        xlabel(tlo, 'Running Speed (cm/s)', 'FontWeight', 'bold');
    end
end
end