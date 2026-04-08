function [response] = processSessionTuningCurves(sessionFileInfo, stimName, useZScoredSignals, plotflag)
% processSessionTuningCurves: 
% 1. Finds data via sessionFileInfo
% 2. Bins every signal (dFF, spks, etc.) using the 7% occupancy rule
% 3. Calculates 2 s.e. error bars
% 4. Performs circular-shift permutation testing for significance (p-values)

if nargin < 3 || isempty(useZScoredSignals), useZScoredSignals = false; end
if nargin < 4, plotflag = true; end

%% 1. Automated Data Loading
% Find stimulus index
stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
if isempty(stimIdx), error('Stimulus %s not found.', stimName); end

% Get Paths
mergedPath = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;
periphPath = sessionFileInfo.stimFiles(stimIdx).processedPeripheralData;

% Load
data2P = load(mergedPath);
dataPeriph = load(periphPath, 'peripheralData');

% Choose Signal Source
sourceName = 'processedSignals';
if useZScoredSignals, sourceName = 'zScoredProcessedSignals'; end
signals = data2P.(sourceName);

% Ensure independent spks are included
if isfield(data2P, 'spks') && ~isempty(data2P.spks), signals.spks = data2P.spks; end

%% 2. Speed Calculation & 7% Binning Setup
fs = 60; 
tickToCm = 3.1415 * 20 / 1024; 
% Calculate speed from wheel ticks and timestamps
wheelSpeed = [0; diff(dataPeriph.peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(dataPeriph.peripheralData.Wheel.sampleTimes)];
wheelSpeed(abs(wheelSpeed) > 150) = NaN; % Glitch removal

mIdx = wheelSpeed > 1 & ~isnan(wheelSpeed); 
sIdx = wheelSpeed <= 1 & ~isnan(wheelSpeed); 

% Calculate number of bins based on 7% occupancy requirement
ptsPerBin = floor(0.07 * length(wheelSpeed)); 
nB = floor(sum(mIdx) / ptsPerBin);

% Define Edges (Quantile for Moving, Fixed for Stationary)
edges = quantile(wheelSpeed(mIdx), linspace(0, 1, max(1, nB) + 1));
edges = unique(edges); 
numMoveBins = length(edges) - 1;

%% permutation
allFields = fieldnames(signals);
response.speedBins = edges;
response.occupancy.stationary = sum(sIdx) / fs;
response.occupancy.moving = zeros(1, numMoveBins);

nPerms = 1000; % Number of shuffles for p-value

for f = 1:length(allFields)
    fn = allFields{f};
    sigData = signals.(fn);
    if size(sigData, 2) ~= length(wheelSpeed), continue; end
    numROIs = size(sigData, 1);
    
    % Initialize Results
    response.(fn).statMean = mean(sigData(:, sIdx), 2, 'omitnan');
    response.(fn).statSEM  = std(sigData(:, sIdx), 0, 2, 'omitnan') / sqrt(sum(sIdx));
    response.(fn).moveMean = zeros(numROIs, numMoveBins);
    response.(fn).moveSEM  = zeros(numROIs, numMoveBins);
    response.(fn).pVal     = zeros(numROIs, 1);

    % Pre-calculate Bin Indices
    binIndices = cell(numMoveBins, 1);
    for b = 1:numMoveBins
        binIndices{b} = wheelSpeed >= edges(b) & wheelSpeed < edges(b+1);
        if f == 1, response.occupancy.moving(b) = sum(binIndices{b}) / fs; end
    end

    % Process each ROI
    fprintf('Processing %s Tuning and Stats...\n', fn);
    for r = 1:numROIs
        roiSig = sigData(r, :);
        
        % Calculate Real Means
        for b = 1:numMoveBins
            response.(fn).moveMean(r, b) = mean(roiSig(binIndices{b}), 'omitnan');
            response.(fn).moveSEM(r, b)  = std(roiSig(binIndices{b}), 0, 'omitnan') / sqrt(sum(binIndices{b}));
        end
        
        % Significance: Circular Shifting
        realCurve = [response.(fn).statMean(r), response.(fn).moveMean(r, :)];
        realVar = nanvar(realCurve);
        permVars = zeros(nPerms, 1);
        
        for iperm = 1:nPerms
            % Shift signal by at least 10 seconds to decouple from behavior
            shuffSig = circshift(roiSig, randi([600, length(roiSig)-600]));
            
            % Re-bin shuffled data
            sMeanShuff = mean(shuffSig(sIdx), 'omitnan');
            mMeanShuff = zeros(1, numMoveBins);
            for b = 1:numMoveBins
                mMeanShuff(b) = mean(shuffSig(binIndices{b}), 'omitnan');
            end
            permVars(iperm) = nanvar([sMeanShuff, mMeanShuff]);
        end
        response.(fn).pVal(r) = sum(permVars > realVar) / nPerms;
    end
end

%% 4. Automated Saving
saveName = sprintf('%s_%s_Tuning_%s.mat', sessionFileInfo.animal_name, sessionFileInfo.session_name, stimName);
savePath = fullfile(sessionFileInfo.Directories.save_folder, saveName);
save(savePath, 'response', '-v7.3');
fprintf('Saved to: %s\n', savePath);

%% 5. Plotting with Occupancy Overlays
if plotflag
    roiIdx = 1; % Default to first ROI for sanity check
    plot_tuning_summary(response, roiIdx, stimName);
end

end

function plot_tuning_summary(response, roiIdx, stimName)
    fields = intersect(fieldnames(response), {'spks', 'spikes', 'dFF', 'dFFNeuropilCorrected'}, 'stable');
    numPlots = length(fields);
    fig = figure('Color', 'w', 'Name', [stimName ' - ROI ' num2str(roiIdx)]);
    tlo = tiledlayout(numPlots, 1, 'TileSpacing', 'compact');
    
    x_vals = [0.5, response.speedBins(1:end-1)];
    occ_vals = [response.occupancy.stationary, response.occupancy.moving];
    
    for p = 1:numPlots
        fn = fields{p};
        ax = nexttile; hold on;
        
        y_m = [response.(fn).statMean(roiIdx); response.(fn).moveMean(roiIdx, :)'];
        y_e = [response.(fn).statSEM(roiIdx);  response.(fn).moveSEM(roiIdx, :)'] * 2; % 2 s.e.
        
        errorbar(ax, x_vals, y_m, y_e, 'o-k', 'MarkerFaceColor', 'k', 'LineWidth', 1.2);
        
        % Add Occupancy Text
        for i = 1:length(x_vals)
            text(ax, x_vals(i), y_m(i) + y_e(i), sprintf(' %.1fs', occ_vals(i)), ...
                'FontSize', 7, 'Horiz', 'center', 'Vert', 'bottom', 'Color', [0.4 0.4 0.4]);
        end
        
        set(ax, 'XScale', 'log'); grid on;
        xticks(ax, [0.5, 1, 5, 10, 20, 40, 60, 100]);
        xticklabels(ax, {'0','1','5','10','20','40','60','100'});
        title(sprintf('%s (p = %.3f)', fn, response.(fn).pVal(roiIdx)));
        ylabel(fn);
    end
    xlabel(tlo, 'Running Speed (cm/s)');
end