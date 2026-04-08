
%% 1. Load Data
% Gray Screen
grayProc = load("Z:\ibn-vision\DATA\SUBJECTS\M26004\Analysis\20260318_oldNeuropilExtraction\M26004_20260318_processed2PData_M26004_GrayScreen_20260318_00001.mat");
grayPeriph = load("Z:\ibn-vision\DATA\SUBJECTS\M26004\Analysis\20260318_oldNeuropilExtraction\M26004_20260318_PeripheralData_M26004_GrayScreen_20260318_00001.mat");

% Darkness
darkProc = load("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_processed2PData_M25132_Darkness_20260226_00001.mat");
darkPeriph = load("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_PeripheralData_M25132_Darkness_20260226_00001.mat");

%% 2. Parameters & PDF Setup
savePath = "Z:\ibn-vision\USERS\Sonali\Figures\DarknessVsGray\M26004\20260318\";
if ~exist(savePath, 'dir'), mkdir(savePath); end
pdfFileName = fullfile(savePath, 'ROI_SpikesAndDFFs.pdf');

roiList = 150:175; 
fs = 60; 
tickToCmConversion = 3.1415 * 20 / 1024; 
stationaryPlotX = 0.5; 
colors = {[0.6 0.6 0.6], [0 0 0]}; % Gray for Gray Screen, Black for Darkness
condNames = {'Gray Screen', 'Darkness'};
procData = {grayProc, darkProc};
periphData = {grayPeriph, darkPeriph};

% Zoom window for the traces (20-second snippet)
zoomStart = 300; 
zoomDuration = 20; 
zoomWin = [zoomStart, zoomStart + zoomDuration];

fprintf('Starting PDF Export. This may take a minute...\n');

%% 2. Loop Through ROIs
for r = 1:length(roiList)
    roiIdx = roiList(r);
    
    % Adjusted figure height to 12 inches to fit the 6th subplot comfortably
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0 0 9 12], 'Visible', 'off');
    
    % 6 rows now: 3 Tuning, 2 Trace Overlays, 1 Occupancy
    tLayout = tiledlayout(6, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
    title(tLayout, sprintf('ROI %d', roiIdx), 'FontSize', 16, 'FontWeight', 'bold');

    %% 
    signalFields = {'spks', 'dFFNeuropilCorrected', 'dFF'};
    signalLabels = {'Average inferred spike amplitudes', 'Average dF/F Neuropil Corrected', 'Average dF/F'};
    
    % 
    allBinWidths = cell(1,2);
    allBinCounts = cell(1,2);
    allBinSpeeds = cell(1,2);

    for s = 1:3
        ax = nexttile; hold(ax, 'on');
        for i = 1:2
            currP = procData{i}; currW = periphData{i};
            v = [0; diff(currW.peripheralData.Wheel.Value * tickToCmConversion)] ./ [1; diff(currW.peripheralData.Wheel.sampleTimes)];
            v(abs(v) > 100) = 0; 
            
            if strcmp(signalFields{s}, 'spks')
                if isfield(currP, 'spks'), sig = currP.spks(roiIdx, :) * fs;
                else, sig = currP.processedSignals.spks(roiIdx, :) * fs; end
            else
                sig = currP.processedSignals.(signalFields{s})(roiIdx, :);
            end
            
            % Binning
            mIdx = v >= 1; sIdx = v < 1;
            ptsPerBin = floor(0.07 * length(v)); 
            nB = floor(sum(mIdx) / ptsPerBin);
            edg = quantile(v(mIdx), linspace(0, 1, max(1, nB) + 1));
            
            mS = []; mA = []; mE = []; mCount = [];
            for b = 1:length(edg)-1
                idx = v >= edg(b) & v < edg(b+1);
                if any(idx)
                    mS(end+1) = mean(v(idx)); 
                    mA(end+1) = mean(sig(idx)); 
                    mE(end+1) = std(sig(idx)) / sqrt(sum(idx));
                    mCount(end+1) = sum(idx); % Save for occupancy
                end
            end
            
            % Save occupancy data only once (during first signal loop)
            if s == 1
                allBinCounts{i} = [sum(sIdx), mCount];
                allBinSpeeds{i} = [stationaryPlotX, mS];
            end
            
            errorbar(ax, [stationaryPlotX; mS'], [mean(sig(sIdx)); mA'], [std(sig(sIdx))/sqrt(sum(sIdx)); mE'], 'o-', ...
                'Color', colors{i}, 'LineWidth', 1.5, 'MarkerFaceColor', colors{i}, 'DisplayName', condNames{i});
        end
        set(ax, 'XScale', 'log'); xlim([0.4, 100]); grid on;
        xticks([0.5, 1, 5, 10, 20, 40, 60, 100]); 
        if s == 3, xticklabels({'0','1','5','10','20','40','60','100'}); else, xticklabels([]); end
        ylabel(signalLabels{s});
    end
    
    %% occupacy 
    axOcc = nexttile; hold on;
    for i = 1:2
        % Convert counts to seconds
        occSeconds = allBinCounts{i} / fs;
        % Plot as a bar chart
        b = bar(allBinSpeeds{i}, occSeconds, 'FaceColor', colors{i}, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.4, 'DisplayName', condNames{i});
        % Adjust bar width for log scale visibility
        if i == 1, b.BarWidth = 0.4; else, b.BarWidth = 0.2; end
    end
    set(axOcc, 'XScale', 'log'); xlim([0.4, 100]);
    xticks([0.5, 1, 5, 10, 20, 40, 60, 100]); xticklabels({'0','1','5','10','20','40','60','100'});
    ylabel('Occupancy (s)'); xlabel('Running Speed (cm/s)');
    title('Time Spent per Speed Bin'); grid on; box off;
    legend('Location', 'northeast', 'FontSize', 8, 'Box', 'off');
    %% line plots 
    for i = 1:2
        axTrace = nexttile; hold on;
        currP = procData{i};
        if isfield(currP, 'spks'), spkDat = currP.spks(roiIdx, :); else, spkDat = currP.processedSignals.spks(roiIdx, :); end
        
        % Try to find F
        if isfield(currP.processedSignals, 'F'), rawF = currP.processedSignals.F(roiIdx, :);
        elseif isfield(currP, 'Fcell'), rawF = currP.Fcell(roiIdx, :);
        else, rawF = currP.processedSignals.dFF(roiIdx, :); end
        
        tv = (0:length(spkDat)-1)/fs;
        zoomIdx = tv >= zoomWin(1) & tv <= zoomWin(2);
        
        yyaxis left; plot(tv, rawF); ylabel('Raw F');
        if any(zoomIdx), ylim([min(rawF(zoomIdx))*0.98, max(rawF(zoomIdx))*1.02]); end
        %set(gca);
        
        yyaxis right; plot(tv, spkDat, 'LineWidth', 1.2); ylabel('Spks');
        if any(zoomIdx), ylim([0, max(spkDat(zoomIdx))*1.2 + 0.01]); end
        %set(gca);
        
        xlim(zoomWin); title(sprintf('%s Overlay', condNames{i}), 'FontSize', 10);
        box off;
    end

    %% 
   

    %% 3. PDF Export
    drawnow;
    if r == 1, exportgraphics(fig, pdfFileName, 'ContentType', 'vector');
    else, exportgraphics(fig, pdfFileName, 'ContentType', 'vector', 'Append', true); end
    close(fig);
    fprintf('ROI %d added to PDF.\n', roiIdx);
end
    close(fig);
    fprintf('ROI %d added to PDF.\n', roiIdx);


fprintf('Full Diagnostic PDF saved to: %s\n', pdfFileName);




useZScoredSignals = false;

% Locate Stimulus
stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
if isempty(stimIdx), error('Stimulus not found!'); end

% Load Data
data2P = load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData);
dataPeriph = load(sessionFileInfo.stimFiles(stimIdx).processedPeripheralData, 'peripheralData');

%% 
% Select Signal (Z-scored or Raw)

sigAll = data2P.spks; 

% Calculate & Clean Speed (CRITICAL for avoiding NaNs)
tickToCm = 3.1415 * 20 / 1024; 
wheelVal = dataPeriph.peripheralData.Wheel.Value;
wheelTime = dataPeriph.peripheralData.Wheel.sampleTimes;

varValues = [0; diff(wheelVal * tickToCm)] ./ [1; diff(wheelTime)];
varValues = abs(varValues); % Edward's logic works best with positive speed
varValues(isnan(varValues)) = 0; % Replace NaNs with 0 to keep varAxis clean

%% 
options = struct;
options.binSize = 0.1;
options.nBins = 10;
options.catBelow = 1;       % Your stationary threshold
options.varLims = [0 80];   % Clip at 60cm/s
options.getSig = true;
options.nPerms = 1000;      % 1000 random shuffles
options.fitCustomFun = true; % Trigger the Global Model fit
options.CVcustomFun = true;  % Trigger the Global-to-Fold R2
options.kfold = 5;

%%
numROIs = 10; % 
sessionResults = cell(numROIs, 1);

fprintf('Running Edward''s Logic for %d ROIs...\n', numROIs);
tic;
for i = 1:numROIs
    fprintf('ROI %d... ', i);
    

    [m, sem, speeds, p, p2, binned, peak, model, cv, vals, finalOpts] = ...
        getTuningCurve_2P_Complete(sigAll(i,:), wheelTime, varValues, options);
    
    sessionResults{i}.meanTuning = m;
    sessionResults{i}.semTuning  = sem;
    sessionResults{i}.speeds     = speeds;
    sessionResults{i}.pVal       = p;      % Full curve significance
    sessionResults{i}.pVal2      = p2;     % Moving-only significance
    sessionResults{i}.model      = model;  % Contains .params and .cv_R2_avg
    
    % Print progress
    if isfield(model, 'cv_R2_avg')
        fprintf('Done (p=%.3f, R2=%.2f)\n', p, model.cv_R2_avg);
    else
        fprintf('Done (p=%.3f, R2=FAILED)\n', p);
    end
end
toc;