% E Horrcks 
% Editting for sonali's twop bouton analyses 
%% specify sessions
% RSP only 
filteredTable = filterMasterTable('MouseID',{'M25132', 'M25133', 'M26003', }, 'HasStimulus', 'DotMotion_SpeedTuning','Suite2PPreprocessing', 1, 'Exclude', 0);
dotFieldRSP_sessionTags = sessionsToProcess(filteredTable);
%% load session
% isession = 1;
% inputSuffix = 'basic.mat';
% tic
% inputFileName = [sessionTags{isession,1},'_', sessionTags{isession,2},'_', inputSuffix];
% load(fullfile(dataDir,inputFileName))
%% get stat/run trials
% change to load response 
if ~exist('response', 'var') || ~exist('bonsaiData', 'var')
    error('Please ensure both "response" and "bonsaiData" structures are loaded in the workspace.');
end

% --- build trials table from 2P structure ---
nGroups = numel(response.wheelData);
trialsSpeed2D = struct('wheelSpeed', {}, 'VelX1', {}, 'Contrast1', {}, 'numDots1', {}, 'runFlag', {}, 'startTime', {});
trialCounter = 1;

for g = 1:nGroups
    grpWheel  = response.wheelData(g);
    grpBonsai = bonsaiData.trialGroups(g); 
    
    speedMatrix  = grpWheel.alignedResponses; 
    tVec         = grpWheel.timeVector;
    stimFramesMask = (tVec >= 0 & tVec <= 2.0); % Evaluate running during active 2s stimulus window
    
    for ti = 1:size(speedMatrix, 2)
        singleTrialTrace = speedMatrix(:, ti);
        if all(isnan(singleTrialTrace)), continue; end
        
        meanSpeed = nanmean(singleTrialTrace(stimFramesMask));
        propRunning = sum(singleTrialTrace(stimFramesMask) > 0.5) / sum(stimFramesMask);
        propStationary = sum(singleTrialTrace(stimFramesMask) < 3) / sum(stimFramesMask);
        
        runFlag = NaN; 
        if propRunning >= 0.75 && meanSpeed > 3
            runFlag = 1;
        elseif propStationary >= 0.75 && meanSpeed < 0.5
            runFlag = 0;
        end
        if isnan(runFlag), continue; end
        
        trialID = grpBonsai.trials(ti); 
        
        trialsSpeed2D(trialCounter).wheelSpeed = singleTrialTrace(stimFramesMask);
        trialsSpeed2D(trialCounter).VelX1      = grpWheel.stimValue;
        
        % edd uses the number of dots to find black trials later in the
        % script 
        if grpWheel.stimValue == 1
            trialsSpeed2D(trialCounter).numDots1   = 0;
        else
            trialsSpeed2D(trialCounter).numDots1   = 573;
        end
        
        trialsSpeed2D(trialCounter).Contrast1  = 1;   
        trialsSpeed2D(trialCounter).runFlag    = runFlag;
        trialsSpeed2D(trialCounter).startTime  = bonsaiData.onARDTimes(trialID);
        
        trialCounter = trialCounter + 1;
    end
end
% --- end conversion ---

trialsSpeed2D(1) = [];
tsd = trialsSpeed2D;

%split trials by state according to wheel data
run_idx = find(cellfun(@(x) (sum(x>0.5)/numel(x))>=0.75 & mean(x)>3, {tsd.wheelSpeed}));
stat_idx = find(cellfun(@(x) (sum(x<3)/numel(x))>=0.75 & mean(x)<0.5, {tsd.wheelSpeed}));
[tsd.runFlag] = deal([nan]);
[tsd(stat_idx).runFlag] = deal([0]);
[tsd(run_idx).runFlag] = deal([1]);
temp_tsd = tsd([tsd.numDots1]==573); % remove blnk trials.
temp_tsd = temp_tsd(~isnan([temp_tsd.runFlag]));
allStimConds = [vertcat(temp_tsd.VelX1), vertcat(temp_tsd.Contrast1), vertcat(temp_tsd.runFlag)];
[uniqueConds, ~, ic] = unique(allStimConds, 'rows');
tally = accumarray(ic, 1);
Result = [uniqueConds tally];
temp_tsd = temp_tsd([temp_tsd.Contrast1]==1); % only full contrast trials

%% get trial-based spike counts
for itrial =1 :numel(temp_tsd)
    temp_tsd(itrial).start_time = temp_tsd(itrial).startTime;
    temp_tsd(itrial).absVel = abs(temp_tsd(itrial).VelX1);
end

% --- Mapping 2P Calcium traces into the units structure ---
% alignedResponses typically includes 2s pre and 4s post for dot fields
nROIs = size(response.psthData(1).alignedResponses, 1);
stimWindowMask = (response.psthData(1).timeVector >= 0 & response.psthData(1).timeVector <= 2.5); % 2s window
uniqueVelocities = unique(abs(allStimConds(:, 1)));
nSpeeds = numel(uniqueVelocities);
units = struct();

for thisROI = 1:nROIs
    units(thisROI).allSpikes = cell(nSpeeds, 2);
    for s = 1:nSpeeds
        for istate = 1:2
            matchingTrials = find(abs([temp_tsd.VelX1]) == uniqueVelocities(s) & [temp_tsd.runFlag] == (istate - 1));
            traceAccumulator = zeros(1, numel(matchingTrials));
            for mt = 1:numel(matchingTrials)
                origGroup = find([response.psthData.stimValue] == temp_tsd(matchingTrials(mt)).VelX1, 1);
                fullTrace = squeeze(response.psthData(origGroup).alignedResponses(thisROI, :, mt));
                traceAccumulator(mt) = nanmean(fullTrace(stimWindowMask));
            end
            units(thisROI).allSpikes{s, istate} = traceAccumulator;
        end
    end
    units(thisROI).tuning = cellfun(@nanmean, units(thisROI).allSpikes);
end
%end mapping

statTsd = temp_tsd([temp_tsd.runFlag]==0);
runTsd = temp_tsd([temp_tsd.runFlag]==1);

% have modified this values to match the stimulus: 2s on and 2s off 
options.intervalStart = 0;
options.intervalEnd = 2.5;
options.binSpacing=2.0;

%% downsample to min # trials available
% set so equal # trials per conditionans
minTrial = min(min(cellfun(@(x) size(x,2), units(1).allSpikes)));
for thisROI = 1:numel(units)
    units(thisROI).allSpikesDownsample = cellfun(@(x) x(:,1:minTrial), units(thisROI).allSpikes,'UniformOutput', false);
end
%% get cross-validated R2 (tuning strength) 
% options
r2opts.nPerms = 10;
r2opts.randFlag = 1;
r2opts.validMeans = 1;
r2opts.kval = 3;
r2opts.nShuffle = 100;
%get the tuning strength and significance
for thisROI = 1:numel(units)
    %stationary trials
    gca = units(thisROI).allSpikesDownsample(:,1)';
    [units(thisROI).statR2, units(thisROI).statR2_pval] = calc_kfold_R2(gca, r2opts.kval, r2opts.nPerms,...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);
    %locomotion trials
    gca = units(thisROI).allSpikesDownsample(:,2)';
    [units(thisROI).runR2, units(thisROI).runR2_pval] = calc_kfold_R2(gca, r2opts.kval, r2opts.nPerms,...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);
end
%% CGF tuning strength (R^2)
figure('Color', 'w', 'Name', 'Population Tuning Quality', 'Position', [200, 200, 500, 400]);
hold on;

%  R2 values across all units
allStatR2 = [units.statR2];
allRunR2  = [units.runR2];
hStat = cdfplot(allStatR2); 
set(hStat, 'Color', 'k', 'LineWidth', 2.5, 'DisplayName', 'Stationary');
hRun = cdfplot(allRunR2); 
set(hRun, 'Color', 'r', 'LineWidth', 2.5, 'DisplayName', 'Locomotion');
title('Cross-validated R^2 distribution', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Cross-validated R^2 score', 'FontSize', 11);
ylabel('Proportion of ROIs', 'FontSize', 11);
legend('Location', 'southeast');
hold off;

%% get fano factor and dynamic range
for thisROI = 1:numel(units)
    istate = 1;
    units(thisROI).dynamicRange_stat = range(units(thisROI).tuning(:,istate))/options.binSpacing;
    units(thisROI).fanoFactor_stat = mean(cellfun(@(x) var(x)/mean(x), units(thisROI).allSpikes(:,istate)),1);
    istate = 2;
    units(thisROI).dynamicRange_run = range(units(thisROI).tuning(:,istate))/options.binSpacing;
    units(thisROI).fanoFactor_run = mean(cellfun(@(x) var(x)/mean(x), units(thisROI).allSpikes(:,istate)),1);
end
%% fit descriptive functions
%clustering analysis reveals 4 main classes of tuning function
%gaussian fit and preferred speed with all trials
%stationary trials
istate = 1;
for thisROI = 1:numel(units)
    %fit gaussians and find best fit
    [units(thisROI).gaussParams_stat(1,:), units(thisROI).gaussChar_stat, units(thisROI).gaussR2_stat] = ...
        fitGaussianTemplates_tuning(units(thisROI).allSpikes(:,istate),0.5,false);
    %get preferred stimulus
    if units(thisROI).gaussChar_stat == 4 % if inverted, pref speed is min fr.
        [~, units(thisROI).prefSpeed_stat] = min(units(thisROI).tuning(:,istate));
    else % max fr
        [~, units(thisROI).prefSpeed_stat] = max(units(thisROI).tuning(:,istate));
    end
end
%locomotion trials
istate = 2;
for thisROI = 1:numel(units)
    %fit gaussians and find best fit
    [units(thisROI).gaussParams_run(1,:), units(thisROI).gaussChar_run, units(thisROI).gaussR2_run] = ...
        fitGaussianTemplates_tuning(units(thisROI).allSpikes(:,istate),0.5,false);
    %get preferred stimulus
    if units(thisROI).gaussChar_run == 4 % if inverted, pref speed is min fr.
        [~, units(thisROI).prefSpeed_run] = min(units(thisROI).tuning(:,istate));
    else % max fr
        [~, units(thisROI).prefSpeed_run] = max(units(thisROI).tuning(:,istate));
    end
end
%% Mutual information analysis
%MI options
optimiseBins = false; nMCSamples = 1000; sigFlag = false; correction = 'MLE'; nBinLim = 6; % nspeeds
for thisROI = 1:numel(units)
    %stationary trials
    gca = units(thisROI).allSpikesDownsample(:,1);
    [units(thisROI).MI_stat, ~, ~, units(thisROI).SSI_stat] =...
        calcMI(gca, correction, optimiseBins, sigFlag, nMCSamples, nBinLim);
    %locomotion trials
    gca = units(thisROI).allSpikesDownsample(:,2);
    [units(thisROI).MI_run, ~, ~, units(thisROI).SSI_run] =...
        calcMI(gca, correction, optimiseBins, sigFlag, nMCSamples, nBinLim);
end
%% plots
r2_thresh = 0.1;
r2p_thresh = 0.05;
validIdx_stat = find(cat(1,units.statR2)>r2_thresh & cat(1,units.statR2_pval)<r2p_thresh);
validIdx_run = find(cat(1,units.runR2)>r2_thresh & cat(1,units.runR2_pval)<r2p_thresh);
validIdx_both = validIdx_stat(ismember(validIdx_stat, validIdx_run));
allTuning = cat(3,units.tuning)/options.binSpacing; % convert to spikes/s
% plot average tuning
figure
if ~isempty(validIdx_stat)
    shadedErrorBar(1:6, mean(allTuning(:,1,validIdx_stat),3),...
        std(allTuning(:,1,validIdx_stat),0,3)./sqrt(numel(validIdx_stat)),'lineProps','k');
end
hold on
if ~isempty(validIdx_run)
    shadedErrorBar(1:6, mean(allTuning(:,2,validIdx_run),3),...
        std(allTuning(:,1,validIdx_run),0,3)./sqrt(numel(validIdx_run)),'lineProps','r');
end
xlabel('Stimulus')
ylabel('Firing Rate (hz)')
title('Mean tuning curves')
figure
plot([units.statR2], [units.runR2],'k.')
hold on
plot([0.1 0.1], [-1 1],'r')
plot([-1, 1], [0.1 0.1],'r')
xlim([-0.6 1])
ylim([-0.6 1])
xlabel('Tuning Strength (stationary)')
ylabel('Tuning Strength (locomotion)')
% MI plots
figure
plot([units.MI_stat], [units.MI_run],'k.')
axis equal
hold on
plot([0 1.5], [0 1.5], 'r')
xlabel('Mutual information (stationary)')
ylabel('Mutual information (locomotion)')
figure
shadedErrorBar(1:6, mean(cat(1,units.SSI_stat),1), std(cat(1,units.SSI_stat),0,1)./sqrt(size(cat(1,units.SSI_stat),1)))
hold on
shadedErrorBar(1:6, mean(cat(1,units.SSI_run),1), std(cat(1,units.SSI_run),0,1)./sqrt(size(cat(1,units.SSI_run),1)),'lineProps','r')
ylabel('SSI')
xlabel('Stimulus')

% preferred speeds joint histogram
if ~isempty(validIdx_both)
    allPrefs = [cat(1,units(validIdx_both).prefSpeed_stat), cat(1,units(validIdx_both).prefSpeed_run)];
    vals = histcounts2(allPrefs(:,1), allPrefs(:,2), 1:7, 1:7);
    vals = vals./sum(vals,2);
    figure
    imagesc(vals')
    hold on
    plot([0.5 6.5],[0.5 6.5],'r')
    axis xy
    colorbar
    xlabel('Preferred stimulus (stationary)')
    ylabel('Preferred stimulus (locomotion)')
end

% plot gauss fit params for tuned bandpass cells
idx = find(cat(1,units.runR2)>r2_thresh & cat(1,units.runR2_pval)<r2p_thresh...
    & cat(1,units.statR2)>r2_thresh & cat(1,units.statR2_pval)<r2p_thresh...
    & cat(1,units.gaussChar_stat)==3 & cat(1,units.gaussChar_run)==3);
if ~isempty(idx)
    allStatParams = cat(1,units(idx).gaussParams_stat);
    allRunParams = cat(1,units(idx).gaussParams_run);
    minVals = min(cat(1,allStatParams,allRunParams));
    maxVals = max(cat(1,allStatParams,allRunParams));
    paramNames = {'baseline','amplitude','mu','sigma'};
    figure
    for iparam = 1:4
        subplot(2,2,iparam)
        plot(allStatParams(:,iparam), allRunParams(:,iparam),'k.')
        hold on
        plot([minVals(iparam), maxVals(iparam)],[minVals(iparam), maxVals(iparam)],'r')
        title(paramNames{iparam})
    end
end
% plot dynamic range and fano factor
if ~isempty(validIdx_both)
    figure
    plot([units(validIdx_both).dynamicRange_stat], [units(validIdx_both).dynamicRange_run],'k.')
    xlabel('Dynamic range (stationary)')
    ylabel('Dynamic range (locomotion)')
    figure
    plot([units(validIdx_both).fanoFactor_stat], [units(validIdx_both).fanoFactor_run],'k.')
    xlabel('Fano Factor (stationary)')
    ylabel('Fano Factor (locomotion)')
end
%% plot example tuning curves
if exist('thal_idx', 'var')
    for i = 1:numel(thal_idx)
        thisROI = thal_idx(i);
        shadedErrorBar(1:6,...
            cellfun(@(x) mean(x,2), units(thisROI).allSpikes(:,1)),...
            cellfun(@(x) (std(x,0,2)/sqrt(numel(x))), units(thisROI).allSpikes(:,1)))
            shadedErrorBar(1:6,...
            cellfun(@(x) mean(x,2), units(thisROI).allSpikes(:,2)),...
            cellfun(@(x) (std(x,0,2)/sqrt(numel(x))), units(thisROI).allSpikes(:,2)),...
            'lineProps','r')
        pause
        close
    end
end



% metrics from the descriptive function 
%% 
if ~exist('units', 'var')
    error('The variable "units" does not exist in the workspace.');
end

% Ensure your cross-validation filters from the previous step are ready
r2_thresh = 0.1;
r2p_thresh = 0.05;
validIdx_stat = find(cat(1,units.statR2) > r2_thresh & cat(1,units.statR2_pval) < r2p_thresh);
validIdx_run  = find(cat(1,units.runR2) > r2_thresh & cat(1,units.runR2_pval) < r2p_thresh);
validIdx_both = validIdx_stat(ismember(validIdx_stat, validIdx_run));

if isempty(validIdx_both)
    warning('No units passed the dual R^2 filter. Plotting all units instead.');
    validIdx_both = 1:numel(units);
end

% Pull the parameters for the filtered cohort
allStatParams = cat(1, units(validIdx_both).gaussParams_stat);
allRunParams  = cat(1, units(validIdx_both).gaussParams_run);

% Label arrays
paramNames = {'Baseline', 'Response amplitude', '\mu (preferred speed location)', '\sigma (tuning width)'};
paramColors = {
    [0, 0.4470, 0.7410],      % Blue
    [0.8500, 0.3250, 0.0980], % Orange
    [0.9290, 0.6940, 0.1250], % Yellow
    [0.4940, 0.1840, 0.5560]  % Purple
};

%% Render 4-Panel Parameter Scatter Grid
figure('Color', 'w', 'Position', [100, 100, 900, 750]);

for iparam = 1:4
    subplot(2, 2, iparam);
    hold on;
    
    x_stat = allStatParams(:, iparam);
    y_run  = allRunParams(:, iparam);
    
    if iparam == 1 || iparam == 2
        x_stat = x_stat / options.binSpacing;
        y_run = y_run / options.binSpacing;
    end
    
    scatter(x_stat, y_run, 45, 'MarkerFaceColor', paramColors{iparam}, ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
    
    % 
    minVal = min([x_stat; y_run]);
    maxVal = max([x_stat; y_run]);
    if minVal == maxVal; minVal = minVal - 0.1; maxVal = maxVal + 0.1; end
    plot([minVal, maxVal], [minVal, maxVal], 'k--', 'LineWidth', 1.2);
    
    % Run a quick statistical test per metric
    p_val = signrank(x_stat, y_run);
    
    % Clean layout
    xlim([minVal, maxVal]); ylim([minVal, maxVal]);
    xlabel('Stationary Value', 'FontSize', 10);
    ylabel('Locomotion Value', 'FontSize', 10);
    title(sprintf('%s', paramNames{iparam}), 'FontSize', 11, 'FontWeight', 'bold');
    hold off;
end

sgtitle('descriptive parameters across behavioral states (filtered boutons)', 'FontSize', 13, 'FontWeight', 'bold');

%% Figure 2: Individual Gaussian Fits 
if ~exist('units', 'var')
    error('The variable "units" does not exist in the workspace.');
end

% highly reliable boutons only
r2_thresh = 0.1;
r2p_thresh = 0.05;

validIdx_stat = find(cat(1, units.statR2) > r2_thresh & cat(1, units.statR2_pval) < r2p_thresh);
validIdx_run  = find(cat(1, units.runR2) > r2_thresh & cat(1, units.runR2_pval) < r2p_thresh);
validIdx_both = validIdx_stat(ismember(validIdx_stat, validIdx_run));

if isempty(validIdx_both)
    error('No boutons passed your cross-validation thresholds. Cannot plot fits.');
end


fprintf('Found %d cross-validated boutons with high-quality fits.\n', numel(validIdx_both));

% We will plot up to the first 6 cross-validated units to see their fit quality
maxPlots = min(6, numel(validIdx_both));
xDense = linspace(1, 6, 100);

figure('Color', 'w', 'Name', 'Cross-Validated Gaussian Fit Overlays', 'Position', [100, 50, 1100, 800]);

for iPlot = 1:maxPlots
    thisROI = validIdx_both(iPlot);
    subplot(2, 3, iPlot); 
    hold on;
    
    % Extract raw binned data points (normalized by binSpacing)
    y_data_stat = units(thisROI).tuning(:, 1) / options.binSpacing;
    y_data_run  = units(thisROI).tuning(:, 2) / options.binSpacing;
    
    % Reconstruct Analytical Stationary Fit Profile
    p_stat = units(thisROI).gaussParams_stat;
    y_fit_stat = p_stat(1) + p_stat(2) * exp(-((xDense - p_stat(3)).^2) / (2 * p_stat(4)^2));
    y_fit_stat = y_fit_stat / options.binSpacing;
    
    % Reconstruct Analytical Locomotion Fit Profile
    p_run = units(thisROI).gaussParams_run;
    y_fit_run = p_run(1) + p_run(2) * exp(-((xDense - p_run(3)).^2) / (2 * p_run(4)^2));
    y_fit_run = y_fit_run / options.binSpacing;
    
   
    % Stationary (Black)
    plot(1:6, y_data_stat, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5.5);
    plot(xDense, y_fit_stat, 'k-', 'LineWidth', 2, 'DisplayName', 'Stationary');
    
    %  Locomotion 
    plot(1:6, y_data_run, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5.5);
    plot(xDense, y_fit_run, 'r-', 'LineWidth', 2, 'DisplayName', 'Locomotion');
    
    % 
    allVals = [y_data_stat; y_data_run; y_fit_stat'; y_fit_run'];
    ylim([min(allVals) - 0.03, max(allVals) + 0.03]);
    xlim([0.5, 6.5]);
    xticks(1:6);
    
    ylabel('\DeltaF/F ');
    xlabel('Visual Speed Index');
    title(sprintf('Bouton ID: %d\nStat Type %d | Run Type %d', ...
        thisROI, units(thisROI).gaussChar_stat, units(thisROI).gaussChar_run), 'FontSize', 10);
    
    if iPlot == 1, legend('Location', 'best'); end
end

sgtitle('Verified Gaussian Model Fits Across States (Cross-Validated Cohort)', 'FontSize', 13, 'FontWeight', 'bold');


sessionNames = keys(sessionMinTrial);
sessionVals  = cell2mat(values(sessionMinTrial));
figure('Color','w');
bar(sessionVals);
set(gca, 'XTick', 1:numel(sessionNames), 'XTickLabel', sessionNames, 'XTickLabelRotation', 45);
ylabel('minTrial (this session)');
title('Per-session trial-count floor used for R^2 downsampling');