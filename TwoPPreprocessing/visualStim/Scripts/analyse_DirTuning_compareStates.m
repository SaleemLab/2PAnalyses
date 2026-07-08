% DirTuning_Pooled_BehaviorSplit.m 
% TODO: (compare 
%
% Pools DirTuning boutons across MULTIPLE mice/sessions, splitting trials
% into RUNNING vs STATIONARY using the same wheel-speed classification
% logic as DotFields_TuningCurveAnalysis_compareStates_2PData.m:
%   propRunning    = fraction of stim-window samples with speed > 0.5
%   propStationary = fraction of stim-window samples with speed < 3
%   runFlag = 1 if propRunning >= 0.75 AND meanSpeed > 3
%   runFlag = 0 if propStationary >= 0.75 AND meanSpeed < 0.5
%   runFlag = NaN (trial excluded) otherwise
%
% ASSUMES response.wheelData(d) exists per direction group, with the
% same structure as in the DotFields script:
%   response.wheelData(d).alignedResponses : [nTimepoints x nTrials] wheel speed
%   response.wheelData(d).timeVector       : [nTimepoints x 1]
% VERIFY this field exists for your DirTuning response.mat files before
% trusting results -- if it doesn't, this script will warn and skip
% behavior classification for that session (all trials excluded).
%
% For each bouton, computes ANOVA+threshold and cross-validated R^2
% SEPARATELY for running trials and stationary trials (mirroring the
% statR2/runR2 pattern in the original DotFields script), then compares
% stationary vs running results.
%
% Requires calc_kfold_R2.m on your MATLAB path.

%% 
mouseList = {'M25132','M25133', 'M26003'};  

baseWin      = [-0.5 0];
respWin_2sOn = [0.2 3];
respWin_1sOn = [0.2 1];   % see note in DirTuning_Pooled_AllTrials.m -- verify this matches your protocol

stimWindowForBehavior = [0 2.0]; % window used to classify running vs stationary (matches DotFields script's 2s eval window)
runSpeedThresh  = 3;    % cm/s (or your wheel's units) -- mean speed above this = candidate "running"
statSpeedThresh = 0.5;  % mean speed below this = candidate "stationary"
propThresh      = 0.75; % fraction of the window that must meet the speed criterion

ALPHA = 0.05;
NSD   = 1;

r2opts.kval       = 3;
r2opts.nPerms     = 10;
r2opts.randFlag   = 1;
r2opts.validMeans = 1;
r2opts.nShuffle   = 100;
%% 

DirTuningTable = filterMasterTable_usingNameSessionPairs('MouseID', mouseList, ...
    'Exclude', 0, 'HasStimulus', {'DirTuning'});
allMice    = DirTuningTable.MouseID;
uniqueMice = unique(allMice, 'stable');

allDirTuning = struct('trialMeanResp_run', {}, 'trialMeanResp_stat', {}, ...
                      'meanDirResponse_run', {}, 'meanDirResponse_stat', {}, ...
                      'stimVariant', {}, 'sessionLabel', {});
haveWarnedShape = false;
haveWarnedWheel = false;

for iMouse = 1:length(uniqueMice)
    thisMouse    = uniqueMice{iMouse};
    mouseSessIdx = find(strcmp(allMice, thisMouse));
    fprintf('MOUSE: %s | %d sessions\n', thisMouse, length(mouseSessIdx));

    for iSess = 1:length(mouseSessIdx)
        tableRow        = DirTuningTable(mouseSessIdx(iSess), :);
        thisSessionName = char(tableRow.Session);
        fprintf('  --- Session: %s ---\n', thisSessionName);

        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath), warning('    sfi missing -- skipping.'); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        DirTuneIdx = find(contains(stimNames, 'DirTuning'));
        if isempty(DirTuneIdx), warning('    No DirTuning files -- skipping.'); continue; end

        try
            load(sessionFileInfo.stimFiles(DirTuneIdx(1)).Response, 'response');
        catch ME
            warning('    Could not load response: %s', ME.message); continue;
        end
        if ~isfield(response, 'psthData') || isempty(response.psthData)
            warning('    psthData empty/missing -- skipping.'); continue;
        end
        if ~isfield(response, 'wheelData') || isempty(response.wheelData)
            if ~haveWarnedWheel
                warning('    No response.wheelData found -- cannot classify behavior for %s. Skipping session.', thisSessionName);
                haveWarnedWheel = true;
            end
            continue;
        end

        switch response.postStimTime
            case 4, respWin = respWin_2sOn; stimVariant = 4;
            case 3, respWin = respWin_1sOn; stimVariant = 3;
            otherwise, warning('    Unexpected postStimTime -- skipping.'); continue;
        end

        nDir      = numel(response.psthData);
        timeVec   = response.psthData(1).timeVector(:)';
        baseIdx   = timeVec >= baseWin(1) & timeVec <= baseWin(2);
        respIdx   = timeVec >= respWin(1) & timeVec <= respWin(2);
        stimVals  = arrayfun(@(x) x.stimValue, response.psthData);

        nBoutons = size(response.psthData(1).alignedResponses, 1);
        if ~haveWarnedShape
            fprintf('    Detected shape: %d boutons x %d timepoints x %d trials (dir 1)\n', ...
                size(response.psthData(1).alignedResponses,1), ...
                size(response.psthData(1).alignedResponses,2), ...
                size(response.psthData(1).alignedResponses,3));
            haveWarnedShape = true;
        end

        %% classify each trial, each direction, as running / stationary / excluded
        runFlagPerDir = cell(nDir, 1); % {nDir}, each [nTrials x 1]: 1=run, 0=stat, NaN=excluded
        for thisDir = 1:nDir
            wheelTVec  = response.wheelData(thisDir).timeVector(:)';
            wheelMask  = wheelTVec >= stimWindowForBehavior(1) & wheelTVec <= stimWindowForBehavior(2);
            speedMat   = response.wheelData(thisDir).alignedResponses; % [nTimepoints x nTrials]
            nTrialsDir = size(speedMat, 2);

            rf = nan(nTrialsDir, 1);
            for ti = 1:nTrialsDir
                trace = speedMat(wheelMask, ti);
                if all(isnan(trace)), continue; end
                meanSpeed      = mean(trace, 'omitnan');
                propRunning    = sum(trace > statSpeedThresh) / numel(trace);
                propStationary = sum(trace < runSpeedThresh)  / numel(trace);
                if propRunning >= propThresh && meanSpeed > runSpeedThresh
                    rf(ti) = 1;
                elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                    rf(ti) = 0;
                end
            end
            runFlagPerDir{thisDir} = rf;
        end

        %% per-bouton extraction, split by behavior state
        for iBouton = 1:nBoutons
            trialMeanResp_run  = cell(nDir, 1);
            trialMeanResp_stat = cell(nDir, 1);
            meanDirResponse_run  = nan(nDir, 1);
            meanDirResponse_stat = nan(nDir, 1);

            for thisDir = 1:nDir
                traceMat = squeeze(response.psthData(thisDir).alignedResponses(iBouton, :, :));
                if isvector(traceMat), traceMat = traceMat(:); end
                traceMat = double(traceMat)'; % [nTrials x nTimepoints]

                perTrialBaseline = mean(traceMat(:, baseIdx), 2, 'omitnan');
                traceMatSub = traceMat - perTrialBaseline;
                trialResp = mean(traceMatSub(:, respIdx), 2, 'omitnan'); % [nTrials x 1]

                rf = runFlagPerDir{thisDir};
                nTrialsHere = min(numel(rf), numel(trialResp)); % safety, should match

                runMask  = rf(1:nTrialsHere) == 1;
                statMask = rf(1:nTrialsHere) == 0;

                trialMeanResp_run{thisDir}  = trialResp(runMask);
                trialMeanResp_stat{thisDir} = trialResp(statMask);
                meanDirResponse_run(thisDir)  = mean(trialResp(runMask), 'omitnan');
                meanDirResponse_stat(thisDir) = mean(trialResp(statMask), 'omitnan');
            end

            allDirTuning(end+1) = struct( ...
                'trialMeanResp_run',    {trialMeanResp_run}, ...
                'trialMeanResp_stat',   {trialMeanResp_stat}, ...
                'meanDirResponse_run',  meanDirResponse_run, ...
                'meanDirResponse_stat', meanDirResponse_stat, ...
                'stimVariant',          stimVariant, ...
                'sessionLabel',         sprintf('%s_%s', thisMouse, thisSessionName));
        end
        fprintf('    Added %d boutons (running total: %d).\n', nBoutons, numel(allDirTuning));
    end
end

nBoutonsTotal = numel(allDirTuning);
fprintf('\nDone pooling. %d boutons total across %d mice.\n', nBoutonsTotal, numel(uniqueMice));

%%  ANOVA + cross-val R^2, separately per state 
states = {'run', 'stat'};
results = struct();
for si = 1:2
    st = states{si};
    fieldResp = sprintf('trialMeanResp_%s', st);
    fieldMean = sprintf('meanDirResponse_%s', st);

    anovaP    = nan(nBoutonsTotal, 1);
    cvR2      = nan(nBoutonsTotal, 1);
    cvPval    = nan(nBoutonsTotal, 1);

    for b = 1:nBoutonsTotal
        s = allDirTuning(b);
        nDir = numel(s.(fieldResp));

        y = []; grp = [];
        for thisDir = 1:nDir
            vals = s.(fieldResp){thisDir}(:);
            vals = vals(~isnan(vals));
            y   = [y; vals];
            grp = [grp; repmat(thisDir, numel(vals), 1)];
        end
        if numel(y) >= nDir * 2
            anovaP(b) = anova1(y, grp, 'off');
        end

        gca_ = cell(1, nDir);
        for thisDir = 1:nDir
            vals = s.(fieldResp){thisDir}(:)';
            vals = vals(~isnan(vals));
            gca_{thisDir} = vals;
        end
        trialCounts = cellfun(@numel, gca_);
        if all(trialCounts > 0)
            minTrial = min(trialCounts);
            if minTrial >= r2opts.kval % need at least kFold trials per direction
                gcaDownsampled = cellfun(@(x) x(1:minTrial), gca_, 'UniformOutput', false);
                [cvR2(b), cvPval(b)] = calc_kfold_R2(gcaDownsampled, r2opts.kval, r2opts.nPerms, ...
                    r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);
            end
        end

        if mod(b, 100) == 0
            fprintf('[%s] Processed %d / %d boutons...\n', st, b, nBoutonsTotal);
        end
    end

    results.(st).anovaP = anovaP;
    results.(st).cvR2   = cvR2;
    results.(st).cvPval = cvPval;
    results.(st).isTunedCV = cvPval < ALPHA;

    fprintf('\n[%s] %d / %d boutons significant cross-validated tuning (p<%.2f).\n', ...
        st, sum(results.(st).isTunedCV, 'omitnan'), nBoutonsTotal, ALPHA);
end

for b = 1:nBoutonsTotal
    allDirTuning(b).anovaP_run    = results.run.anovaP(b);
    allDirTuning(b).cvR2_run      = results.run.cvR2(b);
    allDirTuning(b).cvPval_run    = results.run.cvPval(b);
    allDirTuning(b).isTunedCV_run = results.run.isTunedCV(b);

    allDirTuning(b).anovaP_stat    = results.stat.anovaP(b);
    allDirTuning(b).cvR2_stat      = results.stat.cvR2(b);
    allDirTuning(b).cvPval_stat    = results.stat.cvPval(b);
    allDirTuning(b).isTunedCV_stat = results.stat.isTunedCV(b);
end

%% 
validBoth = ~isnan(results.run.cvR2) & ~isnan(results.stat.cvR2);
fprintf('\n=== Running vs Stationary cross-validated R^2 (n=%d boutons with both states valid) ===\n', sum(validBoth));

figure('Position', [100 100 900 400]);
subplot(1,2,1);
scatter(results.stat.cvR2(validBoth), results.run.cvR2(validBoth), 15, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
minVal = min([results.stat.cvR2(validBoth); results.run.cvR2(validBoth)]);
maxVal = max([results.stat.cvR2(validBoth); results.run.cvR2(validBoth)]);
plot([minVal maxVal], [minVal maxVal], 'k--');
xlabel('Cross-val R^2 (stationary)'); ylabel('Cross-val R^2 (running)');
title('Tuning strength: stationary vs running');

subplot(1,2,2);
nBothTuned  = sum(results.run.isTunedCV(validBoth) & results.stat.isTunedCV(validBoth));
nRunOnly    = sum(results.run.isTunedCV(validBoth) & ~results.stat.isTunedCV(validBoth));
nStatOnly   = sum(~results.run.isTunedCV(validBoth) & results.stat.isTunedCV(validBoth));
nNeither    = sum(~results.run.isTunedCV(validBoth) & ~results.stat.isTunedCV(validBoth));
bar([nBothTuned, nRunOnly, nStatOnly, nNeither]);
xticklabels({'Both states', 'Running only', 'Stationary only', 'Neither'});
ylabel('Number of boutons');
title('Tuning significance by behavioral state');

sgtitle('DirTuning: running vs stationary cross-validated tuning');

fprintf('Both states tuned: %d | Running only: %d | Stationary only: %d | Neither: %d\n', ...
    nBothTuned, nRunOnly, nStatOnly, nNeither);