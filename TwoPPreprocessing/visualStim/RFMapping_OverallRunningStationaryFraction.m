% RFMapping_OverallRunningStationaryFraction.m
%
% Computes the fraction of running vs stationary vs ambiguous trials,
% separately for GRID trials and BLANK trials, pooled across ALL RF
% mapping sessions -- and separately, restricted to just the sessions
% that contributed at least one bouton to your 1SD-responsive
% population. Running/stationary state is a per-trial, session-level
% property (the same for every bouton in a session), so this reports at
% the session level, not per-bouton.
%
% Requires: sessionLabels (cell array from the main pooling loop,
% tracking which session each bouton in allRFMapping came from),
% zScore/anovaSig (or candidateRespIdxList) already computed for the
% responsive-population comparison.

respWin = [0.1 3]; 
runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;
NSD_final = 1;
ALPHA = 0.05;
stimFramesMask_range = [-0.2 2.8];  
isResponsive_1SD = ([allRFMapping.pValANOVA] < ALPHA) & (zScore(:)' > NSD_final);
respIdxList_1SD = find(isResponsive_1SD);
fprintf('%d / %d boutons responsive at 1SD criterion.\n', numel(respIdxList_1SD), numel(allRFMapping));
respSessions = unique(sessionLabels(respIdxList_1SD), 'stable');
allSessions  = unique(sessionLabels, 'stable');
fprintf('%d / %d sessions contain at least one 1SD-responsive bouton.\n', numel(respSessions), numel(allSessions));
sessionStats = struct('sessionLabel', {}, 'nGridRun', {}, 'nGridStat', {}, 'nGridAmbig', {}, ...
    'nBlankRun', {}, 'nBlankStat', {}, 'nBlankAmbig', {});
for s = 1:numel(allSessions)
    thisLabel = allSessions{s};
    us = strsplit(thisLabel, '_');
    thisMouse = us{1};
    thisSessionName = strjoin(us(2:end), '_'); 
    infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
    if ~isfile(infoPath), warning('sfi missing for %s -- skipping.', thisLabel); continue; end
    loadedInfo      = load(infoPath, 'sessionFileInfo');
    sessionFileInfo = loadedInfo.sessionFileInfo;
    stimNamesHere   = {sessionFileInfo.stimFiles.name};
    iStim = find(contains(stimNamesHere, 'RFMapping'), 1);
    if isempty(iStim), warning('No RFMapping file for %s -- skipping.', thisLabel); continue; end
    try
        load(sessionFileInfo.stimFiles(iStim).Response, 'response');
    catch ME
        warning('Could not load response for %s: %s', thisLabel, ME.message); continue;
    end
    if ~isfield(response, 'wheelData') || numel(response.wheelData) ~= numel(response.psthData)
        warning('wheelData missing/mismatched for %s -- skipping.', thisLabel); continue;
    end
    psthDataHere = response.psthData;
    stimVsHere   = vertcat(psthDataHere.stimValue);
    blankIdxHere = find(stimVsHere(:,1) == 200 & stimVsHere(:,2) == 0, 1);
    gridIdxHere  = find(stimVsHere(:,1) ~= 200);
    wheelTimeVecHere = response.wheelData(1).timeVector(:)';
    wheelRespIdxHere = wheelTimeVecHere >= stimFramesMask_range(1) & wheelTimeVecHere <= stimFramesMask_range(2);
    nGridRun = 0; nGridStat = 0; nGridAmbig = 0;
    for g = gridIdxHere(:)'
        targetAz = stimVsHere(g,1);
        targetEl = stimVsHere(g,2);
        wG = [];
        for w = 1:numel(response.wheelData)
            wVal = response.wheelData(w).stimValue;
            if numel(wVal) >= 2 && wVal(1) == targetAz && wVal(2) == targetEl
                wG = w;
                break;
            end
        end
        if isempty(wG)
            continue;
        end
        wheelTrials = response.wheelData(wG).alignedResponses;
        for ti = 1:size(wheelTrials, 2)
            trace = wheelTrials(wheelRespIdxHere, ti);
            if all(isnan(trace)), continue; end
            meanSpeed      = nanmean(trace);
            propRunning    = sum(trace > statSpeedThresh) / sum(wheelRespIdxHere);
            propStationary = sum(trace < runSpeedThresh)  / sum(wheelRespIdxHere);
            if propRunning >= propThresh && meanSpeed > runSpeedThresh
                nGridRun = nGridRun + 1;
            elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                nGridStat = nGridStat + 1;
            else
                nGridAmbig = nGridAmbig + 1;
            end
        end
    end
    nBlankRun = 0; nBlankStat = 0; nBlankAmbig = 0;
    if ~isempty(blankIdxHere)
        targetAz = stimVsHere(blankIdxHere,1);
        targetEl = stimVsHere(blankIdxHere,2);
        wG_blank = [];
        for w = 1:numel(response.wheelData)
            wVal = response.wheelData(w).stimValue;
            if numel(wVal) >= 2 && wVal(1) == targetAz && wVal(2) == targetEl
                wG_blank = w;
                break;
            end
        end
        if ~isempty(wG_blank)
            wheelTrials = response.wheelData(wG_blank).alignedResponses;
            for ti = 1:size(wheelTrials, 2)
                trace = wheelTrials(wheelRespIdxHere, ti);
                if all(isnan(trace)), continue; end
                meanSpeed      = nanmean(trace);
                propRunning    = sum(trace > statSpeedThresh) / sum(wheelRespIdxHere);
                propStationary = sum(trace < runSpeedThresh)  / sum(wheelRespIdxHere);
                if propRunning >= propThresh && meanSpeed > runSpeedThresh
                    nBlankRun = nBlankRun + 1;
                elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                    nBlankStat = nBlankStat + 1;
                else
                    nBlankAmbig = nBlankAmbig + 1;
                end
            end
        end
    end
    sessionStats(end+1) = struct('sessionLabel', thisLabel, 'nGridRun', nGridRun, 'nGridStat', nGridStat, ...
        'nGridAmbig', nGridAmbig, 'nBlankRun', nBlankRun, 'nBlankStat', nBlankStat, 'nBlankAmbig', nBlankAmbig); 
    fprintf('%s: grid run=%d stat=%d ambig=%d | blank run=%d stat=%d ambig=%d\n', ...
        thisLabel, nGridRun, nGridStat, nGridAmbig, nBlankRun, nBlankStat, nBlankAmbig);
end
totalGridRun   = sum([sessionStats.nGridRun]);
totalGridStat  = sum([sessionStats.nGridStat]);
totalGridAmbig = sum([sessionStats.nGridAmbig]);
totalGridAll   = totalGridRun + totalGridStat + totalGridAmbig;
totalBlankRun   = sum([sessionStats.nBlankRun]);
totalBlankStat  = sum([sessionStats.nBlankStat]);
totalBlankAmbig = sum([sessionStats.nBlankAmbig]);
totalBlankAll   = totalBlankRun + totalBlankStat + totalBlankAmbig;
fprintf('\n=== OVERALL (all %d sessions) ===\n', numel(sessionStats));
fprintf('Grid trials:  running=%.1f%% | stationary=%.1f%% | ambiguous=%.1f%% (n=%d)\n', ...
    100*totalGridRun/totalGridAll, 100*totalGridStat/totalGridAll, 100*totalGridAmbig/totalGridAll, totalGridAll);
fprintf('Blank trials: running=%.1f%% | stationary=%.1f%% | ambiguous=%.1f%% (n=%d)\n', ...
    100*totalBlankRun/totalBlankAll, 100*totalBlankStat/totalBlankAll, 100*totalBlankAmbig/totalBlankAll, totalBlankAll);
respMask = ismember({sessionStats.sessionLabel}, respSessions);
respStats = sessionStats(respMask);
totalGridRun_r   = sum([respStats.nGridRun]);
totalGridStat_r  = sum([respStats.nGridStat]);
totalGridAmbig_r = sum([respStats.nGridAmbig]);
totalGridAll_r   = totalGridRun_r + totalGridStat_r + totalGridAmbig_r;
totalBlankRun_r   = sum([respStats.nBlankRun]);
totalBlankStat_r  = sum([respStats.nBlankStat]);
totalBlankAmbig_r = sum([respStats.nBlankAmbig]);
totalBlankAll_r   = totalBlankRun_r + totalBlankStat_r + totalBlankAmbig_r;
fprintf('\n=== RESTRICTED to %d sessions containing 1SD-responsive boutons ===\n', numel(respStats));
fprintf('Grid trials:  running=%.1f%% | stationary=%.1f%% | ambiguous=%.1f%% (n=%d)\n', ...
    100*totalGridRun_r/totalGridAll_r, 100*totalGridStat_r/totalGridAll_r, 100*totalGridAmbig_r/totalGridAll_r, totalGridAll_r);
fprintf('Blank trials: running=%.1f%% | stationary=%.1f%% | ambiguous=%.1f%% (n=%d)\n', ...
    100*totalBlankRun_r/totalBlankAll_r, 100*totalBlankStat_r/totalBlankAll_r, 100*totalBlankAmbig_r/totalBlankAll_r, totalBlankAll_r);
figure('Position', [100 100 700 400]);
subplot(1,2,1);
bar([100*totalGridRun/totalGridAll, 100*totalGridStat/totalGridAll, 100*totalGridAmbig/totalGridAll; ...
     100*totalGridRun_r/totalGridAll_r, 100*totalGridStat_r/totalGridAll_r, 100*totalGridAmbig_r/totalGridAll_r]);
xticklabels({'All sessions', 'Responsive-only sessions'});
legend({'Running', 'Stationary', 'Ambiguous'}, 'Location', 'best');
ylabel('% of grid trials');
title('Grid trials by state');
subplot(1,2,2);
bar([100*totalBlankRun/totalBlankAll, 100*totalBlankStat/totalBlankAll, 100*totalBlankAmbig/totalBlankAll; ...
     100*totalBlankRun_r/totalBlankAll_r, 100*totalBlankStat_r/totalBlankAll_r, 100*totalBlankAmbig_r/totalBlankAll_r]);
xticklabels({'All sessions', 'Responsive-only sessions'});
legend({'Running', 'Stationary', 'Ambiguous'}, 'Location', 'best');
ylabel('% of blank trials');
title('Blank trials by state');
sgtitle('Trial composition by behavioral state: all sessions vs responsive-only sessions');