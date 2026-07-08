% RFMapping_CheckPositionCoverageByState.m
%
% For each hand-picked session, reports how many stationary vs running
% trials exist AT EACH of the 16 grid positions (+ blank). This is a
% SESSION-level property, not bouton-specific -- trial-to-position
% assignment and running/stationary classification are both properties
% of the trial itself, the same for every ROI in that session.
%
% This directly checks the concern that splitting by state can leave
% some grid positions with zero (or very few) trials in one state,
% which would create holes in that state's spatial map and make its
% ANOVA a different, less comparable test than the other state's.
%
% Uses the SAME session list and classification thresholds as
% RFMapping_BehaviorSplit.m -- keep in sync if you change either.

%% ===================== CONFIG (keep in sync with RFMapping_BehaviorSplit.m) =====================
sessionPairs = struct();
sessionPairs.M25132 = {'20260228', '20260306'};
sessionPairs.M26003 = {'20260316', '20260324'};

stimName = 'RFMapping';
respWin  = [0.5 3];

runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;

minTrialWarn = 3; % flag any position x state cell below this count
%% ===================================================

mouseNames = fieldnames(sessionPairs);

for iMouse = 1:numel(mouseNames)
    thisMouse = mouseNames{iMouse};
    sessList  = sessionPairs.(thisMouse);

    for iSess = 1:numel(sessList)
        thisSessionName = sessList{iSess};
        sessionLabel = sprintf('%s_%s', thisMouse, thisSessionName);

        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath), warning('sfi missing for %s -- skipping.', sessionLabel); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        iStim = find(contains(stimNames, stimName), 1);
        if isempty(iStim), warning('No %s file for %s -- skipping.', stimName, sessionLabel); continue; end

        try
            load(sessionFileInfo.stimFiles(iStim).Response, 'response');
        catch ME
            warning('Could not load response for %s: %s', sessionLabel, ME.message); continue;
        end
        if ~isfield(response, 'wheelData') || numel(response.wheelData) ~= numel(response.psthData)
            warning('wheelData missing/mismatched for %s -- skipping.', sessionLabel); continue;
        end

        psthData = response.psthData;
        stimVs   = vertcat(psthData.stimValue);
        blankIdx  = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
        gridMask  = stimVs(:,1) ~= 200;
        gridPSTHIdx = find(gridMask);
        gridStim  = stimVs(gridMask, :);

        wheelTimeVec = response.wheelData(1).timeVector(:)';
        wheelRespIdx = wheelTimeVec >= respWin(1) & wheelTimeVec <= respWin(2);

        fprintf('\n=== %s ===\n', sessionLabel);
        fprintf('%-20s %8s %8s %8s %8s\n', 'Position (Az,El)', 'Stat', 'Run', 'Ambig', 'Total');

        nPosLowStat = 0; nPosLowRun = 0; nPosZeroStat = 0; nPosZeroRun = 0;

        allGroupIdx = [blankIdx, gridPSTHIdx(:)'];
        allLabels = [{'BLANK'}, arrayfun(@(i) sprintf('(%g,%g)', gridStim(i,1), gridStim(i,2)), ...
            1:numel(gridPSTHIdx), 'UniformOutput', false)];

        for gi = 1:numel(allGroupIdx)
            g = allGroupIdx(gi);
            wheelTrials = response.wheelData(g).alignedResponses;
            nTrialsHere = size(wheelTrials, 2);

            nStat = 0; nRun = 0; nAmbig = 0;
            for ti = 1:nTrialsHere
                trace = wheelTrials(wheelRespIdx, ti);
                if all(isnan(trace)), continue; end
                meanSpeed      = nanmean(trace);
                propRunning    = sum(trace > statSpeedThresh) / sum(wheelRespIdx);
                propStationary = sum(trace < runSpeedThresh)  / sum(wheelRespIdx);
                if propRunning >= propThresh && meanSpeed > runSpeedThresh
                    nRun = nRun + 1;
                elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                    nStat = nStat + 1;
                else
                    nAmbig = nAmbig + 1;
                end
            end

            flagStr = '';
            if nStat < minTrialWarn, flagStr = [flagStr ' STAT-LOW']; nPosLowStat = nPosLowStat + 1; end
            if nRun  < minTrialWarn, flagStr = [flagStr ' RUN-LOW'];  nPosLowRun  = nPosLowRun + 1;  end
            if nStat == 0, nPosZeroStat = nPosZeroStat + 1; end
            if nRun == 0,  nPosZeroRun  = nPosZeroRun + 1;  end

            fprintf('%-20s %8d %8d %8d %8d%s\n', allLabels{gi}, nStat, nRun, nAmbig, nTrialsHere, flagStr);
        end

        fprintf('--- Summary for %s ---\n', sessionLabel);
        fprintf('Positions (of %d incl. blank) with < %d stationary trials: %d\n', ...
            numel(allGroupIdx), minTrialWarn, nPosLowStat);
        fprintf('Positions with < %d running trials: %d\n', minTrialWarn, nPosLowRun);
        fprintf('Positions with ZERO stationary trials: %d\n', nPosZeroStat);
        fprintf('Positions with ZERO running trials: %d\n', nPosZeroRun);
    end
end
