% RFMapping_TrialCountMatched_SingleSession.m
%
% Directly tests whether the running/stationary responsiveness gap in
% M25132_20260306 (the one session with no zero-trial position holes)
% survives trial-count matching: for each bouton, at EACH grid position
% (+ blank), downsamples BOTH states to the SAME (smaller) trial count,
% then recomputes the ANOVA+2SD-threshold responsiveness criterion
% separately for each state on this matched data. Repeated with
% different random subsamples (nRepeats) and averaged, since single-digit
% trial counts per position make any one subsample noisy.
%
% Same logic as DotFields_TrialCountMatchedComparison.m, adapted to RF
% mapping's grid-position structure instead of speed conditions.

%%load 
thisMouse       = 'M25132';
thisSessionName = '20260306';
stimName = 'RFMapping';
respWin  = [0.5 3];
baseWin  = [-1.0 0];

runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;

ALPHA = 0.05;
NSD   = 2;
nRepeats = 10;
rng(1);
%% ===================================================

infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
loadedInfo      = load(infoPath, 'sessionFileInfo');
sessionFileInfo = loadedInfo.sessionFileInfo;
stimNames       = {sessionFileInfo.stimFiles.name};
iStim = find(contains(stimNames, stimName), 1);
load(sessionFileInfo.stimFiles(iStim).Response, 'response');

if ~isfield(response, 'wheelData') || numel(response.wheelData) ~= numel(response.psthData)
    error('wheelData missing/mismatched for this session.');
end

psthData = response.psthData;
stimVs   = vertcat(psthData.stimValue);
nROI     = size(psthData(1).alignedResponses, 1);

blankIdx  = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
gridMask  = stimVs(:,1) ~= 200;
gridPSTHIdx = find(gridMask);
gridStim  = stimVs(gridMask, :);
nPositions = numel(gridPSTHIdx);

uAz      = sort(unique(gridStim(:,1)), 'ascend');
uEl_plot = sort(unique(gridStim(:,2)), 'descend');

timeVec = psthData(1).timeVector(:)';
respIdx = timeVec >= respWin(1) & timeVec <= respWin(2);
baseIdx = timeVec >= baseWin(1) & timeVec < baseWin(2);

wheelTimeVec = response.wheelData(1).timeVector(:)';
wheelRespIdx = wheelTimeVec >= respWin(1) & wheelTimeVec <= respWin(2);

%% ===================== classify every trial by state =====================
runFlagByGroup = cell(numel(psthData), 1);
for g = 1:numel(psthData)
    wheelTrials = response.wheelData(g).alignedResponses;
    nTrialsHere = size(wheelTrials, 2);
    rf = nan(nTrialsHere, 1);
    for ti = 1:nTrialsHere
        trace = wheelTrials(wheelRespIdx, ti);
        if all(isnan(trace)), continue; end
        meanSpeed      = nanmean(trace);
        propRunning    = sum(trace > statSpeedThresh) / sum(wheelRespIdx);
        propStationary = sum(trace < runSpeedThresh)  / sum(wheelRespIdx);
        if propRunning >= propThresh && meanSpeed > runSpeedThresh
            rf(ti) = 1;
        elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
            rf(ti) = 0;
        end
    end
    runFlagByGroup{g} = rf;
end

%% ===================== per-position matched trial count =====================
allGroupIdx = [gridPSTHIdx(:)', blankIdx];
nMatchedPerGroup = nan(numel(allGroupIdx), 1);
for gi = 1:numel(allGroupIdx)
    g = allGroupIdx(gi);
    rf = runFlagByGroup{g};
    nStatHere = sum(rf == 0);
    nRunHere  = sum(rf == 1);
    nMatchedPerGroup(gi) = min(nStatHere, nRunHere);
end

fprintf('Matched trial count per position+blank (min across states): %s\n', mat2str(nMatchedPerGroup'));
if any(nMatchedPerGroup < 1)
    warning('At least one position has 0 trials in one state even after matching floor -- those positions will be excluded from the matched ANOVA.');
end

%% ===================== per-bouton matched responsiveness, repeated subsampling =====================
isResp_matched_stat = false(nROI, nRepeats);
isResp_matched_run  = false(nROI, nRepeats);
isResp_unmatched_stat = false(nROI, 1);
isResp_unmatched_run  = false(nROI, 1);

for iROI = 1:nROI
    for si = 1:2
        stateVal = si - 1;

        % --- UNMATCHED (original, all available trials in this state) ---
        bTrials = squeeze(psthData(blankIdx).alignedResponses(iROI, :, :));
        if isvector(bTrials), bTrials = bTrials(:); end
        blankMask = (runFlagByGroup{blankIdx} == stateVal);
        bTrialsState = bTrials(:, blankMask);
        if size(bTrialsState, 2) >= 2
            tB = mean(bTrialsState(baseIdx, :), 1, 'omitnan');
            bCorr = bTrialsState - tB;
            blankMeansUnmatched = mean(bCorr(respIdx, :), 1, 'omitnan')';

            yU = blankMeansUnmatched(:); grpU = repmat(nPositions+1, numel(blankMeansUnmatched), 1);
            gridResp = nan(numel(uEl_plot), numel(uAz));
            for pIdx = 1:nPositions
                g = gridPSTHIdx(pIdx);
                trialsAtPos = squeeze(psthData(g).alignedResponses(iROI, :, :));
                if isvector(trialsAtPos), trialsAtPos = trialsAtPos(:); end
                posMask = (runFlagByGroup{g} == stateVal);
                trialsAtPosState = trialsAtPos(:, posMask);
                if isempty(trialsAtPosState), continue; end
                tB2 = mean(trialsAtPosState(baseIdx, :), 1, 'omitnan');
                corr2 = trialsAtPosState - tB2;
                posMeans = mean(corr2(respIdx, :), 1, 'omitnan')';
                yU = [yU; posMeans(:)]; grpU = [grpU; repmat(pIdx, numel(posMeans), 1)];
                rowIdx = find(uEl_plot == gridStim(pIdx,2),1); colIdx = find(uAz == gridStim(pIdx,1),1);
                gridResp(rowIdx, colIdx) = mean(posMeans, 'omitnan');
            end
            validU = ~isnan(yU) & ~isnan(grpU);
            if numel(unique(grpU(validU))) >= 2
                pU = anova1(yU(validU), grpU(validU), 'off');
                bMean = mean(blankMeansUnmatched,'omitnan'); bStd = std(blankMeansUnmatched,'omitnan');
                [prefValU] = max(gridResp(:), [], 'omitnan');
                respU = (pU < ALPHA) && (prefValU > (bMean + NSD*bStd));
                if si==1, isResp_unmatched_stat(iROI) = respU; else, isResp_unmatched_run(iROI) = respU; end
            end
        end
    end

    % --- MATCHED (downsampled to equal trial count per position+blank) ---
    for rep = 1:nRepeats
        for si = 1:2
            stateVal = si - 1;

            % matched blank
            bTrials = squeeze(psthData(blankIdx).alignedResponses(iROI, :, :));
            if isvector(bTrials), bTrials = bTrials(:); end
            blankMask = find(runFlagByGroup{blankIdx} == stateVal);
            nMatchBlank = nMatchedPerGroup(allGroupIdx == blankIdx);
            if numel(blankMask) < nMatchBlank || nMatchBlank < 2, continue; end
            sampBlank = blankMask(randperm(numel(blankMask), nMatchBlank));
            bTrialsM = bTrials(:, sampBlank);
            tB = mean(bTrialsM(baseIdx, :), 1, 'omitnan');
            bCorr = bTrialsM - tB;
            blankMeansM = mean(bCorr(respIdx, :), 1, 'omitnan')';

            yM = blankMeansM(:); grpM = repmat(nPositions+1, numel(blankMeansM), 1);
            gridResp = nan(numel(uEl_plot), numel(uAz));
            for pIdx = 1:nPositions
                g = gridPSTHIdx(pIdx);
                trialsAtPos = squeeze(psthData(g).alignedResponses(iROI, :, :));
                if isvector(trialsAtPos), trialsAtPos = trialsAtPos(:); end
                posMask = find(runFlagByGroup{g} == stateVal);
                nMatchPos = nMatchedPerGroup(allGroupIdx == g);
                if numel(posMask) < nMatchPos || nMatchPos < 1, continue; end
                sampPos = posMask(randperm(numel(posMask), nMatchPos));
                trialsAtPosM = trialsAtPos(:, sampPos);
                tB2 = mean(trialsAtPosM(baseIdx, :), 1, 'omitnan');
                corr2 = trialsAtPosM - tB2;
                posMeans = mean(corr2(respIdx, :), 1, 'omitnan')';
                yM = [yM; posMeans(:)]; grpM = [grpM; repmat(pIdx, numel(posMeans), 1)];
                rowIdx = find(uEl_plot == gridStim(pIdx,2),1); colIdx = find(uAz == gridStim(pIdx,1),1);
                gridResp(rowIdx, colIdx) = mean(posMeans, 'omitnan');
            end

            validM = ~isnan(yM) & ~isnan(grpM);
            if numel(unique(grpM(validM))) < 2, continue; end
            pM = anova1(yM(validM), grpM(validM), 'off');
            bMean = mean(blankMeansM,'omitnan'); bStd = std(blankMeansM,'omitnan');
            prefValM = max(gridResp(:), [], 'omitnan');
            respM = (pM < ALPHA) && (prefValM > (bMean + NSD*bStd));

            if si == 1
                isResp_matched_stat(iROI, rep) = respM;
            else
                isResp_matched_run(iROI, rep) = respM;
            end
        end
    end

    if mod(iROI, 100) == 0
        fprintf('Processed %d / %d boutons...\n', iROI, nROI);
    end
end

fracResp_matched_stat = mean(isResp_matched_stat, 2);
fracResp_matched_run  = mean(isResp_matched_run, 2);
isResp_matched_stat_final = fracResp_matched_stat >= 0.5;
isResp_matched_run_final  = fracResp_matched_run  >= 0.5;

fprintf('\n=== M25132_20260306: trial-count-matched vs unmatched ===\n');
fprintf('%-30s %10s %10s\n', 'Condition', 'Stationary', 'Running');
fprintf('%-30s %10d %10d\n', 'Unmatched (original)', sum(isResp_unmatched_stat), sum(isResp_unmatched_run));
fprintf('%-30s %10d %10d\n', 'Matched (equal trials/position)', sum(isResp_matched_stat_final), sum(isResp_matched_run_final));

nSplitOnly_unmatched = sum(~isResp_unmatched_stat & isResp_unmatched_run | isResp_unmatched_stat & ~isResp_unmatched_run);
nSplitOnly_matched   = sum(~isResp_matched_stat_final & isResp_matched_run_final | isResp_matched_stat_final & ~isResp_matched_run_final);
fprintf('\nSplit-only count (either state, not both) -- unmatched: %d | matched: %d\n', ...
    nSplitOnly_unmatched, nSplitOnly_matched);
