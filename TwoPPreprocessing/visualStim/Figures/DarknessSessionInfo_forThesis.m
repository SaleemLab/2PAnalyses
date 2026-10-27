%% ===== Running speed tuning in darkness: population summary =====
clear; clc;

%% Settings
useField              = 'dFFNeuropilCorrected';
targetStruct          = 'tuningCurve';
onlyBoutons           = true;
binFraction           = 0.07;
minBins               = 5;
pval_shuffleThreshold = 0.01;
minR2Threshold        = -Inf;     % no R2 threshold: every speed-tuned bouton keeps its shape

%% Sessions (darkness, RSP boutons)
pairs = struct;
pairs.M25132 = {'20260219', '20260223', '20260226', '20260228', '20260303', '20260306', '20260313'};
pairs.M25133 = {'20260219', '20260221', '20260223'};
pairs.M26003 = {'20260316', '20260322', '20260324', '20260325'};
filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'Darkness'});
mouseInfo = sessionsToProcess(filteredTable);

%% Accumulators
tuningTypesPlot    = {'highpass', 'bandpass', 'lowpass', 'trough_inverted', 'unclassified'};
allSessionPcts_fig = [];
allTypes        = {};
allPrefSpeeds   = [];
allR2s          = [];
allBinDurations = [];

shapeNames      = {'highpass', 'lowpass', 'bandpass', 'trough_inverted'};
shapeCountsSess = zeros(0, numel(shapeNames));   % speed-tuned boutons of each shape, per included session
mouseInc        = strings(0, 1);                 % mouse ID per included session
shapeLabelAll   = {};                            % shape of each speed-tuned bouton
sessIDAll       = [];                            % session number of each speed-tuned bouton

logFields = {'Mouse', 'Session', 'Type', 'Status', 'RecDur_s', 'RunTime_s', 'RunFrac', 'nBins', ...
    'BinMean_s', 'BinMin_s', 'BinMax_s', 'Target7pct_s', 'nROIs', ...
    'MedSpeed', 'FastSpeed', 'MaxSpeed', 'nTunedMoving', 'nModFull', 'nBoth'};
sessionLog = cell2struct(cell(numel(logFields), 0), logFields, 1);
newEntry = @(m, s, t, status) cell2struct([{string(m); string(s); string(t); string(status)}; ...
    num2cell(nan(numel(logFields) - 4, 1))], logFields, 1);

%% Load data across sessions
for thisMouse = 1:size(mouseInfo, 1)
    mousenumber  = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};

    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if isempty(infoPath) || ~isfile(infoPath)
            sessionLog(end+1) = newEntry(mousenumber, sessionName, "", "skipped: no session info file"); %#ok<SAGROW>
            continue;
        end

        matchRow = find(strcmp(filteredTable.MouseID, mousenumber) & strcmp(filteredTable.Session, sessionName), 1);
        if isempty(matchRow)
            sessionLog(end+1) = newEntry(mousenumber, sessionName, "", "skipped: not in metadata table"); %#ok<SAGROW>
            continue;
        end
        typeImaged = filteredTable.TypeImaged{matchRow};

        if onlyBoutons && ~strcmpi(typeImaged, 'Boutons')
            sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "skipped: not a bouton session"); %#ok<SAGROW>
            continue;
        end

        try
            loadedInfo      = load(infoPath, 'sessionFileInfo');
            sessionFileInfo = loadedInfo.sessionFileInfo;
            stimNames       = {sessionFileInfo.stimFiles.name};
            targetIdx       = find(contains(stimNames, 'Darkness'));

            if isempty(targetIdx)
                sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "skipped: no darkness stimulus"); %#ok<SAGROW>
                continue;
            end

            for thisStim = 1:length(targetIdx)
                thisStimName = stimNames{targetIdx(thisStim)};
                stimFileName = sprintf('%s_%s_Response_%s.mat', mousenumber, sessionName, thisStimName);
                fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

                if ~isfile(fileFullPath)
                    sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "skipped: no response file"); %#ok<SAGROW>
                    continue;
                end

                loadedData = load(fileFullPath, 'response');
                resp       = loadedData.response;

                if ~isfield(resp, targetStruct) || ~isfield(resp.(targetStruct), useField) || ...
                   ~isfield(resp.(targetStruct).(useField), 'classification')
                    sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "skipped: not classified"); %#ok<SAGROW>
                    continue;
                end

                % --- Behaviour ---
                occ    = resp.(targetStruct).occupancy.moving;       % running time per bin (s)
                statT  = resp.(targetStruct).occupancy.stationary;   % stationary time (s)
                runT   = sum(occ);
                recDur = statT + runT;                               % excludes artefact frames (>150 cm/s)
                nBins  = numel(occ);
                edges  = resp.(targetStruct).speedBins;              % quantiles of running speed

                entry = newEntry(mousenumber, sessionName, typeImaged, "");
                entry.RecDur_s     = recDur;
                entry.RunTime_s    = runT;
                entry.RunFrac      = runT / recDur;
                entry.nBins        = nBins;
                entry.BinMean_s    = mean(occ);
                entry.BinMin_s     = min(occ);
                entry.BinMax_s     = max(occ);
                entry.Target7pct_s = binFraction * recDur;
                entry.MedSpeed     = median(edges);                  % middle edge ~ median running speed
                entry.FastSpeed    = edges(end-1);                   % second-to-last edge ~ 93rd percentile
                entry.MaxSpeed     = edges(end);

                % --- Exclusion before any pooling ---
                if nBins < minBins
                    entry.Status = sprintf("excluded: only %d bin(s)", nBins);
                    sessionLog(end+1) = entry; %#ok<SAGROW>
                    continue;
                end

                % --- Tuning ---
                cls        = resp.(targetStruct).(useField).classification;
                pvalMoving = resp.(targetStruct).(useField).pValMoving;   % tuned to running speed
                pvalFull   = resp.(targetStruct).(useField).pValFull;     % modulated across behavioural states

                classifiedIdx   = find(cls.R2 >= minR2Threshold & pvalMoving <= pval_shuffleThreshold);
                unclassifiedIdx = find(cls.R2 <  minR2Threshold & pvalMoving <= pval_shuffleThreshold);
                untunedIdx      = find(pvalMoving > pval_shuffleThreshold);

                entry.Status       = "included";
                entry.nROIs        = numel(pvalMoving);
                entry.nTunedMoving = sum(pvalMoving <= pval_shuffleThreshold);
                entry.nModFull     = sum(pvalFull   <= pval_shuffleThreshold);
                entry.nBoth        = sum(pvalFull   <= pval_shuffleThreshold & pvalMoving <= pval_shuffleThreshold);
                sessionLog(end+1) = entry; %#ok<SAGROW>
                allBinDurations = [allBinDurations; occ(:)]; %#ok<AGROW>

                % --- Shapes per session (for panel f and the session comparisons) ---
                tunedTypes = cls.tuningType(pvalMoving <= pval_shuffleThreshold);
                shapeCountsSess(end+1, :) = cellfun(@(t) sum(strcmp(tunedTypes, t)), shapeNames); %#ok<AGROW>
                mouseInc(end+1, 1)        = string(mousenumber); %#ok<AGROW>
                shapeLabelAll = [shapeLabelAll; tunedTypes(:)]; %#ok<AGROW>
                sessIDAll     = [sessIDAll; repmat(size(shapeCountsSess, 1), numel(tunedTypes), 1)]; %#ok<AGROW>

                % --- Per-session percentages for figures ---
                totalSig = numel(classifiedIdx) + numel(unclassifiedIdx);
                if totalSig > 0
                    sessionCounts = zeros(1, numel(tuningTypesPlot));
                    for k = 1:numel(tuningTypesPlot)-1
                        sessionCounts(k) = sum(strcmp(cls.tuningType(classifiedIdx), tuningTypesPlot{k}));
                    end
                    sessionCounts(end) = numel(unclassifiedIdx);
                    allSessionPcts_fig = [allSessionPcts_fig; sessionCounts / totalSig * 100]; %#ok<AGROW>
                end

                % --- Pool boutons ---
                allTypes      = [allTypes;      cls.tuningType(classifiedIdx)]; %#ok<AGROW>
                allPrefSpeeds = [allPrefSpeeds; cls.preferredSpeed(classifiedIdx)]; %#ok<AGROW>
                allR2s        = [allR2s;        cls.R2(classifiedIdx)]; %#ok<AGROW>
                allTypes      = [allTypes;      repmat({'unclassified'}, numel(unclassifiedIdx), 1)]; %#ok<AGROW>
                allPrefSpeeds = [allPrefSpeeds; cls.preferredSpeed(unclassifiedIdx)]; %#ok<AGROW>
                allR2s        = [allR2s;        cls.R2(unclassifiedIdx)]; %#ok<AGROW>
                allTypes      = [allTypes;      repmat({'untuned'}, numel(untunedIdx), 1)]; %#ok<AGROW>
                allPrefSpeeds = [allPrefSpeeds; nan(numel(untunedIdx), 1)]; %#ok<AGROW>
                allR2s        = [allR2s;        cls.R2(untunedIdx)]; %#ok<AGROW>
            end
        catch ME
            sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "error: " + string(ME.message)); %#ok<SAGROW>
        end
    end
end

%% ================= SUMMARY =================
S   = sessionLog;
st  = [S.Status];
inc = st == "included";
Si  = S(inc);

fprintf('\n===== SESSIONS =====\n');
fprintf('Listed: %d | included: %d | excluded (<%d bins): %d | skipped/error: %d\n', numel(S), sum(inc), minBins, ...
    sum(startsWith(st, "excluded")), sum(startsWith(st, "skipped") | startsWith(st, "error")));
disp(table([S.Mouse]', [S.Session]', st', 'VariableNames', {'Mouse', 'Session', 'Status'}));
fprintf('Mice included: %d (%s)\n', numel(unique([Si.Mouse])), strjoin(unique([Si.Mouse]), ', '));

% --- Behaviour ---
statTime = [Si.RecDur_s] - [Si.RunTime_s];
fprintf('\n===== BEHAVIOUR (included sessions, mean +/- SD) =====\n');
fprintf('Recording duration:   %.1f +/- %.1f min (range %.1f to %.1f)\n', mean([Si.RecDur_s])/60, std([Si.RecDur_s])/60, min([Si.RecDur_s])/60, max([Si.RecDur_s])/60);
fprintf('Running time:         %.1f +/- %.1f min\n', mean([Si.RunTime_s])/60, std([Si.RunTime_s])/60);
fprintf('Stationary time:      %.1f +/- %.1f min\n', mean(statTime)/60, std(statTime)/60);
fprintf('Fraction running:     %.0f +/- %.0f%% (range %.0f to %.0f%%)\n', 100*mean([Si.RunFrac]), 100*std([Si.RunFrac]), 100*min([Si.RunFrac]), 100*max([Si.RunFrac]));
fprintf('Median running speed: %.1f +/- %.1f cm/s\n', mean([Si.MedSpeed]), std([Si.MedSpeed]));
fprintf('Fast running (~93rd percentile): %.1f +/- %.1f cm/s | maximum %.1f to %.1f cm/s\n', ...
    mean([Si.FastSpeed]), std([Si.FastSpeed]), min([Si.MaxSpeed]), max([Si.MaxSpeed]));
fprintf('Speed bins per session: %d to %d | running time per bin: %.1f s (range %.1f to %.1f)\n', ...
    min([Si.nBins]), max([Si.nBins]), mean(allBinDurations), min(allBinDurations), max(allBinDurations));

% --- Tuning ---
nR = sum([Si.nROIs]);
pctPerFOV = 100 * [Si.nTunedMoving] ./ [Si.nROIs];
fprintf('\n===== TUNING =====\n');
fprintf('Boutons: %d\n', nR);
fprintf('Modulated across behavioural states (pValFull): %d (%.1f%%)\n', sum([Si.nModFull]), 100*sum([Si.nModFull])/nR);
fprintf('Tuned to running speed (pValMoving):           %d (%.1f%%) | per FOV %.1f +/- %.1f%%\n', ...
    sum([Si.nTunedMoving]), 100*sum([Si.nTunedMoving])/nR, mean(pctPerFOV), std(pctPerFOV));
fprintf('Speed-tuned AND modulated across states:       %d\n', sum([Si.nBoth]));
fprintf('Modulated across states but NOT speed-tuned:   %d\n', sum([Si.nModFull]) - sum([Si.nBoth]));

% --- Shapes, pooled ---
sigTypes = allTypes(~strcmp(allTypes, 'untuned'));
fprintf('\n===== SHAPES (pooled speed-tuned boutons) =====\n');
for k = 1:numel(tuningTypesPlot)
    n = sum(strcmp(sigTypes, tuningTypesPlot{k}));
    fprintf('  %-16s %4d (%.1f%%)\n', tuningTypesPlot{k}, n, 100*n/max(numel(sigTypes), 1));
end
for t = {'bandpass', 'trough_inverted'}
    ps = allPrefSpeeds(strcmp(allTypes, t{1}));
    if ~isempty(ps)
        fprintf('  %s preferred speed: median %.1f cm/s | range %.1f to %.1f cm/s\n', t{1}, median(ps), min(ps), max(ps));
    end
end

% --- Shapes per session (panel f) ---
nTunedSess = sum(shapeCountsSess, 2);
pctShape   = 100 * shapeCountsSess ./ nTunedSess;
fprintf('\n===== SHAPES PER SESSION (%% of speed-tuned boutons) =====\n');
disp(array2table([nTunedSess, round(pctShape, 1)], 'VariableNames', [{'nTuned'}, shapeNames], ...
    'RowNames', cellstr([Si.Mouse]' + "_" + [Si.Session]')));
for k = 1:numel(shapeNames)
    fprintf('  %-16s %.1f +/- %.1f%% (mean +/- SD across sessions)\n', shapeNames{k}, mean(pctShape(:,k)), std(pctShape(:,k)));
end

% --- Do shape proportions differ across sessions more than chance? ---
lab = shapeLabelAll;
lab(~ismember(lab, {'highpass', 'lowpass'})) = {'other'};     % rare shapes merged
[~, ~, labID] = unique(lab);
obsChi2 = localChi2(sessIDAll, labID);
nPerm = 10000;  rng(1);  nullChi2 = zeros(nPerm, 1);
for k = 1:nPerm
    nullChi2(k) = localChi2(sessIDAll, labID(randperm(numel(labID))));
end
pPerm = (sum(nullChi2 >= obsChi2) + 1) / (nPerm + 1);
fprintf('\nShape proportions differ across sessions: chi-square = %.1f, permutation p = %.4f\n', obsChi2, pPerm);

% --- Is the high-pass proportion related to behaviour? ---
[rhoRun, pRun] = corr([Si.RunFrac]', pctShape(:,1), 'type', 'Spearman');
[rhoSpd, pSpd] = corr([Si.MedSpeed]', pctShape(:,1), 'type', 'Spearman');
fprintf('\nHigh-pass %% vs fraction of time running: rho = %.2f, p = %.3f (n = %d sessions)\n', rhoRun, pRun, numel(Si));
fprintf('High-pass %% vs median running speed:    rho = %.2f, p = %.3f\n', rhoSpd, pSpd);

% --- By mouse ---
fprintf('\n===== BY MOUSE =====\n');
for m = unique(mouseInc)'
    im = mouseInc == m;
    fprintf('%s (%d sessions): running %.0f%% | high-pass %.0f%% | low-pass %.0f%% | speed-tuned %.1f%%\n', m, sum(im), ...
        100*mean([Si(im).RunFrac]), mean(pctShape(im,1)), mean(pctShape(im,2)), mean(pctPerFOV(im)));
end


for m = unique(mouseInc)'
    im = mouseInc == m;
    fprintf('%s: high-pass %.0f +/- %.0f%% | low-pass %.0f +/- %.0f%%\n', m, ...
        mean(pctShape(im,1)), std(pctShape(im,1)), mean(pctShape(im,2)), std(pctShape(im,2)));
end

%% ---------- local function (must stay at the end of the script) ----------
function c = localChi2(sid, lid)
    T = accumarray([sid(:), lid(:)], 1);
    E = sum(T, 2) * sum(T, 1) / sum(T(:));
    ok = E > 0;
    c = sum((T(ok) - E(ok)).^2 ./ E(ok));
end