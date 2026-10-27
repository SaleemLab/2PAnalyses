function results = summarizeHighlyCorrelatedROIs(sessionMetrics)
    % SUMMARIZEHIGHLYCORRELATEDROIS Reports, across sessions, how many ROIs
    % were pruned as highly correlated duplicates.
    %
    % Works directly off sessionMetrics (e.g. RSPData/VISpData, the output
    % of getTuningDataByCondition), using the same 'uniqueBoutonIdx' field
    % that appendFilteredROIs already reads for its FilterDuplicateBoutons
    % filter -- no separate sessionFileInfo/file path lookup needed.
    %
    % nDiscarded per session = (total ROIs) - length(uniqueBoutonIdx)
    % i.e. the number of duplicate ROIs that were NOT the kept
    % highest-SNR representative of their correlated group.
    %
    % Only applies to sessions where TypeImaged is 'Boutons' (uniqueBoutonIdx
    % is a bouton-specific field); soma sessions are skipped with a note.
    
    numSessions = length(sessionMetrics);
    
    nDiscarded = nan(1, numSessions);
    nTotalROIs = nan(1, numSessions);
    isBoutonSession = false(1, numSessions);
    
    for s = 1:numSessions
        thisSession = sessionMetrics(s);
        
        if ~isfield(thisSession, 'TypeImaged') || ~strcmpi(thisSession.TypeImaged, 'Boutons')
            continue; % not a bouton session, uniqueBoutonIdx doesn't apply
        end
        isBoutonSession(s) = true;
        
        if ~isfield(thisSession, 'ConditionData') || isempty(thisSession.ConditionData)
            warning('Session %d (Mouse %s) has no ConditionData. Skipping.', s, thisSession.MouseID);
            continue;
        end
        condNames = fieldnames(thisSession.ConditionData);
        sampleData = thisSession.ConditionData.(condNames{1}).LapActivity;
        numTotalROIs(s) = size(sampleData, 1); %#ok<AGROW>
        nTotalROIs(s) = numTotalROIs(s);
        
        if isfield(thisSession, 'uniqueBoutonIdx') && ~isempty(thisSession.uniqueBoutonIdx)
            nDiscarded(s) = nTotalROIs(s) - length(thisSession.uniqueBoutonIdx);
        else
            warning('Mouse %s Session %s: TypeImaged is Boutons but uniqueBoutonIdx is missing. Skipping.', ...
                thisSession.MouseID, thisSession.Session);
        end
    end
    
    validSessions = ~isnan(nDiscarded);
    nValid = sum(validSessions);
    nSkippedNonBouton = sum(~isBoutonSession);
    
    if nValid == 0
        error(['No sessions had a usable uniqueBoutonIdx field. Check that TypeImaged is set to ' ...
               '''Boutons'' and that findHightlyCorrelatedROIs / your data assembly step actually ' ...
               'populated uniqueBoutonIdx on sessionMetrics.']);
    end
    
    pctDiscarded = (nDiscarded(validSessions) ./ nTotalROIs(validSessions)) * 100;
    
    meanDiscarded = mean(nDiscarded(validSessions));
    sdDiscarded   = std(nDiscarded(validSessions));
    meanPctDiscarded = mean(pctDiscarded);
    sdPctDiscarded   = std(pctDiscarded);
    
    fprintf('\n=== Highly Correlated ROI Summary (n = %d bouton sessions) ===\n', nValid);
    if nSkippedNonBouton > 0
        fprintf('(%d non-bouton sessions skipped, uniqueBoutonIdx does not apply to somas)\n', nSkippedNonBouton);
    end
    fprintf('Duplicates discarded per session: %.1f +/- %.1f (mean +/- SD)\n', meanDiscarded, sdDiscarded);
    fprintf('  as %% of that session''s ROIs:     %.1f +/- %.1f%%\n', meanPctDiscarded, sdPctDiscarded);
    
    results.nDiscarded       = nDiscarded;
    results.nTotalROIs       = nTotalROIs;
    results.pctDiscarded     = pctDiscarded;
    results.meanDiscarded    = meanDiscarded;
    results.sdDiscarded      = sdDiscarded;
    results.meanPctDiscarded = meanPctDiscarded;
    results.sdPctDiscarded   = sdPctDiscarded;
end