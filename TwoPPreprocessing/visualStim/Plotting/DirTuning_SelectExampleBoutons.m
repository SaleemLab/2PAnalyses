% DirTuning_SelectExampleBoutons.m
%
% Picks example boutons to feature in figures, based on OSI/DSI already
% computed (run computeDirTuningOSI and computeDirTuningDSI first).
%
% IMPORTANT: candidates are filtered to require isTunedCV (cross-
% validated R^2 significance) BEFORE ranking by OSI/DSI. Ranking by raw
% OSI/DSI alone is a selection-bias trap -- taking the single highest
% value across many noisy boutons reliably finds flukes (flat, noisy
% boutons whose noise happens to produce a near-1 OSI/DSI by chance),
% not genuine tuning. Requiring isTunedCV first ensures the tuning shape
% itself was validated on held-out trials before trusting its OSI/DSI.

responsiveIdx = find([allDirTuning.isResponsive_ttest]); % or whichever criterion you're using

reliableIdx = responsiveIdx([allDirTuning(responsiveIdx).isTunedCV] == 1);
fprintf('%d / %d responsive boutons also pass isTunedCV (reliable tuning shape) -- ranking examples from these only.\n', ...
    numel(reliableIdx), numel(responsiveIdx));

if isempty(reliableIdx)
    error('No boutons pass both isResponsive_ttest AND isTunedCV -- cannot select reliable examples. Consider loosening criteria or checking your data.');
end

allOSI = [allDirTuning(reliableIdx).OSI];
allDSI = [allDirTuning(reliableIdx).DSI];

% "Direction selective" example: high DSI, but NOT necessarily high OSI
% (mirrors reference panel E: OSI=0.256, DSI=0.653 -- direction matters
% more than orientation for this cell)
[~, sortedByDSI] = sort(allDSI, 'descend', 'MissingPlacement', 'last');
directionSelectiveExample = reliableIdx(sortedByDSI(1));

% "Orientation selective" example: high OSI, but NOT necessarily high DSI
% (mirrors reference panel H: OSI=0.808, DSI=0.278 -- orientation matters
% more, direction of motion within that orientation doesn't)
[~, sortedByOSI] = sort(allOSI, 'descend', 'MissingPlacement', 'last');
orientationSelectiveExample = reliableIdx(sortedByOSI(1));

fprintf('Direction-selective example: bouton %d (OSI=%.3f, DSI=%.3f, cvR2=%.3f)\n', ...
    directionSelectiveExample, allDirTuning(directionSelectiveExample).OSI, ...
    allDirTuning(directionSelectiveExample).DSI, allDirTuning(directionSelectiveExample).cvR2);
fprintf('Orientation-selective example: bouton %d (OSI=%.3f, DSI=%.3f, cvR2=%.3f)\n', ...
    orientationSelectiveExample, allDirTuning(orientationSelectiveExample).OSI, ...
    allDirTuning(orientationSelectiveExample).DSI, allDirTuning(orientationSelectiveExample).cvR2);

% If you want a cleaner "one but not the other" pick (e.g. high DSI
% specifically paired with LOW OSI, matching the reference figure's
% contrast more precisely), filter first:
lowOSI_highDSI = reliableIdx(allOSI < 0.4 & allDSI > 0.5);
highOSI_lowDSI = reliableIdx(allOSI > 0.6 & allDSI < 0.4);
fprintf('\nBoutons with low OSI + high DSI (cleaner "direction-selective" examples): %s\n', mat2str(lowOSI_highDSI));
fprintf('Boutons with high OSI + low DSI (cleaner "orientation-selective" examples): %s\n', mat2str(highOSI_lowDSI));

% Plot both examples
plotExampleDirTuningBouton(allDirTuning, directionSelectiveExample);
plotExampleDirTuningBouton(allDirTuning, orientationSelectiveExample);
