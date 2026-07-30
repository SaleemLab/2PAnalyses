function allDirTuning = computeDirTuningOSI(allDirTuning, responsivenessField, makePlots)
% COMPUTEDIRTUNINGOSI  Compute orientation selectivity (both vector
% method and simple ratio method) for visually responsive boutons.
%
%   allDirTuning = computeDirTuningOSI(allDirTuning)
%   allDirTuning = computeDirTuningOSI(allDirTuning, responsivenessField)
%   allDirTuning = computeDirTuningOSI(allDirTuning, responsivenessField, makePlots)
%
% INPUTS:
%   allDirTuning        - struct array from your DirTuning pooling script.
%                         Must have fields: meanDirResponse, stimValues,
%                         and whichever responsivenessField you specify.
%   responsivenessField - (optional, default 'isResponsive_ttest') which
%                         field on allDirTuning to use as the
%                         responsiveness gate. Options typically
%                         available: 'isResponsive' (SD-heuristic, NOT
%                         recommended), 'isResponsive_ttest',
%                         'isResponsive_ranksum', 'isTunedCV'.
%   makePlots            - (optional, default true) whether to generate
%                         summary figures.
%
% OUTPUT:
%   allDirTuning - same struct array, with these fields ADDED per bouton:
%       .OSI                 - vector-method OSI (gOSI-style), NaN if
%                               not responsive or not computable
%       .prefOrientationDeg  - preferred orientation angle (0-180 deg)
%                               from the vector method
%       .OSI_simple          - simple ratio OSI, Niell & Stryker-style:
%                               (Rpref - Rortho) / (Rpref + Rortho)
%
% METHOD NOTES:
%   - Vector method: OSI = |sum_i(R_i * exp(1i*2*theta_i))| / sum_i(R_i).
%     Uses ALL directions at once; more robust to noise in any single
%     direction. This is NOT the same formula as classic Niell & Stryker
%     OSI -- see "simple ratio" below for that.
%   - Simple ratio method: collapses opposite-direction pairs (0&180,
%     45&225, 90&270, 135&315) into 4 orientation values, then computes
%     (Rpref - Rortho)/(Rpref + Rortho) using only the preferred
%     orientation and its single 90-degree-orthogonal partner. Requires
%     exactly 4 distinct orientations after folding (i.e. 8 directions,
%     45 degrees apart) -- boutons that don't fit this are skipped
%     (OSI_simple left as NaN), not guessed.
%   - Both methods HALF-WAVE RECTIFY negative (suppressed) responses to
%     0 before computing anything, since a negative "vector length" or
%     ratio term is not well-defined. This is an assumption, not
%     something stated in the original paper -- flag it if you report
%     these numbers.
%   - Neither method changes which trials/response window were used to
%     compute meanDirResponse -- that was already fixed when
%     allDirTuning was built. Make sure respWin matches what you intend
%     to report (e.g. [0.1, 3] to match the paper) BEFORE calling this.

if nargin < 2 || isempty(responsivenessField)
    responsivenessField = 'isTunedCVR2';
end
if nargin < 3 || isempty(makePlots)
    makePlots = true;
end

fprintf('computeDirTuningOSI: using "%s" as the responsiveness criterion.\n', responsivenessField);

nBoutonsTotal = numel(allDirTuning);

%% Timplalexi et al 2025
% VECTOR METHOD -- commented out entirely, keeping only the simple-ratio
% method below active. OSI/prefOrientationDeg are left as all-NaN so
% downstream code that checks for these fields doesn't error, it just
% won't have vector-method values.
OSI = nan(nBoutonsTotal, 1);
prefOrientationDeg = nan(nBoutonsTotal, 1);

% for b = 1:nBoutonsTotal
%     s = allDirTuning(b);
%
%     if ~isfield(s, responsivenessField) || ~s.(responsivenessField)
%         continue;
%     end
%     if ~isfield(s, 'stimValues')
%         continue;
%     end
%
%     thetaDeg = s.stimValues(:)';
%     thetaRad = deg2rad(thetaDeg);
%     R = s.meanDirResponse(:)';
%     R_rect = max(R, 0); % half-wave rectify
%
%     if sum(R_rect) == 0 || all(isnan(R_rect))
%         continue;
%     end
%
%     vectorSum = sum(R_rect .* exp(1i * 2 * thetaRad));
%     OSI(b) = abs(vectorSum) / sum(R_rect);
%
%     prefAngleRad = angle(vectorSum) / 2;
%     prefOrientationDeg(b) = mod(rad2deg(prefAngleRad), 180);
% end

%% simple ratio method (Niell & Stryker) 
OSI_simple = nan(nBoutonsTotal, 1);
prefOrientationDeg_simple = nan(nBoutonsTotal, 1); % preferred orientation (0-180) from the simple-ratio method, since vector-method prefOrientationDeg is no longer computed

for b = 1:nBoutonsTotal
    s = allDirTuning(b);
    if ~isfield(s, responsivenessField) || ~s.(responsivenessField)
        continue;
    end
    if ~isfield(s, 'stimValues')
        continue;
    end

    thetaDeg = mod(round(s.stimValues(:)'), 360);
    R = s.meanDirResponse(:)';
    R_rect = max(R, 0); % half-wave rectify, same convention as vector method above

    orientBins = mod(thetaDeg, 180);
    uniqueOrients = unique(orientBins);
    if numel(uniqueOrients) ~= 4
        continue; % unexpected direction spacing/count -- skip rather than guess
    end

    orientResponse = nan(1, 4);
    for oi = 1:4
        orientResponse(oi) = mean(R_rect(orientBins == uniqueOrients(oi)), 'omitnan');
    end

    [Rpref, prefIdx] = max(orientResponse);
    orthoIdx = mod(prefIdx - 1 + 2, 4) + 1;
    Rortho = orientResponse(orthoIdx);

    if (Rpref + Rortho) == 0
        continue;
    end
    OSI_simple(b) = (Rpref - Rortho) / (Rpref + Rortho);
    prefOrientationDeg_simple(b) = uniqueOrients(prefIdx);
end

%%  attach results 
for b = 1:nBoutonsTotal
    allDirTuning(b).OSI                       = OSI(b);
    allDirTuning(b).prefOrientationDeg        = prefOrientationDeg(b);        % vector method -- all-NaN now, kept for compatibility
    allDirTuning(b).prefOrientationDeg_simple = prefOrientationDeg_simple(b); % simple-ratio method -- use this instead
    allDirTuning(b).OSI_simple                = OSI_simple(b);
end

nResponsive  = sum(arrayfun(@(s) isfield(s, responsivenessField) && s.(responsivenessField), allDirTuning));
nValidOSI    = sum(~isnan(OSI));
nValidSimple = sum(~isnan(OSI_simple));
fprintf('%d / %d boutons passed the "%s" responsiveness criterion.\n', nResponsive, nBoutonsTotal, responsivenessField);
fprintf('%d boutons have a valid vector-method OSI.\n', nValidOSI);
fprintf('%d boutons have a valid simple-ratio OSI.\n', nValidSimple);

%% plots 
if makePlots
    figure('Position', [100 100 500 400]);
    histogram(OSI_simple(~isnan(OSI_simple)), 20);
    xlabel('OSI (simple ratio, N&S-style)'); ylabel('Number of boutons');
    title(sprintf('Simple-ratio OSI (n=%d, criterion: %s)', nValidSimple, responsivenessField), ...
        'Interpreter', 'none');

    figure;
    polarhistogram(deg2rad(prefOrientationDeg_simple(~isnan(prefOrientationDeg_simple))) * 2, 16);
    title('Preferred orientation distribution (simple-ratio method, doubled for display)');

    % vector-method plots (histogram, method-agreement scatter,
    % preferred-orientation polar histogram) -- commented out along with
    % the vector method itself, above
%     figure('Position', [100 100 1300 400]);
%     subplot(1,3,1);
%     histogram(OSI(~isnan(OSI)), 20);
%     xlabel('OSI (vector method / gOSI)'); ylabel('Number of boutons');
%     title(sprintf('Vector-method OSI (n=%d)', nValidOSI));
%
%     subplot(1,3,2);
%     histogram(OSI_simple(~isnan(OSI_simple)), 20);
%     xlabel('OSI (simple ratio, N&S-style)'); ylabel('Number of boutons');
%     title(sprintf('Simple-ratio OSI (n=%d)', nValidSimple));
%
%     subplot(1,3,3);
%     bothValid = ~isnan(OSI) & ~isnan(OSI_simple);
%     scatter(OSI_simple(bothValid), OSI(bothValid), 15, 'filled', 'MarkerFaceAlpha', 0.5);
%     hold on;
%     plot([0 1], [0 1], 'k--');
%     xlabel('OSI (simple ratio)'); ylabel('OSI (vector method)');
%     title(sprintf('Method agreement (n=%d)', sum(bothValid)));
%     axis square; xlim([0 1]); ylim([0 1]);
%
%     sgtitle(sprintf('OSI: vector method vs simple ratio (criterion: %s)', responsivenessField), ...
%         'Interpreter', 'none');
%
%     figure;
%     polarhistogram(deg2rad(prefOrientationDeg(~isnan(prefOrientationDeg))) * 2, 16);
%     title('Preferred orientation distribution (vector method, doubled for display)');
end

end