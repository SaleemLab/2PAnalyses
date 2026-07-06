function allDirTuning = computeDirTuningDSI(allDirTuning, responsivenessField, makePlots)
% COMPUTEDIRTUNINGDSI  Compute direction selectivity (both vector method
% and simple ratio method) for visually responsive boutons.
%
%   allDirTuning = computeDirTuningDSI(allDirTuning)
%   allDirTuning = computeDirTuningDSI(allDirTuning, responsivenessField)
%   allDirTuning = computeDirTuningDSI(allDirTuning, responsivenessField, makePlots)
%
% Same interface and conventions as computeDirTuningOSI.m -- see that
% file for full documentation of the shared assumptions (half-wave
% rectification, responsivenessField options, etc.). This function
% measures DIRECTION selectivity (does the bouton prefer one direction
% of motion over its exact opposite), which is a DIFFERENT question from
% orientation selectivity (does it prefer an axis, regardless of which
% way it drifts). Run both if you want the full picture.
%
% INPUTS:
%   allDirTuning        - struct array from your DirTuning pooling script.
%                         Must have fields: meanDirResponse, stimValues,
%                         and whichever responsivenessField you specify.
%   responsivenessField - (optional, default 'isResponsive_ttest')
%   makePlots            - (optional, default true)
%
% OUTPUT:
%   allDirTuning - same struct array, with these fields ADDED per bouton:
%       .DSI                - vector-method DSI
%       .prefDirectionDeg   - preferred direction angle (0-360 deg) from
%                             the vector method
%       .DSI_simple         - simple ratio DSI: (Rpref - Ropp)/(Rpref + Ropp)
%                             where Ropp is the response 180 degrees
%                             from the preferred DIRECTION (no folding
%                             needed here, unlike OSI -- 180 degrees away
%                             is always one of your 8 measured directions
%                             if they're 45 degrees apart)
%
% METHOD NOTES:
%   - Vector method: DSI = |sum_i(R_i * exp(1i*theta_i))| / sum_i(R_i).
%     Uses theta_i directly (NOT 2*theta_i like OSI), since direction is
%     periodic over 360 degrees, not 180.
%   - Simple ratio method: DSI = (Rpref - Ropp)/(Rpref + Ropp), using the
%     raw preferred direction and the single direction opposite it.
%   - Both HALF-WAVE RECTIFY negative responses to 0 before computing,
%     same convention as the OSI function -- flag this if reporting.
%   - Response window (respWin) used to compute meanDirResponse was
%     already fixed when allDirTuning was built; this function doesn't
%     change it.

if nargin < 2 || isempty(responsivenessField)
    responsivenessField = 'isResponsive_ttest';
end
if nargin < 3 || isempty(makePlots)
    makePlots = true;
end

fprintf('computeDirTuningDSI: using "%s" as the responsiveness criterion.\n', responsivenessField);

nBoutonsTotal = numel(allDirTuning);

%% ===================== vector method =====================
DSI = nan(nBoutonsTotal, 1);
prefDirectionDeg = nan(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s = allDirTuning(b);

    if ~isfield(s, responsivenessField) || ~s.(responsivenessField)
        continue;
    end
    if ~isfield(s, 'stimValues')
        continue;
    end

    thetaDeg = s.stimValues(:)';
    thetaRad = deg2rad(thetaDeg);
    R = s.meanDirResponse(:)';
    R_rect = max(R, 0); % half-wave rectify

    if sum(R_rect) == 0 || all(isnan(R_rect))
        continue;
    end

    vectorSum = sum(R_rect .* exp(1i * thetaRad)); % NOTE: theta, not 2*theta
    DSI(b) = abs(vectorSum) / sum(R_rect);

    prefDirectionDeg(b) = mod(rad2deg(angle(vectorSum)), 360);
end

%% ===================== simple ratio method =====================
DSI_simple = nan(nBoutonsTotal, 1);

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
    R_rect = max(R, 0);

    nDir = numel(thetaDeg);
    [~, sortIdx] = sort(thetaDeg);
    thetaSorted = thetaDeg(sortIdx);
    Rsorted = R_rect(sortIdx);

    % check directions are evenly spaced (e.g. 45 deg apart for 8 directions)
    spacing = mode(diff(thetaSorted));
    stepsFor180 = round(180 / spacing);
    if abs(stepsFor180 * spacing - 180) > 1e-6 || mod(nDir, 2) ~= 0
        continue; % directions don't evenly support a true 180-deg opposite -- skip rather than guess
    end

    [Rpref, prefIdx] = max(Rsorted);
    oppIdx = mod(prefIdx - 1 + stepsFor180, nDir) + 1;
    Ropp = Rsorted(oppIdx);

    if (Rpref + Ropp) == 0
        continue;
    end
    DSI_simple(b) = (Rpref - Ropp) / (Rpref + Ropp);
end

%% ===================== attach results =====================
for b = 1:nBoutonsTotal
    allDirTuning(b).DSI              = DSI(b);
    allDirTuning(b).prefDirectionDeg = prefDirectionDeg(b);
    allDirTuning(b).DSI_simple       = DSI_simple(b);
end

nResponsive  = sum(arrayfun(@(s) isfield(s, responsivenessField) && s.(responsivenessField), allDirTuning));
nValidDSI    = sum(~isnan(DSI));
nValidSimple = sum(~isnan(DSI_simple));
fprintf('%d / %d boutons passed the "%s" responsiveness criterion.\n', nResponsive, nBoutonsTotal, responsivenessField);
fprintf('%d boutons have a valid vector-method DSI.\n', nValidDSI);
fprintf('%d boutons have a valid simple-ratio DSI.\n', nValidSimple);

%% ===================== plots =====================
if makePlots
    figure('Position', [100 100 1300 400]);
    subplot(1,3,1);
    histogram(DSI(~isnan(DSI)), 20);
    xlabel('DSI (vector method)'); ylabel('Number of boutons');
    title(sprintf('Vector-method DSI (n=%d)', nValidDSI));

    subplot(1,3,2);
    histogram(DSI_simple(~isnan(DSI_simple)), 20);
    xlabel('DSI (simple ratio)'); ylabel('Number of boutons');
    title(sprintf('Simple-ratio DSI (n=%d)', nValidSimple));

    subplot(1,3,3);
    bothValid = ~isnan(DSI) & ~isnan(DSI_simple);
    scatter(DSI_simple(bothValid), DSI(bothValid), 15, 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    plot([0 1], [0 1], 'k--');
    xlabel('DSI (simple ratio)'); ylabel('DSI (vector method)');
    title(sprintf('Method agreement (n=%d)', sum(bothValid)));
    axis square; xlim([0 1]); ylim([0 1]);

    sgtitle(sprintf('DSI: vector method vs simple ratio (criterion: %s)', responsivenessField), ...
        'Interpreter', 'none');

    figure;
    polarhistogram(deg2rad(prefDirectionDeg(~isnan(prefDirectionDeg))), 16); % NOT doubled -- direction spans full 360
    title('Preferred direction distribution (vector method)');
end

end