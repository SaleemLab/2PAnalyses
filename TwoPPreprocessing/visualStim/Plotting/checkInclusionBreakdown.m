% printClearSummary.m
%
% Run AFTER the main fitting script has populated gaussFitResults. No refitting needed.
%
% Prints, in plain text:
%   - Total number of fits (above screening floor)
%   - How many pass ALL criteria safely (R^2 + bootstrap-CI + non-degenerate-sigma = isTrusted)
%   - Of the ones that do NOT pass, how many fail ONLY because of broad sigma
%   - Of the ones that do NOT pass, how many fail ONLY because of a beyond-range center
%   - Mean +/- SD sigma for each of these groups

%% params
r2Thresh = 0.2;
ciLevel  = 0.95;
azBroadCutoff = 35;
elBroadCutoff = 20;
lowerPct = 100 * (1 - ciLevel);

%% gather
R2vals   = [gaussFitResults.R2];
sigAzAll = [gaussFitResults.sigmaX];
sigElAll = [gaussFitResults.sigmaY];
isBeyondAll = [gaussFitResults.isBeyondRange];
isDegenAll  = [gaussFitResults.isDegenerateSigma];
isLargeAll  = (sigAzAll > azBroadCutoff) | (sigElAll > elBroadCutoff);

isRobust = false(size(gaussFitResults));
for gi = 1:numel(gaussFitResults)
    validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));
    if numel(validBoot) >= 20
        isRobust(gi) = prctile(validBoot, lowerPct) >= r2Thresh;
    end
end

passR2CI  = (R2vals >= r2Thresh) & isRobust;
isTrusted = passR2CI & ~isDegenAll;   % your actual criterion (R^2 + CI + non-degenerate-sigma)

nTotal      = numel(gaussFitResults);
nTrusted    = sum(isTrusted);
nNotTrusted = nTotal - nTrusted;

% Among the ones NOT trusted, how many are held back ONLY by broad sigma, ONLY by beyond-range,
% both, or neither (failed for some other reason, e.g. low R^2/CI or degenerate-sigma itself).
notTrustedIdx  = ~isTrusted;
failOnlyBroad  = notTrustedIdx & isLargeAll  & ~isBeyondAll;
failOnlyBeyond = notTrustedIdx & isBeyondAll & ~isLargeAll;
failBoth       = notTrustedIdx & isLargeAll  & isBeyondAll;
failOther      = notTrustedIdx & ~isLargeAll & ~isBeyondAll;

%% ===================== PRINT =====================
fprintf('\n========================================\n');
fprintf('TOTAL FITS (above screening floor)         : n = %d\n', nTotal);
fprintf('  sigAz = %.2f +/- %.2f | sigEl = %.2f +/- %.2f\n\n', ...
    mean(sigAzAll), std(sigAzAll), mean(sigElAll), std(sigElAll));

fprintf('PASS ALL CRITERIA SAFELY (isTrusted)        : n = %d\n', nTrusted);
fprintf('  sigAz = %.2f +/- %.2f | sigEl = %.2f +/- %.2f\n\n', ...
    mean(sigAzAll(isTrusted)), std(sigAzAll(isTrusted)), mean(sigElAll(isTrusted)), std(sigElAll(isTrusted)));

fprintf('NOT TRUSTED (fail R^2/CI/degenerate-sigma)  : n = %d\n\n', nNotTrusted);

fprintf('  Of the %d not trusted:\n', nNotTrusted);
fprintf('    ONLY broad sigma (az>%d or el>%d) is the issue   : n = %d\n', azBroadCutoff, elBroadCutoff, sum(failOnlyBroad));
if sum(failOnlyBroad) > 0
    fprintf('      sigAz = %.2f +/- %.2f | sigEl = %.2f +/- %.2f\n', ...
        mean(sigAzAll(failOnlyBroad)), std(sigAzAll(failOnlyBroad)), mean(sigElAll(failOnlyBroad)), std(sigElAll(failOnlyBroad)));
end
fprintf('    ONLY beyond-range center is the issue           : n = %d\n', sum(failOnlyBeyond));
if sum(failOnlyBeyond) > 0
    fprintf('      sigAz = %.2f +/- %.2f | sigEl = %.2f +/- %.2f\n', ...
        mean(sigAzAll(failOnlyBeyond)), std(sigAzAll(failOnlyBeyond)), mean(sigElAll(failOnlyBeyond)), std(sigElAll(failOnlyBeyond)));
end
fprintf('    BOTH broad sigma AND beyond-range               : n = %d\n', sum(failBoth));
fprintf('    Neither (failed R^2/CI/degenerate for other reasons) : n = %d\n', sum(failOther));
fprintf('========================================\n');