function S = sk_load_session(iSess, T, mode)
% Load the binned export of one session as firing-rate arrays (Hz), smoothed
% along position, and attach the matching rows of the cell table T.
%   S.D, S.V, S.F : [nCells x nTrials x nBins] data, vision+speed prediction,
%                   full-model prediction (rates, smoothed, NaN where unvisited)
%   S.speed       : [nTrials x nBins] mean running speed
%   S.cond        : [nTrials x 1] condition code
%   S.T           : table rows for these cells, in Spk column order
% mode (optional) replaces the data with simulated spikes, for calibration:
%   'real'    recorded spikes (default)
%   'simAdd'  Poisson spikes from the stored full GLM (vision+speed + additive position)
%   'simGain' Poisson spikes from vision+speed x a base-condition gain g(x)
if nargin < 3, mode = 'real'; end

C = sk_config();
mouse = C.sessions{iSess,1}; date = C.sessions{iSess,2}; sid = C.sessions{iSess,3};
f = fullfile(C.outRoot, sprintf('%s_%d_binned.mat', mouse, date));
L = load(f);
occ = reshape(L.occ, [1 size(L.occ)]);
toRate = @(M) double(M) ./ occ * L.meta.sampleRate;      % counts per bin -> Hz
S.V = sk_smooth(toRate(L.predVS));
S.F = sk_smooth(toRate(L.predFull));
switch mode
    case 'real'
        S.D = sk_smooth(toRate(L.spk));
    case 'simAdd'
        lam = max(double(L.predFull), 0);
        S.D = sk_smooth(toRate(poissrnd(lam)));
    case 'simGain'
        Dreal = sk_smooth(toRate(L.spk));
        b = L.trialCond == 1;
        Db = mean(Dreal(:,b,:),2,'omitnan'); Vb = mean(S.V(:,b,:),2,'omitnan');
        [~, ~, gp] = sk_gain(Db, Vb);
        rate = max(gp(S.V), 0);                          % Hz, every trial
        rate(isnan(rate)) = 0;
        lam = rate .* occ / L.meta.sampleRate;
        S.D = sk_smooth(toRate(poissrnd(lam)));
end
S.mode = mode;
S.speed = L.speedTB;
S.cond  = L.trialCond;
S.trialSpeed = L.trialSpeed;
S.tag = sprintf('%s_%d', mouse, date);
S.sid = sid;

rows = T(T.sessionid == sid, :);
rows = sortrows(rows, 'cellNum');
assert(height(rows) == size(S.D,1) && all(rows.cellNum == (1:height(rows))'), ...
    'table rows do not match Spk cell order for %s', S.tag);
S.T = rows;
end

function Y = sk_smooth(X)
% Gaussian smoothing along dim 3 (position) ignoring NaNs; sigma = 1.5 bins (3 cm)
sig = 1.5; h = exp(-0.5*((-5:5)/sig).^2); h = reshape(h/sum(h), 1, 1, []);
m = ~isnan(X); X(~m) = 0;
num = convn(X, h, 'same'); den = convn(double(m), h, 'same');
Y = num ./ den; Y(den < 0.2) = NaN;
end
