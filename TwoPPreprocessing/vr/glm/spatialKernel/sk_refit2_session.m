function [Rf, Ex] = sk_refit2_session(iSess, cellSel, opts)
% Refit v2 (Poisson / log link only). Changes from sk_refit_session:
%   1. Omit is NOT part of the vision+speed baseline. Position is tested
%      against vision(no omit)+speed; the omit block is added afterwards,
%      on top of the fitted model (offset fixed), so it can only explain what
%      vision, speed and position could not.
%   2. Corridor onset/offset are long temporal kernels: raised-cosine basis
%      (log-stretched, opts.nOnOffBasis functions over 0-opts.onOffSpan s, default 8 over 1.5 s)
%      on time since corridor onset and time since offset (blank), replacing
%      the original 25-lag (~0.4 s) XvisualONOFF block.
%   3. Two penalties: ridge and elastic net (alpha = opts.alpha, glmnet-style
%      standardised penalty, active-set solver). Lambda chosen by CV
%      separately for each model; the elastic-net path is stopped early once
%      held-out LL has fallen for opts.patience lambdas in a row.
%   opts.models picks which of VS, VSP, VSI, VSIP to fit (default [3 4]:
%   only the vision x speed models); the others are left NaN.
% Same sample subset, original CV folds and per-cell visual delay as before.
%
% Blocks: 1 OnOff (temporal basis), 2 RelSpd, 3 Tex, 4 BG, 5 End, 6 Omit,
%         7 Space, 8 vision x speed (OnOff, Tex, BG, End columns x speed z)
% Models: VS [1-5], VSP [1-5 7], VSI [1-5 8], VSIP [1-5 7 8];
%         VSP+O, VSIP+O: omit block fitted after VSP / VSIP.
%
% Rf columns per penalty pen in {ridge, en}:
%   LLHi_<pen>_<model>   CV log-likelihood gain over constant rate, bits per unit activity
%   lam_<pen>_<model>
%   LLHrel_<pen> (VSP/VS), LLHrelI_<pen> (VSIP/VSI)
%   p_<pen>_pos, p_<pen>_posI   one-sided paired t over the 10 folds
%   LR_<pen>_pos(I), pLR_*      likelihood-ratio test as in Julien's RunGLMs
%                               (2 x summed held-out LL gain, chi2, dof = nPosBins)
%   phi_<pen>_pos(I), pLRq_*    same divided by the Pearson dispersion (quasi-
%                               likelihood; needed because activity is deconvolved)
%   LLHi_<pen>_VSPO, LLHi_<pen>_VSIPO, p_<pen>_omitI (omit adds to VSIP)
% Ex: foldLL [cells x K x models x pens], kernels.<pen>_VSP/_VSIP [cells x 100]
%     (log gain), onoff.<pen>_VSIP [cells x nBasis] coefficients and the basis.

if nargin < 2, cellSel = []; end
if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'kfold', 10);
opts = setDefault(opts, 'lamRidge', [1e2 1e3 1e4]);
opts = setDefault(opts, 'lamOmit', [1e1 1e2 1e3 1e4 1e5 1e6]);
opts = setDefault(opts, 'alpha', 0.9);
opts = setDefault(opts, 'nLamEN', 10);
opts = setDefault(opts, 'lamRatioEN', 1e-3);
opts = setDefault(opts, 'onOffSpan', 1.5);     % s (covers visual + indicator delays)
opts = setDefault(opts, 'nOnOffBasis', 8);
opts = setDefault(opts, 'useParfor', true);
opts = setDefault(opts, 'models', [3 4]);       % VSI, VSIP; [1 2 3 4] adds the additive VS, VSP
opts = setDefault(opts, 'patience', 2);         % stop the lambda path after 2 worse lambdas
pens = {'ridge','en'};

C = sk_config();
mouse = C.sessions{iSess,1}; date = C.sessions{iSess,2}; sid = C.sessions{iSess,3};
tag = sprintf('%s_%d', mouse, date);
L = load(fullfile(C.outRoot, [tag '_samples.mat']));
fit = L.fit; keep = fit.inGLM; bi0 = fit.blockIdx; fs = L.meta.sampleRate;
Xfull = L.X; Yall = double(L.Y(keep,:));
trial = fit.trialID(keep);
nC = size(Yall,2); if isempty(cellSel), cellSel = 1:nC; end
cellSel = cellSel(:)';

delayMs = sk_read_best_delay(L.meta.file);
delaySmp = round(delayMs(cellSel)' / 1000 * fs);
visAll0 = ismember(bi0, [1 3 4 5 6]);
onoffCols = find(bi0 == 1);
assert(numel(onoffCols) == 50, 'expected 25 onset + 25 offset lag columns');

[foldOfTrial, ut, exact] = sk_recover_folds(iSess);
if ~exact, warning('%s: original folds not recovered exactly', tag); end
[~, ti] = ismember(trial, ut);
fold = foldOfTrial(ti);

% raised-cosine basis on lag (samples), log-stretched (Pillow et al. 2008)
maxLag = round(opts.onOffSpan * fs);
lagAx = (0:maxLag)';
Bas = rcBasis(lagAx, opts.nOnOffBasis, 3);      % [maxLag+1 x nBasis]

models = struct('name', {'VS','VSP','VSI','VSIP'}, 'blocks', {1:5, [1:5 7], [1:5 8], [1:5 7 8]});
nM = numel(models); nP = numel(pens); K = opts.kfold;
nSel = numel(cellSel);
nSpk = sum(Yall(:,cellSel),1)';
nPosBins = nnz(bi0 == 7);

dlys = unique(delaySmp);
GD = cell(numel(dlys),1);
grpOf = zeros(nSel,1); colOf = zeros(nSel,1);
for gi = 1:numel(dlys)
    dly = dlys(gi);
    grp = find(delaySmp == dly);
    grpOf(grp) = gi; colOf(grp) = 1:numel(grp);
    X = Xfull;
    if dly > 0, X(:,visAll0) = [zeros(dly, nnz(visAll0), 'uint8'); Xfull(1:end-dly, visAll0)]; end

    % onset / offset events (lag-0 columns of the original block, after the delay shift)
    onLag  = lagSince(X(:, onoffCols(1))  > 0, maxLag);
    offLag = lagSince(X(:, onoffCols(26)) > 0, maxLag);
    XonOff = [basisAt(onLag, Bas), basisAt(offLag, Bas)];
    X = X(keep,:); XonOff = XonOff(keep,:);

    % rebuild blocks: 1 OnOff basis, 2-7 original
    keepOrig = bi0 ~= 1;
    Xb = [XonOff, double(X(:, keepOrig))];
    bb = [ones(size(XonOff,2),1); bi0(keepOrig)];

    % vision x speed (not omit): OnOff, Tex, BG, End columns x speed-bin z
    spdCols = find(bb == 2);
    [~, spdBin] = max(Xb(:,spdCols), [], 2); spdBin(~any(Xb(:,spdCols),2)) = NaN;
    sz = (spdBin - mean(spdBin,'omitnan')) / std(spdBin,'omitnan'); sz(isnan(sz)) = 0;
    visCols = find(ismember(bb, [1 3 4 5]));
    Xall = [Xb, Xb(:,visCols) .* sz]; bAll = [bb; 8*ones(numel(visCols),1)];
    keepCol = any(Xall ~= 0, 1)';
    posBin = find(keepCol(bb == 7));
    Xall = Xall(:, keepCol); bAll = bAll(keepCol);

    [U, ~, g] = unique([fold, Xall], 'rows');
    foldU = U(:,1); XU = sparse(U(:,2:end)); nU = accumarray(g, 1);
    fprintf('[%s] delay %d samples: %d cells, %d rows -> %d aggregated, %d columns\n', ...
        tag, dly, numel(grp), numel(g), size(U,1), size(XU,2));
    Y = Yall(:, cellSel(grp));
    SY = zeros(size(U,1), numel(grp)); SY2 = SY;
    for c = 1:numel(grp)
        SY(:,c)  = accumarray(g, Y(:,c));
        SY2(:,c) = accumarray(g, Y(:,c).^2);     % for the Pearson dispersion
    end
    inM = cell(nM,1); posIdx = cell(nM,1);
    for m = 1:nM
        inM{m} = ismember(bAll, models(m).blocks);
        posIdx{m} = 1 + find(bAll(inM{m}) == 7);
    end
    onIdx = 1 + find(bAll(inM{4}) == 1);
    GD{gi} = struct('XU',XU,'inM',{inM},'inO',bAll == 6,'posIdx',{posIdx},'onIdx',onIdx, ...
        'posBin',posBin,'nU',nU,'foldU',foldU,'SY',SY,'SY2',SY2);
end

nOut = nM + 2;                                   % + VSP+O, VSIP+O
resAll = nan(nSel, nOut, nP); lbAll = resAll; foldAll = nan(nSel, K, nOut, nP);
kP = nan(nSel, nPosBins, 2, nP); kOn = nan(nSel, 2*opts.nOnOffBasis, nP);
t0 = tic; tCell = nan(nSel,1); phiAll = nan(nSel, nM, nP);
if opts.useParfor
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(v) fprintf('[%s] cell %d done in %.0f s (%s)\n', tag, v(1), v(2), datestr(now, 'HH:MM:SS')));
    parfor c = 1:nSel
        [resAll(c,:,:), lbAll(c,:,:), foldAll(c,:,:,:), kP(c,:,:,:), kOn(c,:,:), tCell(c), phiAll(c,:,:)] = ...
            fitCellAll(GD{grpOf(c)}, colOf(c), K, opts, nM, nPosBins); %#ok<PFBNS>
        send(dq, [cellSel(c) tCell(c)]);
    end
else
    for c = 1:nSel
        [resAll(c,:,:), lbAll(c,:,:), foldAll(c,:,:,:), kP(c,:,:,:), kOn(c,:,:), tCell(c), phiAll(c,:,:)] = ...
            fitCellAll(GD{grpOf(c)}, colOf(c), K, opts, nM, nPosBins);
        fprintf('[%s] cell %d done in %.0f s\n', tag, cellSel(c), tCell(c));
    end
end
fprintf('[%s] %d cells fitted in %.0f s (median %.0f s per cell per worker)\n', tag, nSel, toc(t0), median(tCell));
LLHi = resAll ./ nSpk / log(2);

outNames = [{models.name}, {'VSPO','VSIPO'}];
nPosFit = cellfun(@(g) numel(g.posBin), GD); nPosFit = nPosFit(grpOf);   % dof of the LR test
Rf = table(repmat(sid,nSel,1), cellSel(:), nSpk, delaySmp(:), nPosFit, ...
    'VariableNames', {'sessionid','cellNum','nSpikes','visDelaySamples','nPosBins'});
for p = 1:nP
    for m = 1:nOut
        Rf.(sprintf('LLHi_%s_%s', pens{p}, outNames{m})) = LLHi(:,m,p);
        Rf.(sprintf('lam_%s_%s',  pens{p}, outNames{m})) = lbAll(:,m,p);
    end
    Rf.(sprintf('LLHrel_%s', pens{p}))  = LLHi(:,2,p) ./ LLHi(:,1,p);
    Rf.(sprintf('LLHrelI_%s', pens{p})) = LLHi(:,4,p) ./ LLHi(:,3,p);
    Rf.(sprintf('p_%s_pos', pens{p}))   = foldTest(foldAll(:,:,2,p) - foldAll(:,:,1,p));
    Rf.(sprintf('p_%s_posI', pens{p}))  = foldTest(foldAll(:,:,4,p) - foldAll(:,:,3,p));
    Rf.(sprintf('p_%s_omitI', pens{p})) = foldTest(foldAll(:,:,6,p) - foldAll(:,:,4,p));
    % Likelihood-ratio test as in Julien's RunGLMs: summed held-out LL of the
    % model with vs without position, chi2 with dof = number of position
    % bins. Raw (as for true Poisson counts) and divided by the Pearson
    % dispersion of the model with position (quasi-likelihood; deconvolved
    % activity is not in count units, so the raw test is not calibrated).
    for q = 1:2                                   % q=1: VSP vs VS, q=2: VSIP vs VSI
        sfx = ternary(q == 1, 'pos', 'posI'); mB = 2*q - 1; mP = 2*q;
        st = 2 * sum(foldAll(:,:,mP,p) - foldAll(:,:,mB,p), 2);
        Rf.(sprintf('LR_%s_%s', pens{p}, sfx))     = st;
        Rf.(sprintf('phi_%s_%s', pens{p}, sfx))    = phiAll(:,mP,p);
        Rf.(sprintf('pLR_%s_%s', pens{p}, sfx))    = chi2cdf(max(st,0), nPosFit, 'upper');
        Rf.(sprintf('pLRq_%s_%s', pens{p}, sfx))   = chi2cdf(max(st,0) ./ phiAll(:,mP,p), nPosFit, 'upper');
    end
end
Ex = struct('sessionid', sid, 'cellNum', cellSel(:), 'foldLL', foldAll, 'models', {outNames}, ...
    'pens', {pens}, 'binCentres', C.binCentres, 'onOffBasis', Bas, 'onOffLagS', lagAx/fs);
for p = 1:nP
    Ex.kernels.([pens{p} '_VSP'])  = kP(:,:,1,p);
    Ex.kernels.([pens{p} '_VSIP']) = kP(:,:,2,p);
    Ex.onoff.([pens{p} '_VSIP'])   = kOn(:,:,p);
end
end

% ------------------------------------------------------------------------
function lag = lagSince(ev, maxLag)
% samples since the most recent event (0 at the event); NaN if > maxLag or none yet
n = numel(ev); lag = nan(n,1); last = -Inf;
for i = 1:n
    if ev(i), last = i; end
    d = i - last; if d <= maxLag, lag(i) = d; end
end
end

function B = basisAt(lag, Bas)
B = zeros(numel(lag), size(Bas,2));
ok = ~isnan(lag); B(ok,:) = Bas(lag(ok)+1, :);
B = round(B, 6);                                 % exact repeats -> rows still aggregate
end

function B = rcBasis(t, nB, c)
% log-stretched raised cosines tiling t (samples) from 0 to max(t)
y = log(t + c); yr = [log(c), log(max(t) + c)];
w = diff(yr) / (nB - 1);                         % centre spacing
ctr = yr(1) + (0:nB-1) * w;
B = zeros(numel(t), nB);
for k = 1:nB
    a = (y - ctr(k)) * pi / (2*w);
    B(:,k) = (abs(a) < pi) .* 0.5 .* (1 + cos(a));
end
end

function p = foldTest(dF)
n = sum(~isnan(dF),2); m = mean(dF,2,'omitnan'); s = std(dF,0,2,'omitnan');
p = 1 - tcdf(m ./ (s ./ sqrt(n)), n-1);
end

function [res, lb, foldLL, kP, kOn, tCell, phi] = fitCellAll(G, j, K, opts, nM, nPosBins)
t0 = tic;
nOut = nM + 2; res = nan(1, nOut, 2); lb = res; foldLL = nan(1, K, nOut, 2);
kP = nan(1, nPosBins, 2, 2); kOn = nan(1, numel(G.onIdx), 2); phi = nan(1, nM, 2);
sy = G.SY(:,j); sy2 = G.SY2(:,j); N = sum(G.nU);
Xm = cell(nM,1);
for m = opts.models, Xm{m} = [ones(size(G.XU,1),1), G.XU(:, G.inM{m})]; end
XO = G.XU(:, G.inO);
for p = 1:2
    pen = struct('type', ternary(p == 1, 'ridge', 'en'), 'alpha', opts.alpha, ...
        'nLam', opts.nLamEN, 'ratio', opts.lamRatioEN);
    B = cell(nM,1); bAll = cell(nM,1);
    for m = opts.models     % lambda chosen by CV separately for each model
        if p == 1, lams = opts.lamRidge; else, lams = []; end
        [res(1,m,p), lb(1,m,p), foldLL(1,:,m,p), B{m}, bAll{m}] = ...
            cvPath(Xm{m}, G.nU, sy, G.foldU, K, lams, pen, opts.patience);
        % Pearson dispersion of the all-data fit: sum over samples of
        % (y - mu)^2 / mu, from the per-row sums of y and y^2, over N - df
        mu = exp(Xm{m} * bAll{m});
        pear = sum((sy2 - 2*mu.*sy + G.nU.*mu.^2) ./ mu);
        phi(1,m,p) = pear / (N - nnz(bAll{m}));
    end
    % omit fitted after VSP (model 2) and VSIP (model 4): fold betas fixed as offset
    for q = 1:2
        m = 2*q; o = nM + q;
        if ~ismember(m, opts.models), continue; end
        [res(1,o,p), lb(1,o,p), foldLL(1,:,o,p)] = cvOmitAfter(Xm{m}, XO, B{m}, G.nU, sy, G.foldU, K, opts.lamOmit);
        % kernels: all-data fit at the chosen lambda (fitted alongside the folds)
        kP(1, G.posBin, q, p) = bAll{m}(G.posIdx{m});
        if m == 4, kOn(1,:,p) = bAll{m}(G.onIdx); end
    end
end
tCell = toc(t0);
end

function o = ternary(c, a, b), if c, o = a; else, o = b; end, end

function [dLL, lamBest, foldLL, Bbest, bAllBest] = cvPath(Xm, n, sy, fold, K, lams, pen, patience)
% CV log-likelihood gain over a constant-rate model at the best lambda.
% Lambda is the outer loop: all K folds and the all-data fit step down the
% path together (each warm-started from its previous lambda), and the path is
% abandoned once the summed held-out LL has fallen for `patience` lambdas in a
% row. Returns per-fold betas at the best lambda (for fitting omit afterwards)
% and the all-data betas at the best lambda (kernels).
if isempty(lams), lams = enLambdaPath(Xm, n, sy, pen); end
lams = sort(lams, 'descend'); nL = numel(lams);
% fold 0 = all data (no test set)
S = cell(K+1,1);
for k = 0:K
    if k == 0, tr = true(size(fold)); else, tr = fold ~= k; end
    te = ~tr;
    if k > 0 && ~any(te), continue; end
    s.X = Xm(tr,:); s.n = n(tr); s.sy = sy(tr);
    s.Xte = Xm(te,:); s.nte = n(te); s.syte = sy(te);
    s.b = zeros(size(Xm,2),1); s.b(1) = log(max(sum(s.sy)/sum(s.n), 1e-12));
    if strcmp(pen.type, 'en'), [s.w1, s.w2] = stdWeights(s.X, s.n); end
    s.ll0 = poisLL(repmat(s.b(1), nnz(te), 1), s.nte, s.syte);
    S{k+1} = s;
end
ll = nan(K, nL); Bpath = cell(K+1, nL);
best = -Inf; ib = 1; worse = 0;
for il = 1:nL
    for k = 0:K
        s = S{k+1}; if isempty(s), continue; end
        switch pen.type
            case 'ridge', s.b = fitRidge(s.X, s.n, s.sy, lams(il), s.b);
            case 'en',    s.b = fitEN(s.X, s.n, s.sy, lams(il), pen.alpha, s.w1, s.w2, s.b);
        end
        S{k+1}.b = s.b; Bpath{k+1, il} = s.b;
        if k > 0, ll(k,il) = poisLL(s.Xte*s.b, s.nte, s.syte); end
    end
    tot = sum(ll(:,il), 'omitnan');
    if tot > best, best = tot; ib = il; worse = 0; else, worse = worse + 1; end
    if worse >= patience, break; end
end
ll0 = nan(K,1);
for k = 1:K, if ~isempty(S{k+1}), ll0(k) = S{k+1}.ll0; end, end
dLL = best - sum(ll0, 'omitnan'); lamBest = lams(ib);
foldLL = (ll(:,ib) - ll0)';
Bbest = Bpath(2:end, ib);
bAllBest = Bpath{1, ib};
end

function [dLL, lamBest, foldLL] = cvOmitAfter(Xm, XO, Bf, n, sy, fold, K, lams)
% Omit block fitted on top of the fold's base fit (fixed offset), ridge.
nL = numel(lams); ll = nan(K, nL); ll0 = nan(K,1);
for k = 1:K
    tr = fold ~= k; te = ~tr;
    if ~any(te) || isempty(Bf{k}), continue; end
    offTr = Xm(tr,:) * Bf{k}; offTe = Xm(te,:) * Bf{k};
    Otr = XO(tr,:); Ote = XO(te,:);
    act = any(Otr ~= 0, 1);                      % omit columns seen in training
    bo = zeros(nnz(act),1);
    for il = 1:nL
        bo = fitOffsetRidge(Otr(:,act), n(tr), sy(tr), offTr, lams(il), bo);
        ll(k,il) = poisLL(offTe + Ote(:,act)*bo, n(te), sy(te));
    end
    b0 = log(max(sum(sy(tr))/sum(n(tr)), 1e-12));
    ll0(k) = poisLL(repmat(b0, nnz(te), 1), n(te), sy(te));
end
[best, ib] = max(sum(ll,1,'omitnan'));
dLL = best - sum(ll0,'omitnan'); lamBest = lams(ib);
foldLL = (ll(:,ib) - ll0)';
end

function b = fitOffsetRidge(X, n, sy, off, lam, b)
% ridge Poisson with offset, no intercept (the base model has one)
if isempty(X), b = zeros(0,1); return; end
p = size(X,2);
f = @(bb) -(sy' * (off + X*bb) - n' * exp(off + X*bb)) + 0.5*lam*sum(bb.^2);
fb = f(b);
for it = 1:50
    mu = n .* exp(off + X*b);
    gr = full(X' * (mu - sy)) + lam*b;
    H = full(X' * (X .* mu)); H(1:p+1:end) = H(1:p+1:end) + lam + 1e-8;
    step = cholSolve(H, gr); t = 1;
    while t > 1e-4
        bn = b - t*step; fn = f(bn);
        if fn <= fb, break; end
        t = t/2;
    end
    if fb - fn < 1e-6 * max(1, abs(fb)), b = bn; break; end
    b = bn; fb = fn;
end
end

function ll = poisLL(eta, n, sy)
ll = sy' * eta - n' * exp(eta);                  % log(y!) omitted: same for every model
end

% ---- fitting: ridge (Newton) or elastic net (proximal Newton + coordinate descent)
function lams = enLambdaPath(X, n, sy, pen)
N = sum(n); [w1, ~] = stdWeights(X, n);
b0 = log(max(sum(sy)/N, 1e-12));
g = full(X' * (n .* exp(b0) - sy)) / N;
ok = w1 > 0;
lmax = max(abs(g(ok)) ./ (pen.alpha * w1(ok)));
lams = lmax * logspace(0, log10(pen.ratio), pen.nLam);
end

function [w1, w2] = stdWeights(X, n)
% glmnet-style standardisation as penalty weights: |b_j| sd_j and b_j^2 sd_j^2
N = sum(n);
m1 = full(X' * n) / N; m2 = full((X.^2)' * n) / N;
sd = sqrt(max(m2 - m1.^2, 0)); sd(1) = 0;       % intercept unpenalised
w1 = sd; w2 = sd.^2;
end

function b = fitRidge(X, n, sy, lam, b)
p = size(X,2); pen = lam*ones(p,1); pen(1) = 0;
f = @(bb) -(sy' * (X*bb) - n' * exp(X*bb)) + 0.5*lam*sum(bb(2:end).^2);
fb = f(b);
for it = 1:50
    mu = n .* exp(X*b);
    gr = full(X' * (mu - sy)) + pen .* b;
    H = full(X' * (X .* mu)); H(1:p+1:end) = H(1:p+1:end) + pen' + 1e-8;
    step = cholSolve(H, gr); t = 1;
    while t > 1e-4
        bn = b - t*step; fn = f(bn);
        if fn <= fb, break; end
        t = t/2;
    end
    if fb - fn < 1e-5 * max(1, abs(fb)), b = bn; break; end
    b = bn; fb = fn;
end
end

function b = fitEN(X, n, sy, lam, alpha, w1, w2, b)
% minimise (1/N)[n'exp(Xb) - sy'Xb] + lam*sum(alpha*w1|b| + (1-alpha)/2*w2*b.^2)
% Active set (as glmnet): solve over the nonzero coefficients plus those that
% violate the optimality conditions, then re-check all columns (KKT) and add
% violators until none remain. Same optimum, much smaller Hessians when the
% penalty is strong.
N = sum(n);
l1 = lam*alpha*w1; l2 = lam*(1-alpha)*w2;
g = full(X' * (n .* exp(X*b) - sy)) / N;
A = b ~= 0 | abs(g) > l1 | l1 == 0;
for outer = 1:50
    b(~A) = 0;
    b(A) = fitENactive(X(:,A), n, sy, l1(A), l2(A), b(A));
    g = full(X' * (n .* exp(X*b) - sy)) / N;
    viol = ~A & abs(g) > l1 * (1 + 1e-4) + 1e-10;
    if ~any(viol), break; end
    A = A | viol;
end
end

function b = fitENactive(X, n, sy, l1, l2, b)
N = sum(n); p = size(X,2);
F = @(bb) (n' * exp(X*bb) - sy' * (X*bb)) / N + sum(l1.*abs(bb)) + 0.5*sum(l2.*bb.^2);
fb = F(b);
% The Hessian is the costly part (sparse X'WX). It is formed once per call and
% reused (chord Newton) while full steps keep being accepted; it is refreshed
% only when a step has to be shortened or fails. Warm starts along the lambda
% path keep the stale Hessian accurate.
needH = true;
for it = 1:30
    mu = n .* exp(X*b);
    g = full(X' * (mu - sy)) / N;
    if needH, H = full(X' * (X .* mu)) / N; hjj = diag(H); fresh = true; else, fresh = false; end
    needH = false;
    d = zeros(p,1); Hd = zeros(p,1);
    for sweep = 1:200                             % coordinate descent on the quadratic model
        maxCh = 0;
        for jj = 1:p
            if hjj(jj) + l2(jj) <= 0, continue; end
            u = hjj(jj)*(b(jj) + d(jj)) - (g(jj) + Hd(jj));
            z = sign(u) * max(abs(u) - l1(jj), 0) / (hjj(jj) + l2(jj));
            dz = z - (b(jj) + d(jj));
            if dz ~= 0
                d(jj) = d(jj) + dz; Hd = Hd + H(:,jj)*dz;
                maxCh = max(maxCh, abs(dz) * sqrt(hjj(jj)));
            end
        end
        if maxCh < 1e-7, break; end
    end
    t = 1;
    while t > 1e-4
        bn = b + t*d; fn = F(bn);
        if fn <= fb, break; end
        t = t/2;
    end
    if fn > fb
        if fresh, break; end                      % no descent even with a fresh Hessian
        needH = true; continue;                   % stale Hessian: refresh and retry
    end
    conv = fb - fn < 1e-7 * max(1, abs(fb));
    b = bn; fb = fn;
    if conv, break; end
    needH = t < 1;
end
end

function x = cholSolve(A, y)
[R, flag] = chol(A);
if flag == 0, x = R \ (R' \ y); else, x = A \ y; end
end

function o = setDefault(o, f, v), if ~isfield(o,f), o.(f) = v; end, end
