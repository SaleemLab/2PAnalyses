function [Rf, Ex] = sk_refit_session(iSess, cellSel, opts)
% Refit the GLM for one session with Gaussian (identity) and Poisson (log)
% likelihoods, with and without position and with and without an explicit
% vision x speed interaction. Matches the original GLM1 where it can: same
% sample subset, the original 10 CV folds (recovered exactly, sk_recover_folds),
% and each cell's original visual delay (sk_read_best_delay). Differences:
% ridge instead of lasso, and the RF-position option is not reproduced.
%
% Models (blocks as in GLM1: 1 ONOFF, 2 RelSpd, 3 Tex, 4 BG, 5 End, 6 Omit, 7 Space):
%   VS      blocks 1-6
%   VSP     blocks 1-7                      (the original "full" model)
%   VSI     blocks 1-6 + vision x speed     (each visual column times centred speed-bin index)
%   VSIP    blocks 1-7 + vision x speed
% Fitting: ridge-penalised maximum likelihood (Newton/IRLS). Gaussian: lambda
% chosen by CV per model. Poisson: lambda chosen by CV on VSP and shared by all
% four models. Samples are aggregated over identical design rows within each
% fold (exact for both likelihoods, since all predictors are discrete).
% opts.Y replaces the recorded activity (e.g. simulated cells from sk_simulate_null).
%
% Rf (table, one row per selected cell, Spk column order):
%   LLHi_<dist>_<model> : CV log-likelihood gain over a constant-rate model, bits/spike
%   LLHrel_<dist>       : LLHi_VSP / LLHi_VS   (same definition as the table's LLHrel)
%   LLHrelI_<dist>      : LLHi_VSIP / LLHi_VSI (position's gain once vision x speed is in)
%   p_<dist>_pos, p_<dist>_posI : one-sided paired t-test over the 10 CV folds that
%                         position improves held-out likelihood (VSP>VS, VSIP>VSI)
% Ex (struct):
%   foldLL   [nCells x K x nModels x nDists] held-out log-likelihood minus null, nats
%   kernels  struct of [nCells x 100] position kernels fitted on all data at the
%            chosen lambda: gauss_VSP (additive, activity units), pois_VSP and
%            pois_VSIP (log multiplicative gain); NaN for never-visited bins

if nargin < 2, cellSel = []; end
if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'kfold', 10);
opts = setDefault(opts, 'lamGauss', [1e2 1e3 1e4 1e5]);
opts = setDefault(opts, 'lamPois', [1e2 1e3 1e4]);
opts = setDefault(opts, 'useParfor', true);
dists = {'gauss','pois'};

C = sk_config();
mouse = C.sessions{iSess,1}; date = C.sessions{iSess,2}; sid = C.sessions{iSess,3};
tag = sprintf('%s_%d', mouse, date);
L = load(fullfile(C.outRoot, [tag '_samples.mat']));
if isfield(opts, 'Y'), L.Y = opts.Y; end         % simulated activity (sk_simulate_null), same size as L.Y
fit = L.fit; keep = fit.inGLM; bi = fit.blockIdx;
Xfull = L.X; Yall = double(L.Y(keep,:));
trial = fit.trialID(keep);
nC = size(Yall,2); if isempty(cellSel), cellSel = 1:nC; end
cellSel = cellSel(:)';

% per-cell visual delay from the original fit (GLMs{1}.Perf.bestVisdelay, ms):
% visual columns are shifted so the input d samples earlier drives the response now
delayMs = sk_read_best_delay(L.meta.file);
delaySmp = round(delayMs(cellSel)' / 1000 * L.meta.sampleRate);   % row vector
visAll = ismember(bi, [1 3 4 5 6]);

% folds: the original GLM1 partition (contiguous blocks of whole trials),
% recovered from the stored per-fold activity totals
[foldOfTrial, ut, exact] = sk_recover_folds(iSess);
if ~exact, warning('%s: original folds not recovered exactly', tag); end
[~, ti] = ismember(trial, ut);
fold = foldOfTrial(ti);

models = struct('name', {'VS','VSP','VSI','VSIP'}, 'blocks', {1:6, 1:7, [1:6 8], 1:8});
nM = numel(models); nD = numel(dists); K = opts.kfold;
nSel = numel(cellSel);
nSpk = sum(Yall(:,cellSel),1)';
nPosBins = nnz(bi == 7);

% Build the aggregated design once per delay group, then fit all cells in one
% parfor (good load balance regardless of group sizes).
dlys = unique(delaySmp);
GD = cell(numel(dlys),1);
grpOf = zeros(nSel,1); colOf = zeros(nSel,1);
for gi = 1:numel(dlys)
    dly = dlys(gi);
    grp = find(delaySmp == dly);
    grpOf(grp) = gi; colOf(grp) = 1:numel(grp);
    X = Xfull;
    if dly > 0, X(:,visAll) = [zeros(dly, nnz(visAll), 'uint8'); Xfull(1:end-dly, visAll)]; end
    X = X(keep,:);

    % vision x speed interaction: visual columns times centred speed-bin index
    spdCols = find(bi == 2);
    [~, spdBin] = max(X(:,spdCols), [], 2); spdBin(~any(X(:,spdCols),2)) = NaN;
    sz = (spdBin - mean(spdBin,'omitnan')) / std(spdBin,'omitnan'); sz(isnan(sz)) = 0;
    visCols = find(visAll);
    XI = double(X(:,visCols)) .* sz;             % discrete values -> still aggregatable
    Xall = [double(X), XI]; bAll = [bi; 8*ones(numel(visCols),1)];
    keepCol = any(Xall ~= 0, 1)';
    posBin = find(keepCol(bi == 7));             % corridor bin of each kept position column
    Xall = Xall(:, keepCol); bAll = bAll(keepCol);

    % aggregate identical rows within fold
    [U, ~, g] = unique([fold, Xall], 'rows');
    foldU = U(:,1); XU = sparse(U(:,2:end)); nU = accumarray(g, 1);
    Y = Yall(:, cellSel(grp));
    SY = zeros(size(U,1), numel(grp)); SY2 = SY;
    for c = 1:numel(grp)
        SY(:,c) = accumarray(g, Y(:,c)); SY2(:,c) = accumarray(g, Y(:,c).^2);
    end
    fprintf('[%s] delay %d samples: %d cells, %d rows -> %d aggregated\n', tag, dly, numel(grp), numel(g), size(U,1));
    % store the design once; each worker selects a model's columns itself
    % (storing one matrix per model made the broadcast too large for 16 workers)
    inM = cell(nM,1); posIdx = cell(nM,1);
    for m = 1:nM
        inM{m} = ismember(bAll, models(m).blocks);
        posIdx{m} = 1 + find(bAll(inM{m}) == 7); % position coefficients within that model's design
    end
    GD{gi} = struct('XU',XU,'inM',{inM},'posIdx',{posIdx},'posBin',posBin,'nU',nU,'foldU',foldU,'SY',SY,'SY2',SY2);
end

lamSets = {opts.lamGauss, opts.lamPois};
resAll = nan(nSel, nM, nD); lbAll = resAll; foldAll = nan(nSel, K, nM, nD);
kGV = nan(nSel, nPosBins); kPV = kGV; kPVI = kGV;
t0 = tic;
if opts.useParfor
    parfor c = 1:nSel
        [resAll(c,:,:), lbAll(c,:,:), foldAll(c,:,:,:), kGV(c,:), kPV(c,:), kPVI(c,:)] = ...
            fitCellAll(GD{grpOf(c)}, colOf(c), K, lamSets, dists, nM, nPosBins); %#ok<PFBNS>
    end
else
    for c = 1:nSel
        [resAll(c,:,:), lbAll(c,:,:), foldAll(c,:,:,:), kGV(c,:), kPV(c,:), kPVI(c,:)] = ...
            fitCellAll(GD{grpOf(c)}, colOf(c), K, lamSets, dists, nM, nPosBins);
    end
end
fprintf('[%s] %d cells fitted in %.0f s\n', tag, nSel, toc(t0));
LLHi = resAll ./ nSpk / log(2);

Rf = table(repmat(sid,nSel,1), cellSel(:), nSpk, delaySmp(:), 'VariableNames', {'sessionid','cellNum','nSpikes','visDelaySamples'});
for d = 1:nD
    for m = 1:nM
        Rf.(sprintf('LLHi_%s_%s', dists{d}, models(m).name)) = LLHi(:,m,d);
        Rf.(sprintf('lam_%s_%s', dists{d}, models(m).name)) = lbAll(:,m,d);
    end
    Rf.(sprintf('LLHrel_%s', dists{d}))  = LLHi(:,2,d) ./ LLHi(:,1,d);
    Rf.(sprintf('LLHrelI_%s', dists{d})) = LLHi(:,4,d) ./ LLHi(:,3,d);
    Rf.(sprintf('p_%s_pos', dists{d}))  = foldTest(foldAll(:,:,2,d) - foldAll(:,:,1,d));
    Rf.(sprintf('p_%s_posI', dists{d})) = foldTest(foldAll(:,:,4,d) - foldAll(:,:,3,d));
end
Ex = struct('sessionid', sid, 'cellNum', cellSel(:), 'foldLL', foldAll, ...
    'models', {{models.name}}, 'dists', {dists}, ...
    'kernels', struct('gauss_VSP', kGV, 'pois_VSP', kPV, 'pois_VSIP', kPVI), ...
    'binCentres', C.binCentres);
end

% ------------------------------------------------------------------------
function p = foldTest(dF)
% one-sided paired t-test over folds that the difference is > 0
n = sum(~isnan(dF),2); m = mean(dF,2,'omitnan'); s = std(dF,0,2,'omitnan');
t = m ./ (s ./ sqrt(n));
p = 1 - tcdf(t, n-1);
end

function [res, lb, foldLL, kGV, kPV, kPVI] = fitCellAll(G, j, K, lamSets, dists, nM, nPosBins)
% Gaussian: lambda chosen by CV separately for every model (closed form, cheap).
% Poisson: lambda chosen by CV on the VSP model (model 2) and reused for the
% other models, so all four Poisson models of a cell share one penalty.
% Then position kernels are fitted on all data at the chosen lambda.
res = nan(1, nM, numel(dists)); lb = res; foldLL = nan(1, K, nM, numel(dists));
sy = G.SY(:,j); sy2 = G.SY2(:,j);
G.Xm = cell(nM,1);
for m = 1:nM, G.Xm{m} = [ones(size(G.XU,1),1), G.XU(:, G.inM{m})]; end
for d = 1:numel(dists)
    if strcmp(dists{d}, 'pois')
        [res(1,2,d), lb(1,2,d), foldLL(1,:,2,d)] = cvCell(G.Xm{2}, G.nU, sy, sy2, G.foldU, K, lamSets{d}, dists{d});
        for m = [1 3:nM]
            [res(1,m,d), lb(1,m,d), foldLL(1,:,m,d)] = cvCell(G.Xm{m}, G.nU, sy, sy2, G.foldU, K, lb(1,2,d), dists{d});
        end
    else
        for m = 1:nM
            [res(1,m,d), lb(1,m,d), foldLL(1,:,m,d)] = cvCell(G.Xm{m}, G.nU, sy, sy2, G.foldU, K, lamSets{d}, dists{d});
        end
    end
end
kern = @(m, d) posKernel(G, m, fitOne(G.Xm{m}, G.nU, sy, sy2, lb(1,m,d), dists{d}, []), nPosBins);
kGV  = kern(2, 1);
kPV  = kern(2, 2);
kPVI = kern(4, 2);
end

function k = posKernel(G, m, b, nPosBins)
k = nan(1, nPosBins);
k(G.posBin) = b(G.posIdx{m});
end

function [dLL, lamBest, foldLL] = cvCell(Xm, n, sy, sy2, fold, K, lams, dist)
% CV log-likelihood gain over a constant-rate model, summed over folds,
% at the lambda with the best total CV log-likelihood; also per fold.
lams = sort(lams, 'descend');                    % strong -> weak penalty, warm-started
ll = zeros(K, numel(lams)); ll0 = zeros(K,1);
for k = 1:K
    tr = fold ~= k; te = ~tr;
    if ~any(te), ll(k,:) = NaN; ll0(k) = NaN; continue; end
    Xtr = Xm(tr,:); Xte = Xm(te,:);              % slice the sparse design once per fold
    ntr = n(tr); sytr = sy(tr); sy2tr = sy2(tr);
    b = []; XtX = [];
    if strcmp(dist, 'gauss'), XtX = full(Xtr' * (Xtr .* ntr)); end   % shared by all lambdas
    for il = 1:numel(lams)
        b = fitOne(Xtr, ntr, sytr, sy2tr, lams(il), dist, b, XtX);
        ll(k,il) = testLL(Xte, n(te), sy(te), sy2(te), b, dist, Xtr, ntr, sytr, sy2tr);
    end
    b0 = fitOne(ones(nnz(tr),1), n(tr), sy(tr), sy2(tr), 0, dist, []);
    ll0(k) = testLL(ones(nnz(te),1), n(te), sy(te), sy2(te), b0, dist, ones(nnz(tr),1), n(tr), sy(tr), sy2(tr));
end
[best, ib] = max(sum(ll,1,'omitnan'));
dLL = best - sum(ll0,'omitnan'); lamBest = lams(ib);
foldLL = (ll(:,ib) - ll0)';
end

function b = fitOne(X, n, sy, sy2, lam, dist, b0, XtX) %#ok<INUSL>
% X may be sparse (one-hot design); the normal/Hessian matrices are formed
% dense, since they are nearly full and a dense Cholesky is much faster than
% a sparse solve on them.
p = size(X,2); pen = lam*ones(p,1); pen(1) = 0;  % intercept unpenalised
if nargin < 8, XtX = []; end
switch dist
    case 'gauss'      % weighted ridge on row means
        if isempty(XtX), XtX = full(X' * (X .* n)); end
        A = XtX; A(1:p+1:end) = A(1:p+1:end) + pen';
        b = cholSolve(A, full(X' * sy));
    case 'pois'       % Newton with exposure n
        if isempty(b0), b = zeros(p,1); b(1) = log(max(sum(sy)/sum(n), 1e-6)); else, b = b0; end
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
end

function x = cholSolve(A, y)
[R, flag] = chol(A);
if flag == 0, x = R \ (R' \ y); else, x = A \ y; end
end

function ll = testLL(X, n, sy, sy2, b, dist, Xtr, ntr, sytr, sy2tr)
eta = X*b;
switch dist
    case 'gauss'      % sigma^2 from training residuals
        etr = Xtr*b;
        ss_tr = sum(sy2tr - 2*etr.*sytr + ntr.*etr.^2);
        s2 = max(ss_tr / sum(ntr), 1e-12);
        ss = sum(sy2 - 2*eta.*sy + n.*eta.^2);
        ll = -0.5*sum(n)*log(2*pi*s2) - ss/(2*s2);
    case 'pois'       % log(y!) omitted: identical for every model on the same data
        ll = sy' * eta - n' * exp(eta);
end
end

function o = setDefault(o, f, v), if ~isfield(o,f), o.(f) = v; end, end
