function [foldOfTrial, trials, exact] = sk_recover_folds(iSess)
% Recover the original GLM1 trial -> CV-fold assignment for one session.
% The stored per-cell, per-fold activity totals (GLMs{1}.Tuning(1).nSpikes_fold)
% pin the partition down: with hundreds of cells, only one assignment of whole
% trials to folds reproduces them. Solved as a small integer program on the 60
% most active cells, then checked against all cells. Cached in outRoot.
C = sk_config();
mouse = C.sessions{iSess,1}; date = C.sessions{iSess,2};
tag = sprintf('%s_%d', mouse, date);
cache = fullfile(C.outRoot, [tag '_recoveredFolds.mat']);
if exist(cache,'file')
    S = load(cache);
    if isfield(S,'exact'), foldOfTrial = S.foldOfTrial; trials = S.ut; exact = S.exact; return; end
end
L = load(fullfile(C.outRoot, [tag '_samples.mat']), 'Y', 'fit', 'meta');
G = load(L.meta.file, 'GLMs'); nsf = G.GLMs{1}.Tuning(1).nSpikes_fold;   % [cells x K]
k = L.fit.inGLM; Y = double(L.Y(k,:)); tr = L.fit.trialID(k);
ut = unique(tr); nT = numel(ut); K = size(nsf,2); nC = size(Y,2);
V = zeros(nC, nT); for t = 1:nT, V(:,t) = sum(Y(tr==ut(t),:),1)'; end
[~,ord] = sort(sum(V,2),'descend'); cs = ord(1:min(60,nC));
nx = nT*K; ns = numel(cs)*K;
Aeq = [kron(ones(1,K), speye(nT)), sparse(nT,2*ns); ...
       kron(speye(K), sparse(V(cs,:))), speye(ns), -speye(ns)];
beq = [ones(nT,1); reshape(nsf(cs,:),[],1)];
fobj = [zeros(nx,1); ones(2*ns,1)];
o = optimoptions('intlinprog','Display','off','MaxTime',600);
z = intlinprog(fobj, 1:nx, [], [], Aeq, beq, zeros(nx+2*ns,1), [ones(nx,1); inf(2*ns,1)], o);
X = reshape(round(z(1:nx)), nT, K);
[~, foldOfTrial] = max(X, [], 2);
M = V * X;
relErr = max(abs(M(:)-nsf(:))) / max(nsf(:));
exact = relErr < 1e-6;
fprintf('[%s] recovered folds: max rel. error over all cells %.2g (%s); trials per fold %s\n', ...
    tag, relErr, ternary(exact,'exact','NOT exact'), mat2str(sum(X,1)));
save(cache, 'foldOfTrial', 'ut', 'exact', 'relErr');
trials = ut;
end

function r = ternary(c,a,b), if c, r = a; else, r = b; end, end
