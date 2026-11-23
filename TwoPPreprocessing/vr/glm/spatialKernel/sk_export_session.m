function sk_export_session(iSess)
% Export one SpkLin session into two compact files, reading only the fields
% needed (h5read on the v7.3 file, so the 10 GB structs are never loaded):
%
%   <mouse>_<date>_binned.mat   trial x 2cm-bin arrays for every cell
%       spk, predFull, predVS, predPos : [nCells x nTrials x nBins] single,
%           summed over samples (counts per bin per trial)
%       occ   : [nTrials x nBins] number of 60 Hz samples
%       speed : [nTrials x nBins] mean smthRunSpeed
%       trialCond : [nTrials x 1] condition code (1 base,2 swap23,3 swap34,4 omit2,5 omit3,6 omit4)
%       condDesign : struct, per condition x bin mean of each visual design block
%   <mouse>_<date>_samples.mat  sample-level data for refitting the GLM
%       X (uint8 one-hot design, all 7 blocks), blockIdx, Y (single spikeTrain),
%       inGLM (exact original GLM1 sample subset), inFit (inGLM with a valid
%       position), trialID, condition, speed, pos, blanks
%
% NOTE: these are 2-photon recordings. Spk.spikeTrain is deconvolved activity
% (non-integer, arbitrary scale), not spike counts; Spk.spikeTimes only lists
% the sample indices where it is non-zero.
%
% Cell order is the Spk column order, i.e. T.cellNum within the session.

C = sk_config();
mouse = C.sessions{iSess,1}; date = C.sessions{iSess,2}; sid = C.sessions{iSess,3};
d = dir(fullfile(C.rawRoot, mouse, num2str(date), 'SpkLin_EXP_*.mat'));
if isempty(d), fprintf('session %s %d not found, skipping\n', mouse, date); return; end
f = fullfile(d(1).folder, d(1).name);
if ~exist(C.outRoot,'dir'), mkdir(C.outRoot); end
tag = sprintf('%s_%d', mouse, date);
fprintf('[%s] reading %s\n', tag, f); t0 = tic;

rd = @(p) h5read(f, p);
trialID   = double(rd('/Nav/trialID'));
condition = double(rd('/Nav/condition'));
pos       = double(rd('/Nav/X'));
speed     = double(rd('/Nav/smthRunSpeed'));
tooShort  = double(rd('/Nav/TrialTooShort'));
blanks    = double(rd('/Nav/blanks'));

% Original GLM1 subset: condition 1-6, smthRunSpeed > 1, blanks 0/1, TrialTooShort == 0.
% inGLM is exactly that (it includes blank-screen samples, where X is NaN);
% inFit additionally requires a valid corridor position, for the binned maps.
inGLM = ismember(condition,1:6) & speed > 1 & ismember(blanks,[0 1]) & tooShort == 0;
inFit = inGLM & trialID > 0 & ~isnan(pos);

Y        = single(rd('/Spk/spikeTrain'));
predFullS = single(rd('/Spk/GLMprediction1'));
predVSS   = single(rd('/Spk/GLMprediction_noSpace1'));
predPosS  = single(rd('/Spk/GLMprediction_Space1'));
nS = numel(trialID); nC = size(Y,2);
assert(size(Y,1)==nS && size(predVSS,1)==nS, 'sample count mismatch');

% ---------- binned trial x position arrays ----------
trials = unique(trialID(inFit));
nT = numel(trials); nB = C.nBins;
bin = discretize(pos, C.binEdges); bin(pos>=200) = nB;
ok = inFit & ~isnan(bin);
[~, tIdx] = ismember(trialID, trials);
lin = sub2ind([nT nB], tIdx(ok), bin(ok));     % trial-bin index per kept sample

occ   = reshape(accumarray(lin, 1, [nT*nB 1]), nT, nB);
speedTB = reshape(accumarray(lin, speed(ok), [nT*nB 1]), nT, nB) ./ occ;
binsum = @(M) permute(reshape(sparse(lin, 1:numel(lin), 1, nT*nB, numel(lin)) * double(M(ok,:)), nT, nB, nC), [3 1 2]);
spk      = single(binsum(Y));
predFull = single(binsum(predFullS));
predVS   = single(binsum(predVSS));
predPos  = single(binsum(predPosS));
clear predFullS predVSS predPosS

trialCond = zeros(nT,1);
for t = 1:nT, trialCond(t) = mode(condition(trialID==trials(t) & inFit)); end
trialSpeed = zeros(nT,1);
for t = 1:nT, trialSpeed(t) = mean(speed(trialID==trials(t) & inFit)); end

% ---------- design (one-hot blocks) ----------
Xb = cell(1,numel(C.blocks)); blockIdx = zeros(0,1);
for k = 1:numel(C.blocks)
    B = rd(['/Nav/' C.blocks{k}]);
    B = reshape(B, nS, []);                 % flatten e.g. [nS 19 9] -> [nS 171]
    Xb{k} = uint8(B);
    blockIdx = [blockIdx; k*ones(size(B,2),1)]; %#ok<AGROW>
end
X = [Xb{:}]; clear Xb

% per-condition x bin mean of each visual design column (to locate landmarks)
condDesign = struct();
visCols = find(blockIdx ~= 7 & blockIdx ~= 2);
for c = 1:6
    m = ok & condition == c;
    if ~any(m), condDesign.(C.condNames{c}) = []; continue; end
    A = accumarray(bin(m), 1, [nB 1]);
    S = zeros(nB, numel(visCols));
    for j = 1:numel(visCols), S(:,j) = accumarray(bin(m), double(X(m,visCols(j))), [nB 1]); end
    condDesign.(C.condNames{c}) = S ./ A;
end

cellNum = (1:nC)';
meta = struct('mouse',mouse,'date',date,'sessionid',sid,'file',f,'sampleRate',60, ...
    'binEdges',C.binEdges,'blocks',{C.blocks},'visCols',visCols);

save(fullfile(C.outRoot,[tag '_binned.mat']), 'spk','predFull','predVS','predPos','occ', ...
    'speedTB','trialCond','trialSpeed','trials','condDesign','cellNum','meta','-v7.3');

fit = struct('trialID',trialID,'condition',condition,'speed',single(speed),'pos',single(pos), ...
    'inFit',inFit,'inGLM',inGLM,'blanks',blanks,'blockIdx',blockIdx);
save(fullfile(C.outRoot,[tag '_samples.mat']), 'X','Y','fit','cellNum','meta','-v7.3');
fprintf('[%s] done: %d cells, %d trials, %.0f s\n', tag, nC, nT, toc(t0));
end
