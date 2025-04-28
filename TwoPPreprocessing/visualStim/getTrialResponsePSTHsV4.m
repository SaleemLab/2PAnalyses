function [response, sessionFileInfo] = getTrialResponsePSTHsV4( ...
    sessionFileInfo, stimName, method, interpRate, frameRate, ...
    applyDeltaFoverF, applyNeuropilCorrection)

% Defaults
if nargin<3, method = 1; end
if nargin<4, interpRate = 60; end
if nargin<5, frameRate = 7.28; end
if nargin<6, applyDeltaFoverF = true; end
if nargin<7, applyNeuropilCorrection = true; end

% Locate stimulus
iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}),1);
assert(~isempty(iStim), 'Stimulus "%s" not found.', stimName);

% Load data
s = load(sessionFileInfo.stimFiles(iStim).Response, 'response');
response = s.response;
b = load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
bonsaiData = b.bonsaiData;
t = load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData');
twoPData = t.twoPData;

% Build combined (F_c or dF/F) signal matrix: [nNeurons x nFrames]
signalMatrix = [];
planeIdx      = [];
for p = 1:numel(twoPData)
    F    = twoPData(p).F;       % [nROI x nFrames]
    Fneu = twoPData(p).Fneu;    % [nROI x nFrames]
    if applyNeuropilCorrection
        [Fc,~,~,~] = correct_neuropil(F',Fneu',frameRate);
        Fc = Fc';  % back to [nROI x nFrames]
    else
        Fc = F;
    end
    if applyDeltaFoverF
        F0  = get_F0(Fc', frameRate)';        % [nROI x nFrames] -> transpose
        dFF = get_delta_F_over_F(Fc, F0);      % [nROI x nFrames]
        M   = dFF;
    else
        M = Fc;
    end
    signalMatrix = [signalMatrix; M];
    planeIdx      = [planeIdx; repmat(p, size(M,1),1)];
end
nNeurons = size(signalMatrix,1);

% Prepare output struct
nGroups = numel(bonsaiData.trialGroups);
pd = struct( ...
    'stimValue',         [], ...
    'alignedResponses',  [], ...
    'meanResponse',      [], ...
    'stdResponse',       [], ...
    'semResponse',       [], ...
    'timeVector',        [], ...
    'responseType',      [] );
response.psthData = repmat(pd, nGroups,1);

% If method 4: build a global interpolated timeseries
if method == 4
    % collect ALL frame times across every trial & plane:
    allTimes = [];
    for p = 1:numel(response)
        for tr = bonsaiData.trialGroups(:)'  % use trials from bonsaiData
            idxs = tr.trials;
            rt   = response(p).responseFrameRelTimes(idxs);
            allTimes = [allTimes; vertcat(rt{:})];
        end
    end
    tmin = min(allTimes);  tmax = max(allTimes);
    globalTime = linspace(tmin, tmax, round((tmax-tmin)*interpRate));
    interpSignal = nan(nNeurons, numel(globalTime));
    for n = 1:nNeurons
        p = planeIdx(n);
        % pick the first trial of that plane as a proxy
        idx0 = bonsaiData.trialGroups(1).trials(1);
        RT   = response(p).responseFrameRelTimes{idx0};
        FM   = response(p).responseFrameIdx{idx0};
        rawF = signalMatrix(n, FM);
        interpSignal(n,:) = interp1(RT, rawF, globalTime, 'linear', NaN);
    end
end

% Loop through each stimulus group
for g = 1:nGroups
    grp    = bonsaiData.trialGroups(g);
    trIdxs = grp.trials;
    nTrials= numel(trIdxs);
    
    % Build common time vector for methods 1–3
    if method ~= 4
        timesAll = [];
        for p = 1:numel(response)
            rt = response(p).responseFrameRelTimes(trIdxs);
            timesAll = [timesAll; vertcat(rt{:})];
        end
        tmin = min(timesAll); tmax = max(timesAll);
        timeVector = linspace(tmin, tmax, round((tmax-tmin)*interpRate));
    else
        timeVector = globalTime;
    end
    
    aligned = nan(nNeurons, numel(timeVector), nTrials);
    
    for p = 1:numel(twoPData)
        neuronIdxs = find(planeIdx == p);
        nROI       = numel(neuronIdxs);
        
        for ti = 1:nTrials
            trialID = trIdxs(ti);
            frameMask   = response(p).responseFrameIdx{trialID};
            relTimes    = response(p).responseFrameRelTimes{trialID};
            
            switch method
              case 1  % per-trial interpolation
                D = signalMatrix(neuronIdxs, frameMask);
                for ni = 1:nROI
                    aligned(neuronIdxs(ni),:,ti) = ...
                      interp1(relTimes, D(ni,:), timeVector, 'linear', NaN);
                end
                
              case 2  % smoothing then interp
                % five‐point Gaussian smoothing
                w = gausswin(5); w = w/sum(w);
                Dsm = filtfilt(w,1,signalMatrix(neuronIdxs,:)')';
                D   = Dsm(:, frameMask);
                for ni = 1:nROI
                    aligned(neuronIdxs(ni),:,ti) = ...
                      interp1(relTimes, D(ni,:), timeVector, 'linear', NaN);
                end
                
              case 3  % direct mean (no interp)
                Dfull = signalMatrix(neuronIdxs, frameMask);
                nT = size(Dfull,2);
                aligned(neuronIdxs,1:nT,ti) = Dfull;
                
              case 4  % global interpolation then split
                aligned(neuronIdxs,:,ti) = interpSignal(neuronIdxs,:);
                
              otherwise
                error('Invalid method %d', method);
            end
        end
    end
    
    % Compute mean/std/sem across trials
    mResp = nanmean(aligned,3);
    sResp = nanstd(aligned,0,3);
    semR  = sResp ./ sqrt(sum(~isnan(aligned),3));
    
    % Store in response
    response.psthData(g).stimValue        = grp.value;
    response.psthData(g).alignedResponses = aligned;
    response.psthData(g).meanResponse     = mResp;
    response.psthData(g).stdResponse      = sResp;
    response.psthData(g).semResponse      = semR;
    response.psthData(g).timeVector       = timeVector;
    response.psthData(g).responseType     = method;
end

% save back out
save(sessionFileInfo.stimFiles(iStim).Response, 'response');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end
