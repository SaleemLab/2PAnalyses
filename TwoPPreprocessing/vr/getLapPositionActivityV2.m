function [response, sessionFileInfo] = getLapPositionActivityV2(sessionFileInfo, signalField, applySmoothing, VRStimName, onlyIncludeROIs, applyDeltaFoverF, applyNeuropilCorrection)
%% Split the delta f/f into a different function and save this as a new session; Using median as the metric. 
% Handle optional inputs
if nargin < 2 || isempty(signalField), signalField = 'F'; end
if nargin < 3, applySmoothing = false; end
if nargin < 4, error('VRStimName must be provided'); end
if nargin < 5, onlyIncludeROIs = true; end
if nargin < 6, applyDeltaFoverF = true; end
if nargin < 7, applyNeuropilCorrection = true; end

%% Find VR stimulus
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName, sessionFileInfo.stimFiles(iStim).name);
end
iStim = find(bonsaiData.isVRstim == 1);
if isempty(iStim), error('No VRCorr stimulus found in sessionFileInfo.'); end

%% Load data
disp('Loading processedTwoPData and response...');
load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'processedTwoPData')
load(sessionFileInfo.stimFiles(iStim).Response, 'response')

%% Set signal source
F = processedTwoPData.F;        % ROI x time
Fneu = processedTwoPData.Fneu;  % ROI x time
fs = processedTwoPData.ops{1}.fs;  % sampling rate

% Get cell ROIs if needed
if onlyIncludeROIs
    isCell = logical(processedTwoPData.iscell(:, 1));
    cellROIs = find(isCell);
    numCells = length(cellROIs);
else
    numCells = size(F, 1);
    cellROIs = 1:numCells;
end

% Optional: neuropil correction
if applyNeuropilCorrection
    disp('Applying Nueropil corrections ...');
    [Fc, ~, ~, ~] = correct_neuropil(F', Fneu', fs);  % outputs time x ROI
    Fc = Fc';  % back to ROI x time
else
    Fc = F;
end

% Optional: delta F over F
if applyDeltaFoverF
    disp('Calculating delta F over F ...');
    F0 = get_F0(Fc', fs);    % outputs time x ROI
    F0 = F0';
    signalMatrix = get_delta_F_over_F(Fc, F0);  % ROI x time
    signalField = 'deltaF_over_F';
else
    signalMatrix = Fc;
end

%% Bin parameters
binCentres = 0.5:1:139.5;
numBins = length(binCentres);
nLaps = length(response.completedStartTimes);
lapPositionActivity = nan(numCells, nLaps, numBins);

%% Extract binned mean signals
% Compute mean signal per bin per lap
for thisCell = 1:numCells
    if onlyIncludeROIs
        roiIdx = cellROIs(thisCell);
    else
        roiIdx = thisCell;
    end
    for lapIdx = 1:nLaps
        for binIdx = 1:numBins
            frameIdx = response.lapPosition2PFrameIdx{roiIdx, lapIdx, binIdx};
            if ~isempty(frameIdx)
                lapPositionActivity(thisCell, lapIdx, binIdx) = ...
                    median(signalMatrix(roiIdx, frameIdx));
            end
        end
    end
end

%% Optional smoothing
if applySmoothing
    w = gausswin(9); w = w / sum(w);
    for thisCell = 1:numCells
        for lapIdx = 1:nLaps
            signal = squeeze(lapPositionActivity(thisCell, lapIdx, :));
            if all(isnan(signal)), continue; end
            nanMask = isnan(signal);
            signal(nanMask) = 0;
            smoothed = filtfilt(w, 1, signal);
            smoothed(nanMask) = NaN;
            lapPositionActivity(thisCell, lapIdx, :) = smoothed;
        end
    end
end

%% Save and return
response.lapPositionActivity = lapPositionActivity;
response.signalUsed = signalField;
response.smoothingApplied = applySmoothing;
response.deltaFoverFApplied = applyDeltaFoverF;
response.neuropilCorrected = applyNeuropilCorrection;
if onlyIncludeROIs
    response.cellROIs = cellROIs;
end

disp('Saving updated response with lapPositionActivity...');
save(sessionFileInfo.stimFiles(iStim).Response, 'response', '-v7.3');

end
