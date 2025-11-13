function [response, sessionFileInfo] = getLapPositionActivityV3(sessionFileInfo, VRStimName, applyTemporalSmoothing, onlyIncludeROIs)
% Calculates and saves binned lap activity for four signal types:
% F_raw, dFF_raw, F_neuropil_corrected, and dFF_neuropil_corrected.
% MODIFIED: Temporal smoothing is now applied to F and Fneu at the start.
%
% Aman and Sonali October 2025

%% Handle optional inputs
if nargin < 3, applyTemporalSmoothing = true; end
if nargin < 4, onlyIncludeROIs = false; end

%% Find VR stimulus
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end

%% Load data
disp('Loading processedTwoPData and response...');
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');
load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');

%% Pre-processing: Optional Temporal Smoothing

F = processedTwoPData.F;    % Original F (ROI x time)
Fneu = processedTwoPData.Fneu; % Original Fneu (ROI x time)

if applyTemporalSmoothing
    disp('Applying temporal smoothing to F and Fneu time-series...');
    w = gausswin(9);
    w = w / sum(w);
    
    numROIs = size(F, 1);
    % Using fSmoothed and fneuSmoothed >>>
    fSmoothed = zeros(size(F));
    fneuSmoothed = zeros(size(Fneu));
    
    % Loop through each ROI to apply the filter along the time dimension
    for i = 1:numROIs
        fSmoothed(i, :) = filtfilt(w, 1, F(i, :));
        fneuSmoothed(i, :) = filtfilt(w, 1, Fneu(i, :));
    end
else
    disp('Skipping temporal smoothing.');
    % Might be good to change variable name here to something more neutral 
    fSmoothed = F;
    fneuSmoothed = Fneu;
end

%% Prepare all four signal matrices
disp('Preparing all four signal types from data: fRaw, fNeuropilCorrected, dFF and dFFNeuropilCorrected');
fs = processedTwoPData.ops{1}.fs;  % sampling rate

% Raw F
signals.fRaw = fSmoothed;

% Neuropil-Corrected F (Fc)
disp('Computing neuropil correction...');
[Fc, ~, ~, ~] = correct_neuropil(fSmoothed', fneuSmoothed', fs);
signals.fNeuropilCorrected = Fc'; % Transpose back to ROI x time

% Delta F/F on Raw F
disp('Computing delta f/f...');
f0Raw = get_F0(fSmoothed', fs)';
signals.dFF = get_delta_F_over_F(fSmoothed, f0Raw);

% Delta F/F on Neuropil-Corrected F
disp('Computing delta f/f on neuropil corrected f...');
F0_c = get_F0(Fc, fs)'; 
signals.dFFNeuropilCorrected = get_delta_F_over_F(signals.fNeuropilCorrected, F0_c);

signalNames = fieldnames(signals);
disp('Signal matrices created: fRaw, fNeuropilCorrected, dFF & dFFNeuropilCorrected');

%% Get cell ROIs if needed
if onlyIncludeROIs
    isCell = logical(processedTwoPData.iscell(:, 1));
    cellROIs = find(isCell);
else
    cellROIs = 1:size(F, 1); 
end
numCells = length(cellROIs);

%% Binning parameters
numBins = 140; % 1cm bins for a 140cm track
nLaps = length(response.completedStartTimes);
LabPositionActivity = struct(); 

%% Extract binned mean signals for all four types
for thisSignal = 1:length(signalNames)
    currentSignalName = signalNames{thisSignal};
    currentSignalMatrix = signals.(currentSignalName);
    disp(['Binning lab-position-activity for: ' currentSignalName '...']);
    
    tempActivity = nan(numCells, nLaps, numBins);
    
    for thisCell = 1:numCells
        roiIdx = cellROIs(thisCell);
        for thisLap = 1:nLaps
            for thisBin = 1:numBins
                frameIdx = response.lapPosition2PFrameIdx{roiIdx, thisLap, thisBin};
                if ~isempty(frameIdx)
                    tempActivity(thisCell, thisLap, thisBin) = mean(currentSignalMatrix(roiIdx, frameIdx));
                end
            end
        end
    end
    LabPositionActivity.(currentSignalName) = tempActivity;
end

%% Save all four activity matrices to the response struct
if isfield(response, 'LabPositionActivities'), response = rmfield(response, 'LabPositionActivities'); end

response.lapPositionActivities = LabPositionActivity;
response.smoothingApplied = applyTemporalSmoothing; 
if onlyIncludeROIs
    response.cellROIs = cellROIs;
end

disp('Saving response with all lapPositionActivities...');
save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');
end