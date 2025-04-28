% Load mergedBonsaiSuite2PData from M25041 20250328 SparseNoise 
isRoiIdx = find(twoPData.iscell(:, 1) == 1);
rawFRois = twoPData.F(isRoiIdx, :);
rawFneuRois = twoPData.Fneu(isRoiIdx,:);
frameRate = 8.31;

%Calculates the corrected neuropil traces and the specific values that
% were used to determine the correction factor (intercept and slope of
% linear fits, F traces bin values, N traces bin values). Refer to function
% for further details.
[Fc, regPars, F_binValues, N_binValues] = correct_neuropil( ...
    rawFRois', ...
    rawFneuRois', ...
    frameRate);

% Calculates the baseline fluorescence F0 used to calculate delta F over F.
F0 = get_F0( ...
    Fc, ...
    frameRate);

% Calculates delta F oer F given the corrected neuropil traces and the
% baseline fluorescence.
dF = get_delta_F_over_F(Fc, F0)';

% Temporal smoothning:  
w = gausswin(10); w = w / sum(w);
smoothed = filtfilt(w, 1, dF(roiIdx,:));
roiIdx = 70;

figure; 
plot(twoPData.ArduinoTime, dF(roiIdx,:))
hold on 
plot(twoPData.ArduinoTime,smoothed, LineWidth=3)

% Interpolation 
