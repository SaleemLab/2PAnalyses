function [ratioVarToTuningRange, ratioVarToTuningVar, meanTuning]= computeVarianceAcrossPositionBins(sessionFileInfo, response, signalToUse)

% calculates metrics of spatial selectivity and lap-to-lap reliability 
% for each ROI by comparing the variance of the average
% spatial tuning curve against the average activity 
% variance within individual position bins .

%% Handle optional arguments 
if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end

%% Get data 
if response.signalsZScored
    disp('Using zScored-lapPositionActivity for computations..')
else
    disp('Using lapPositionActivity (without zscoring) for computations..')
end 

lapPositionActivity = response.lapPositionActivity.(signalToUse);
% Mean tuning across laps 
meanTuning   = squeeze(nanmean(lapPositionActivity,2));
% Measure variance across lap-averaged position bins
tuningVar = var(meanTuning,0,2);

% Tuning range 
% Max - Min Mean Activity. This is the "signal magnitude." 
% It measures the difference between the peak and trough of the 
% average spatial tuning curve. It represents the maximum possible 
% change in the signal due to position.
tuningRange = (max(meanTuning')-min(meanTuning'))';

% Measure mean variance across all laps
% lap-to-lap variability?
% It measures the average variance of the activity within a single position
% bin across all laps. High posAllVar means the cell's activity at any 
% given location is inconsistent from one lap to the next.
posAllVar = nanmean(nanvar(lapPositionActivity,0,2),3);


%% How noisy is the roi's activity compared to the total span/magnitude of its signal
% LOW = the cell's activity across laps is stable and precise relative 
% to the magnitude of its spatial change. The position-dependent change
%  (signal) is much larger than the lap-to-lap fluctuation (noise).

% HIGH = he cell's activity is highly variable and that the lap-to-lap 
% inconsistency (noise) is a substantial fraction of the cell's total 
% firing range (signal).
ratioVarToTuningRange = posAllVar./tuningRange;

%% tuning quality or spatial selectivity relative to noise. 
% It compares the inherent noise to the actual structure of the spatial tuning curve.
% LOW = strongly-tuned neuron. 
% It means the variance of the average spatial tuning curve 
% is much larger than the lap-to-lap variability (low noise). 
% The cell is selective for a position and reliably fires at that position.

% HIGH = lap-to-lap variability is high.
% suggests that the lap-to-lap variability is high relative to the 
% strength of the spatial tuning. The cell might have a weak or
% poor place field structure that is easily obscured by the noise.
ratioVarToTuningVar   = posAllVar./tuningVar;

%% 
tuningCurveVariance.ratioVarToTuningVar = ratioVarToTuningVar;
tuningCurveVariance.ratioVarToTuningRange = ratioVarToTuningRange; 


%% Check if sessionROIData exists and append 
if exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
       "tuningCurveVariance", ...
         '-append')
else
    warning('sessionROIData file not found at: %s. Cannot append variance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
end

end 


% threshold = 10;
% selective_ROI_indices = find(ratioVarToTuningVar <= threshold);
% selective_tuning_curves = meanTuning(selective_ROI_indices, :);
% figure;
% num_selective_rois = length(selective_ROI_indices);
% for i = 1:num_selective_rois
%     current_roi_number = selective_ROI_indices(i);
%     clf; 
%     current_tuning_curve = meanTuning(current_roi_number, :);
%     plot(current_tuning_curve, 'LineWidth', 2, 'Color', 'b');
%     xlabel('Position Bin');
%     ylabel('Mean Activity');
%     title_text = ['Selective ROI #', num2str(current_roi_number), ...
%                   ' (', num2str(i), ' of ', num2str(num_selective_rois), ')'];
%     title(title_text);
%     pause(0.5); 
% 
% end
