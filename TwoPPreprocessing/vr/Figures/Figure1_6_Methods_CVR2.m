% example figure to show cvr2 (methhods)
sessionROIData = load("Z:\ibn-vision\DATA\SUBJECTS\M26004\Analysis\20260318\M26004_20260318_sessionROIData.mat")
response = load("Z:\ibn-vision\DATA\SUBJECTS\M26004\Analysis\20260318\M26004_20260318_Response_M26004_LandManipCorridor_20260318_00002.mat", "lapPositionActivity", "trialIndicesByCondition")


baseIdx = response.trialIndicesByCondition.Baseline;

lapPositionActivity = response.lapPositionActivity.dFFNeuropilCorrected(:,baseIdx, :);


meanExpVar = sessionROIData.crossValExpVar.dFFNeuropilCorrected.meanExpVar;
cvExpVar = sessionROIData.crossValExpVar.dFFNeuropilCorrected.cvExpVar;
pVals = sessionROIData.crossValExpVar.dFFNeuropilCorrected.pValues;
[~, sortIdx] = sort(meanExpVar, 'descend');
highIdx = sortIdx(1);
lowIdx  = sortIdx(round(length(sortIdx)*0.97)); % a middling/lower example

figABCD = plotTuningCurveExamples(lapPositionActivity, cvExpVar, highIdx, lowIdx);
fprintf('p (high) = %.3f, p (low) = %.3f\n', pVals(highIdx), pVals(lowIdx));

set(figABCD, 'Visible', 'on');
saveFigureFormats(figABCD, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\Methods\Fig2_7_CVR2_examples', 'tuned_untuned_eg'));