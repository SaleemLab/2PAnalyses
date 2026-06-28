response = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_BaselineCorridor_20260322_00002.mat")
response = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_BaselineCorridor_20260226_CombinedRuns.mat")

responseDark = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_Darkness_20260226_00001.mat")


roiIdx = 335;
speedData = response.lapPositionRunningSpeed; 
rawROIData = response.lapPositionActivity.dFFNeuropilCorrected;

landmarkBin = 40;
window = landmarkBin-10 : landmarkBin+10;  % 

speedAtLandmark = mean(speedData(:, window), 2); 

roiActivity = squeeze(rawROIData(roiIdx, :, window));  
peakResponse = max(roiActivity, [], 2);        

sum(isnan(speedAtLandmark))
sum(isnan(peakResponse))

validIdx = ~isnan(speedAtLandmark) & ~isnan(peakResponse);
[r, p]=corr(speedAtLandmark(validIdx), peakResponse(validIdx))

scatter(speedAtLandmark, peakResponse)
%%

numROIs = size(rawROIData, 1);
rValues = nan(numROIs, 1);
pValues = nan(numROIs, 1);

for roi = 1:numROIs
    peakResponse = max(squeeze(rawROIData(roi, :, window)), [], 2);
    validIdx = ~isnan(speedAtLandmark) & ~isnan(peakResponse);
    [r, p] = corr(speedAtLandmark(validIdx), peakResponse(validIdx));
    rValues(roi) = r;
    pValues(roi) = p;
end
figure; 
histogram(rValues)

[h, p] = ttest(rValues)
mean(rValues)

%%

landmarks = [40, 80, 120, 160];
numROIs = size(rawROIData, 1);

rValues = nan(numROIs, length(landmarks));
pValues = nan(numROIs, length(landmarks));

for lm = 1:length(landmarks)
    landmarkBin = landmarks(lm);
    window = landmarkBin-15 : landmarkBin+15;
    speedAtLandmark = mean(speedData(:, window), 2);
    
    for roi = 1:numROIs
        peakResponse = max(squeeze(rawROIData(roi, :, window)), [], 2);
        validIdx = ~isnan(speedAtLandmark) & ~isnan(peakResponse);
        [r, p] = corr(speedAtLandmark(validIdx), peakResponse(validIdx));
        rValues(roi, lm) = r;
        pValues(roi, lm) = p;
    end
end

for lm = 1:length(landmarks)
    figure;
    histogram(rValues(:, lm))
    title(['Landmark at ' num2str(landmarks(lm)) ' cm'])
    [~, p] = ttest(rValues(:, lm));
    xlabel(['mean r = ' num2str(mean(rValues(:, lm), 'omitnan'), '%.3f') ', p = ' num2str(p, '%.3f')])
end

%% 
landmarks = [40, 80, 120, 160];
numROIs = size(rawROIData, 1);

meanLowAll  = nan(numROIs, length(landmarks));
meanHighAll = nan(numROIs, length(landmarks));

for lm = 1:length(landmarks)
    window = landmarks(lm)-15 : landmarks(lm)+15;
    speedAtLandmark = mean(speedData(:, window), 2);
    medianSpeed = median(speedAtLandmark, 'omitnan');

    for roi = 1:numROIs
        peakResponse = max(squeeze(rawROIData(roi, :, window)), [], 2);
        validIdx = ~isnan(speedAtLandmark) & ~isnan(peakResponse);
        speedValid   = speedAtLandmark(validIdx);
        responseValid = peakResponse(validIdx);

        lowIdx  = speedValid <= medianSpeed;
        highIdx = speedValid > medianSpeed;

        meanLowAll(roi, lm)  = mean(responseValid(lowIdx));
        meanHighAll(roi, lm) = mean(responseValid(highIdx));
    end
end

figure;
scatter(meanLowAll(:), meanHighAll(:), 10, 'filled')
hold on
plot([0 7], [0 7], 'k--')
xlabel('Low speed response')
ylabel('High speed response')

[h, p] = ttest(meanLowAll(:), meanHighAll(:));
fprintf('Low: %.3f, High: %.3f, p = %.4f\n', mean(meanLowAll(:), 'omitnan'), mean(meanHighAll(:), 'omitnan'), p)

%
edges = quantile(speedValid, [0 1/3 2/3 1]);
lowIdx    = speedValid <= edges(2);
medIdx    = speedValid > edges(2) & speedValid <= edges(3);
highIdx   = speedValid > edges(3);

meanLow    = mean(responseValid(lowIdx));
meanMed    = mean(responseValid(medIdx));
meanHigh   = mean(responseValid(highIdx));

%%

landmarks = [40, 80, 120, 160];
numROIs = size(rawROIData, 1);

meanLowAll  = nan(numROIs, length(landmarks));
meanMedAll  = nan(numROIs, length(landmarks));
meanHighAll = nan(numROIs, length(landmarks));

for lm = 1:length(landmarks)
    window = landmarks(lm)-5 : landmarks(lm)+5;
    speedAtLandmark = mean(speedData(:, window), 2);

    for roi = 1:numROIs
        peakResponse = mean(squeeze(rawROIData(roi, :, window)), 2, 'omitnan');
        validIdx = ~isnan(speedAtLandmark) & ~isnan(peakResponse);
        speedValid    = speedAtLandmark(validIdx);
        responseValid = peakResponse(validIdx);

        edges = quantile(speedValid, [0 1/3 2/3 1]);
        lowIdx  = speedValid <= edges(2);
        medIdx  = speedValid > edges(2) & speedValid <= edges(3);
        highIdx = speedValid > edges(3);

        meanLowAll(roi, lm)  = mean(responseValid(lowIdx),  'omitnan');
        meanMedAll(roi, lm)  = mean(responseValid(medIdx),  'omitnan');
        meanHighAll(roi, lm) = mean(responseValid(highIdx), 'omitnan');
    end
end

% population means and SEM
popLow  = mean(meanLowAll(:),  'omitnan');
popMed  = mean(meanMedAll(:),  'omitnan');
popHigh = mean(meanHighAll(:), 'omitnan');

semLow  = std(meanLowAll(:),  'omitnan') / sqrt(sum(~isnan(meanLowAll(:))));
semMed  = std(meanMedAll(:),  'omitnan') / sqrt(sum(~isnan(meanMedAll(:))));
semHigh = std(meanHighAll(:), 'omitnan') / sqrt(sum(~isnan(meanHighAll(:))));

% bar plot
figure;
bar([popLow popMed popHigh], 'FaceColor', [0.4 0.6 0.8])
hold on
errorbar([1 2 3], [popLow popMed popHigh], [semLow semMed semHigh], 'k.', 'LineWidth', 1.5)
xticks([1 2 3])
xticklabels({'Low', 'Medium', 'High'})
xlabel('Running speed')
ylabel('Peak \DeltaF/F')
title('Landmark response by speed group')
% ylim([1.5 1.9])

% stats
[~, p_lm] = ttest(meanLowAll(:), meanMedAll(:));
[~, p_mh] = ttest(meanMedAll(:), meanHighAll(:));
[~, p_lh] = ttest(meanLowAll(:), meanHighAll(:));
fprintf('Low vs Med:  p = %.4f\n', p_lm)
fprintf('Med vs High: p = %.4f\n', p_mh)
fprintf('Low vs High: p = %.4f\n', p_lh)

%%
% get tuning types
tuningTypes = responseDark.response.tuningCurve.dFFNeuropilCorrected.classification.tuningType;

% find indices for each type
highpassIdx = strcmp(tuningTypes, 'highpass');
lowpassIdx  = strcmp(tuningTypes, 'lowpass');
bandpassIdx = strcmp(tuningTypes, 'bandpass');
untunedIdx  = strcmp(tuningTypes, 'untuned');

% compare speed modulation across groups
% meanLowAll and meanHighAll are [numROIs x numLandmarks] from before
diffAll = mean(meanLowAll - meanHighAll, 2, 'omitnan'); % one value per ROI

figure;
boxplot(diffAll, tuningTypes)
xlabel('Darkness tuning type')
ylabel('Low - High speed response in VR')
yline(0, 'k--')

% within your landmark window loop
framesLow  = sum(lowIdx);
framesMed  = sum(medIdx);
framesHigh = sum(highIdx);