twop0001 = load("Z:\ibn-vision\DATA\SUBJECTS\M25126\Analysis\20260123\M25126_20260123_processed2PData_M25126_VRCorr_20260123_00001.mat");
twop0002 = load("Z:\ibn-vision\DATA\SUBJECTS\M25126\Analysis\20260123\M25126_20260123_processed2PData_M25126_BaselineCorridor_20260123_00002.mat");

response0001 = load("Z:\ibn-vision\DATA\SUBJECTS\M25126\Analysis\20260123\M25126_20260123_Response_M25126_VRCorr_20260123_00001.mat");
response0002 = load("Z:\ibn-vision\DATA\SUBJECTS\M25126\Analysis\20260123\M25126_20260123_Response_M25126_BaselineCorridor_20260123_00002.mat");



lapActivity = response0001.lapPositionActivity.dFF; 
w = gausswin(5); w = w / sum(w);

for iCell = 1:size(lapActivity, 1)
    for iLap = 1:size(lapActivity, 2)
        trace = squeeze(lapActivity(iCell, iLap, :));
        if all(isnan(trace)), continue; end
        nanMask = isnan(trace);
        trace(nanMask) = 0;
        smoothed = filtfilt(w, 1, trace);
        smoothed(nanMask) = NaN;
        lapActivity(iCell, iLap, :) = smoothed;
    end
end

oddLaps = lapActivity(:, 1:2:end, :);
evenLaps = lapActivity(:, 2:2:end, :);

% Average across laps
meanOdd = squeeze(mean(oddLaps, 2, 'omitnan'));
meanEven = squeeze(mean(evenLaps, 2, 'omitnan'));

% Normalize across position bins
normOdd = normalize(meanOdd, 2, 'range');
normEven = normalize(meanEven, 2, 'range');

% Sort cells by peak location in odd lap average
[~, peakIdx] = max(normOdd, [], 2);
[~, sortIdx] = sort(peakIdx);


roisPeakEarly = sortIdx(1:20);
roisPeakLate = sortIdx(115:135);



lapActivity2 = response0002.lapPositionActivity.dFF; 
w = gausswin(5); w = w / sum(w);

for iCell = 1:size(lapActivity2, 1)
    for iLap = 1:size(lapActivity2, 2)
        trace = squeeze(lapActivity2(iCell, iLap, :));
        if all(isnan(trace)), continue; end
        nanMask = isnan(trace);
        trace(nanMask) = 0;
        smoothed = filtfilt(w, 1, trace);
        smoothed(nanMask) = NaN;
        lapActivity2(iCell, iLap, :) = smoothed;
    end
end

oddLaps = lapActivity2(:, 1:2:end, :);
evenLaps = lapActivity2(:, 2:2:end, :);

% Average across laps
meanOdd2 = squeeze(mean(oddLaps, 2, 'omitnan'));
meanEven2 = squeeze(mean(evenLaps, 2, 'omitnan'));

% Normalize across position bins
normOdd2 = normalize(meanOdd2, 2, 'range');
normEven2 = normalize(meanEven2, 2, 'range');


figure; imagesc(normOdd2(roisPeakLate, :))
figure; imagesc(normOdd(roisPeakLate, :))

n=6;
figure; 
plot(twop0002.TwoPFrameTime,twop0002.processedSignals.dFF(roisPeakEarly(n),:)*1000)
hold on;
plot(twop0002.TwoPFrameTime, response0002.mouseVirtualPosition);
hold on;
xline(response0002.completedStartTimes);

figure;
imagesc(squeeze(response0002.lapPositionActivity.dFF(roisPeakEarly(n),:,:)))



figure; 
plot(twop0001.TwoPFrameTime,twop0001.processedSignals.dFF(roisPeakEarly(n),:)*1000)
hold on;
plot(twop0001.TwoPFrameTime, response0001.mouseVirtualPosition);
hold on;
xline(response0001.completedStartTimes);

figure;
imagesc(squeeze(response0001.lapPositionActivity.dFF(roisPeakEarly(n),:,:)))



subplot(121)
imagesc(normOdd(roisPeakEarly, 50:140));
subplot(122)
imagesc(normOdd2(roisPeakEarly, 50:140));