%% Bin-wise speed stratification 
target = 11; 

ROIData = squeeze(response2.lapPositionActivity.dFFNeuropilCorrected(target, :, :)); 
speedData  = response2.lapPositionRunningSpeed; % Speed per position bin

% Initialize result vectors
lowCurve  = nan(1, 200);
medCurve  = nan(1, 200);
highCurve = nan(1, 200);

allSpeeds = speedData(:);
activeSpeeds = allSpeeds(allSpeeds > 1);

T = quantile(activeSpeeds, [0.33, 0.66]);
lowThresh  = T(1); 
highThresh = T(2); 

fprintf('Data-Driven Thresholds: Low < %.1f, Med < %.1f, High > %.1f\n', ...
         lowThresh, highThresh, highThresh);



for thisbin = 1:200
    % Get all 57 speed samples for THIS position bin
    binSpeeds = speedData(:, thisbin);
    binActivity = ROIData(:, thisbin);
    
    % Create masks for this specific bin
    isLow  = binSpeeds > 1 & binSpeeds < lowThresh;
    isMed  = binSpeeds >= lowThresh & binSpeeds < highThresh;
    isHigh = binSpeeds >= highThresh;
    
    % Average the activity only if there are enough samples (e.g., > 2)
    if sum(isLow) > 2,  lowCurve(thisbin)  = mean(binActivity(isLow));  end
    if sum(isMed) > 2,  medCurve(thisbin)  = mean(binActivity(isMed));  end
    if sum(isHigh) > 2, highCurve(thisbin) = mean(binActivity(isHigh)); end
end

figure; hold on;

plot(smoothdata(lowCurve, 'gaussian', 9), 'c', 'LineWidth', 2); 
plot(smoothdata(medCurve, 'gaussian', 9), 'k', 'LineWidth', 2); 
plot(smoothdata(highCurve, 'gaussian', 9), 'm', 'LineWidth', 2);

legLow  = sprintf('Low < %.1f cm/s)', lowThresh);
legMed  = sprintf('Med < %.1f cm/s)', highThresh);
legHigh = sprintf('High > %.1f cm/s)', highThresh);


title(sprintf('Bouton %d - Bin-wise running stratification', target));
xlabel('Position (cm)', 'FontSize',12); 
ylabel('\DeltaF/F (Neu)', 'FontSize',12);
xline([40, 80, 120, 160], '--k', 'Alpha', 0.4, 'LineWidth', 1.5);
legend({legLow, legMed, legHigh}, 'Location', 'best');
defaultAxesProperties(gca, true);

%%  Using discrite bins -- version this 

currentROIData = ROIData;   
currentSpeedData = speedData; 
numPosBins = 200; 

%% 10 LINEAR Speed Bins
minSpeed = min(currentSpeedData(:), [], 'omitnan');; % They often group everything <= 1 into the first bin
maxSpeed = max(currentSpeedData(:), [], 'omitnan');

% Create LINEAR bins (Linear spacing like Figure 1e/f)
% This creates 10 equal-sized bins (e.g., 1-5, 5-10, 10-15...)
speedBins = linspace(minSpeed, maxSpeed, 11); 

% Calculate Centers
nSpeedBins = length(speedBins) - 1;
speedCenters = (speedBins(1:end-1) + speedBins(2:end)) / 2;


tuningSurfaceRaw = nan(nSpeedBins, numPosBins);
%%
%%
for b = 1:numPosBins
    % Extract speed and activity for this 1cm chunk
    binSpeeds = currentSpeedData(:, b);
    binActivity = currentROIData(:, b);
    
    for s = 1:nSpeedBins
        %  Find the logical index for this Speed/Position box
        idx = binSpeeds >= speedBins(s) & binSpeeds < speedBins(s+1);
        
        % Count samples
        sampleCount = sum(idx);
        
        % Apply Noise Filter; Increase sample to 5? Or use 7
        if sampleCount >= 3
            tuningSurfaceRaw(s, b) = mean(binActivity(idx), 'omitnan');
        else
            tuningSurfaceRaw(s, b) = nan; 
        end
    end
end

% Apply a 2D Gaussian filter to the surface
% [0.5, 1.2] means: smooth slightly in speed (0.5) and more in position (1.2)
tuningSurfaceRaw = imgaussfilt(tuningSurfaceRaw, [0.1, 0.1], 'Padding', 'replicate');

%% 
figure('Name', sprintf('Bouton %d - 10 Bin Raw 1cm', target), 'Position', [100 100 650 400]);


imagesc(1:200, speedCenters, tuningSurfaceRaw);
set(gca, 'YDir', 'normal'); 


% 
colormap(parula);
c = colorbar;
c.Label.String = '\DeltaF/F [NeuC]';

activeData = tuningSurfaceRaw(tuningSurfaceRaw > 0);
if ~isempty(activeData)
    maxVal = prctile(activeData, 98.5);
else
    maxVal = 1;
end
set(gca, 'CLim', [0, maxVal]);


hold on;
xline([40, 80, 120, 160, 200], '--w', 'Alpha', 0.4, 'LineWidth', 1.5);
xticks([1 40 80 120 160 200]);
xticklabels({'1', '40', '80', '120', '160', '200'});
set(gca, 'YScale', 'log'); 
yticks([2, 5, 10, 20, 30]); % Match the paper's tick style
yticklabels({'2', '5', '10', '20', '30'});

title(['Bouton ' num2str(target) ': 10 Speeds (1cm Bins)']);
xlabel('Position (cm)');
ylabel('Running Speed (cm/s)');
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);



%% - OCCUPANCY 

debugOccupancy = zeros(nSpeedBins, 100);

% Calculate Sample Counts 
for b = 1:200
    binSpeeds = currentSpeedData(:, b);
    for s = 1:nSpeedBins
        idx = binSpeeds >= speedBins(s) & binSpeeds < speedBins(s+1);
        debugOccupancy(s, b) = sum(idx); % Count of samples
    end
end


figure('Name', 'Occupancy');
imagesc(2:2:200, speedCenters, debugOccupancy);
set(gca, 'YDir', 'normal');
colormap(hot); % Hot colormap: Black = 0 samples (Gaps), Yellow/White = Lots of data
c = colorbar;
c.Label.String = 'Number of Samples (Time points)';

hold on; 
xline([40, 80, 120, 160], '--w', 'Alpha', 0.4, 'LineWidth', 1.5);
title('Occupancy Map (Dark Areas = No Data Recorded)');




%% Same as above but 1cm bin  --- > use this!!!! 
currentROIData = ROIData;   
currentSpeedData = speedData; 


numPosBins = 200; 

%% estimate speed bins
allSpeeds = currentSpeedData(:);
runningIdx = allSpeeds > 1 & ~isnan(allSpeeds);
runningSpeeds = allSpeeds(runningIdx);

% Since we are using 1cm bins, we have fewer samples per spatial bin.
% We adjust the target samples per box to 5 to help fill the gaps.
targetSamplesPerBox = 5; 
nSpeedBins = floor(length(runningSpeeds) / (numPosBins * targetSamplesPerBox));
nSpeedBins = max(3, min(nSpeedBins, 10)); 

fprintf('Using %d speed bins for 1cm resolution.\n', nSpeedBins);

% calculate qunatile edges 
speedBins = quantile(runningSpeeds, linspace(0, 1, nSpeedBins + 1));
speedCenters = (speedBins(1:end-1) + speedBins(2:end)) / 2;

%% 1cm
tuningSurfaceRaw = nan(nSpeedBins, numPosBins);

for b = 1:numPosBins
    binSpeeds = currentSpeedData(:, b);
    binActivity = currentROIData(:, b);
    for s = 1:nSpeedBins
        idx = binSpeeds >= speedBins(s) & binSpeeds < speedBins(s+1);
        
        % Using threshold of 2 to ensure we see the pillar structure in sparse data
        if sum(idx) >= 2 
            tuningSurfaceRaw(s, b) = mean(binActivity(idx), 'omitnan');
        end
    end
end

mask = ~isnan(tuningSurfaceRaw);

% Replace NaNs with 0s so the math works
dataZeroed = tuningSurfaceRaw;
dataZeroed(isnan(tuningSurfaceRaw)) = 0;

% Filter both the data and the mask
blurredData = imgaussfilt(dataZeroed, [1, 1.2], 'Padding', 'replicate');
blurredMask = imgaussfilt(double(mask), [1, 1.2], 'Padding', 'replicate');

%Divide data by mask (this re-normalizes the average)
tuningSurfaceRaw = blurredData ./ blurredMask;

%%  Visualization
figure('Name', 'Bouton 11 - 1cm Log Quantile', 'Position', [100 100 750 500]);

% xCoords now corresponds to the 1cm bins (1 to 200)
xCoords = 1:numPosBins; 
[X, Y] = meshgrid(xCoords, speedCenters);

h = pcolor(X, Y, tuningSurfaceRaw);
shading flat; 

% set log scale and direction
set(gca, 'YScale', 'log', 'YDir', 'normal');

% include minor ticks 
set(gca, 'YMinorTick', 'on', 'TickDir', 'out');

colormap(parula); 
set(gca, 'Color', [0.15 0.15 0.15]); % Gray for empty areas

% Brighten the contrast
activeData = tuningSurfaceRaw(~isnan(tuningSurfaceRaw));
if ~isempty(activeData)
    set(gca, 'CLim', [0, prctile(activeData, 99)]);
end

% 
tickVals = [1, 2, 5, 10, 20, 40, 80];
yticks(tickVals); 
yticklabels(cellstr(num2str(tickVals')));

%
xlim([1, 200]);
ylim([min(speedCenters), max(speedCenters)]);

% Formatting
hold on;
xline([40, 80, 120, 160], '--w');
xticks(0:40:200);

ylabel('Running Speed (cm/s)'); 
xlabel('Position (cm)');
title(sprintf('Bouton %s: %d Speed Bins / 1cm Position', num2str(target), nSpeedBins));
grid off;
box off;
c = colorbar;
c.Label.String = '\DeltaF/F';

%%
%% occupancy matrix [sanity check]
% We use the same bins and centers from  previous quantile calculation
debugOccupancy = zeros(nSpeedBins, numPosBins);

for b = 1:numPosBins
    binSpeeds = currentSpeedData(:, b);
    for s = 1:nSpeedBins
        % Logic: Count how many frames fall into this Speed x Position box
        idx = binSpeeds >= speedBins(s) & binSpeeds < speedBins(s+1);
        debugOccupancy(s, b) = sum(idx); 
    end
end

%% 
figure('Name', 'Quantile Normalization Check', 'Position', [150 150 900 450]);


subplot(1, 4, 1:3);
imagesc(1:numPosBins, speedCenters, debugOccupancy);
set(gca, 'YScale', 'log', 'YDir', 'normal');
colormap(hot);
c = colorbar('Location', 'southoutside');
c.Label.String = 'Samples per Bin (Frames)';
xlabel('Position (cm)'); ylabel('Speed (cm/s)');
title('Speed-Position Occupancy (Equal-Sample Bins)');


subplot(1, 4, 4);
totalSamplesPerSpeedBin = sum(debugOccupancy, 2);
barh(speedCenters, totalSamplesPerSpeedBin, 'FaceColor', [0.8 0.2 0.2]);
set(gca, 'YScale', 'log', 'YDir', 'normal');
xlabel('Total Samples');
title('Total Data per Row');
grid on;