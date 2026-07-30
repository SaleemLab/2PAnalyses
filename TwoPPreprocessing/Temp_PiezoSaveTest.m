filePath = 'Z:\ibn-vision\USERS\Sonali\SS131-EDDBoutons\2025-02-04 A\2\NiDaqInput0.bin';

fid = fopen(filePath, 'r');
if fid == -1, error('File not found'); end
rawBuffer = fread(fid, Inf, '*int16'); 
fclose(fid);

numChannels = 5;
% each row is a channel, each column is a value in time?
dataMatrix = reshape(rawBuffer, numChannels, []);


niDaq.photodiode = double(dataMatrix(1, :));
niDaq.frameclock = double(dataMatrix(2, :));
niDaq.pockel     = double(dataMatrix(3, :));
niDaq.piezo      = double(dataMatrix(4, :));
niDaq.sync       = double(dataMatrix(5, :));

figure('Name', 'NiDaq Channel Check');
t = 1:min(100, length(niDaq.photodiode)); % Plot first 10k samples

subplot(3,1,1); plot(niDaq.photodiode(t)); title('Photodiode');
subplot(3,1,2); plot(niDaq.frameclock(t)); title('Frame Clock');
subplot(3,1,3); plot(niDaq.piezo(t));      title('Piezo');
xlabel('Samples');



%% load piezo data and two-p frame times
piezoData = readtable("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\PiezoReadTest\PRT-piezo\EyeTracking\20260203\PRT_PRT8P2S_20260203_00001_Piezo2026-02-03T18_46_45.csv");
frameData = readtable("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\PiezoReadTest\PRT-frametimes\Bonsai\20260203\PRT_PRT8P2S_20260203_00001_2P2026-02-03T18_47_10.csv");

%%  convert ADC to microns 
piezoData.Microns = (double(piezoData.Piezo) / 1023) * 400;

figure('Name', 'Piezo Signal in Microns');
plot(piezoData.ArduinoTime / 1000, piezoData.Microns);
xlabel('Arduino Time (seconds)');
ylabel('Position (microns)');
title('Full Piezo Trace Converted to Microns');


%% align piezo to twopframe times 
piezoPulses = unique(piezoData.LastSyncPulsetime);
framePulses = unique(frameData.LastSyncPulseTime);
pIntervals = diff(piezoPulses);
fIntervals = diff(framePulses);

numToMatch = 10;
minError = inf;
bestShift = 0;
for s = 1:(length(pIntervals) - numToMatch)
    currentError = sum(abs(pIntervals(s:s+numToMatch-1) - fIntervals(1:numToMatch)));
    if currentError < minError
        minError = currentError;
        bestShift = s;
    end
end

globalOffset = piezoPulses(bestShift) - framePulses(1);
frameData.AlignedTime = frameData.TwoPFrameTime + globalOffset;

%% extract planes
numPlanes = 8;
frameDurationMs = 16; %
numFrames = height(frameData);
numVolumes = floor(numFrames / numPlanes);


frameData.PlaneIdx = mod(0:numFrames-1, numPlanes)';

figure('Name', 'Average Piezo Shape per Plane');
hold on;
colors = lines(numPlanes);
timeX = (0:0.1:frameDurationMs); % Sub-millisecond time axis for plotting

for p = 0:numPlanes-1
    % all frames assigned to this plane
    planeFrames = find(frameData.PlaneIdx == p);
    allSegments = [];
    
    for i = 1:length(planeFrames)
        fIdx = planeFrames(i);
        tStart = frameData.AlignedTime(fIdx);
        tEnd = tStart + frameDurationMs;
        
        %piezo samples during this specific frame
        mask = piezoData.ArduinoTime >= tStart & piezoData.ArduinoTime < tEnd;
        seg = piezoData.Microns(mask);
        
        % Normalize length for averaging (resample to 100 points)
        if ~isempty(seg)
            allSegments = [allSegments, interp1(linspace(0, frameDurationMs, length(seg)), seg, timeX)'];
        end
    end
    
    % Plot the average shape for this plane
    avgShape = mean(allSegments, 2, 'omitnan');
    plot(timeX, avgShape, 'Color', colors(p+1,:), 'LineWidth', 2, 'DisplayName', ['Plane ' num2str(p)]);
end

xlabel('Time within Frame (ms)');
ylabel('Position (microns)');
title('Average piezo trajectory for each plane');
legend('show', 'Location', 'eastoutside');



%%
function [alignedEyeData] = alignEyeData(EyeData,EyeTimestamps,twoPLog)
% SGS SDL 02/2026
% Function to align the Eyedata collected in 2p rig with the rest of the
% data
%
% EyeData = EyeCamLog file
% EyeTimestamps = EyeCamTimeStamps
% TwoPLog = 2P bonsai file
intEyeData = struct();
alignedEyeData = struct();
% Remove duplicates of samples in EyeTimeStamps
[tt,tp] = unique(EyeTimestamps.ArduinoTime);
EyeTimestamps = EyeTimestamps(tp,:);
% Interpolate to the EyeTimestamps X all the EyeData
intEyeData.Centroid_X = interp1(EyeData.eyeMsSinceStartOfDay,EyeData.Centroid_X,EyeTimestamps.EyeCamTime,'linear','extrap');
intEyeData.Centroid_Y = interp1(EyeData.eyeMsSinceStartOfDay,EyeData.Centroid_Y,EyeTimestamps.EyeCamTime,'linear','extrap');
intEyeData.Area = interp1(EyeData.eyeMsSinceStartOfDay,EyeData.Area,EyeTimestamps.EyeCamTime,'linear','extrap');
intEyeData.MajorAxisLength = interp1(EyeData.eyeMsSinceStartOfDay,EyeData.MajorAxisLength,EyeTimestamps.EyeCamTime,'linear','extrap');
intEyeData.MinorAxisLength = interp1(EyeData.eyeMsSinceStartOfDay,EyeData.MinorAxisLength,EyeTimestamps.EyeCamTime,'linear','extrap');
% Interpolate to twoPLog.TwoPFrameTime all the EyeTimestamps dependent data
uSyncEye = unique(EyeTimestamps.LastSyncPulseTime);
uSyncTwoP = unique(twoPLog.LastSyncPulseTime);
EyeTimestamps.newArduinoTime = align2PSyncPulses(uSyncEye,uSyncTwoP,EyeTimestamps.ArduinoTime);
alignedEyeData.Centroid_X = interp1(EyeTimestamps.newArduinoTime,intEyeData.Centroid_X,twoPLog.TwoPFrameTime,'linear','extrap');
alignedEyeData.Centroid_Y = interp1(EyeTimestamps.newArduinoTime,intEyeData.Centroid_Y,twoPLog.TwoPFrameTime,'linear','extrap');
alignedEyeData.Area = interp1(EyeTimestamps.newArduinoTime,intEyeData.Area,twoPLog.TwoPFrameTime,'linear','extrap');
alignedEyeData.MajorAxisLength = interp1(EyeTimestamps.newArduinoTime,intEyeData.MajorAxisLength,twoPLog.TwoPFrameTime,'linear','extrap');
alignedEyeData.MinorAxisLength = interp1(EyeTimestamps.newArduinoTime,intEyeData.MinorAxisLength,twoPLog.TwoPFrameTime,'linear','extrap');
alignedEyeData.TimeBase = twoPLog.TwoPFrameTime;

end 


