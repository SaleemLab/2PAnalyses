function [aligned] = alignPupilData(PupilData, TwoP_LastSyncPulseTime, TwoP_TwoPFrameTime)
% ALIGNPUPILDATA Aligns pupil camera data to 2P imaging timebase
%
% PupilData: The 'peripheralData.Pupil' structure containing .raw and .int fields
% TwoP_LastSyncPulseTime: Sync pulse times from the 2P microscope log
% TwoP_TwoPFrameTime: The master timebase (timestamps for each 2P frame)

aligned = struct();

%% Sync arduino pupil time clock with twop clock 

uSyncEye = unique(PupilData.raw.LastSyncPulseTime);
uSyncTwoP = unique(TwoP_LastSyncPulseTime);

% Calculate the new timebase for the eye data aligned to the 2P clock
% Note: Using .raw.ArduinoTime as the reference for the EyeCam timestamps
newArduinoTime = align2PSyncPulses(uSyncEye, uSyncTwoP, PupilData.raw.ArduinoTime);

%% 2. Interpolate Pupil Features to 2P Frame Times
% We take the data already interpolated to Arduino time (.int) 
% and move it to the 2P frame timebase (.TwoPFrameTime)

aligned.Centroid_X = interp1(newArduinoTime, PupilData.int.CentroidX, TwoP_TwoPFrameTime, 'linear', 'extrap');
aligned.Centroid_Y = interp1(newArduinoTime, PupilData.int.CentroidY, TwoP_TwoPFrameTime, 'linear', 'extrap');
aligned.Area       = interp1(newArduinoTime, PupilData.int.Area,      TwoP_TwoPFrameTime, 'linear', 'extrap');

aligned.MajorAxisLength = interp1(newArduinoTime, PupilData.int.MajorAxisLength, TwoP_TwoPFrameTime, 'linear', 'extrap');
aligned.MinorAxisLength = interp1(newArduinoTime, PupilData.int.MinorAxisLength, TwoP_TwoPFrameTime, 'linear', 'extrap');

% Store the timebase used for reference
aligned.TimeBase = TwoP_TwoPFrameTime;

end