function [plane_data, plane_data_new] = get_bonsai_twopframetimes_by_planes(filepath, nplanes)
% [plane_data] = get_bonsai_twopframetimes_by_planes(filepath, nplanes)
% Function to read csv file of frame times and split into frametimes per plane.


% Read CSV file generated by Bonsai containing all the TwoP frame timestamps 
twop_csv = readtable(filepath);

% Temp clean up; remove rows where TwoPFrameTime is 0;
if twop_csv.TwoPFrameTime(1) == 0
    twop_csv = twop_csv(2:end, :);
end
twop_csv.TwoPFrameTime = twop_csv.TwoPFrameTime./1000;
twop_csv.LastSyncPulseTime = twop_csv.LastSyncPulseTime./1000;
twop_csv.ArduinoTime = twop_csv.ArduinoTime./1000;

% Initialize empty tables for each plane's data with the same variable names
for p = 1:nplanes
    plane_data.(['plane' num2str(p-1)]) = twop_csv(p:nplanes:end, :); 
    num_frames_per_plane(p) = size(plane_data.(['plane' num2str(p-1)]),1);
    plane_data_new{p} = twop_csv(p:nplanes:end, :);
end
% limit_frames = min(num_frames_per_plane);
% for p = 1:nplanes
%     plane_data.(['plane' num2str(p-1)]) = plane_data.(['plane' num2str(p-1)])(1:limit_frames,:); 
%     plane_data_new{p} = plane_data_new{p}(1:limit_frames,:);
% end
