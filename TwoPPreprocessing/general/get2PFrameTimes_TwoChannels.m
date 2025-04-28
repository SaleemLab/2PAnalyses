function [sessionFileInfo] = get2PFrameTimes_TwoChannels(sessionFileInfo, nPlanesPreZCorrection, nChannels)
% [sessionFileInfo] = get2PFrameTimes_TwoChannels(sessionFileInfo, nPlanesPreZCorrection, nChannels)
%
% This function calculates how many frames belong to each imaging plane
% for each stimulus, taking into account the number of imaging channels.
%
% Arguments:
%   - sessionFileInfo: struct containing session and stimulus information
%   - nPlanesPreZCorrection: number of planes (optional)
%   - nChannels: number of imaging channels (optional, default = 1)

if nargin > 1
    nPlanes = nPlanesPreZCorrection;
else
    nPlanes = sessionFileInfo.numPlanes;
end

if nargin < 3
    nChannels = 1; % default to 1 channel
end

rootDir = sessionFileInfo.Directories.rootDir;
save_folder = fullfile(rootDir, 'Analysis', sessionFileInfo.session_name);
save_fileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_sessionFileInfo.mat'];

%% This gets all the tiff lengths
tiff_struct = get_tiff_lengths(sessionFileInfo);

%% ----- Initialize a structure to hold the plane frame counts -----
plane_names = arrayfun(@(x) ['plane' num2str(x)], 0:nPlanes-1, 'UniformOutput', false);
nStimuli = length(sessionFileInfo.stimFiles);

% Loop through each stimulus and calculate frame counts for each plane
for iStim = 1:nStimuli
    stimuli{iStim} = sessionFileInfo.stimFiles(iStim).name;
    stim_name = sessionFileInfo.stimFiles(iStim).name;
    frame_range = tiff_struct.(stim_name);
    totalFrames = diff(frame_range) + 1;

    framesPerPlane = floor(totalFrames / (nPlanes * nChannels));
    extraFrames = rem(totalFrames, nPlanes * nChannels);

    counts = framesPerPlane .* ones(1, nPlanes);  % Initialize counts for the current stimulus

    for iPlane = 1:nPlanes
        if iPlane <= floor(extraFrames / nChannels)
            counts(iPlane) = counts(iPlane) + 1;
        end
        perPlaneCounts{iPlane}(iStim) = counts(iPlane);
    end
end

for iPlane = 1:nPlanes
    end_idxs{iPlane} = cumsum(perPlaneCounts{iPlane});
    start_counts{iPlane} = [1 end_idxs{iPlane}(1:end-1)+1];
end

for iStim = 1:nStimuli
    stim_name = sessionFileInfo.stimFiles(iStim).name;
    for iPlane = 1:nPlanes
        framerun.(plane_names{iPlane}).(stim_name) = [start_counts{iPlane}(iStim) end_idxs{iPlane}(iStim)];
    end
end

disp('Frame Counts by Plane for Each Stimulus:');

%% Restructure Struct to Cell
framerun_table_data = cell(nPlanes, nStimuli);
for kPlane = 1:nPlanes
    plane_name = plane_names{kPlane};
    for iStim = 1:nStimuli
        stim_name = sessionFileInfo.stimFiles(iStim).name;
        frame_range = framerun.(plane_name).(stim_name);
        framerun_table_data{kPlane, iStim} = {frame_range};  % Store as cell
    end
end

stim_frameruns = cell2table(framerun_table_data, 'VariableNames', stimuli, 'RowNames', plane_names);
disp(stim_frameruns);

sessionFileInfo.stim_framerun = stim_frameruns;
sessionFileInfo.Directories.save_folder = save_folder;
save([save_folder filesep save_fileName], 'sessionFileInfo', '-append');

%% Nested function: Count TIFF frames for each stimulus
function tiff_struct = get_tiff_lengths(sessionFileInfo)
    cumulative_frames = 0;
    tiff_struct = struct();
    for iStimu = 1:length(sessionFileInfo.stimFiles)
        stimName = sessionFileInfo.stimFiles(iStimu).name;
        disp(['Counting frames for stimulus: ', stimName])
        num_frames_in_stim = 0;
        tif_files = sessionFileInfo.stimFiles(iStimu).tif_filepaths;
        for iTiff = 1:length(tif_files)
            disp(['Counting frames for tiff file ', num2str(iTiff)])
            tic;
            info = imfinfo(tif_files{iTiff});
            disp(['Reading Tiff took ', num2str(toc)])
            num_frames_in_stim = num_frames_in_stim + numel(info);
        end
        start_range = cumulative_frames + 1;
        end_range = cumulative_frames + num_frames_in_stim;
        field_name = stimName;
        tiff_struct.(field_name) = [start_range, end_range];
        cumulative_frames = end_range;
    end
    disp('Tiff Ranges by Recording Name:');
    disp(tiff_struct);
end
end
