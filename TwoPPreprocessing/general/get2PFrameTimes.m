function [sessionFileInfo] = get2PFrameTimes(sessionFileInfo, nPlanesPreZCorrection)
% [sessionFileInfo] = get2PFrameTimes(sessionFileInfo)
% 
% animal_name = 'M24049';
% session_number = '20240816';
if nargin>1
    nPlanes = nPlanesPreZCorrection; 
else 
    nPlanes = sessionFileInfo.numPlanes;
end

rootDir = sessionFileInfo.Directories.rootDir;
save_folder = fullfile(rootDir, 'Analysis', sessionFileInfo.session_name);
save_fileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_sessionFileInfo.mat'];


%% This gets all the tiff lengths
tiff_struct = get_tiff_lengths(sessionFileInfo);

% %%% debugging mode %%%
% tiff_struct = struct();
% load('Z:\ibn-vision\DATA\SUBJECTS\M24048\Analysis\20240816\tiff_struct_temp.mat');


%% ----- Initialize a structure to hold the plane frame counts -----
% Specify the number of planes
% save('tiff_struct_temp','tiff_struct');

plane_names = arrayfun(@(x) ['plane' num2str(x)], 0:nPlanes-1, 'UniformOutput', false);
% stimuli = fieldnames(sessionFileInfo.stimFiles.name);
nStimuli = length(sessionFileInfo.stimFiles);

% Loop through each stimulus and calculate frame counts for each plane
for iStim = 1:nStimuli
    stimuli{iStim} = sessionFileInfo.stimFiles(iStim).name; 
    stim_name = sessionFileInfo.stimFiles(iStim).name; 
    frame_range = tiff_struct.(stim_name);
    base_count = floor((diff(frame_range)+1)/nPlanes);
    extra_frames = rem((diff(frame_range)+1),nPlanes);
    counts = base_count.*ones(1, nPlanes);  % Initialize counts for the current stimulus
    % Calculate frame indices for each plane
    for iPlane = 1:nPlanes
        if iPlane<=extra_frames
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

% Initialize a cell array to store the data for the table
framerun_table_data = cell(nPlanes, nStimuli);

% Loop through each plane and each stimulus to extract the frame ranges
for kPlane = 1:nPlanes
    plane_name = plane_names{kPlane};  % 'plane0', 'plane1', etc.

    for iStim = 1:nStimuli
        stim_name = sessionFileInfo.stimFiles(iStim).name;% stimuli{iStim};  % 'VRCorr', 'SparseNoise', 'DotMotion_SpeedTuning', 'DirTuning'

        % Extract the frame range for the current plane and stimulus
        frame_range = framerun.(plane_name).(stim_name);

        % Store the frame range in the cell array
        framerun_table_data{kPlane, iStim} = {frame_range};  % Store as a cell to match the desired format
    end
end

% Create the table with row and column names
stim_frameruns = cell2table(framerun_table_data, 'VariableNames', stimuli, 'RowNames', plane_names);

% Display the final table
disp(stim_frameruns);
sessionFileInfo.stim_framerun = stim_frameruns;

% % % Save the table and the save directory.. 
sessionFileInfo.Directories.save_folder = save_folder; 
save([save_folder filesep save_fileName], 'sessionFileInfo',  '-append');


%% Use metadata from .tif file to count the number of frames per stimulus
    function tiff_struct = get_tiff_lengths(sessionFileInfo)

        cumulative_frames = 0;
        tiff_struct = struct();
        for iStimu = 1:length(sessionFileInfo.stimFiles)
            stimName = sessionFileInfo.stimFiles(iStimu).name;
            disp(['Counting frames for stimulus:', stimName])
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
            field_name = stimName; % Use the corresponding name from sorted_stim_names
            tiff_struct.(field_name) = [start_range, end_range];
            cumulative_frames = end_range;
        end
        disp('Tiff Ranges by Recording Name:');
        disp(tiff_struct);
    end
end