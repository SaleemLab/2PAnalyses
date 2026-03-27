function [sessionFileInfo] = get2PFrameTimes_SpeedyVersion(sessionFileInfo, isZcorrected, nPlanes)
% Calculates stimulus frame ranges from Suite2p output for multiple planes.
%
% This function performs a direct 1-to-1 mapping. It assumes that the number
% of stimuli in 'sessionFileInfo.stimFiles' is equal to the number of TIFF
% files processed by Suite2p, as reported in 'ops.frames_per_folder'.


if nargin < 2
    error('A second argument, isZcorrected (true/false), is required.');
end
if nargin < 3
    twoPMetaData = load(sessionFileInfo.stimFiles(1).TwoPMetaData);
    nPlanes = twoPMetaData.twopMetadata.numSlices;
end 
rootDir = sessionFileInfo.Directories.rootDir;
save_folder = fullfile(rootDir, 'Analysis', sessionFileInfo.session_name);
save_fileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_sessionFileInfo.mat'];
suite2pDir = fullfile(rootDir, 'Processed', sessionFileInfo.session_name, 'suite2p');
nStimuli = length(sessionFileInfo.stimFiles);
stimuli_names = {sessionFileInfo.stimFiles.name};
names_to_remove = contains(stimuli_names, 'CombinedRuns') & contains(stimuli_names, 'VRCorr');

% Condition in case
if any(names_to_remove)
    disp('Found stimuli containing both CombinedRuns and VRCorr. Removing...')
    stimuli_names = stimuli_names(~names_to_remove);
    disp('Removal complete.')
else
    disp('No stimuli found containing both "CombinedRuns" and "VRCorr". List unchanged.')
end
framerun = struct();

%% --- Calculate Frame Runs Based on Processing Pipeline ---

if isZcorrected
    %  Z-motion correction applied (One 'plane' folder) ---
    disp('Z-corrected mode: Calculating frame runs once and replicating for all planes...');
    plane_folder_path = fullfile(suite2pDir, 'plane_z');
    if ~exist(plane_folder_path, 'dir')
        error('Z-corrected mode failed: Expected folder not found at %s', plane_folder_path);
    end
    fallMatPath = fullfile(plane_folder_path, 'fall.mat');
    s2p_data = load(fallMatPath, 'ops');
    frames_per_tiff = s2p_data.ops.frames_per_folder;

    % Calculate the frame runs for the single dataset
    stim_ranges = calculate_stim_ranges_direct_map(sessionFileInfo, frames_per_tiff);
    
    % Replicate these ranges for every plane
    plane_names_out = arrayfun(@(x) ['plane' num2str(x)], 0:nPlanes-1, 'UniformOutput', false);
    for iPlane = 1:nPlanes
        for iStim = 1:nStimuli
            framerun.(plane_names_out{iPlane}).(stimuli_names{iStim}) = stim_ranges.(stimuli_names{iStim});
        end
    end

else
    % No Z-motion correction (Multiple 'planeN' folders) ---
    disp('Multi-plane mode: Calculating frame runs for each plane individually...');
    planeItems = dir(fullfile(suite2pDir, 'plane*'));
    planeFolders = planeItems([planeItems.isdir]);
    if isempty(planeFolders)
        error('Multi-plane mode failed: No ''plane*'' folders found in %s', suite2pDir);
    end
    [~, sortIdx] = sort( cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), {planeFolders.name}) );
    planeFolders = planeFolders(sortIdx);

    if length(planeFolders) ~= nPlanes
        warning('Mismatch: Expected %d planes, but found %d folders. Adjusting to use %d.', nPlanes, length(planeFolders), length(planeFolders));
        nPlanes = length(planeFolders);
    end

    % Loop through each plane folder to get its unique frame counts
    for iPlane = 1:nPlanes
        plane_name = planeFolders(iPlane).name;
        fprintf('Processing %s...\n', plane_name);
        fallMatPath = fullfile(suite2pDir, plane_name, 'fall.mat');
        s2p_data = load(fallMatPath, 'ops');
        frames_per_tiff = s2p_data.ops.frames_per_folder;
        stim_ranges = calculate_stim_ranges_direct_map(sessionFileInfo, frames_per_tiff);
        framerun.(plane_name) = stim_ranges;
    end
    plane_names_out = {planeFolders.name};
end

%% --- Restructure Struct to Table for display and saving ---
disp('Final Frame Counts by Plane for Each Stimulus:');
framerun_table_data = cell(nPlanes, nStimuli);
for kPlane = 1:nPlanes
    plane_name = plane_names_out{kPlane};
    for iStim = 1:nStimuli
        stim_name = stimuli_names{iStim};
        frame_range = framerun.(plane_name).(stim_name);
        framerun_table_data{kPlane, iStim} = {frame_range};
    end
end
stim_frameruns = cell2table(framerun_table_data, 'VariableNames', stimuli_names, 'RowNames', plane_names_out);
disp(stim_frameruns);

%% --- Save the results back to the sessionFileInfo struct ---
sessionFileInfo.stim_framerun = stim_frameruns;
sessionFileInfo.Directories.save_folder = save_folder;
if ~exist(save_folder, 'dir'), mkdir(save_folder); end
save(fullfile(save_folder, save_fileName), 'sessionFileInfo', '-append');
fprintf('\nSuccessfully updated and saved sessionFileInfo.\n');

end

% --- REWRITTEN HELPER FUNCTION FOR DIRECT MAPPING ---
function stim_ranges = calculate_stim_ranges_direct_map(sessionFileInfo, frames_per_tiff)
    stim_ranges = struct();
    nStimuli = length(sessionFileInfo.stimFiles);

    % --- The NEW, CORRECT Safety Check ---
    % It compares the number of STIMULI to the number of frame counts.
    if nStimuli ~= length(frames_per_tiff)
        err_msg = sprintf(['Fatal Mismatch: The number of stimuli in sessionFileInfo (%d) does not match the number of TIFFs processed by Suite2p (%d).\n' ...
        'This function assumes that each stimulus entry corresponds to exactly one TIFF file in the Suite2p run.'], ...
        nStimuli, length(frames_per_tiff));
        error('get2PFrameTimes:mismatch', err_msg);
    end

    end_points = cumsum(frames_per_tiff);
    start_points = [1, end_points(1:end-1) + 1];

    for iStim = 1:nStimuli
        stimName = sessionFileInfo.stimFiles(iStim).name;
        stim_ranges.(stimName) = [start_points(iStim), end_points(iStim)];
    end
end