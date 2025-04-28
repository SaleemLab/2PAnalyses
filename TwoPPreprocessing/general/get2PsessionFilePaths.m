function [sessionFileInfo] = get2PsessionFilePaths(animal_name, session_name, stim_list, rerun_process, fileNameAddition)
% Function to collect all the file paths for Bonsai, ophys(tiff) and Suite2p
% files.
%
% Usage: [sessionFileInfo] = get2PsessionFilePaths(animal_name, session_name, (optional) stim_list)
% Example: sessionFileInfo = get2PsessionFilePaths('M24048', '20240816')
% Example 'stim_list': stim_list = {'VRCorr', 'SparseNoise', 'DotMotion_SpeedTuning', 'DirTuning'};
% stim_list needs to have the list in order of presentation
%
% Sonali and Aman - Oct 2024

if nargin<5
    fileNameAddition = [];
end

if nargin<4
    rerun_process = 0;
end

if nargin<3
    % Stim list
    stim_list = {'VRCorr', 'SparseNoise', 'DotMotion_SpeedTuning', 'DirTuning'};
end

sessionFileInfo.animal_name = animal_name;
sessionFileInfo.session_name = session_name;


% Directories
rootDir = ['Z:' filesep fullfile('ibn-vision','DATA','SUBJECTS',animal_name)];
bonsai_folder = fullfile(rootDir, 'Bonsai', session_name);
ophys_folder = fullfile(rootDir, 'Ophys', session_name);
processed_folder = fullfile(rootDir, 'Processed', session_name);
suite2p_folder = fullfile(processed_folder, 'suite2p');

save_folder = fullfile(rootDir, 'Analysis', session_name);
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
save_fileName = [animal_name '_' session_name '_sessionFileInfo' fileNameAddition '.mat'];

if exist([save_folder filesep save_fileName])~=2 || rerun_process == 1
    planes = dir(fullfile(suite2p_folder, 'plane*'));
    sessionFileInfo.numPlanes = length(planes);

    sessionFileInfo.Directories.rootDir = rootDir;
    sessionFileInfo.Directories.bonsai = bonsai_folder;
    sessionFileInfo.Directories.ophys = ophys_folder;
    sessionFileInfo.Directories.suite2p = suite2p_folder;
    sessionFileInfo.Directories.save_folder = fullfile(rootDir, 'Analysis', sessionFileInfo.session_name);

    % Stim list
%     stim_list = {'FOV3_SurfaceExpressionTest_FOV1_00001'}; %VRCorridor
    % Define files to exclude
    exclude = {'excluded_file1.csv', 'excluded_file2.tif'};

    %% Process 'Bonsai' files (CSV and BIN files)
    disp('Collecting Bonsai')
    for iStim = 1:length(stim_list)
        sessionFileInfo.stimFiles(iStim).name = stim_list{iStim};

        % Process CSV files
        csvFiles = dir(fullfile(bonsai_folder, '*.csv'));
        sessionFileInfo.stimFiles(iStim).bonsai_filepaths = {};
        for iPlane = 1:length(csvFiles)
            csvFilePath = fullfile(bonsai_folder, csvFiles(iPlane).name);
            if ~ismember(csvFiles(iPlane).name, exclude) && contains(csvFiles(iPlane).name, stim_list{iStim})
                sessionFileInfo.stimFiles(iStim).bonsai_filepaths{end+1} = csvFilePath;
            end
        end

        % Process BIN files
        binFiles = dir(fullfile(bonsai_folder, '*.bin'));
        for iPlane = 1:length(binFiles)
            binFilePath = fullfile(bonsai_folder, binFiles(iPlane).name);
            if ~ismember(binFiles(iPlane).name, exclude) && contains(binFiles(iPlane).name, stim_list{iStim})
                sessionFileInfo.stimFiles(iStim).bonsai_filepaths{end+1} = binFilePath;
            end
        end

        %% Process 'Ophys' TIFF files
        disp('Collecting Ophys Tiff')
        % Look one more level down for TIFF files

        stimDirs = dir(fullfile(ophys_folder, ['*', sessionFileInfo.stimFiles(iStim).name, '*']));
        for iStimDirs = 1:length(stimDirs)
            stimDirPath = fullfile(ophys_folder, stimDirs(iStimDirs).name);
            sessionFileInfo.stimFiles(iStim).tif_filepaths = {};
            if isfolder(stimDirPath)
                tifFiles = dir(fullfile(stimDirPath, '*.tif'));
                for l = 1:length(tifFiles)
                    tifFilePath = fullfile(stimDirPath, tifFiles(l).name);
                    sessionFileInfo.stimFiles(iStim).tif_filepaths{end+1} = tifFilePath;
%                     if ~ismember(tifFiles(l).name, exclude) && contains(tifFiles(l).name, session_name)
%                         sessionFileInfo.stimFiles(iStim).tif_filepaths{end+1} = tifFilePath;
%                     end
                end
            end
        end
    end

    %% Process 'Processed' Suite2p files (plane0, plane1, etc.)
    disp('Collecting suite2p processed files')
    planes = dir(fullfile(suite2p_folder, 'plane*'));
    for iPlane = 1:length(planes)
        disp(['Plane: ' num2str(iPlane)])
        planePath = fullfile(suite2p_folder, planes(iPlane).name);
        sessionFileInfo.suite2pFiles(iPlane).planeName = planes(iPlane).name;

        % Initialize a column for suite2p file paths
        suite2p_column = {};

        % Get all files in the plane directory
        allFiles = dir(fullfile(planePath, '*'));
        allFiles = allFiles(~ismember({allFiles.name}, {'.', '..'}));
        for l = 1:length(allFiles)
            filePath = fullfile(planePath, allFiles(l).name);
            suite2p_column{end+1} = filePath; % Append file path
        end
        sessionFileInfo.suite2pFiles(iPlane).planes = suite2p_column;
    end

    % Display the structured session_planes data
    disp(sessionFileInfo);
    % save([save_folder filesep animal_name '_vr_session_planes_info.mat'], 'sessionFileInfo');
    sessionFileInfo.sessionFileInfo_filepath = [save_folder filesep save_fileName];
    save([save_folder filesep save_fileName], 'sessionFileInfo');
   
else
    load([save_folder filesep save_fileName], 'sessionFileInfo');
end