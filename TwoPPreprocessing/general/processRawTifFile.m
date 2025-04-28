function [tempsessionFileInfo] = processRawTifFile(animal_name, session_name)
% Usage: [tempsessionFileInfo] = processRawTifFile('M25012', '20250210')
% Aman and Sonali - Feb 2025
%% Create a temp sessionFileInfo struct (that is currently not saved) to find the tif paths 
rootDir = ['Z:' filesep fullfile('ibn-vision','DATA','SUBJECTS',animal_name)];
ophysFolder = fullfile(rootDir, 'Ophys', session_name);
entries = dir(ophysFolder);
% Filter out only directories (excluding '.' and '..')
StimDirs = {entries([entries.isdir] & ~ismember({entries.name}, {'.', '..'})).name};
% Remove subfolders that contain 'zstack'
StimDirs = StimDirs(~contains(StimDirs, 'zstack', 'IgnoreCase', true));

for iStimDirs = 1:length(StimDirs)
    stimDirPath = fullfile(ophysFolder, StimDirs{iStimDirs});
    tempsessionFileInfo.stimFiles(iStimDirs).tif_filepaths = {};
    if isfolder(stimDirPath)
        tifFiles = dir(fullfile(stimDirPath, '*.tif'));
        for l = 1:length(tifFiles)
            tifFilePath = fullfile(stimDirPath, tifFiles(l).name);
            tempsessionFileInfo.stimFiles(iStimDirs).tif_filepaths{end+1} = tifFilePath;
            if contains(tifFiles(l).name, session_name)
               tempsessionFileInfo.stimFiles(iStimDirs).tif_filepaths{end+1} = tifFilePath;
            end
        end
    end
end
   
%% Loop through all stimulus files

for iStim = 1:length(tempsessionFileInfo.stimFiles) 
    % Pick the last TIFF file in this stimulus acquisition
    lastTiffilePath = tempsessionFileInfo.stimFiles(iStim).tif_filepaths{end};
    fprintf('Processing tif files for stimulus: %s\n', StimDirs{iStim})

    % ----- Find total number of frames for this stimulus -------
    % Read current tiff file to extract the index of the last frame
    lastTifFile = Tiff(lastTiffilePath, 'r'); 
%     lastFrameNum = lastTifFile.numberOfDirectories();
    lastFrameNum = numel(imfinfo(lastTiffilePath));

    % Set directory to last frame and extract metadata
    lastTifFile.setDirectory(lastFrameNum);
    lastFrameDescription = lastTifFile.getTag('ImageDescription');
    totalFramesRecorded = parseFrameDescription(lastFrameDescription);
    
    % Get twop metadata from the open tiff file 
    twopMetadata = [];
    metadataStr = lastTifFile.getTag('Software'); 
    % Parse metadata and store in struct
    twopMetadata = parseMetadata(metadataStr);

    lastTifFile.close();
     % ---- Check if last grab was complete -----

%     expectedFramesPerPlanePerChannel = totalFramesRecorded / (twopMetadata.numSlices * twopMetadata.channelSave);
    % 1 'grab' = twopMetadata.numSlices * twopMetadata.channelSave
    % The frame numbers for both channels are saved the same. eg. last grab
    % Ch1(last grab) 35002 
    % length(str2num(twopMetadata.channelSave))) 
    framesToRemove = mod(totalFramesRecorded*length(str2num(twopMetadata.channelSave)), twopMetadata.numSlices *length(str2num(twopMetadata.channelSave))); 
    
    if framesToRemove ~= 0 
        fprintf('Last acquisition grab was incomplete! %d missing frames detected. Last grab will be removed.\n ', ...
            twopMetadata.numSlices * length(twopMetadata.channelSave) - framesToRemove);
        
        % Trim the incomplete grab
        trimmedTifPath = trimLastGrab(lastTiffilePath, framesToRemove, lastFrameNum);
        fprintf('Updating sessionFileInfo with the new tif file path, \n');
        %tempsessionFileInfo.stimFiles(iStim).tif_filepaths{end} = trimmedTifPath;
        %save(tempsessionFileInfo.sessionFileInfo_filepath, "tempsessionFileInfo", '-append')
    else
        fprintf('Frame acquisition complete for stimulus %d\n. No frames were trimmed and you can proceed. \n', iStim);
    end

end



%% Nested function to parse frame description 
function totalFramesRecorded = parseFrameDescription(lastFrameDescription)
    match = regexp(lastFrameDescription, 'frameNumbers\s*=\s*(\d+)', 'tokens');
    if ~isempty(match)
        totalFramesRecorded = str2double(match{1}{1}); % Convert to integer 
    else
        totalFramesRecorded = NaN; 
        warning('frameNumbers not found in metadata. \n');
    end
end

%% Nested function to trim the last incomplete grab and remove frames in place (modifies the TIFF in place)
function trimmedTifPath = trimLastGrab(lastTiffilePath, framesToRemove, lastFrameNum)
    % Define archive folder
    archiveFolder = fullfile(fileparts(lastTiffilePath), 'rawTiffArchived');

    % Create archive folder if it doesn't exist
    if ~exist(archiveFolder, 'dir')
        mkdir(archiveFolder);
    end

    % Extract file name and extension
    [folderPath, fileName, ext] = fileparts(lastTiffilePath);

    % Define paths for archived files
    archivedTiffPath = fullfile(archiveFolder, [fileName,ext]);  % Original moved file
    archivedCopyPath = fullfile(archiveFolder, [fileName,'_copy',ext]);  % Duplicate copy
    trimmedTifPath = fullfile(folderPath, [fileName,'_trimmed',ext]);  % Trimmed file
   
    if exist(trimmedTifPath, 'file')
        fprintf('Trimmed TIFF already exists: %s. Skipping...\n', trimmedTifPath);
        return;  % Skip this iteration and move to the next
    end

    % Move the original TIFF file to the archive folder
    if ~exist(archivedTiffPath, 'file') 
        movefile(lastTiffilePath, archivedTiffPath, 'f');  % Move instead of copying
        fprintf('Original TIFF moved to archive: %s\n', archivedTiffPath);
    else
        fprintf('TIFF file already archived. Proceeding with modification.\n');
    end

    % Create a duplicate of the archived TIFF file
    if ~exist(archivedCopyPath, 'file')
        copyfile(archivedTiffPath, archivedCopyPath);
        fprintf('Backup copy created at: %s\n', archivedCopyPath);
    else
        fprintf('Backup copy already exists.\n');
    end

    % Calculate frames to keep
    framesToKeep = lastFrameNum - framesToRemove;
    fprintf('Total frames: %d, removing last %d frame(s)...\n', lastFrameNum, framesToRemove);

    % Start processing timer
    tic;

    % Open the archived TIFF file (this is now the working copy); 
    TifToRead = Tiff(archivedTiffPath, 'r');

    % Create the trimmed TIFF file in the original location
    newTif = Tiff(trimmedTifPath, 'w');

    % Copy frames that we want to keep from the original TIFF to the new trimmed file
    for i = 1:framesToKeep 
        TifToRead.setDirectory(i);
        imgData = TifToRead.read();

        % Copy tags
        newTif.setTag('Photometric', TifToRead.getTag('Photometric'));
        newTif.setTag('ImageWidth', TifToRead.getTag('ImageWidth'));
        newTif.setTag('ImageLength', TifToRead.getTag('ImageLength'));
        newTif.setTag('Compression', Tiff.Compression.None);
        newTif.setTag('BitsPerSample', TifToRead.getTag('BitsPerSample'));
        newTif.setTag('SampleFormat', TifToRead.getTag('SampleFormat'));
        newTif.setTag('PlanarConfiguration', TifToRead.getTag('PlanarConfiguration'));
        newTif.setTag('RowsPerStrip', TifToRead.getTag('RowsPerStrip'));
        newTif.setTag('SamplesPerPixel', TifToRead.getTag('SamplesPerPixel'));
        newTif.setTag('Orientation', TifToRead.getTag('Orientation'));

        % Preserve metadata
        try
            newTif.setTag('ImageDescription', TifToRead.getTag('ImageDescription'));
        catch
            warning('ImageDescription tag not found.');
        end

        % Write frame
        newTif.write(imgData);

        % Prepare for next frame
        if i < framesToKeep
            newTif.writeDirectory();
        end
    end

    % Close TIFF files
    TifToRead.close();
    newTif.close();

    % Release file handles
    clear TifToRead newTif;
    fclose('all');
    pause(1);  % Allow MATLAB to release file locks

    fprintf('Trimmed TIFF saved at: %s\n', trimmedTifPath);

    % Stop timer
    toc;

    % Return the number of frames after trimming
    trimmedFrames = framesToKeep;
end
end

%% Helper function to parse metadata - entirely ChatGPT :P
function metadata = parseMetadata(metadataStr)
% Parses ScanImage metadata from a string and extracts key-value pairs.

    % Initialize struct
    metadata = struct();
    
    % Define regex pattern to match key-value pairs - what is regex? 
    pattern = '(\S+)\s*=\s*(.+?)(?=(\s*\S+\s*=)|$)';

    % Extract matches
    tokens = regexp(metadataStr, pattern, 'tokens');

    % Iterate over matches and store them in struct
    for i = 1:length(tokens)
        fullKey = tokens{i}{1};  % Extract full key (e.g., SI.hStackManager.numSlices)
        valueStr = strtrim(tokens{i}{2}); % Extract value string and trim whitespace

        % Convert value only if it is purely numeric
        value = convertValue(valueStr);

        % Extract only the last part of the key (e.g., 'numSlices')
        keyParts = strsplit(fullKey, '.');
        fieldName = keyParts{end}; 

        % Convert field name to a valid MATLAB struct field name
        fieldName = matlab.lang.makeValidName(fieldName);

        % Assign value to metadata struct
        metadata.(fieldName) = value;
    end
end

%% Helper function to convert values 
function value = convertValue(valueStr)
% Converts a metadata value string into a numeric or keeps it as a string.

    % Try converting to a number
    numericValue = str2num(valueStr); 
    
    if isempty(numericValue) || contains(valueStr, ' ') || contains(valueStr, '{') || contains(valueStr, '[')
        % Keep as string if conversion fails or if it's an array-like structure
        value = valueStr;
    else
        % Store as numeric if conversion succeeds
        value = numericValue;
    end
end