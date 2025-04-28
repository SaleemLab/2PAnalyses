function [sessionFileInfo] = get2PMetadata(sessionFileInfo)
% Extracts and parses metadata from the first TIFF file of each stimulus. 

for iStim = 1:length(sessionFileInfo.stimFiles)
    % Save 2p metadata filepath to sessionFileInfo for future reference
    stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name ...
        '_2pMetaData' '_' sessionFileInfo.stimFiles(iStim).name '.mat'];
    sessionFileInfo.stimFiles(iStim).TwoPMetaData = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
    
    % Load metadata from the first TIFF file; All tiffs will have the same
    % 'software' metadata unless trimmed! 
    if isempty(sessionFileInfo.stimFiles(iStim).tif_filepaths)
        continue       
    end
    
    % add an if statement to deal with trimmed tiff files with only one tif
    % saved:
    if length(sessionFileInfo.stimFiles(iStim).tif_filepaths) == 1 && contains(sessionFileInfo.stimFiles(iStim).tif_filepaths(1), '_trimmed')
        % load the original tif (which contains software) in the
        trimmedPath = sessionFileInfo.stimFiles(iStim).tif_filepaths{1};
        % Remove '_trimmed' from the filename
        [folderPath, fileName, ext] = fileparts(trimmedPath);
        fileNameRaw = strrep(fileName, '_trimmed', '');
        % Go one level up from the folderPath
        parentDir = fileparts(trimmedPath);
        % Create new folder path by appending 'rawTiffArchived'
        newFolderPath = fullfile(parentDir, 'rawTiffArchived');
        % Create the new full path
        TiffPathToRead = fullfile(newFolderPath, [fileNameRaw ext]);
        % archivedTiff 
        disp('Reading Metadata from Original rather than trimmed tif file ')
    else
        TiffPathToRead = sessionFileInfo.stimFiles(iStim).tif_filepaths{1};
    end 

    tiffFile = Tiff(TiffPathToRead, 'r');
    metadataStr = tiffFile.getTag('Software'); 
    % Parse metadata and store in struct
    twopMetadata = parseMetadata(metadataStr);

    if iStim==1
        sessionFileInfo.origNPlanes = twopMetadata.numSlices; 
    elseif twopMetadata.numSlices ~= sessionFileInfo.origNPlanes
        WARNING('The number of slices do not match between the different stimuli!!!')
    end
    
    save(sessionFileInfo.stimFiles(iStim).TwoPMetaData, 'twopMetadata');
end
% Append the sessionFileIngo with the filepaths to the twopMetadata
% files
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

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