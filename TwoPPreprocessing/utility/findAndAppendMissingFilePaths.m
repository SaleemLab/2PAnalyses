function [sessionFileInfo, allFilesFound] = findAndAppendMissingFilePaths(sessionFileInfo)
%FINDANDAPPENDMISSINGFILEPATHS Checks for missing files using a correct field-to-filename map.
%
%   This version uses a containers.Map to correctly generate the expected
%   filenames, accounting for discrepancies between struct field names and
%   the strings used in the actual .mat filenames.
%
%
%   OUTPUTS:
%       sessionFileInfo - The updated struct with missing paths filled in.
%       allFilesFound   - A logical scalar, true if all expected files were found.

    % --- Initialisation ---
    allFilesFound = true;
    rootDir = 'Z:\ibn-vision\DATA\SUBJECTS';   
    mouseID = sessionFileInfo.animal_name;
    session = sessionFileInfo.session_name;
    searchDirectory = fullfile(rootDir, mouseID, 'Analysis', session);


    fprintf('Checking session for Mouse: %s, Date: %s\n', mouseID, session);
    if ~isfolder(searchDirectory)
        error('Analysis directory does not exist: %s', searchDirectory);
    end

    % *** KEY CHANGE: Define the mapping from struct field name to filename component ***
    fieldToFilenameMap = containers.Map(...
        {'TwoPMetaData', 'processedPeripheralData', 'BonsaiData', 'mergedBonsai2PSuite2pData', 'processedMergedBonsaiSuite2pData', 'Response'}, ...
        {'2pMetaData',   'PeripheralData',          'BonsaiData', '2pData',                     'processed2PData',                  'Response'} ...
    );
    
    fieldsToCheck = keys(fieldToFilenameMap);

    % --- Loop through each stimulus presentation ---
    for i = 1:length(sessionFileInfo.stimFiles)
        stimName = sessionFileInfo.stimFiles(i).name;
        fprintf('Verifying files for stimulus: %s\n', stimName);

        % --- Loop through the specific fields we need to check ---
        for j = 1:length(fieldsToCheck)
            fieldName = fieldsToCheck{j};
            
            if ~isfield(sessionFileInfo.stimFiles(i), fieldName) || isempty(sessionFileInfo.stimFiles(i).(fieldName))
                
                % Use the map to get the correct string for the filename
                fileNameComponent = fieldToFilenameMap(fieldName);
                
                % Construct the expected filename using the correct component
                expectedFileName = sprintf('%s_%s_%s_%s.mat', mouseID, session, fileNameComponent, stimName);
                expectedFilePath = fullfile(searchDirectory, expectedFileName);

                if isfile(expectedFilePath)
                    sessionFileInfo.stimFiles(i).(fieldName) = expectedFilePath;
                    fprintf('Found and appended path for field: %s\n', fieldName);
                else
                    warning('Could not find file: %s', expectedFileName);
                    allFilesFound = false;
                end
            end
        end
    end

    % --- Final Summary ---
    if allFilesFound
        fprintf('\nFile verification complete. All paths are up to date.\n');
    else
        fprintf('\nFile verification failed. One or more files were not found.\n');
    end
end