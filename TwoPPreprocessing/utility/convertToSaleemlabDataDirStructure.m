function convertToSaleemlabDataDirStructure(mouseID, sessionID)

    uclBaseDir = 'C:\UCLOpenDATA\'; 
    saleemBaseDir = 'C:\DATA\Subjects';
    
    % Construct source and destination path
    sourcePath = fullfile(uclBaseDir, ['sub-' mouseID], sessionID); 
    destPath = fullfile(saleemBaseDir, mouseID, 'Bonsai' ,sessionID);
    
    if ~exist(sourcePath, 'dir')
        error('Source directory not found: %s', sourcePath);
    end
    
    if ~exist(destPath, 'dir')
        mkdir(destPath);
        fprintf('Created destination folder: %s\n', destPath);
    end
    
    % Identify all files recursively
    fileList = dir(fullfile(sourcePath, '**', '*.*'));
    fileList = fileList(~[fileList.isdir]);
    
    fprintf('Found %d files. Starting copy...\n', length(fileList));
    copyCount = 0;
    skipCount = 0;
    
    for thisFile = 1:length(fileList)
        sourceFile = fullfile(fileList(thisFile).folder, fileList(thisFile).name);
        
        % parent folder path
        [parentPath, fileName, fileExt] = fileparts(sourceFile);
        
        % grandparent folder path -- is this what it is called?
        [grandParentPath, ~] = fileparts(parentPath);
        
        % grandparent folder name :P
        [~, grandParentFolderName] = fileparts(grandParentPath);
        
        % create new name 
        newName = [grandParentFolderName, '_', fileName, fileExt];
        
        destFile = fullfile(destPath, newName);
  
        % Safety check for collisions
        if exist(destFile, 'file')
            warning('Collision: %s already exists. Skipping.', newName);
            skipCount = skipCount + 1;
            continue;
        end
        
        try
            copyfile(sourceFile, destFile);
            copyCount = copyCount + 1;
        catch ME
            fprintf('Error copying %s: %s\n', newName, ME.message);
        end
    end
    
    fprintf('\n--- Process Complete ---\n');
    fprintf('Files Copied: %d\n', copyCount);
    fprintf('Files Skipped: %d\n', skipCount);
    fprintf('Destination: %s\n', destPath);
end