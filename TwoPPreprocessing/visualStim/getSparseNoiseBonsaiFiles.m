function [bonsaiData] = getSparseNoiseBonsaiFiles(sessionFileInfo, gridSize)
    % Import and save stimulus matrix (position on screen) from Sparse Noise Bin file 
    % 
    % Inputs:
    %   - sessionFileInfo (.mat) : File info .mat file
    %   - gridSize (list) : Grid size x and grid size y (e.g., [12 8])  
    %  
    % Output:
    %   - bonsaiData (struct) : Collated stim matrix showing position of grids
    %   
    % Sam, Aman and Sonali 

    if nargin < 2
        gridSize = [6 4];
    end 
    for iStim = 1:length(sessionFileInfo.stimFiles)
        bonsaiData.isSparseNoise(iStim) = strcmp('SparseNoise',sessionFileInfo.stimFiles(iStim).name);
    end
    
    iStim = find(bonsaiData.isSparseNoise==1);
    
    % Construct file path for saving Bonsai data
    stimFileName = sprintf('%s_%s_BonsaiData_%s.mat', ...
            sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
    % Save filepath to sessionFileInfo     
    sessionFileInfo.stimFiles(iStim).BonsaiData = ...
        fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
    
    % Locate and read stimulus events table
    binFilePath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'Log');
    
    fileID=fopen(binFilePath);
    thisBinFile=fread(fileID);
    fclose(fileID);
    
    % Translate stimulus into -1:1 scale
    stimMatrix = zeros(1,length(thisBinFile));
    stimMatrix(thisBinFile==0)=-1;
    stimMatrix(thisBinFile==255)=1;
    stimMatrix(thisBinFile==128)=0;
    
    % Make a NxM grid from the stimulus log
    stimMatrix = reshape(stimMatrix, [gridSize(1), gridSize(2), length(thisBinFile)/gridSize(1)/gridSize(2)]);
    stimMatrix = stimMatrix(:,:,1:end-1); % The last 'stimulus'
    
    for thisTrial = 1:size(stimMatrix,3)
        bonsaiData.stimMatrix{thisTrial,1} = squeeze(stimMatrix(:,:,thisTrial));
    end
    
    % Save the extracted Bonsai data
    save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo')
end