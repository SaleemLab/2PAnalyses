function analyseSparseNoise(sessionFileInfo, gridSize, plotflag, framesToShow)
    %
    % TODO:
    % (1) Fix the number of rois to run svd on
    % (2) Add the option to save the plots in a pdf 
    %
    % Inputs:
    %   - sessionFileInfo (struct): Structure containing stimFiles and associated data file paths.
    %   - gridSize (array): Sparse Noise Grid Sizes (e.g., [x y]); Defult is [6 4] 
    %   - plotflag (int): To plot the SVD outputs; Default is to plot (1)
    %   - framesToShow (int:int): Number of frames to show in the plot; Default is 1:3. 
    %

    % Set default values if not provided
    if nargin < 2
        gridSize = [6 4]; % Default grid size [grid_x grid_y] 
    end 

    if nargin < 3
        plotflag = 1; % Defualt if to plot 
    end 
    
    if nargin < 4
        framesToShow = 1:3; 
    end  

    %% Load session file information and relevant data files
    isStim = false(1, length(sessionFileInfo.stimFiles));
    for iStim = 1:length(sessionFileInfo.stimFiles)
        isStim(iStim) = strcmp('SparseNoise', sessionFileInfo.stimFiles(iStim).name);
    end
    iStim = find(isStim, 1);

    if isempty(iStim)
        error('No SparseNoise stimulus file found in sessionFileInfo.');
    end

    % Load data files
    if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
       exist(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'file') && ...
       exist(sessionFileInfo.stimFiles(iStim).Response, 'file') && ...
       exist(sessionFileInfo.stimFiles(iStim).TwoPMetaData, 'file')

        % Load stimulus data
        load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData'); 
        load(sessionFileInfo.stimFiles(iStim).Response, 'response'); 
        load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData'); 
        load(sessionFileInfo.stimFiles(iStim).TwoPMetaData, 'twopMetadata');
    else 
        error('Required data files (BonsaiData.mat, Response.mat, TwoPData.mat) are missing.');
    end

    %% Extract stimulus timing and matrix
    onARDStimTimes = sort([bonsaiData.onARDTimes; bonsaiData.offARDTimes]);
    stimMatrix = bonsaiData.stimMatrix;
    stimMatrix = stimMatrix(1:length(onARDStimTimes));

    %% Loop through each imaging plane
    numPlanes = length(twoPData);
    for thisPlane = 1:numPlanes
        disp("Processing plane " + thisPlane);
        
        % Extract fluorescence data for this plane
        F = twoPData(thisPlane).F;
        validROIs = twoPData(thisPlane).iscell(:,1) == 1; 
        rois = F(validROIs, :); 
        
        % Number of neurons and trials
        numRois = size(rois, 1);
        numTrials = length(response(thisPlane).responseFrameIdx);

        % Process responseFramesIdx for this plane
        framesToAnalyse = cellfun(@find, response(thisPlane).responseFrameIdx, 'UniformOutput', false);
        maxFrames = max(cellfun(@numel, framesToAnalyse));
        twopIndices = nan(numTrials, maxFrames);

        for trial = 1:numTrials
            twopIndices(trial, 1:numel(framesToAnalyse{trial})) = framesToAnalyse{trial};
        end

        %% Extract responses for each neuron
        roiStimResponses = zeros(numRois, numTrials, maxFrames);
        validMask = ~isnan(twopIndices);

        for neuron = 1:numRois
            tempF = rois(neuron, :);
            tempIndices = twopIndices(validMask);
            roiStimResponses(neuron, validMask) = tempF(tempIndices);
        end

        %% Reshape stimulus matrix
        stimulusMatrixCells = cellfun(@(x) x(:)', stimMatrix, 'UniformOutput', false);
        stimulusMatrix = cat(1, stimulusMatrixCells{:});
        stimMatrix = reshape(permute(stimulusMatrix, [2 1]), [gridSize(2), gridSize(1), size(stimulusMatrix, 1)]);

        %% Define options
        sn_options.grid_size = [size(stimMatrix, 1), size(stimMatrix, 2)];
        sn_options.mapSampleRate = twopMetadata.scanFrameRate / twopMetadata.numSlices; 
        sn_options.mapsToShow = {'linear', 'black', 'white', 'contrast'};
        sn_options.mapMethod = 'fitlm';
        sn_options.framesToShow = framesToShow;  
        sn_options.plotflag = plotflag;

        %% Analyse and plot
        initMap = cell(numRois, 1); 

        for iN = 1:numRois
            roiRespTmp = squeeze(roiStimResponses(iN, :, :)); 
            initMap{iN} = sparseNoiseAnalysis(stimMatrix, roiRespTmp, [], [], sn_options);
        end
    end
end
