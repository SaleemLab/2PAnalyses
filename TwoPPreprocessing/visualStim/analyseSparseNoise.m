function [initMap, sessionFileInfo] = analyseSparseNoise(sessionFileInfo,signalToUse, plotflag, framesToShow, saveData)
% analyseSparseNoise: Extracts RFs, saves to sessionROIData, and exports PDF.
%
% Inputs:
%   - sessionFileInfo: struct containing stimFiles and data file paths.
%   - plotflag: 1 to plot/save PDF, 0 to just compute (Default: 1)
%   - framesToShow: indices for the SVD/time-bin plots (Default: [1 2 4 6 8 10 12 14])
%   - saveData: 1 to save 'sparseNoiseRF' struct to sessionROIData (Default: 1)

%% Default parameters
if nargin < 2, signalToUse = 'dFFNeuropilCorrected'; end 
if nargin < 3, plotflag = 1; end
if nargin < 4, framesToShow = [1 2 4 6 8 10 12 14]; end
%if nargin < 4, framesToShow = [1 5 10 15 20 25 30]; end
if nargin < 5, saveData = 1; end 



%% Identify Sparse Noise Stimulus
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isSparseNoise(iStim) = contains(sessionFileInfo.stimFiles(iStim).name, 'SparseNoiseTexture');
end
iStim = find(bonsaiData.isSparseNoise==1);
if isempty(iStim)
    error('No SparseNoise stimulus file found in sessionFileInfo.');
end

%% Load Data Files
if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
   exist(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'file') && ...
   exist(sessionFileInfo.stimFiles(iStim).Response, 'file')
    
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).Response, 'response');
    processedtwoPData = load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData);
else
    error('Required Sparse Noise data files are missing.');
end

%% Setup Directory and PDF
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
pdfFileName = fullfile(figSaveDir, [sessionFileInfo.stimFiles(iStim).name '_RFs.pdf']);

% Delete existing PDF to avoid appending to old results
if plotflag && exist(pdfFileName, 'file')
    delete(pdfFileName);
end

%% Extract Fluorescence and Timing
if contains(signalToUse, 'spks')
    signal = processedtwoPData.spks; 
else   
    signal = processedtwoPData.processedSignals.(signalToUse);

end 
numRois = size(signal, 1);
numTrials = length(response.responseFrameIdx);

% Map 2P frame indices for each stimulus trial
framesToAnalyse = cellfun(@find, response.responseFrameIdx, 'UniformOutput', false);
maxFrames = max(cellfun(@numel, framesToAnalyse));
twopIndices = nan(numTrials, maxFrames);
for trial = 1:numTrials
    twopIndices(trial, 1:numel(framesToAnalyse{trial})) = framesToAnalyse{trial};
end

% Extract responses per neuron
roiStimResponses = zeros(numRois, numTrials, maxFrames);
validMask = ~isnan(twopIndices);
for neuron = 1:numRois
    tempF = signal(neuron, :);
    roiStimResponses(neuron, validMask) = tempF(twopIndices(validMask));
end

%% Format Stimulus Matrix
stimulusMatrixCells = cellfun(@(x) x(:)', bonsaiData.stimMatrix(1:numTrials), 'UniformOutput', false);
stimulusMatrix = cat(1, stimulusMatrixCells{:});
stimMatrix = reshape(permute(stimulusMatrix, [2 1]), bonsaiData.gridSize(1), bonsaiData.gridSize(2), size(stimulusMatrix, 1));

%% Analysis Options
sn_options.grid_size = [bonsaiData.gridSize(1), bonsaiData.gridSize(2)];
sn_options.mapSampleRate = 60; 
sn_options.mapsToShow = {'linear', 'black', 'white', 'contrast'};
sn_options.mapMethod = 'fitlm';
sn_options.framesToShow = framesToShow;
sn_options.plotflag = plotflag;

%% Main Analysis Loop
initMap = cell(numRois, 1);
initPMap = cell(numRois, 1);

for iN = 1:numRois
    roiRespTmp = squeeze(roiStimResponses(iN, :, :));
    
    % Capture current figures to distinguish from new ones
    existingFigs = findobj('Type', 'figure');
    
    trialResponse = mean(roiRespTmp, 2, 'omitnan');   
    [chi2stat(iN), chi2pVal(iN)] = computeRFChiSquareSignificance(stimMatrix, trialResponse, 1000);

    % Run Sparse Noise Analysis
    [initMap{iN}, initPMap{iN}, initR2{iN} ~] = sparseNoiseAnalysis(stimMatrix, roiRespTmp, [], [], sn_options);
    
    if plotflag
        % Identify new figures created for this ROI
        allFigs = findobj('Type', 'figure');
        newFigs = setdiff(allFigs, existingFigs);
        
        if ~isempty(newFigs)
            for iF = 1:length(newFigs)
                fig = newFigs(iF);
                ax = findall(fig, 'type', 'axes');
                
                if ~isempty(ax)
                    title(ax(1), sprintf('ROI %d - Map %d', iN, iF), 'FontSize', 12);
                    % Export and append to PDF
                    exportgraphics(fig, pdfFileName, 'Append', true, 'ContentType', 'vector');
                end
            end
            % Explicitly close new figures to free memory and prevent empty pages
            close(newFigs);
        end
    end
    
    if mod(iN, 20) == 0 || iN == numRois
        fprintf('Processed ROI %d of %d\n', iN, numRois);
    end
end

%% Saving Data to sessionROIData 
if saveData && isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file')
    sparseNoiseRF.initPMap = initPMap; 
    sparseNoiseRF.initMap = initMap;
    sparseNoiseRF.gridSize = bonsaiData.gridSize;
    sparseNoiseRF.options = sn_options;
    sparseNoiseRF.stimName = sessionFileInfo.stimFiles(iStim).name;
    sparseNoiseRF.chi2stat = chi2stat;
    sparseNoiseRF.chi2pVal = chi2pVal;
    
    disp(['Appending RF data to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, 'sparseNoiseRF', '-append');
    
elseif saveData
    warning('sessionROIData path not found. Data not appended.');
end

disp('Sparse Noise Analysis Complete.');
end