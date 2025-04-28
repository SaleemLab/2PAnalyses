function [processedTwoPData] = resampleMergedBonsaiSuite2PDataV0(sessionFileInfo, samplingRate, mainTimeToUse, plotFlag)

    % Resample Bonsai + Suite2p Merged Data and Standardise Time Across Planes
    %
    % This function concatenates and resamples all relevant time and calcium 
    % imaging data fields from the merged Bonsai-Suite2P data, aligning them 
    % to a unified sampling rate and time base.
    %
    % Inputs:
    %   - sessionFileInfo (struct): 
    %       Contains paths and metadata related to the session, including
    %       stimulus files, save directories, and session identifiers.
    %
    %   - samplingRate (float, optional): 
    %       The desired output sampling rate in Hz. 
    %       Default is 60 Hz if not provided.
    %
    %   - mainTimeToUse (char, optional): 
    %       The time vector to use as reference for resampling ('TwoPFrameTime', 
    %       'ArduinoTime', etc.). Default is 'TwoPFrameTime'.
    %
    %   - plotFlag (logical or char, optional): 
    %       Whether to plot a summary after resampling. Can be
    %       logical true/false or the string 'true'/'false'. Default is 'true'.
    %
    % Outputs:
    %   - processedTwoPData (struct):
    %       Contains resampled and merged imaging and time fields, including:
    %         - F, Fneu, spks           - Interpolated calcium traces
    %         - TwoPFrameTime, etc.     - Resampled time fields
    %         - iscell, redcell, stat   - ROI metadata from Suite2p
    %         - ops, planeName          - Per-plane metadata
    %         - roiPlaneIdentity        - ROI-to-plane mapping
    %         - resamplingInfo          - Struct with:
    %             - samplingRate
    %             - mainTimeToUse
    %             - mainTime
    %             - uniqueTime
    %             - uniqueIdx
    %             - resampledTime 
    %
    %
    % Usage Example:
    %   [processedTwoPData] = resampleMergedBonsaiSuite2PData(sessionFileInfo, 60, 'TwoPFrameTime', true)
    %
    % Dependencies:
    %  - processPeripheralFiles.m
    %  - getVRBonsaiFiles.m
    %  - mergeBonsaiSuite2pFiles.m
    % 
    % Sonali and Aman - March 2025
     

    if nargin < 2
        samplingRate = 60; 
    end 
    
    if nargin < 3
        mainTimeToUse = 'TwoPFrameTime';
    end 
    
    if nargin < 4
        plotFlag = 'true';
    end 
    
    for iStim = 1:length(sessionFileInfo.stimFiles)
        bonsaiData.isVRstim(iStim) = strcmp('VRCorr',sessionFileInfo.stimFiles(iStim).name);
    end
    
    iStim = find(bonsaiData.isVRstim==1);
    
    % Load all the relevent data files to save a copy of the resampledTime. 
    if exist(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'file') && ...
            exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
            exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file')
        load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData')
        load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
        load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
        
    else
        error('MergedBonsaiSuite2P Data missing for stimulus "%s".', stimName);
    end
    
    % Create a new .mat file to save the interpolated twoPdata 
    stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_processedBonsai2PData' '_' sessionFileInfo.stimFiles(iStim).name '.mat'];
    sessionFileInfo.stimFiles(iStim).processedBonsai2PSuite2pData = fullfile(sessionFileInfo.Directories.save_folder,stimFileName);
    
    %% Concatenate time varaibles across planes and resample time vectors 
    bonsaiFieldsToResample = {'TwoPFrameTime', 'BonsaiTime', 'ArduinoTime', 'RenderFrameCount', 'LastSyncPulseTime'};
    
    % All fields share the same length
    mainTime = vertcat(twoPData.(mainTimeToUse));  % or use ArduinoTime as the reference
    % Finds the first occurrence of each unique value to remove repated measures 
    [uniqueTime, uniqueIdx] = unique(mainTime);
    % Interpolate each field onto a common time base
    resampledTime = uniqueTime(1):1/samplingRate:uniqueTime(end);
    
    for thisField = 1:numel(bonsaiFieldsToResample)
        fieldNameResample = bonsaiFieldsToResample{thisField};
        concatenatedTime = vertcat(twoPData.(fieldNameResample));
        concatenatedTime = concatenatedTime(uniqueIdx); % Apply the same filtering to remove duplicates
        processedTwoPData.(fieldNameResample) = interp1(uniqueTime,concatenatedTime,resampledTime,'linear')'; % Linear for all time variables 
    end
    
    
    %% Resample suite2p data and append other plane-related information
    
    suite2PFieldsToResample = {'F', 'Fneu', 'spks'};
    interpType = {'linear', 'linear', 'nearest'};  % Interpolation method per signal type
    
    % Initialize fields
    for idx = 1:numel(suite2PFieldsToResample)
        processedTwoPData.(suite2PFieldsToResample{idx}) = [];
    end
    processedTwoPData.roiPlaneIdentity = [];
    
    % Append fields
    processedTwoPData.iscell = [];
    processedTwoPData.redcell = [];
    processedTwoPData.stat = {};
    processedTwoPData.ops = {};
    processedTwoPData.planeName = {};
    
    for thisPlane = 1:numel(twoPData)
        originalTime = double(twoPData(thisPlane).(mainTimeToUse));
    
        % Resample ROI signals
        for idx = 1:numel(suite2PFieldsToResample)
            field = suite2PFieldsToResample{idx};
            signal = double(twoPData(thisPlane).(field));
            interpolatedSignal = interp1(originalTime, signal', resampledTime, interpType{idx})';
            processedTwoPData.(field) = [processedTwoPData.(field); interpolatedSignal];
        end
    
        % Append ROI-wise
        processedTwoPData.iscell = [processedTwoPData.iscell; twoPData(thisPlane).iscell];
        processedTwoPData.redcell = [processedTwoPData.redcell; twoPData(thisPlane).redcell];
        
        % ROI identity
        nROIs = size(twoPData(thisPlane).F, 1);
        processedTwoPData.roiPlaneIdentity = [processedTwoPData.roiPlaneIdentity; repmat(thisPlane-1, nROIs, 1)];
    
        % Append cell-wise
        processedTwoPData.stat = [processedTwoPData.stat, twoPData(thisPlane).stat];
    
        % Append per-plane
        processedTwoPData.ops{end+1} = twoPData(thisPlane).ops;
        processedTwoPData.planeName{end+1} = twoPData(thisPlane).planeName;
       
    end
    
    %% Save the resamping information as a copy in all data steams
    processedTwoPData.resamplingInfo.samplingRate = samplingRate;
    processedTwoPData.resamplingInfo.resampledTime = resampledTime;
    processedTwoPData.resamplingInfo.mainTimeToUse = mainTimeToUse;
    processedTwoPData.resamplingInfo.mainTimeRaw = mainTime;
    processedTwoPData.resamplingInfo.uniqueTimeRaw = uniqueTime;
    processedTwoPData.resamplingInfo.uniqueIdx = uniqueIdx;
    
    % Also add to bonsaiData and peripheralData
    bonsaiData.resamplingInfo = processedTwoPData.resamplingInfo;
    peripheralData.resamplingInfo = processedTwoPData.resamplingInfo;

    
    %% Optional plotting (sanity check)
    if islogical(plotFlag) && plotFlag || ischar(plotFlag) && strcmpi(plotFlag, 'true')
    
        figure;
        hold on;
        histogram(diff(mainTime), 'BinWidth', 0.001, 'DisplayName', 'Original');
        histogram(diff(uniqueTime), 'BinWidth', 0.001, 'DisplayName', 'Unique');
        histogram(diff(processedTwoPData.(mainTimeToUse)), 'BinWidth', 0.001, 'DisplayName', 'Resampled');
        xlabel('Time Difference (s)');
        ylabel('Count');
        title('Distribution of Time Differences');
        legend('Location', 'best');
        xlim([0 0.2]);
    
    
        % Neuron trace comparison
        planeIndex = 1;
        fOrig = double(twoPData(planeIndex).F);
        originalTime = double(twoPData(planeIndex).(mainTimeToUse));
        nROIs = size(fOrig, 1);
        roiIndices = randperm(nROIs, 5);
        roiMask = processedTwoPData.roiPlaneIdentity == (planeIndex - 1);
        fResampled = processedTwoPData.F(roiMask, :);
        
        figure('Name', 'Calcium Trace Comparison');
        for idx = 1:numel(roiIndices)
            roi = roiIndices(idx);
            subplot(numel(roiIndices), 1, idx);
            hold on;
            plot(originalTime, fOrig(roi, :), 'k.', 'DisplayName', 'Original');
            plot(processedTwoPData.(mainTimeToUse), fResampled(roi, :), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Resampled');
            title(sprintf('ROI %d (Plane %d)', roi, planeIndex));
            ylabel('F');
            if idx == 1
                legend();
            end
        end
        xlabel('Time (s)');
        sgtitle('Neuron Trace Resampling: Original vs Interpolated');

    end
    

fprintf('Saving processed data files...\n');
tic; 
save(sessionFileInfo.stimFiles(iStim).processedBonsai2PSuite2pData, 'processedTwoPData', '-v7.3');
fprintf('Saved processedTwoPData\n');
save(sessionFileInfo.stimFiles(iStim).BonsaiData, "bonsaiData");
fprintf('Saved bonsaiData\n');
save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData");
fprintf('Saved peripheralData\n');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
fprintf('Saved sessionFileInfo\n');
elapsedTime = toc; 
fprintf('✅ All files saved in %.2f seconds.\n', elapsedTime);

end

