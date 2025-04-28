function [bonsaiData, sessionFileInfo] = getStimTimes(sessionFileInfo, stimName, thresholdOn, thresholdOff, useQuadState, plotFlag)
    % Extracts stimulus ON and OFF times using PD (and QuadState as a reference, if required)
    % and removes stimulus events after two-photon imaging has stopped.
    % Loads previously generated bonsaiData (from getTuningStimEventsFile.m for eg) 
    % if available or creates a new one for the stimulus if not previously created.  
    %
    % Inputs:
    %   sessionFileInfo (struct): Structure with file and session information. Used to load the peripheralData and twoPData             
    %   stimName (str): Name of the stimulus to extract times for
    %   thresholdOn (int/float): Threshold for detecting PD ON events
    %   thresholdOff (int/float): (Optional) Threshold for detecting PD OFF events (default: same as thresholdOn)
    %   useQuadState (bool): (Optional) Boolean flag; if true, use QuadState to refine PD events (default: false)
    %   plotFlag (bool): (Optional) Boolean flag; if true, plot results (default: true)
    %
    % Outputs:
    %   bonsaiData         - Structure with fields:
    %                           .onARDTimes  - Arduino times for detected ON events
    %                           .offARDTimes - Arduino times for detected OFF events
    %                             
    %
    % Aman and Sonali - Feb 2025

    % Set default values if not provided
    if nargin < 4
        thresholdOff = thresholdOn;
    end 
    if nargin < 5
        useQuadState = false;
    end 
    if nargin < 6
        plotFlag = true;
    end

    %% Load relevant paths and data
    
    % Find stimulus index: Locate the current stimulus in sessionFileInfo
    isStim = false(1, length(sessionFileInfo.stimFiles));
    for iStim = 1:length(sessionFileInfo.stimFiles)
        isStim(iStim) = strcmp(stimName, sessionFileInfo.stimFiles(iStim).name);
    end
    iStim = find(isStim, 1);  % take the first match
    if isempty(iStim)
        error('Stimulus name not found in sessionFileInfo.');
    end

    % BonsaiData: Check if the BonsaiData file exists; if yes, load it; if not, create a new file path.
    if isfield(sessionFileInfo.stimFiles(iStim), 'BonsaiData') && ...
            exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file')
        load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    else
        stimFileName = sprintf('%s_%s_BonsaiData_%s.mat', ...
            sessionFileInfo.animal_name, sessionFileInfo.session_name, sessionFileInfo.stimFiles(iStim).name);
        sessionFileInfo.stimFiles(iStim).BonsaiData = ...
            fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
    end

    % PeripheralData & twoPData
    if exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file')
        load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
        load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData');
    else
        error('Missing PeripheralData and/or twoPData.');
    end

    %% Process PD

    % Extract PD signals
    pdArduinoTime = peripheralData.Photodiode.rawArduinoTime;
    pdValue       = peripheralData.Photodiode.rawValue;
    % Smoothning 
    pdValue = movmedian(peripheralData.Photodiode.rawValue,20);

    % Compute ON and OFF transitions
    pdON  = pdValue >= thresholdOn;
    pdOFF = pdValue <= thresholdOff;

    pdONDiff  = [0; diff(pdON)];
    pdOFFDiff = [0; diff(pdOFF)];

    pdOnsetIndex  = pdONDiff > 0;
    pdOffsetIndex = pdOFFDiff > 0;

    % Store PD ON and OFF events
    onARDTimes  = pdArduinoTime(pdOnsetIndex);
    offARDTimes = pdArduinoTime(pdOffsetIndex);

    % Compute average ON interval for refining OFF events (if possible)
    avgOnInterval = 0;
    if numel(onARDTimes) > 1
        avgOnInterval = mean(diff(onARDTimes));
    end

    % Refine OFF events by keeping only the first OFF event within the expected interval
    refinedOffTimes = [];
    for i = 1:length(onARDTimes)
        validOffTimes = offARDTimes(offARDTimes > onARDTimes(i) & offARDTimes <= onARDTimes(i) + avgOnInterval);
        if ~isempty(validOffTimes)
            refinedOffTimes = [refinedOffTimes; validOffTimes(1)];
        end
    end
    offARDTimes = refinedOffTimes;

    %% If QuadState exists, use it to refine both ON and OFF detection
    if useQuadState && isfield(peripheralData, 'quadstate') && isfield(peripheralData.Quadstate, 'Value')
        quadTimes  = peripheralData.Quadstate.rawArduinoTime;
        quadValues = peripheralData.Quadstate.rawValue;

        % Determine QuadState ON (value 1) and OFF (value 0) times
        quadOnArduinoTime  = quadTimes(quadValues == 1);
        quadOffArduinoTime = quadTimes(quadValues == 0);

        stimOnPd  = zeros(size(quadOnArduinoTime));
        stimOffPd = zeros(size(quadOffArduinoTime));

        % Find the first PD ON after each QuadState ON event
        for i = 1:length(quadOnArduinoTime)
            tempInd = pdArduinoTime >= quadOnArduinoTime(i);
            tempTimeWindow = pdArduinoTime(tempInd);
            transientOnInd = find(pdValue(tempInd) >= thresholdOn, 1);
            if ~isempty(transientOnInd)
                stimOnPd(i) = tempTimeWindow(transientOnInd);
            else
                stimOnPd(i) = NaN;
            end
        end

        % Find the first PD OFF after each QuadState OFF event
        for i = 1:length(quadOffArduinoTime)
            tempInd = pdArduinoTime >= quadOffArduinoTime(i);
            tempTimeWindow = pdArduinoTime(tempInd);
            transientOffInd = find(pdValue(tempInd) <= thresholdOff, 1);
            if ~isempty(transientOffInd)
                stimOffPd(i) = tempTimeWindow(transientOffInd);
            else
                stimOffPd(i) = NaN;
            end
        end

        % Replace PD-based times with those aligned to QuadState
        onARDTimes  = stimOnPd(~isnan(stimOnPd));
        offARDTimes = stimOffPd(~isnan(stimOffPd));
    else
        warning('Only using PD to detect stim events.');
    end

    %% Remove stimulus events that occur after two-photon imaging has stopped 
    % Can move out??

    idxToDrop = find(onARDTimes >= max(twoPData(1).TwoPFrameTime));

    if ~isempty(idxToDrop) && max(idxToDrop) <= length(offARDTimes)
        onARDTimes(idxToDrop)  = [];
        offARDTimes(idxToDrop) = [];
    else
        warning('Skipping deletion: No valid indices or out-of-bounds index.');
    end


    %% Plot the results if requested
    if plotFlag
        figure;
        plot(pdArduinoTime, pdValue, 'k'); hold on;
        scatter(onARDTimes, repmat(10, size(onARDTimes)), 'g', 'filled');
        scatter(offARDTimes, repmat(10, size(offARDTimes)), 'r', 'filled');
        legend('Raw PD', 'PD ON', 'PD OFF');
        
        % If QuadState data was used, plot its timing as well
        if exist('quadOffArduinoTime', 'var')
            scatter(quadOffArduinoTime, repmat(6, size(quadOffArduinoTime)), 'b', 'filled');
            scatter(quadOnArduinoTime, repmat(6, size(quadOnArduinoTime)), 'y', 'filled');
            legend('Raw PD', 'PD ON', 'PD OFF', 'Quad OFF', 'Quad ON');
        end

        xlabel('Arduino Time (s)');
        ylabel('Photodiode Value');
        title('Stimulus Onset Detection');
        hold off;
    end

    %% Build the output structure and save it
    % Should we save in responseData??
    bonsaiData.onARDTimes  = onARDTimes;
    bonsaiData.offARDTimes = offARDTimes;
    bonsaiData.isStim = isStim;
    save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo')
end