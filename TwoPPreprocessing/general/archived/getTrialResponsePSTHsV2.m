function psthData = getTrialResponsePSTHsV2(sessionFileInfo, stimName, method, interpRate)
    % getTrialResponsePSTHs - Extract and process trial-aligned responses using real time values.
    %
    % Inputs:
    %   sessionFileInfo (struct) - Structure containing file paths and session details.
    %   stimName (char) - Name of the stimulus condition.
    %   method (int) - Choice of PSTH calculation method:
    %                  1: Interpolate & average
    %                  2: Collapse bins & smooth & SEM 
    %                  3: Direct mean response (without interpolation) 
    %   interpRate (float, optional) - Interpolation rate in Hz (default: 60 Hz)
    %
    % Output:
    %   psthData (struct) - Contains mean, std, and SEM of trial responses per stimulus.

    if nargin < 4 
        interpRate = 60; 
    end

    % Locate the requested stimulus in sessionFileInfo
    isStim = strcmp(stimName, {sessionFileInfo.stimFiles.name});
    iStim = find(isStim, 1);
    
    % Check if the requested stimulus exists
    if isempty(iStim)
        error('Stimulus name "%s" not found in sessionFileInfo.', stimName);
    end
    
    % Load Response Data
    if exist(sessionFileInfo.stimFiles(iStim).Response, 'file') && ...
            exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file')
        load(sessionFileInfo.stimFiles(iStim).Response, 'response');
        load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData')
    else
        error('Missing response file for stimulus "%s".', stimName);
    end
    
    % Load Two-Photon Data
    if exist(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'file')
        load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData');
    else
        error('Missing two-photon data file for stimulus "%s".', stimName);
    end
    
    % Initialize output structure
    psthData = struct();
    numStimuli = length(bonsaiData.trialGroups);
    numNeurons = size(twoPData.F, 1);
    
    for i = 1:numStimuli
        % Get stimulus condition details
        stimValue = bonsaiData.trialGroups(i).value;
        trialIndices = bonsaiData.trialGroups(i).trials;
        trialIndices = trialIndices(~isnan(trialIndices)); % Remove NaNs
        numTrials = length(trialIndices);
        
        % Determine max trial duration
        maxTime = max(cellfun(@max, response.responseFrameRelTimes(trialIndices)));
        minTime = min(cellfun(@min, response.responseFrameRelTimes(trialIndices)));
        timeVector = linspace(minTime, maxTime, interpRate);
        
        % Preallocate response matrix (Neurons x Time x Trials)
        alignedResponses = NaN(numNeurons, length(timeVector), numTrials);
        
        % Extract trial-aligned responses
        for t = 1:numTrials
            trialIdx = trialIndices(t);
            selectedFrames = response.responseFrameIdx{trialIdx};
            trialData = twoPData.F(:, selectedFrames);
            trialTimes = response.responseFrameRelTimes{trialIdx};
            
            % Compute F0 using pre-stimulus frames
            preStimDuration = response.preStimTime; % in seconds
            preStimFrames = find(trialTimes <= 0 & trialTimes >= -preStimDuration);
            if isempty(preStimFrames)
                error('No frames found for F0 computation in pre-stimulus period.');
            end
            prctl_F = 8; % Default percentile for F0 calculation
            
            % Compute actual sampling rate from raw data
            fs = 1 / mean(diff(trialTimes)); 
            
            % Set window size based on pre-stimulus period using actual sampling rate
            window_size = round(preStimDuration * fs); 
            
            % Compute F0 using pre-stim frames and actual frame rate
            F0 = get_F0(trialData(:, preStimFrames), fs, prctl_F, window_size);
%             F0 = mean(F0, 2, 'omitnan'); % Compute mean across pre-stim frames for each neuron separately
%             F0 = repmat(F0, 1, size(trialData, 2)); % Expand to match all time points
%             
            % Compute ΔF/F
            trialData = get_delta_F_over_F(trialData, F0);
            
            % Interpolate onto common time vector
            for neuronIdx = 1:numNeurons
                alignedResponses(neuronIdx, :, t) = interp1(trialTimes, trialData(neuronIdx, :), timeVector, 'linear', NaN);
            end
        end
        
        % Compute PSTH based on selected method
        switch method
            case 1 % Interpolate & average
                psthMean = nanmean(alignedResponses, 3);
                psthStd = nanstd(alignedResponses, 0, 3);
            
            case 2 % Temporal smoothing
                w = gausswin(4); % 9 for VR 
                w = w / sum(w);
                psthMean = nanmean(alignedResponses, 3);
                nan_values = isnan(psthMean);
                psthMean(nan_values) = nanmean(psthMean(:));              
                psthMean = filtfilt(w, 1, psthMean')';
                psthStd = nanstd(alignedResponses, 0, 3);
                psthStd = filtfilt(w, 1, psthStd')';
            
            case 3 % Direct mean without interpolation
                maxFrames = max(cellfun(@(x) sum(x), response.responseFrameIdx(trialIndices)));
                directMeanResponse = NaN(numNeurons, maxFrames, numTrials);
                
                for t = 1:numTrials
                    trialIdx = trialIndices(t);
                    selectedFrames = response.responseFrameIdx{trialIdx};
                    trialData = twoPData.F(:, selectedFrames);
                    directMeanResponse(:, 1:size(trialData, 2), t) = trialData;
                end
                
                psthMean = nanmean(directMeanResponse, 3);
                psthStd = nanstd(directMeanResponse, 0, 3);
            
            otherwise
                error('Invalid method. Choose 1 (Interpolate), 2 (Smooth), or 3 (Mean)');
        end
        
        % Compute Standard Error of Mean (SEM)
        psthSEM = psthStd ./ sqrt(sum(~isnan(alignedResponses), 3));
        
        % Store results in struct
        psthData(i).stimValue = stimValue;
        psthData(i).alignedResponses = alignedResponses;
        psthData(i).meanResponse = psthMean;
        psthData(i).stdResponse = psthStd;
        psthData(i).semResponse = psthSEM;
        psthData(i).timeVector = timeVector; % Store actual time axis
        psthData(i).responseType = method; % Store method used for future reference
    end
end


