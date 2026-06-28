function response = getRunningAndPupilDataByTrials(sessionFileInfo, stimName)
% Aligns wheel speed and pupil area to trial onset.
%
% Inputs:
%   sessionFileInfo - Struct containing path information for the session
%   stimName        - String, name of the stimulus file to process

%% Find and load stimulus 
iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
if isempty(iStim)
    error('Stimulus name "%s" not found in sessionFileInfo.', stimName);
end

load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
load(sessionFileInfo.stimFiles(iStim).Response, 'response');

%% time vector
pre = response(1).preStimTime;
post = response(1).postStimTime;
Time = -pre : (1/60) : post;  % 60 is the interpolation rate
nWheelValue = length(Time);
nGroups = numel(bonsaiData.trialGroups);

%%  compute continuous wheel speed
tickToCmConversion = 3.1415 * 20 / 1024;  % Wheel radius 20 cm, 1024 ticks per revolution

%  displacement 
wheelVals = peripheralData.Wheel.Value;
displacement = [0; diff(wheelVals(:) * tickToCmConversion)];
 
% continuous speed (in cm/s)
wheelTimes = peripheralData.Wheel.sampleTimes(:);
timeDiffs = [0; diff(wheelTimes)];
% If two samples happen within less than 1 ms, use a safe nominal 
% frame rate delta (e.g., 1/60s), or drop the sample to prevent division explosion.
timeDiffs(timeDiffs < 0.001) = 0.0167;
% timeDiffs(timeDiffs <= 0) = 0.001; % Prevent division by zero

wheelSpeedContinuous = displacement ./ timeDiffs;

wheelSpeedContinuous(wheelSpeedContinuous > 100) = 0;
wheelSpeedContinuous(wheelSpeedContinuous < -100) = 0;

%% align data to trial onsets (based on trial groups; eg speed 16 group for dot fields or 180 deg grating for direction tuning) 
for thisGroup = 1:nGroups
    grp = bonsaiData.trialGroups(thisGroup);
    trIdxs = grp.trials;
    nTrialsInGrp = numel(trIdxs);
    
    % 2D matrices for 1D behavioral data (Frames x Trials)
    alignedWheel = nan(nWheelValue, nTrialsInGrp);
    alignedPupil = nan(nWheelValue, nTrialsInGrp);
    
    for ti = 1:nTrialsInGrp
        trialID = trIdxs(ti);
        
        % Skip bad trials
        if response(1).badTrialMask(trialID)
            continue; 
        end
        
        % Check if response frames exist for this trial
        if trialID > numel(response(1).responseFrameIdx) || isempty(response(1).responseFrameIdx{trialID})
            continue;
        end
        
        fMask = response(1).responseFrameIdx{trialID};
        rawRelTimes = response(1).responseFrameRelTimes{trialID}; % Jittered times
        
        % First align wheel speed 
        rawWheel = wheelSpeedContinuous(fMask);
        validWheel = ~isnan(rawRelTimes) & ~isnan(rawWheel');
        
        if sum(validWheel) > 10
            % Re-align to the time base 
            traceWheel = interp1(rawRelTimes(validWheel), rawWheel(validWheel), Time, 'linear', 'extrap');
            % Running is usually not baseline-subtracted, absolute cm/s is kept
            alignedWheel(:, ti) = traceWheel;
        end
        
       % Then align the pupil area if present; only handing pupil here 
        if isfield(peripheralData, 'Pupil') && isfield(peripheralData.Pupil, 'Value')
            rawPupil = peripheralData.Pupil.Value.Area(fMask);
            validPupil = ~isnan(rawRelTimes) & ~isnan(rawPupil');

            if sum(validPupil) > 10
                tracePupil = interp1(rawRelTimes(validPupil), rawPupil(validPupil), Time, 'linear', 'extrap');

                % Baseline subtract pupil if this makes sense? (relative dilation; currently unsure)
                % baselineVal = nanmean(tracePupil(Time < 0));
                % alignedPupil(:, ti) = tracePupil - baselineVal;

                alignedPupil(:, ti) = tracePupil; % Storing raw area for now
            end
        end
    end
    
    % compute and store for wheeldata 
    response(1).wheelData(thisGroup).stimValue        = grp.value;
    response(1).wheelData(thisGroup).alignedResponses = alignedWheel;
    % response(1).wheelData(thisGroup).meanResponse     = mean(alignedWheel, 2, 'omitnan');
    % response(1).wheelData(thisGroup).stdResponse      = std(alignedWheel, 0, 2, 'omitnan');
    % response(1).wheelData(thisGroup).semResponse      = std(alignedWheel, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(alignedWheel), 2));
    response(1).wheelData(thisGroup).timeVector       = Time;
    
    % compute for pupil data 
    if isfield(peripheralData, 'Pupil')
        response(1).pupilData(thisGroup).stimValue        = grp.value;
        response(1).pupilData(thisGroup).alignedResponses = alignedPupil;
        % response(1).pupilData(thisGroup).meanResponse     = mean(alignedPupil, 2, 'omitnan');
        % response(1).pupilData(thisGroup).stdResponse      = std(alignedPupil, 0, 2, 'omitnan');
        % response(1).pupilData(thisGroup).semResponse      = std(alignedPupil, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(alignedPupil), 2));
        response(1).pupilData(thisGroup).timeVector       = Time;
    end
end

%% Save Results
save(sessionFileInfo.stimFiles(iStim).Response, 'response');
fprintf('Wheel and Pupil data saved to response for %s.\n', stimName);

end