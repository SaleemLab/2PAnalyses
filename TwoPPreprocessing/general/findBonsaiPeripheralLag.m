function [bonsaiData] = findBonsaiPeripheralLag(sessionFileInfo, method, samplingRate)
% Find the lag between Bonsai-Peripheral data streams.
%
% Inputs:
%   - sessionFileInfo (struct): Contains paths and information about the session.
%   - method (int): Lag calculation method:
%       1: Use cross correlation of the Quad and PD signals (value)
%       across the full length of the recording and use the
%       maximum/best lag (r)
%       2: For each PDOn-Off 'block' align the quad start time with the
%       PD on Edges and interpolates time per block. This is a good technique
%       when the lag is inconsistent across the recording.
%       This is done for the quad and will be used to correct across all other bonsai
%       variables in the next function. @Does not work well when
%       Arduino is sampled at 100Hz.
%
%       - samplingRate (float, optional): Interpolation rate in Hz
%       (default: Sampling rate taken from TwoP resampled data);
%       Note: Keep sampingRate consistent across 2P,
%       Bonsai and Peripheral resamping scripts.
%
% Outputs:
%   - bonsaiData.LagInfo (depending on method used):
%         .method         - Method used to find lag.
%         .xcorrBestLag   - Frames to shift.
%         .samplingRate   - Use as input in this function.
%         .lagCorrBonsaiArduinoTime - Lag corrected Quad time.
%
% Usage Example:
%   [bonsaiData] = findBonsaiPeripheralLag(sessionFileInfo, 1, 60)
%
% Edited Masa and Diao's code - Sonali March 2025

if nargin < 2; method = 1; end % Set default
if nargin < 3; samplingRate = 60; end % Set default

for iStim = 1:length(sessionFileInfo.stimFiles)
    % Load response and stimulus data
    if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
            exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file')
        load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
        load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
    else
        warning('Missing Peripheral or Bonsai data for stimulus "%s".', sessionFileInfo.stimFiles(iStim).name);
        continue
    end

    % Pick out raw values
    rawPDValue              = peripheralData.Photodiode.rawValue;
    rawPDArduinoTime        = peripheralData.Photodiode.rawArduinoTime;
    rawQuadValue            = bonsaiData.Quadstate.rawValue;
    rawQuadArduinoTime      = bonsaiData.Quadstate.rawArduinoTime;
    bonsaiData.LagInfo.method = method;

    switch method
        case 1 % Cross correlation of the two signals to find the best lag.
            % Interpolate quad and pd to 60hz using the quadTime
            tempQuadValueInterpolated = interp1(rawQuadArduinoTime, rawQuadValue, rawQuadArduinoTime(1):1/samplingRate:rawQuadArduinoTime(end), 'previous')';
            tempPDValueInterpolated = interp1(rawPDArduinoTime, rawPDValue, rawQuadArduinoTime(1):1/samplingRate:rawQuadArduinoTime(end), 'linear')';
            [r, lags] = xcorr(normalize(tempQuadValueInterpolated), normalize(tempPDValueInterpolated));
            [~, jointIdx] = max(r);
            bestLag = lags(jointIdx);

            % Save lag info
            bonsaiData.LagInfo.xcorrBestLag = bestLag;
            bonsaiData.LagInfo.samplingRate = samplingRate;
            % For plotting
            lagCorrArdinoTime = rawQuadArduinoTime-(bestLag*(1/samplingRate));

            figure;
            hold on;
            plot(rawPDArduinoTime(1:1000), rawPDValue(1:1000), 'DisplayName', 'Original Photodiode Signal');
            plot(rawQuadArduinoTime(1:1000), 10*rawQuadValue(1:1000), 'r', 'DisplayName', 'Uncorrected Quad');
            plot(lagCorrArdinoTime(1:1000), 12*rawQuadValue(1:1000), 'g', 'DisplayName', 'Corrected Quad');
            title('Using (xcorr) method 1: Photodiode vs Corrected/Uncorrected Quad Signal');
            legend();
            hold off;

        case 2
            rawPDValueSmoothed = smoothdata(rawPDValue,'movmedian',10);
            pdONOFF = rawPDValueSmoothed >= mean(rawPDValue);  % Photodiode ON when output is above or equal to 4
            pdON = find(diff(pdONOFF) == 1) + 1;  % Photodiode turns ON
            pdOFF = find(diff(pdONOFF) == -1) + 1; % Photodiode turns OFF

            blockLength=[];
            % Calculate the duration of each ON block (ON to OFF period)
            if length(pdON) == length(pdOFF)
                blockLength = abs(pdON - pdOFF); % Block length in terms of index differences
            elseif length(pdON) > length(pdOFF)
                for i = 1:length(pdOFF)
                    blockLength(i) = abs(pdON(i) - pdOFF(i)); % Handle mismatch in length
                end
            else
                for i = 1:length(pdON)
                    blockLength(i) = abs(pdON(i) - pdOFF(i)); % Handle mismatch in length
                end
            end

            blocksInd1 = find(blockLength > 1);
            % Identify trial start and end points based on quad state changes
            idxTrialStart = find(diff(rawQuadValue) == 1) + 1; % Index of quad state change from 0 to 1 (start)
            idxTrialEnd = find(diff(rawQuadValue) == -1) + 1;  % Index of quad state change from 1 to 0 (end)
            quadTrialStartTime = rawQuadArduinoTime(idxTrialStart)'; % Get start times
            quadTrialEndTime = rawQuadArduinoTime(idxTrialEnd)';     % Get end times

            % Clean up: Remove trials with invalid start or end times (where Arduino inputs equals 0)
            quadTrialStartTime(isnan(quadTrialStartTime)) = [];
            idxTrialStart(isnan(quadTrialStartTime)) = [];
            quadTrialEndTime(isnan(quadTrialEndTime)) = [];
            idxTrialEnd(isnan(quadTrialEndTime)) = [];

            nTrials = length(quadTrialStartTime); % Number of trials


            if abs(length(blocksInd1) - length(idxTrialStart)) / length(idxTrialStart) < 1
                disp('Aligning Arduino timestamps based on the delay between photodiode and quad state');
                pdstart = []; % Initialize for photodiode start times
                tic
                H = waitbar(0, 'Finding photodiode ON edges'); % Progress bar for long-running operation
                % Loop through each trial to align quad start times with photodiode ON edges
                for itrial = 1:length(quadTrialStartTime)
                    waitbar(itrial / length(quadTrialStartTime), H); % Update progress bar

                    % Find the next photodiode ON time that aligns with the quad trial start
                    tempIdx = find(rawPDArduinoTime(pdON(blocksInd1)) >= quadTrialStartTime(itrial), 1, 'first');
                    idxStart = pdON(blocksInd1(tempIdx)) - 1;

                    if ~isempty(idxStart)
                        pdstart(itrial) = rawPDArduinoTime(idxStart); % Store aligned photodiode start time
                    else
                        pdstart(itrial) = nan; % Assign NaN if no alignment found
                    end
                end
                toc
                pdend = []; % Initialize for photodiode end times
                tic
                H = waitbar(0, 'Finding photodiode OFF edges'); % Progress bar for long-running operation
                % Loop through each trial to align quad end times with photodiode OFF edges
                for itrial = 1:length(quadTrialEndTime)
                    waitbar(itrial / length(quadTrialEndTime), H); % Update progress bar

                    % Find the next photodiode OFF time that aligns with the quad trial end
                    tempIdx = find(rawPDArduinoTime(pdOFF(blocksInd1)) >= quadTrialEndTime(itrial), 1, 'first');
                    idxEnd = pdOFF(blocksInd1(tempIdx)) - 1;

                    if ~isempty(idxEnd)
                        pdend(itrial) = rawPDArduinoTime(idxEnd); % Store aligned photodiode end time
                    else
                        pdend(itrial) = nan; % Assign NaN if no alignment found
                    end
                end

                toc
                peripheralData.LagInfo.corrPDOn = pdstart';
                peripheralData.LagInfo.corrPDOff = pdend';


                %             meanLag = mean(quadTrialEndTime-pdend');
                %             CorrQuadArduinoTime = rawQuadArduinoTime-meanLag;


                syncTimesQuad = quadTrialEndTime;
                syncTimesPhotodiode = pdend';

                [r, lags] = xcorr(normalize(diff(syncTimesQuad)), normalize(diff(syncTimesPhotodiode)));
                [~, jointIdx] = max(r);
                bestLag = lags(jointIdx);

                if bestLag > 0
                    syncTimesQuad = syncTimesQuad(syncTimesQuad>syncTimesQuad(bestLag));
                end
                % ...and rerun the xcorr
                [r, lags] = xcorr(diff(syncTimesQuad), diff(syncTimesPhotodiode));
                [~, jointIdx] = max(r);
                bestLag = lags(jointIdx);

                if bestLag < 0
                    nSyncOffset = -bestLag+1;
                    t_npix = syncTimesPhotodiode(nSyncOffset:end);
                    if length(t_npix) > length(syncTimesQuad)
                        t_npix = t_npix(1:length(syncTimesQuad));
                    end
                    t_bonsai = syncTimesQuad(1:numel(t_npix));
                else
                    nSyncOffset = bestLag+1;
                    t_npix = syncTimesPhotodiode(1:end);
                    if length(t_npix) > length(syncTimesQuad)
                        t_npix = t_npix(1:length(syncTimesQuad));
                    end
                    t_bonsai = syncTimesQuad(nSyncOffset:numel(t_npix)+nSyncOffset-1);
                end
                tt = interp1(t_bonsai, t_npix, rawQuadArduinoTime,'linear','extrap');
                %     tt = interp1(unique(bonsai_data.corrected_sglxTime),unique(bonsai_data.corrected_sglxTime),linspace(bonsai_data.corrected_sglxTime(1),bonsai_data.corrected_sglxTime(end),length(bonsai_data.corrected_sglxTime)),'linear')';
                [unique_t,index,~] = unique(tt);
                %     xx= interp1(index,unique_t,1:length(tt),'linear');
                lagCorrQuadArduinoTime = interp1(index,unique_t,1:length(tt),'linear','extrap');

                figure;
                hold on;
                plot(rawPDArduinoTime(1:1000), rawPDValue(1:1000), 'DisplayName', 'Original Photodiode Signal');
                plot(rawQuadArduinoTime(1:1000), 10*rawQuadValue(1:1000), 'r', 'DisplayName', 'Uncorrected Quad');
                plot(lagCorrQuadArduinoTime(1:1000), 12*rawQuadValue(1:1000), 'g', 'DisplayName', 'Corrected Quad');
                scatter(pdstart(1:40), 7*ones(1, 40), 'DisplayName', 'PD On Time'); % Scatter plot photodiode start times
                scatter(pdend(1:40), 6*ones(1, 40), 'DisplayName', 'PD off Time');   % Scatter plot photodiode end times
                title('Method 2: Photodiode vs Corrected Quad Signal');
                legend();
                hold off;

                bonsaiData.LagInfo.lagCorrBonsaiArduinoTime = lagCorrQuadArduinoTime';
            else
                disp('Mismatch in the number of photodiode ON events and quad trial start times. ReRun with Method 1');
            end
        otherwise
            error('Invalid method. Choose 1 (xCorr across the full length), 2 (Align quad start with PD On edge)');
    end

   
    save( sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData', '-append' );
    save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

end
end

