function plotBandpassTroughCheck(pairs, targetStruct, useField, minR2Threshold, pval_shuffleThreshold)
% PLOTBANDPASSTROUGHCHECK
% Loops over all sessions and plots all band-pass and trough-inverted cells
% that pass R2 and p-value thresholds, for visual inspection.
%
% USAGE:
%   pairs.M25132 = {'20260226', '20260228', '20260313'};
%   pairs.M25133 = {'20260224'};
%   pairs.M26003 = {'20260322', '20260324', '20260325'};
%   plotBandpassTroughCheck(pairs, 'tuningCurve', 'dFFNeuropilCorrected', 0.1, 0.01)

    if nargin < 3, useField           = 'dFFNeuropilCorrected'; end
    if nargin < 4, minR2Threshold     = 0.1; end
    if nargin < 5, pval_shuffleThreshold = 0.01; end

    filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'Darkness', 'GrayScreen'});
    mouseInfo     = sessionsToProcess(filteredTable);

    for thisMouse = 1:size(mouseInfo, 1)
        mousenumber  = mouseInfo{thisMouse, 1};
        sessionNames = mouseInfo{thisMouse, 2};

        for thisSession = 1:length(sessionNames)
            sessionName = sessionNames{thisSession};

            infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
            if isempty(infoPath) || ~isfile(infoPath), continue; end

            loadedInfo      = load(infoPath, 'sessionFileInfo');
            sessionFileInfo = loadedInfo.sessionFileInfo;
            stimNames       = {sessionFileInfo.stimFiles.name};
            targetIdx       = find(contains(stimNames, {'Darkness', 'GrayScreen'}));

            for thisStim = 1:length(targetIdx)
                thisStimName = stimNames{targetIdx(thisStim)};
                stimFileName = sprintf('%s_%s_Response_%s.mat', mousenumber, sessionName, thisStimName);
                fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

                if ~isfile(fileFullPath), continue; end

                loadedData = load(fileFullPath, 'response');
                resp       = loadedData.response;

                if ~isfield(resp, targetStruct) || ...
                   ~isfield(resp.(targetStruct), useField) || ...
                   ~isfield(resp.(targetStruct).(useField), 'classification')
                    fprintf('Skipping %s %s %s — no classification\n', mousenumber, sessionName, thisStimName);
                    continue;
                end

                cls        = resp.(targetStruct).(useField).classification;
                pvalMoving = resp.(targetStruct).(useField).pValMoving;

                bpIdx = find(strcmp(cls.tuningType, 'bandpass') & ...
                             cls.R2 >= minR2Threshold & pvalMoving <= pval_shuffleThreshold);
                trIdx = find(strcmp(cls.tuningType, 'trough_inverted') & ...
                             cls.R2 >= minR2Threshold & pvalMoving <= pval_shuffleThreshold);

                fprintf('%s %s %s — BP: %d, TR: %d\n', ...
                    mousenumber, sessionName, thisStimName, numel(bpIdx), numel(trIdx));

                if isempty(bpIdx) && isempty(trIdx), continue; end

                edges         = resp.(targetStruct).speedBins;
                movingCenters = (edges(1:end-1) + diff(edges)/2)';
                yMean         = resp.(targetStruct).(useField).moveMean;
                ySEM  = resp.(targetStruct).(useField).moveSEM;

                plotIdx    = [bpIdx(:); trIdx(:)];
                typeLabels = [repmat({'BP'}, numel(bpIdx), 1); ...
                              repmat({'TR'}, numel(trIdx), 1)];
                typeColors = [repmat({[0.96 0.70 0.50]}, numel(bpIdx), 1); ...
                              repmat({[0.53 0.78 0.92]}, numel(trIdx), 1)];

                nPlot = numel(plotIdx);
                nCols = min(4, nPlot);
                nRows = ceil(nPlot / nCols);

                figure('Color', 'w', 'Position', [50 50 300*nCols 250*nRows], ...
                    'Name', sprintf('%s %s %s', mousenumber, sessionName, thisStimName));

                for k = 1:nPlot
                    r = plotIdx(k);
                    subplot(nRows, nCols, k);
                    hold on;

                    % Raw data
                    %                     plot(movingCenters, yMean(r,:), 'ko-', ...
                    %                         'MarkerSize', 4, 'MarkerFaceColor', [0.2 0.2 0.2], 'LineWidth', 0.8);

                    % Raw data

         
                   
                    errorbar(movingCenters, yMean(r,:), ySEM(r,:), 'ko', ...
                        'MarkerSize', 4, 'MarkerFaceColor', [0.2 0.2 0.2], 'LineWidth', 0.8, ...
                        'LineStyle', 'none');
                    % Fitted curve
                    p = cls.fitParams(r,:);
                    if all(isfinite(p))
                        gaussFun = @(params, xdata) params(1) + params(2) .* ...
                            exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
                        xDense = linspace(min(movingCenters), max(movingCenters), 200);
                        yFit   = gaussFun(p, xDense);
                        plot(xDense, yFit, '-', 'Color', typeColors{k}, 'LineWidth', 2);
                    end

                    title(sprintf('%s ROI %d | R2=%.2f | p=%.3f', ...
                        typeLabels{k}, r, cls.R2(r), pvalMoving(r)), 'FontSize', 7);
                    xlabel('Speed (cm/s)', 'FontSize', 7);
                    box off;
                    set(gca, 'TickDir', 'out', 'FontSize', 7);
                end

                sgtitle(sprintf('%s  %s  %s  |  BP=%d  TR=%d', ...
                    mousenumber, sessionName, thisStimName, numel(bpIdx), numel(trIdx)), ...
                    'FontSize', 9, 'FontWeight', 'bold');

            end % thisStim
        end % thisSession
    end % thisMouse
end
