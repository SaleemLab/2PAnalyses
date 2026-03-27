function plotPopulationConditionAverage(response, applyNormalisation)
% plotPopulationConditionAverage: Averages activity across all 
% ROIs for each condition, correctly mapping trial IDs to matrix rows.

if nargin < 2, applyNormalisation = true; end

% Data extraction 
if isfield(response.lapPositionActivity, 'dFF')
    data = response.lapPositionActivity.dFF;
else
    error('dFF field not found in response.lapPositionActivity.');
end

[nROIs, nRows, nBins] = size(data); % nRows is the number of completed laps 
conds = fieldnames(response.trialIndicesByCondition);
colors = lines(length(conds));
compLaps = response.completedLaps; % The "Translation Key" to match Trial IDs to Matrix Rows 

% Initialize ROI-level means: [ROIs x Conditions x Bins]
roiCondMeans = nan(nROIs, length(conds), nBins);

fprintf('Calculating population average for %d ROIs across %d matrix rows.\n', nROIs, nRows);

for iC = 1:length(conds)
    cName = conds{iC};
    absIDs = response.trialIndicesByCondition.(cName); % These are absolute Trial IDs 
    
    % MAP TRIAL IDs TO MATRIX ROWS
    % Find which indices in 'compLaps' match the 'absIDs' for this condition 
    rowIdx = find(ismember(compLaps, absIDs));
    
    if ~isempty(rowIdx)
        for iN = 1:nROIs
            % Extract specific matrix rows for this ROI and condition
            roiTraceLaps = squeeze(data(iN, rowIdx, :));
            
            % Handle single vs multiple laps for mean calculation 
            if length(rowIdx) == 1
                meanTrace = roiTraceLaps'; % Ensure it is a row vector
            else
                meanTrace = mean(roiTraceLaps, 1, 'omitnan');
            end
            
            % Normalise ROI to its own peak if requested 
            if applyNormalisation
                peakVal = max(meanTrace);
                if peakVal > 0
                    meanTrace = meanTrace / peakVal;
                end
            end
            
            roiCondMeans(iN, iC, :) = meanTrace;
        end
    end
end

% Average across the ROI dimension (dim 1) 
grandMean = squeeze(mean(roiCondMeans, 1, 'omitnan'));
% Calculate SEM across ROIs 
grandSEM  = squeeze(std(roiCondMeans, 0, 1, 'omitnan')) ./ sqrt(nROIs);

% --- Plotting ---
figure('Color', 'w', 'Name', 'Population Condition Average', 'Position', [100 100 900 600]);
hold on;
legendEntries = {};

for iC = 1:length(conds)
    % Skip if no data for this condition
    if all(isnan(grandMean(iC, :))), continue; end
    
    x = 1:nBins;
    mu = grandMean(iC, :);
    err = grandSEM(iC, :);
    
    % SEM Shading 
    fill([x fliplr(x)], [mu+err, fliplr(mu-err)], colors(iC,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
    
    % Grand Mean Line 
    plot(x, mu, 'Color', colors(iC,:), 'LineWidth', 3);
    
    % Count contributing ROIs for the legend
    nActiveROIs = sum(~isnan(roiCondMeans(:, iC, 1)));
    legendEntries{end+1} = [strrep(conds{iC}, '_', ' ') ' (n=' num2str(nActiveROIs) ' ROIs)'];
end

% Visual Landmark References 
for pos = [40 80 120 160]
    xline(pos, 'k--', 'Alpha', 0.3, 'HandleVisibility', 'off');
end

grid on; box off;
xlabel('Position (cm)'); 

% FIXED SYNTAX ERROR HERE:
if applyNormalisation
    ylabel('Normalized Pop. Activity');
else
    ylabel('Average dFF');
end

title(['Population Response (n=' num2str(nROIs) ' Total ROIs)']);
xticks([1 40 80 120 160 200]);
xticklabels({'0', '40', '80', '120', '160', '200'});
legend(legendEntries, 'Location', 'northeastoutside');

end