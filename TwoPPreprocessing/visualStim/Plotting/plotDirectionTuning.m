function plotDirectionTuning(sessionFileInfo, response, doSmooth, onlyTuned)
% plotDirectionTuning: Linear PSTHs (Left) and Polar Tuning Plot (Right)
%
% Inputs:
%   sessionFileInfo : Struct containing animal_name, session_name, and Directories
%   response        : Struct containing psthData and results
%   doSmooth        : Boolean to apply gaussian smoothing (default: true)
%   onlyTuned       : Boolean to filter by responsiveness (default: false)

if nargin < 3; doSmooth = false; end
if nargin < 4; onlyTuned = false; end

%% --- Setup Directory and PDF Path ---
saveFolder = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~isfolder(saveFolder); mkdir(saveFolder); end

pdfPath = fullfile(saveFolder, [sessionFileInfo.animal_name, '_' ...
    sessionFileInfo.session_name '_DirTuning.pdf']);

if exist(pdfPath, 'file'), delete(pdfPath); end

%% --- Data Preparation ---
psthData = response.psthData; 
timeVec  = psthData(1).timeVector;
nROI     = size(psthData(1).alignedResponses, 1);
allDirs  = arrayfun(@(x) x.stimValue, psthData);

[uDirs, ~, ~] = unique(allDirs);
[sortedDirs, sortIdx] = sort(uDirs); 

% Analysis Window
respWin = [0.5 3]; 
respIdx = timeVec >= respWin(1) & timeVec <= respWin(2);
tWin = 5; % Smoothing window

%% --- Figure Loop ---
hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 500]);
fprintf('Generating Direction Tuning PDF for %d ROIs...\n', nROI);

for iROI = 1:nROI
    % Optional filter: check if ROI is responsive (if field exists in response)
    if onlyTuned && isfield(response, 'isResponsive')
        if ~response.isResponsive(iROI); continue; end
    end

    clf(hFig);
    means = zeros(length(uDirs), 1);
    sems  = zeros(length(uDirs), 1);
    
    % --- Panel A: PSTH (Linear) ---
    subplot(1,2,1); hold on;
    cmap = copper(length(uDirs)); 
    pHandles = gobjects(length(uDirs), 1);
    
    for d = 1:length(uDirs)
        trials = squeeze(psthData(d).alignedResponses(iROI, :, :)); 
        meanTrace = mean(trials, 2, 'omitnan');
        semTrace  = std(trials, 0, 2, 'omitnan') ./ sqrt(size(trials, 2));
        
        if doSmooth
            meanTrace = smoothdata(meanTrace, 'gaussian', tWin);
        end
        
        % Calculate average response for polar plot
        trialAverages = mean(trials(respIdx, :), 1, 'omitnan');
        means(d) = mean(trialAverages);
        sems(d)  = std(trialAverages) / sqrt(length(trialAverages));
        
        % Shaded Error
        t_col = timeVec(:);
        m_col = meanTrace(:);
        s_col = semTrace(:);
        fill([t_col; flipud(t_col)], [m_col-s_col; flipud(m_col+s_col)], ...
            cmap(d,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        pHandles(d) = plot(t_col, m_col, 'Color', cmap(d,:), 'LineWidth', 2);
    end
    
    xlim([-1, 3]);
    xline(0, 'k--', 'LineWidth', 1); 
    xlabel('Time (s)'); 
    ylabel('\DeltaF/F');
    set(gca, 'Box', 'off', 'TickDir', 'out');
   
    legLabels = string(uDirs) + "°";
    leg = legend(pHandles, legLabels, 'Location', 'bestoutside');
    title(leg, 'Directions');
    
    % --- Panel B: Polar Plot ---
    subplot(1,2,2);
    theta = deg2rad(sortedDirs);
    rho = means(sortIdx);
    
    % Close the loop for the polar plot
    theta_closed = [theta; theta(1)];
    rho_closed = [rho; rho(1)];
    
    polarplot(theta_closed, rho_closed, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'w');
    hold on;
    
    % Add error bars to polar plot
    for k = 1:length(theta)
        polarplot([theta(k) theta(k)], ...
            [means(sortIdx(k))-sems(sortIdx(k)), means(sortIdx(k))+sems(sortIdx(k))], 'k-');
    end
    
    ax = gca;
    ax.ThetaTick = 0:45:315;
    set(ax, 'TickDir', 'out');
    
    % --- Title ---
    sgtitle(sprintf('%s | %s | Bouton %d: Direction Tuning', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, iROI), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    exportgraphics(hFig, pdfPath, 'Append', true);
end

close(hFig);
fprintf('PDF saved to: %s\n', pdfPath);

end