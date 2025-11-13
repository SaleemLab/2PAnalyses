function figHandle = plotBoutonSomaComparison(allData, mouseID, varargin)
% plotBoutonSomaComparison Plots Day 1-5 (Boutons) and Pooled Somas for ONE mouse.

    p = inputParser;
    addParameter(p, 'SavePath', '', @ischar);
    parse(p, varargin{:});

    figHandle = figure('Position', [100 100 1800 400]);
    tiledlayout(figHandle, 1, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 1. Plot Boutons (Days 1-5)
    for day = 1:5
        nexttile;
        % Find data for this mouse, this day, AND Type='Boutons'
        session = allData(strcmp({allData.MouseID}, mouseID) & ...
                          [allData.Day] == day & ...
                          strcmpi({allData.Type}, 'Boutons'));
        
        if isempty(session)
            title(['Day ' num2str(day) ' (No Data)']); axis off; continue;
        end
        
        plotTuningInAxes(gca, session(1).OddMean, session(1).EvenMean);
        title(['Day ' num2str(day) ' Boutons']);
        if day==1, ylabel(mouseID, 'Interpreter','none'); end
    end

    % 2. Plot Pooled Somas (All Days combined)
    nexttile;
    % Find ALL Soma sessions for this mouse
    somaSessions = allData(strcmp({allData.MouseID}, mouseID) & ...
                           strcmpi({allData.Type}, 'Somas'));
                           
    if isempty(somaSessions)
        title('Somas (No Data)'); axis off;
    else
        % Concatenate all soma sessions together!
        pooledOdd = vertcat(somaSessions.OddMean);
        pooledEven = vertcat(somaSessions.EvenMean);
        
        plotTuningInAxes(gca, pooledOdd, pooledEven);
        title('Pooled Somas');
        colorbar; % Only on last plot
    end
    
    if ~isempty(p.Results.SavePath), saveas(figHandle, p.Results.SavePath); end
end