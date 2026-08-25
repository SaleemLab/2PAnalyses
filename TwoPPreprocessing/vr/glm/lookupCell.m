function lookupCell(T, cellID)
% LOOKUPCELL  Print category membership for a single cell.
%   T       CellCategoryTable (from 'CellCategory-table' figname), or
%           call with [] to auto-load 'CellCategoryTable' from base workspace.
%   cellID  Either a cell-name string (e.g. 'M25131_20260318_cell#5')
%           or a numeric Row index.
%
% USAGE
%   lookupCell([], 'M25131_20260318_cell#5')
%   lookupCell([], 4856)   % by Row index

    if isempty(T)
        T = evalin('base', 'CellCategoryTable');
    end

    if ischar(cellID) || isstring(cellID)
        row = find(strcmp(T.CellName, cellID));
    else
        row = find(T.Row == cellID);
    end

    if isempty(row)
        fprintf('Cell not found (not in goodcells, or name/ID typo): %s\n', string(cellID));
        return;
    end

    r = T(row,:);
    fprintf('\n--- %s ---\n', r.CellName{1});
    fprintf('Row: %d | Animal: %d | Session: %d\n', r.Row, r.Animal, r.Session);
    fprintf('  Spatial:         %d  (LLHrel = %.4f)\n', r.IsSpatial, r.LLHrel);
    fprintf('  Significant pos: %d\n', r.IsSignif);
    fprintf('  Omission:        %d  (LLHrel_omit = %.4f, p = %.3g)\n', r.IsOmit, r.LLHrel_omit, r.pval_omit);
    fprintf('  Landmark-driven: %d  (LLHrel_L1L2 = %.4f)\n', r.IsLandmarkDriven, r.LLHrel_L1L2);
    fprintf('  BG-driven:       %d  (LLHrel_BG = %.4f)\n', r.IsBGDriven, r.LLHrel_BG);
    fprintf('  Speed-driven:    %d\n', r.IsSpeedDriven);
end