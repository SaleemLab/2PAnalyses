function [cellIdx, cellNames, T] = lookupCell(csvPath, criteria)
%FINDCELLS_SONALI_JF  Select cells from the exported metrics CSV
%   (cell_metrics_summary.csv, from plotFinalFigures_Sonali_JF's
%   'ExportCellMetrics' option) that satisfy a named set of criteria.
%
% EXAMPLE USAGE
%   [idx, names] = findCells_Sonali_JF('cell_metrics_summary.csv', 'visual')
%   [idx, names] = findCells_Sonali_JF('cell_metrics_summary.csv', 'spatial')
%   [idx, names] = findCells_Sonali_JF('cell_metrics_summary.csv', 'nonspatial')
%   [idx, names] = findCells_Sonali_JF('cell_metrics_summary.csv', 'omission')
%   [idx, names, T] = findCells_Sonali_JF('cell_metrics_summary.csv', 'flatnonspatial')
%
% INPUTS
%   csvPath    Path to the CSV exported by plotFinalFigures_Sonali_JF's
%              'ExportCellMetrics' option.
%   criteria   String selecting the cell-selection rule:
%       'good'          - goodcells == 1
%       'visual'        - goodcells == 1 & LLHi_vis is finite and > 0
%                          (i.e. vision explains this cell better than
%                          the null model; there's no separate "visual
%                          significance" column in the CSV, so this uses
%                          LLHi_vis as the closest available proxy - see
%                          note below)
%       'spatial'       - spatialcells == 1
%       'nonspatial'    - goodcells == 1 & spatialcells == 0
%       'flatnonspatial'- nonspatial cells with an exactly-flat position
%                         kernel (posKernel_maxabs < 1e-6), i.e. clean
%                         examples like cell 57/289 where VS and VSP are
%                         identical
%       'omission'      - goodOmitcells == 1
%
% OUTPUTS
%   cellIdx    Row indices into T (not EXP.Spk row numbers, unless your
%              CSV row order matches EXP.Spk order 1:1, which it does if
%              exported straight from ExportCellMetrics).
%   cellNames  Cell array of CellName strings for the matching rows.
%   T          Filtered table (all columns from the CSV) for those cells.
%
% NOTE ON 'visual': the exported CSV doesn't currently include a
% standalone "vision-only model p-value" column, only LLHi_vis (the
% log-likelihood improvement of the vision-only model over the null). If
% you want a true significance-based 'visual' criterion instead of this
% effect-size proxy, add EXP.GLMs{1}.Tuning(iVismodel).pval as its own
% column in the ExportCellMetrics block and reference it here.

T_full = readtable(csvPath);

switch criteria
    case 'good'
        mask = T_full.goodcells == 1;
    case 'visual'
        mask = T_full.goodcells == 1 & isfinite(T_full.LLHi_vis) & T_full.LLHi_vis > 0;
    case 'spatial'
        mask = T_full.spatialcells == 1;
    case 'nonspatial'
        mask = T_full.goodcells == 1 & T_full.spatialcells == 0;
    case 'flatnonspatial'
        mask = T_full.goodcells == 1 & T_full.spatialcells == 0 & T_full.posKernel_maxabs < 1e-6;
    case 'omission'
        mask = T_full.goodOmitcells == 1;
    otherwise
        error('Unrecognized criteria: %s', criteria);
end

T = T_full(mask, :);
cellNames = T.CellName;
cellIdx = find(mask);

end