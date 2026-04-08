function dF_F = get_delta_F_over_F(Fc, F0)
% Based on Sylvia Schroeder's Python function [get_delta_F_over_F] in repository 'depth-for-2p'
% Calculates delta F over F.
%
% Parameters
% ----------
% Fc : [t x nROIs] matrix
%     Calcium traces (measured signal) of ROIs.
% F0 : [t x nROIs] matrix
%     The baseline fluorescence (F0) traces of ROIs.
%
% Returns
% -------
% dF_F : [t x nROIs] matrix
%     Change in fluorescence (dF/F) of ROIs.

% Sylvia Schröder's function from 'depth-for-2p' repo
% Translated from python to matlab


% dF_F = (Fc - F0) ./ max(1, nanmean(F0, 1));
dF_F = (Fc - F0) ./ max(1, F0);
end
