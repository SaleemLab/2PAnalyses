function dF_F = get_delta_F_over_F(Fc, F0)
    % Function converted to Matlab from Schroeder Lab  
    % Calculates delta F over F.
    % Instead of simply dividing (F-F0) by F0,
    % the mean of F0 is used and only values above 1 are taken.
    % This prevents incorrectly increasing F values if F0 is smaller than 1.
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
    % 

    dF_F = (Fc - F0) ./ max(1, nanmean(F0, 1));
    % Alternative: dF_F = (Fc - F0) ./ max(1, F0);
end
