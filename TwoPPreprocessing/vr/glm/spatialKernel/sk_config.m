function C = sk_config()
% Shared paths and constants for the spatial-kernel (position vs gain) analysis.

C.rawRoot   = 'D:\DATA\Sonali\GLM\ProcessedData\V1\analyzed';
C.tablePath = 'D:\DATA\Sonali\GLM\ForAS\V1_GLMData_table.csv';
C.outRoot   = 'D:\DATA\Sonali\GLM\TrialLevel';
C.figRoot   = fullfile(C.outRoot, 'figures');

% Sessions (mouse, date, sessionid as in V1_GLMData_table.csv)
C.sessions = {
    'M26004', 20260318, 11
    'M26004', 20260321, 12
    'M26004', 20260322, 13
    'M26005', 20260318, 21
    'M26005', 20260321, 22
    'M26005', 20260322, 23
    'M25126', 20260311, 31
    'M25126', 20260312, 32
    'M25126', 20260313, 33
    'M25131', 20260318, 41
    'M25131', 20260321, 42
    'M25131', 20260322, 43 };

C.condNames = {'base','swap23','swap34','omit2','omit3','omit4'};
C.nBins     = 100;          % 2 cm bins over 0-200 cm, same as Nav.Xspace
C.binEdges  = linspace(0, 200, C.nBins + 1);
C.binCentres = C.binEdges(1:end-1) + diff(C.binEdges)/2;

% Design blocks of the original GLM1, in order
C.blocks = {'XvisualONOFF','XrelativeSpd','XvisualTex','XvisualBG','XvisualEnd','XvisualOmit','Xspace'};
end
