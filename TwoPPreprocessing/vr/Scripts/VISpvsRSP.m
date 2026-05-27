% plotting function after the VISp and RSP data have been loaded into the
% workspace 

% Plot smi cumulative probablity and histogram showing distribution 
compareAndPlot_SMI_RSP_vs_VISp(RSPData, VISpData)

% plot peak smis for boutons and somas 
plotSMI_PeakDistributions_VISP_RSP(RSPData, VISpData)

% finds preferred peak position to either gratings or plaid 
% plots bar chart with percentage of selective rois 
plotSMI_LandmarkIdentityPreference_Both(RSPData, VISpData)
