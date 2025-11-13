function plot_neuropil_QC(F, N, signal, regPars, F_binValues, N_binValues, iROI)

figure;

% --- Scatter plot + bins + fit
subplot(2,1,1);
scatter(N(:,iROI), F(:,iROI), 5, 'k', 'filled'); hold on;
plot(N_binValues(:,iROI), F_binValues(:,iROI), 'ro-', 'LineWidth', 2);
a = regPars(1,iROI); b = regPars(2,iROI);
xline = linspace(min(N(:,iROI)), max(N(:,iROI)), 100);
yfit = a + b*xline;
plot(xline, yfit, 'b-', 'LineWidth', 2);

xlabel('Neuropil (N)'); ylabel('ROI (F)');
legend('All points','5th percentile/bin','Linear fit');
title(sprintf('ROI %d: slope = %.2f', iROI, b));

% --- Time traces
subplot(2,1,2);
plot(F(:,iROI), 'k-', 'DisplayName','Raw ROI');
hold on;
plot(signal(:,iROI), 'r-', 'LineWidth',1.5,'DisplayName','Corrected ROI');
plot(N(:,iROI), 'b--', 'DisplayName','Neuropil');
xlabel('Frame'); ylabel('Fluorescence');
legend; title('Time traces');
end
