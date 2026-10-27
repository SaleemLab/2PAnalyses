%% checkGLMSpeedKernel_corrected.m
% Same as checkGLMSpeedKernel.m, but with validFit properly excluding
% all-zero placeholder kernel rows (cells never actually given a speed
% fit), on top of the all-NaN exclusion. Also drops the extreme-tail
% cells identified via diagnoseGLMSpeedOutliers.m pending your visual
% check of their raw kernels.
%
% Run after speedKernel, speedSE, speedPval, speedAxis, glmMod have
% already been computed (from checkGLMSpeedKernel.m).

%% Corrected valid-fit mask
isAllZero = false(size(speedKernel,1), 1);
for c = 1:size(speedKernel,1)
    isAllZero(c) = all(speedKernel(c,:) == 0);
end

validFit = ~all(isnan(speedKernel), 2) & ~isAllZero;

extremeThresh = -50; % from diagnoseGLMSpeedOutliers.m — adjust once you've
                      % confirmed these are genuine fit artifacts by eye
validFit = validFit & (glmMod >= extremeThresh | isnan(glmMod));

sigSpeed = validFit & speedPval <= 0.05;

fprintf('Valid-fit cells (excluding placeholders + extreme outliers): %d / %d\n', ...
        sum(validFit), size(speedKernel,1));
fprintf('Significant speed term among these: %d / %d (%.1f%%)\n', ...
        sum(sigSpeed), sum(validFit), sum(sigSpeed)/sum(validFit)*100);

%% Corrected summary stats
fprintf('\nAll valid-fit cells (corrected):\n');
fprintf('  n = %d\n', sum(validFit));
fprintf('  mean modulation   = %.4f\n', mean(glmMod(validFit), 'omitnan'));
fprintf('  median modulation = %.4f\n', median(glmMod(validFit), 'omitnan'));
fprintf('  low>high: %d (%.1f%%)\n', sum(glmMod(validFit) > 0), mean(glmMod(validFit) > 0)*100);
fprintf('  high>low: %d (%.1f%%)\n', sum(glmMod(validFit) < 0), mean(glmMod(validFit) < 0)*100);
fprintf('  exactly zero (should be ~0 now): %d\n', sum(glmMod(validFit) == 0));

fprintf('\nSignificant-speed-term cells only (corrected):\n');
fprintf('  n = %d\n', sum(sigSpeed));
fprintf('  mean modulation   = %.4f\n', mean(glmMod(sigSpeed), 'omitnan'));
fprintf('  median modulation = %.4f\n', median(glmMod(sigSpeed), 'omitnan'));
fprintf('  low>high: %d (%.1f%%)\n', sum(glmMod(sigSpeed) > 0), mean(glmMod(sigSpeed) > 0)*100);
fprintf('  high>low: %d (%.1f%%)\n', sum(glmMod(sigSpeed) < 0), mean(glmMod(sigSpeed) < 0)*100);

%% Corrected histogram
figure('Color', 'w', 'Position', [100 100 500 350]);
histogram(glmMod(validFit & ~sigSpeed), 40, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on
histogram(glmMod(sigSpeed & glmMod > 0), 40, 'FaceColor', [0.0 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
histogram(glmMod(sigSpeed & glmMod < 0), 40, 'FaceColor', [0.8 0.0 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
xline(0, '--k', 'LineWidth', 1.5);
xlabel('GLM speed modulation (25th - 75th percentile speed)');
ylabel('Number of cells');
legend('Non-sig', 'Sig, low>high', 'Sig, high>low', 'Box', 'off');
title(sprintf('GLM speed modulation, corrected (n=%d)', sum(validFit)));
set(gca, 'Box', 'off', 'TickDir', 'out');%% checkGLMSpeedKernel_corrected.m
% Same as checkGLMSpeedKernel.m, but with validFit properly excluding
% all-zero placeholder kernel rows (cells never actually given a speed
% fit), on top of the all-NaN exclusion. Also drops the extreme-tail
% cells identified via diagnoseGLMSpeedOutliers.m pending your visual
% check of their raw kernels.
%
% Run after speedKernel, speedSE, speedPval, speedAxis, glmMod have
% already been computed (from checkGLMSpeedKernel.m).

%% Corrected valid-fit mask
isAllZero = false(size(speedKernel,1), 1);
for c = 1:size(speedKernel,1)
    isAllZero(c) = all(speedKernel(c,:) == 0);
end

validFit = ~all(isnan(speedKernel), 2) & ~isAllZero;

extremeThresh = -50; % from diagnoseGLMSpeedOutliers.m — adjust once you've
                      % confirmed these are genuine fit artifacts by eye
validFit = validFit & (glmMod >= extremeThresh | isnan(glmMod));

sigSpeed = validFit & speedPval <= 0.05;

fprintf('Valid-fit cells (excluding placeholders + extreme outliers): %d / %d\n', ...
        sum(validFit), size(speedKernel,1));
fprintf('Significant speed term among these: %d / %d (%.1f%%)\n', ...
        sum(sigSpeed), sum(validFit), sum(sigSpeed)/sum(validFit)*100);

%% Corrected summary stats
fprintf('\nAll valid-fit cells (corrected):\n');
fprintf('  n = %d\n', sum(validFit));
fprintf('  mean modulation   = %.4f\n', mean(glmMod(validFit), 'omitnan'));
fprintf('  median modulation = %.4f\n', median(glmMod(validFit), 'omitnan'));
fprintf('  low>high: %d (%.1f%%)\n', sum(glmMod(validFit) > 0), mean(glmMod(validFit) > 0)*100);
fprintf('  high>low: %d (%.1f%%)\n', sum(glmMod(validFit) < 0), mean(glmMod(validFit) < 0)*100);
fprintf('  exactly zero (should be ~0 now): %d\n', sum(glmMod(validFit) == 0));

fprintf('\nSignificant-speed-term cells only (corrected):\n');
fprintf('  n = %d\n', sum(sigSpeed));
fprintf('  mean modulation   = %.4f\n', mean(glmMod(sigSpeed), 'omitnan'));
fprintf('  median modulation = %.4f\n', median(glmMod(sigSpeed), 'omitnan'));
fprintf('  low>high: %d (%.1f%%)\n', sum(glmMod(sigSpeed) > 0), mean(glmMod(sigSpeed) > 0)*100);
fprintf('  high>low: %d (%.1f%%)\n', sum(glmMod(sigSpeed) < 0), mean(glmMod(sigSpeed) < 0)*100);

%% Corrected histogram
figure('Color', 'w', 'Position', [100 100 500 350]);
histogram(glmMod(validFit & ~sigSpeed), 40, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on
histogram(glmMod(sigSpeed & glmMod > 0), 40, 'FaceColor', [0.0 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
histogram(glmMod(sigSpeed & glmMod < 0), 40, 'FaceColor', [0.8 0.0 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
xline(0, '--k', 'LineWidth', 1.5);
xlabel('GLM speed modulation (25th - 75th percentile speed)');
ylabel('Number of cells');
legend('Non-sig', 'Sig, low>high', 'Sig, high>low', 'Box', 'off');
title(sprintf('GLM speed modulation, corrected (n=%d)', sum(validFit)));
set(gca, 'Box', 'off', 'TickDir', 'out');