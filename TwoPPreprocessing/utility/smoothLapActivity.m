function smoothed = smoothLapActivity(activity)
% Spatial smoothning
% smoothLapActivity Applies 6-point Gaussian smoothing to lap traces.
    w = gausswin(6); w = w / sum(w);
    smoothed = activity;
    for c = 1:size(activity, 1)
        for l = 1:size(activity, 2)
            trace = squeeze(activity(c, l, :));
            mask = isnan(trace);
            if all(mask), continue; end
            trace(mask) = 0;
            filt = filtfilt(w, 1, trace);
            filt(mask) = NaN;
            smoothed(c, l, :) = filt;
        end
    end
end