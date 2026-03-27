% Program to align two sync pulses from logs (or from log and npix)
%   In this function, the timebase in log1 (log1.Time)
%   is shifted to the timebase in log2 (log2.Time),
%   and the resultant timebase is appended to log1 as log1.realignedTime
%
%   If you would like to align log timebase to npix timebase, you should therefore
%   include the npix as the second input
% SGS 03/04/2022
function [outTimes2] = align2PSyncPulses(syncTimes1,syncTimes2,inTimes1)
minDetects = min(length(syncTimes1),length(syncTimes2));
syncTimes1 = syncTimes1(1:minDetects);
syncTimes2 = syncTimes2(1:minDetects);
% Find correlation between two
[r, lags] = xcorr(diff(syncTimes1), diff(syncTimes2),100,'unbiased');
[~, joint_idx] = max(r);
best_lag = lags(joint_idx);
% check which async pulse sequence is available first (ie. which was turned
% on first)
if best_lag < 0
    nSyncOffset = -best_lag+1;
    t2 = syncTimes2(nSyncOffset:end); % sync sglx times
    if length(t2) > length(syncTimes1)
        t2 = t2(1:length(syncTimes1));
    end
    t1 = syncTimes1(1:numel(t2)); %sync bonsai times
else
    nSyncOffset = best_lag+1;
    t2 = syncTimes2(1:end-nSyncOffset+1);
    if length(t2) > length(syncTimes1)
        t2 = t2(1:length(syncTimes1));
    end
    t1 = syncTimes1(nSyncOffset:end);
end
outTimes2 = interp1(t1, t2, inTimes1,'linear','extrap');