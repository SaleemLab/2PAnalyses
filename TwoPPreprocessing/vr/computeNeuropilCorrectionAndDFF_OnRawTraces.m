function [sessionFileInfo] = computeNeuropilCorrectionAndDFF_OnRawTraces(sessionFileInfo, subtractAbsosuteZero, applyTemporalSmoothing, prctl_F0, prctl_F, windowSize, smoothW, numN, minNp, maxNp)
% Computes dff using raw traces for all stimuli before data streams are
% interpolated. 

%% Set dafaults
if nargin < 2, subtractAbsosuteZero = true; end
if nargin < 3, applyTemporalSmoothing = false; end 
if nargin < 4 || isempty(prctl_F0), prctl_F0 = 8; end
if nargin < 5 || isempty(prctl_F), prctl_F = 5; end
if nargin < 6 || isempty(windowSize), windowSize = 60; end
if nargin < 7 || isempty(smoothW), smoothW = 15; end
if nargin < 8 || isempty(numN), numN = 20; end
if nargin < 9 || isempty(minNp), minNp = 10; end
if nargin < 10 || isempty(maxNp), maxNp = 90; end

for iStim = 1:length(sessionFileInfo.stimFiles)
    stimName = sessionFileInfo.stimFiles(iStim).name;
    stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
    if isempty(stimIdx), error('Specified stimName not found.'); end
    
    try
        fprintf('Performing Neuropil Correction and computing Delta FF for %s\n', stimName);
        % PMT Offset logic
        % This value represents the absolute zero signal and was obtained by averaging the darkest frame over many imaging sessions.
        % It is important to note that the absolute zero value is arbitrary and depends on the voltage range of the PMTs.
        % Below is the absolute zero value for our B-scope PMTs; Sylvia's scope
        % absolute zero value is 19520 [and thats probably why she implements this].
        absZero = -23;

        %% Load raw F,Neu and ops
        disp('Loading F, FNeu, ops...');
        TwoPDataStruct = load(sessionFileInfo.stimFiles(stimIdx).mergedBonsai2PSuite2pData);
        twoPData = TwoPDataStruct.twoPData;
        for thisPlane=1:length(twoPData)
             F = twoPData(thisPlane).F;
             Fneu = twoPData(thisPlane).Fneu;
             planeRate = twoPData(thisPlane).ops.fs; % [not using the volume rate which would be 60Hz or 30Hz] 

            %% Absolute zero subtraction [a copy of this is not saved back in F, FNeu FYI]
            if subtractAbsosuteZero
               disp('Subtracting absolute zero from F and Fneu traces')
               F = F - absZero;
               Fneu = Fneu - absZero;
            end 

            %% Optional smoothning [usually skipped]
            if applyTemporalSmoothing
                fprintf('Applying smoothing (gausswin %d)...\n', smoothW);
                w = gausswin(smoothW); 
                w = w / sum(w);
                fSmoothed = filtfilt(w, 1, F')';
                fneuSmoothed = filtfilt(w, 1, Fneu')';
            else
                fSmoothed = F; fneuSmoothed = Fneu;
            end
            
            %% Baseline Calculation & Centering
            disp('Calculating slow-drift baselines...');
            % We calculate these once and reuse them for dF/F calculation to save time
            % get_F0 expects [frames x ROIs] and transpose back to [ROIs x
            % frames] 
            f0_F = get_F0(fSmoothed', prctl_F0, windowSize, planeRate)';
            % f0_N = get_F0(fneuSmoothed', prctl_F0, windowSize, fs)';
            %% Neuropil correction 
            
            % Computes F0_F and F0_Neu internally, subtracts these traces
            % from F and FNeu respectively [i.e., centers these traces]
            % before using polyfit to estimate the correction factor [r].
            % F0 is added back to the corrected signal [called Fc] 
                    

            % % Centering signals for estimateNeuropil regression
            % % Sylvia's pipline subtracts F0 before estimating neuropil: 
            % % This removes slow drifts and aligns 
            % % both signals at a zero-baseline. It might be important for neuropil correction to 
            % % ensure the correction factor (r) is estimated based on common 
            % % high-frequency fluctuations rather than divergent absolute offsets.
            % fCentered = fSmoothed - f0_F;
            % fneuCentered = fneuSmoothed - f0_N;
            %disp('Computing neuropil correction using estimateNeuropil..code from +preproc CortexLab');
            % opt.numN = numN;
            % opt.minNp = minNp;
            % opt.maxNp = maxNp;
            % opt.pCell = prctl_F;
            % opt.noNeg = 1;
            % opt.constrainedFit = 0;
            % opt.window = Inf;
            % 
            % % estimateNeuropil output (signalTrace) is centered at 0
            % [signalTraceCentered, neuropCorrPars] = estimateNeuropil(fCentered, fneuCentered, opt);
            % 
            % % RESTORE BASELINE: Add f0_F back to the corrected centered signal
            % % This is the "Corrected Raw Fluorescence" (Fc) required for dF/F
            % % From Sylvia's depth from 2p @ signal[:, iROI] = iF - (b * iN + a) +
            % % F0[:, iROI];
            % Fc = signalTraceCentered + f0_F;
            % 
            



            %% Delta F/F Calculation
            disp('Computing delta F/F signals...');

            % Raw dF/F (using  F and the F0 we already computed)
            dFF = get_delta_F_over_F(fSmoothed', f0_F')';

            % Neuropil Corrected dF/F
            % Per Sylvia: Normalise corrected signal (Fc) by the original Raw F0 (f0_F)
            dFFNeuropilCorrected = get_delta_F_over_F(Fc', f0_F')';
            
            
            twoPData(thisPlane).dFF = dFF; 
            twoPData(thisPlane).dFFNeuropilCorrected = dFFNeuropilCorrected;

            %% Z-Scoring
            % if zScoreProcessedSignals
            %     disp('Z-scoring signals...');
            %     twoPData(thisPlane).zScored_dFFNeuropilCorrected = zscore(dFFNeuropilCorrected, 0, 2);
            %     twoPData(thisPlane).zScored_dFF = zscore(dFF, 0, 2);
            % end
         
        end
        
        %% Saving 
        disp('Saving to 2p data...');
        save(sessionFileInfo.stimFiles(stimIdx).mergedBonsai2PSuite2pData, "twoPData", '-append');

    catch ME
        warning('Failed to process stimulus "%s". Error: %s', stimName, ME.message);
        sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData = [];
    end
end
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end