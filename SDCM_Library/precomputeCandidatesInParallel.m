%ABSTRACT
% Library function for SDCM. Parallel pre-computation of the next candidates for 
% initial representatives of signatures (in step 1) or for additional representatives
% towards correlation maximization and convergence (in step 2). 
% Precomputes all signature axes, correlations of all genes and sample to them, 
% the signature focus (i.e. weights), etc. Stores all results in the global 
% eDState structure. They are collected from there by getSignatureAxesAndScores.

function precomputeCandidatesInParallel(ImJ, inInfo, BsPerformanceSubspace4Correlations)
  %% Initialize:
    global eDState;
    %Handle clearFromCache command:
      if(ischar(inInfo) && strcmp(inInfo, 'clearFromCache'))
        assert(all(ismember(ImJ, eDState.current.precompute.ImJ)), 'code validation: cannot clear ImJ that are not part of the current cache');
        BDelete4ImJ = ismember(eDState.current.precompute.ImJ, ImJ);
        BDelete4I = ismember(eDState.current.precompute.I, ImJ);
        BDelete4J = ismember(eDState.current.precompute.J, -ImJ);
        for fn=fieldnames(eDState.current.precompute)'; fn=fn{1};
          switch(fn)
            %Combined fields for genes and samples:
              case {'ImJ','BAlreadyTestedAndDisqualifiedInStep1','candidateAccus','BcandidateAccusUpToDate'}
                eDState.current.precompute.(fn)(BDelete4ImJ) = [];
            %Reset temporary fields:
              case {'ii','jj'}
                eDState.current.precompute.(fn) = []; %clear locally set fieldnames by getSignatureAxesAndScores
            %Scalars:
              case {'I','norm4twinGeneAxes4I_sdUnits','norm4twinGeneAxes4I_origUnits','norm4sampleAxes4I_sdUnits','norm4sampleAxes4I_origUnits','signatureAbsMean2D4I','sampleSize4signatureAbsMean2D4I','signatureAbsSD2D4I','signatureCorrInExtendedFocus4I','log10_p4I','log10_p4Correlations4I','log10_p4SignalStrength4I','signatureSizeByCorrSum2D4I','signatureSizeByCorrSum4G4I','signatureSizeByCorrSum4P4I'}
                eDState.current.precompute.(fn)(BDelete4I) = [];
              case {'J','norm4geneAxes4J_sdUnits','norm4geneAxes4J_origUnits','norm4twinSampleAxes4J_sdUnits','norm4twinSampleAxes4J_origUnits','signatureAbsMean2D4J','sampleSize4signatureAbsMean2D4J','signatureAbsSD2D4J','signatureCorrInExtendedFocus4J','log10_p4J','log10_p4Correlations4J','log10_p4SignalStrength4J','signatureSizeByCorrSum2D4J','signatureSizeByCorrSum4G4J','signatureSizeByCorrSum4P4J'}
                eDState.current.precompute.(fn)(BDelete4J) = [];
            %Column vectors:
              case {'twinGeneAxes4I_sdUnits','twinGeneAxes4I_origUnits','W4sampleAxes4I','R4Gs4I','sampleSizes4R4Gs4I','P4R4Gs4I','signedExtendedW4Gs4Itwin','signedFocusedW4Gs4I','signedFocusedW4Gs4I_withPerpendicularSpace'}
                eDState.current.precompute.(fn)(:,BDelete4I) = [];
              case {'geneAxes4J_sdUnits','geneAxes4J_origUnits','W4twinSampleAxes4J','R4Gs4Jtwin','sampleSizes4R4Gs4Jtwin','P4R4Gs4J','signedExtendedW4Gs4J','signedFocusedW4Gs4J','signedFocusedW4Gs4J_withPerpendicularSpace'}
                eDState.current.precompute.(fn)(:,BDelete4J) = [];
            %Row vectors:
              case {'sampleAxes4I_sdUnits','sampleAxes4I_origUnits','W4twinGeneAxes4I','R4Ps4Itwin','sampleSizes4R4Ps4Itwin','P4R4Ps4I','signedExtendedW4Ps4I','signedFocusedW4Ps4I','signedFocusedW4Ps4I_withPerpendicularSpace'}
                eDState.current.precompute.(fn)(BDelete4I,:) = [];
              case {'twinSampleAxes4J_sdUnits','twinSampleAxes4J_origUnits','W4geneAxes4J','R4Ps4J','sampleSizes4R4Ps4J','P4R4Ps4J','signedExtendedW4Ps4Jtwin','signedFocusedW4Ps4J','signedFocusedW4Ps4J_withPerpendicularSpace'}
                eDState.current.precompute.(fn)(BDelete4J,:) = [];
            otherwise
              error('code validation: fieldname %s unhandled', fn);
          end
        end
        %Assert consistency:
          assert(length(eDState.current.precompute.ImJ) == length(eDState.current.precompute.I)+length(eDState.current.precompute.J), 'code validation: cache is inconsistent; check elements in .ImJ and compare with .I and .J');
          for fn=fieldnames(eDState.current.precompute)'; fn=fn{1};
            if(length(fn)>2 && strcmp(fn(end-1:end),'4I'))
              assert(...
                 any(size(eDState.current.precompute.(fn)) == length(eDState.current.precompute.I)) ...
                ,'code validation: cache is inconsistent; .precompute.%s has size %s, but .precompute.I has length %d', fn, mat2str(size(eDState.current.precompute.(fn))), length(eDState.current.precompute.I) ...
              );
            elseif(length(fn)>2 && strcmp(fn(end-1:end),'4J'))
              assert(...
                 any(size(eDState.current.precompute.(fn)) == length(eDState.current.precompute.J)) ...
                ,'code validation: cache is inconsistent; .precompute.%s has size %s, but .precompute.I has length %d', fn, mat2str(size(eDState.current.precompute.(fn))), length(eDState.current.precompute.I) ...
              );
            end
          end
        return;
      end
    %reset on empty call:
      if(isempty(ImJ))
        eDState.current.precompute = struct();
      end
    %Shortcuts:
      nG = size(eDState.current.L2Rs,1);
      nP = size(eDState.current.L2Rs,2);            
    %Performance subspace logic:
      bRestrictToPerformanceSubspace = nargin>=3;
      if(~bRestrictToPerformanceSubspace)
        BsPerformanceSubspace4Correlations = {true(nG,1),true(1,nP)}; 
      end
      bPerformanceSubspaceTooSmall = sum(BsPerformanceSubspace4Correlations{1})<2 || sum(BsPerformanceSubspace4Correlations{2})<2;
      if(bPerformanceSubspaceTooSmall)
        warning('only %d/%d genes and %d/%d samples were in the performance subspace; this is not enough for computation; now increasing to the full space; consider less tight .searchStrategy.performance.alpha4Signal4inclusionInComputation4G/P respectively .focusing.performance.minAbsR4inclusionInCorrComputation configuration'...
          ,sum(BsPerformanceSubspace4Correlations{1}),nG, sum(BsPerformanceSubspace4Correlations{2}),nP ...
        );
        BsPerformanceSubspace4Correlations = {true(nG,1),true(1,nP)}; 
      end
  %% Calc in parallel or single threaded, depending on ImJ:
    precomputedChunks = {};
    %Cache logic:
      if(isempty(ImJ))
        eDState.current.precompute.ImJ = [];
        BNew = false(1,0);
      else
        BNew = ~ismember(ImJ, eDState.current.precompute.ImJ);
          assert(all(BNew) || sum(~BNew)==length(eDState.current.precompute.ImJ) && all(ImJ(~BNew)==eDState.current.precompute.ImJ), 'code validation: when precomputing new ImJ while keeping previous ones, ImJ(~BNew) must be identical to eDState.current.precompute.ImJ');
      end
      newImJToPrecomputeNow = ImJ(BNew);
      newImJToPrecomputeNow = sort(newImJToPrecomputeNow); %performance: make the chunks of genes respectively samples as large as possible per worker (with the number of signatures of the other dimension being zero in the optimum case)
    %To minimize the data amount sent to workers, make local copies of only the needed state data:
      lc_currentData = struct();
      lc_currentData.current.L2Rs = eDState.current.L2Rs;
      lc_currentData.current.sdL2Rs = eDState.current.sdL2Rs;
      %Select the base signal for the search according to the preprocessing configuration:
        lc_currentData.current.multiPassesLevel = eDState.current.multiPassesLevel;
        sL = min(lc_currentData.current.multiPassesLevel, inInfo.searchStrategy.nDefinedPasses);
        lc_currentData.current.sdL2Rs = eDState.current.sdL2Rs;
        lc_currentData.current.L2Rs = eDState.current.L2Rs;
      %Pass on the current noise estimation:
        if(eDState.current.k>1)
          lc_currentData.noiseEstimation.halfNormalAbsNoiseMedian = eDState.noiseEstimation.halfNormalAbsNoiseMedian(eDState.current.k);
          lc_currentData.noiseEstimation.preciseNoiseSD = eDState.noiseEstimation.preciseNoiseSD(eDState.current.k);
          lc_currentData.noiseEstimation.absMean2D = eDState.noiseEstimation.absMean2D(eDState.current.k);
          lc_currentData.noiseEstimation.nsBehindAbsMean2D = eDState.noiseEstimation.nsBehindAbsMean2D(eDState.current.k);
          lc_currentData.noiseEstimation.SD4absMean2D = eDState.noiseEstimation.SD4absMean2D(eDState.current.k);
        end
        lc_currentData.noiseEstimation.PSignalStrength4G = 10.^eDState.noiseEstimation.log10PNoise4G;
        lc_currentData.noiseEstimation.PSignalStrength4P = 10.^eDState.noiseEstimation.log10PNoise4P;
    %Decide on parallelization and precompute locally or on workers:
      nWorkers = ilv(gcp('nocreate'),@(pool)iif(isempty(pool),0,pool.NumWorkers)); %matlabpool('size');
      nParallelizationThreshold = max(1,nWorkers)/2;
      bSendToWorkers = length(newImJToPrecomputeNow) > nParallelizationThreshold && nWorkers>0;
      if(~bSendToWorkers)
        precomputedChunks{1} = precomputeNextSignatureCandidatesInParallel_chunk(...
           newImJToPrecomputeNow, inInfo, BsPerformanceSubspace4Correlations, lc_currentData ...
        );
      else %parallel:
        nChunks = max(1, nWorkers);
        nChunkSize = ceil(length(newImJToPrecomputeNow)/nChunks);
        chunkedImJ = cell(1,nChunks);
          for c=1:nChunks
            chunkedImJ{c} = newImJToPrecomputeNow((c-1)*nChunkSize+1 : min(end,c*nChunkSize));
          end
        bFinished = false;
        if(~bFinished && nWorkers>0)
          for retryRuns = 1:2;
            try
              spmd
                precomputedChunks_composite = precomputeNextSignatureCandidatesInParallel_chunk(...
                   chunkedImJ{labindex}, inInfo, BsPerformanceSubspace4Correlations, lc_currentData ...
                );
              end
              bFinished = true;
              break;
            catch ex
              warning('ERROR during parallel precomputation of signature candidates => force-restarting the MATLAB pool and retrying; error details: %s; use >>dbcont to retry', ex.message);
              keyboard;
              %Restart pool:
                try
                  delete(gcp);
                catch
                end
                parpool(nWorkers);
            end
          end
            if(~bFinished)
              warning('ERROR: repeated parallelization errors occurred; use >>dbcont to retry single-threaded or inspect.');
              keyboard;
            end
        end
        %retry single-threaded on error:
          if(~bFinished) 
            for l = 1:nChunks 
              precomputedChunks_composite{l} = precomputeNextSignatureCandidatesInParallel_chunk(...
                 chunkedImJ{l}, inInfo, BsPerformanceSubspace4Correlations, lc_currentData ...
              );
            end
            bFinished = true;
          end
            if(isstruct(precomputedChunks_composite)) precomputedChunks_composite = {precomputedChunks_composite}; end 
        precomputedChunks = cell(1,length(precomputedChunks_composite));
          for c=1:length(precomputedChunks_composite)
            precomputedChunks{c} = precomputedChunks_composite{c};
          end
      end
    %Append all precomputed candidates to the cache:
      if(isempty(ImJ))
        eDState.current.precompute = precomputedChunks{1};
      end
      for c=1:length(precomputedChunks)
        for fn=fieldnames(eDState.current.precompute)'; fn=fn{1};
          switch(fn)
            %Combined fields for genes and samples:
              case {'ImJ','BAlreadyTestedAndDisqualifiedInStep1'}
                eDState.current.precompute.(fn) = [eDState.current.precompute.(fn), precomputedChunks{c}.(fn)];
              case {'candidateAccus','BcandidateAccusUpToDate'}
                %let step 2 fill in the candidateAccus.
            %Reset temporary fields:
              case {'ii','jj'}
                eDState.current.precompute.(fn) = []; %clear locally set fieldnames by getSignatureAxesAndScores
            %Scalars:
              case {'I','norm4twinGeneAxes4I_sdUnits','norm4twinGeneAxes4I_origUnits','norm4sampleAxes4I_sdUnits','norm4sampleAxes4I_origUnits','signatureAbsMean2D4I','sampleSize4signatureAbsMean2D4I','signatureAbsSD2D4I','signatureCorrInExtendedFocus4I','log10_p4I','log10_p4Correlations4I','log10_p4SignalStrength4I','signatureSizeByCorrSum2D4I','signatureSizeByCorrSum4G4I','signatureSizeByCorrSum4P4I'}
                eDState.current.precompute.(fn) = [eDState.current.precompute.(fn), precomputedChunks{c}.(fn)];
              case {'J','norm4geneAxes4J_sdUnits','norm4geneAxes4J_origUnits','norm4twinSampleAxes4J_sdUnits','norm4twinSampleAxes4J_origUnits','signatureAbsMean2D4J','sampleSize4signatureAbsMean2D4J','signatureAbsSD2D4J','signatureCorrInExtendedFocus4J','log10_p4J','log10_p4Correlations4J','log10_p4SignalStrength4J','signatureSizeByCorrSum2D4J','signatureSizeByCorrSum4G4J','signatureSizeByCorrSum4P4J'}
                eDState.current.precompute.(fn) = [eDState.current.precompute.(fn), precomputedChunks{c}.(fn)];
            %Column vectors:
              case {'twinGeneAxes4I_sdUnits','twinGeneAxes4I_origUnits','W4sampleAxes4I','R4Gs4I','sampleSizes4R4Gs4I','P4R4Gs4I','signedExtendedW4Gs4Itwin','signedFocusedW4Gs4I','signedFocusedW4Gs4I_withPerpendicularSpace'}
                eDState.current.precompute.(fn) = [eDState.current.precompute.(fn), precomputedChunks{c}.(fn)];
              case {'geneAxes4J_sdUnits','geneAxes4J_origUnits','W4twinSampleAxes4J','R4Gs4Jtwin','sampleSizes4R4Gs4Jtwin','P4R4Gs4J','signedExtendedW4Gs4J','signedFocusedW4Gs4J','signedFocusedW4Gs4J_withPerpendicularSpace'}
                eDState.current.precompute.(fn) = [eDState.current.precompute.(fn), precomputedChunks{c}.(fn)];
            %Row vectors:
              case {'sampleAxes4I_sdUnits','sampleAxes4I_origUnits','W4twinGeneAxes4I','R4Ps4Itwin','sampleSizes4R4Ps4Itwin','P4R4Ps4I','signedExtendedW4Ps4I','signedFocusedW4Ps4I','signedFocusedW4Ps4I_withPerpendicularSpace'}
                eDState.current.precompute.(fn) = [eDState.current.precompute.(fn); precomputedChunks{c}.(fn)];
              case {'twinSampleAxes4J_sdUnits','twinSampleAxes4J_origUnits','W4geneAxes4J','R4Ps4J','sampleSizes4R4Ps4J','P4R4Ps4J','signedExtendedW4Ps4Jtwin','signedFocusedW4Ps4J','signedFocusedW4Ps4J_withPerpendicularSpace'}
                eDState.current.precompute.(fn) = [eDState.current.precompute.(fn); precomputedChunks{c}.(fn)];
            otherwise
              error('code validation: fieldname %s unhandled', fn);
          end
        end
        %Assert consistency:
          assert(length(eDState.current.precompute.ImJ) == length(eDState.current.precompute.I)+length(eDState.current.precompute.J), 'code validation: cache is inconsistent; check elements in .ImJ and compare with .I and .J');
          for fn=fieldnames(eDState.current.precompute)'; fn=fn{1};
            if(length(fn)>2 && strcmp(fn(end-1:end),'4I'))
              assert(...
                 any(size(eDState.current.precompute.(fn)) == length(eDState.current.precompute.I)) ...
                ,'code validation: cache is inconsistent; .precompute.%s has size %s, but .precompute.I has length %d', fn, mat2str(size(eDState.current.precompute.(fn))), length(eDState.current.precompute.I) ...
              );
            elseif(length(fn)>2 && strcmp(fn(end-1:end),'4J'))
              assert(...
                 any(size(eDState.current.precompute.(fn)) == length(eDState.current.precompute.J)) ...
                ,'code validation: cache is inconsistent; .precompute.%s has size %s, but .precompute.J has length %d', fn, mat2str(size(eDState.current.precompute.(fn))), length(eDState.current.precompute.J) ...
              );
            end
          end
      end
      if(inInfo.internal.bDevEnableInteractiveBreaks)
        interactiveCheckPoint();
      end
end

  function precompute = precomputeNextSignatureCandidatesInParallel_chunk(ImJ, inInfo, BsPerformanceSubspace4Correlations, lc_currentData)
    %% Initialize:
      precompute = struct();
      precompute.ImJ = ImJ;
      precompute.BAlreadyTestedAndDisqualifiedInStep1 = false(size(ImJ));
      bStatusOutput = ~isempty(ImJ);
      if(inInfo.export.nStatusOutputLevel<=2 && exist('labindex')==5)
        bStatusOutput = bStatusOutput && labindex==1; %for .nStatusOutputLevel<=2 only output status messages from the first worker to not pollute the log.
        if(labindex>1)
          SDCM_printStatus(5,'     (-)Not reporting status messaged from worker %d, as they are the same as reported from worker (labindex) 1.',labindex);
        end
      end
      %Shortcuts:
        nG = size(lc_currentData.current.L2Rs,1);
        nP = size(lc_currentData.current.L2Rs,2);            
        sL = min(lc_currentData.current.multiPassesLevel, inInfo.searchStrategy.nDefinedPasses);
        numericTargetPrecision = class(lc_currentData.current.L2Rs);
      %Select base signal for signature search:
        fn4L2Rs_origUnits = 'L2Rs';
        fn4L2Rs_fullG_origUnits = 'L2Rs_fullG';
        fn4L2Rs_fullP_origUnits = 'L2Rs_fullP';
        fn4L2Rs_sdUnits = 'sdL2Rs';
        fn4L2Rs_fullG_sdUnits = 'sdL2Rs_fullG';
        fn4L2Rs_fullP_sdUnits = 'sdL2Rs_fullP';
      %Performance subspace logic:
        bRestrictToPerformanceSubspace = ~all(BsPerformanceSubspace4Correlations{1}) || ~all(BsPerformanceSubspace4Correlations{2});
        bOnlyOneSidedPerfSubspace4Corr = true; 
        reduced_currentData.L2Rs = lc_currentData.current.(fn4L2Rs_origUnits);
        reduced_currentData.sdL2Rs = lc_currentData.current.(fn4L2Rs_sdUnits);
          if(bRestrictToPerformanceSubspace)
            reduced_currentData.L2Rs = reduced_currentData.L2Rs(BsPerformanceSubspace4Correlations{1}, BsPerformanceSubspace4Correlations{2});
            reduced_currentData.sdL2Rs = reduced_currentData.sdL2Rs(BsPerformanceSubspace4Correlations{1}, BsPerformanceSubspace4Correlations{2});
          end
        reduced_currentData.(fn4L2Rs_fullG_origUnits) = lc_currentData.current.(fn4L2Rs_origUnits);
        reduced_currentData.(fn4L2Rs_fullG_sdUnits) = lc_currentData.current.(fn4L2Rs_sdUnits);
          if(bRestrictToPerformanceSubspace)
            reduced_currentData.(fn4L2Rs_fullG_origUnits) = lc_currentData.current.(fn4L2Rs_origUnits)(:, BsPerformanceSubspace4Correlations{2});
            reduced_currentData.(fn4L2Rs_fullG_sdUnits) = lc_currentData.current.(fn4L2Rs_sdUnits)(:, BsPerformanceSubspace4Correlations{2});
          end
        reduced_currentData.(fn4L2Rs_fullP_origUnits) = lc_currentData.current.(fn4L2Rs_origUnits);
        reduced_currentData.(fn4L2Rs_fullP_sdUnits) = lc_currentData.current.(fn4L2Rs_sdUnits);
          if(bRestrictToPerformanceSubspace)
            reduced_currentData.(fn4L2Rs_fullP_origUnits) = lc_currentData.current.(fn4L2Rs_origUnits)(BsPerformanceSubspace4Correlations{1}, :);
            reduced_currentData.(fn4L2Rs_fullP_sdUnits) = lc_currentData.current.(fn4L2Rs_sdUnits)(BsPerformanceSubspace4Correlations{1}, :);
          end
      %Get basic math function handles according to bDataContainsNaNs:
        BM = getBasicMathFunctions(inInfo.preprocessing.bDataContainsNaNs);

    %% Signatures directly based on first representatives (i.e. based on single genes or samples in the signal):
      if(true)
        %% Primary axes (i.e. first representatives, i.e. single genes/samples form the current signal):
          %Next gene initial representative candidates:
            precompute.I = precompute.ImJ(precompute.ImJ>0);
            precompute.sampleAxes4I_origUnits = lc_currentData.current.(fn4L2Rs_origUnits)(precompute.I,:);
            precompute.sampleAxes4I_sdUnits = lc_currentData.current.(fn4L2Rs_sdUnits)(precompute.I,:);
            %Compute vector norms:
              precompute.norm4sampleAxes4I_origUnits = BM.euclidW(precompute.sampleAxes4I_origUnits, 1, 2)';
              precompute.norm4sampleAxes4I_sdUnits = BM.euclidW(precompute.sampleAxes4I_sdUnits, 1, 2)';
          %Next sample initial representative candidates:
            precompute.J = -precompute.ImJ(precompute.ImJ<0);
            precompute.geneAxes4J_origUnits = lc_currentData.current.(fn4L2Rs_origUnits)(:,precompute.J);
            precompute.geneAxes4J_sdUnits = lc_currentData.current.(fn4L2Rs_sdUnits)(:,precompute.J);
            %Compute vector norms:
              precompute.norm4geneAxes4J_origUnits = BM.euclidW(precompute.geneAxes4J_origUnits, 1, 1);
              precompute.norm4geneAxes4J_sdUnits = BM.euclidW(precompute.geneAxes4J_sdUnits, 1, 1);          

        %% Initial signature focus for primary axes:
          precompute.signedFocusedW4Ps4I = calcSignatureFocus(precompute.sampleAxes4I_sdUnits, lc_currentData.noiseEstimation.PSignalStrength4P, [], [], 2, inInfo.searchStrategy.signatureFocus, sL);
          precompute.signedFocusedW4Gs4J = calcSignatureFocus(precompute.geneAxes4J_sdUnits,   lc_currentData.noiseEstimation.PSignalStrength4G, [], [], 1, inInfo.searchStrategy.signatureFocus, sL);
          %Add .full2focusWeightRatio weights from the current and the twin space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the twin space:
            precompute.signedFocusedW4Ps4I_withPerpendicularSpace = addFlatWeights(nan(nG,0), precompute.signedFocusedW4Ps4I, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
            precompute.signedFocusedW4Gs4J_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4J, nan(0,nP), 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
          %<-Note: For specific axes, it is important to restrict to only the points in signature. Just using signal strengths is only a very rough initial signature focus that needs to be refined below with correlations.
        %% Twin axes based on the initial signature focus:
          if(bStatusOutput) SDCM_printStatus(5,'      - aggregating twin axes...\n'); end
          %Twin gene axes for .sampleAxes4I:
            if(true)
              sampleAxes4I_origUnits = iif(bRestrictToPerformanceSubspace, precompute.sampleAxes4I_origUnits(:,BsPerformanceSubspace4Correlations{2}), precompute.sampleAxes4I_origUnits);
              sampleAxes4I_sdUnits = iif(bRestrictToPerformanceSubspace, precompute.sampleAxes4I_sdUnits(:,BsPerformanceSubspace4Correlations{2}), precompute.sampleAxes4I_sdUnits);
              signedFocusedW4Ps4I = iif(bRestrictToPerformanceSubspace, precompute.signedFocusedW4Ps4I(:,BsPerformanceSubspace4Correlations{2}), precompute.signedFocusedW4Ps4I);
              precompute.twinGeneAxes4I_origUnits = zeros(nG,size(signedFocusedW4Ps4I,1));
              precompute.twinGeneAxes4I_sdUnits = zeros(nG,size(signedFocusedW4Ps4I,1));
              for ii=1:size(signedFocusedW4Ps4I,1)
                precompute.twinGeneAxes4I_origUnits(BsPerformanceSubspace4Correlations{1},ii) = projectW(...
                  reduced_currentData.L2Rs ... vectorsInSourceSpace...
                 ,sampleAxes4I_origUnits(ii,:) ... axesInSourceSpace
                 ,abs(signedFocusedW4Ps4I(ii,:)) ... weights4sourceDims
                 ,2 ... sourceSpaceMatrixDim
                ,BM.meanW, BM.euclidW) / BM.euclidW(ones(1,size(signedFocusedW4Ps4I,2)),abs(signedFocusedW4Ps4I(ii,:)),2); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
                precompute.twinGeneAxes4I_sdUnits(BsPerformanceSubspace4Correlations{1},ii) = projectW(...
                  reduced_currentData.sdL2Rs ... vectorsInSourceSpace...
                 ,sampleAxes4I_sdUnits(ii,:) ... axesInSourceSpace
                 ,abs(signedFocusedW4Ps4I(ii,:)) ... weights4sourceDims
                 ,2 ... sourceSpaceMatrixDim
                ,BM.meanW, BM.euclidW) / BM.euclidW(ones(1,size(signedFocusedW4Ps4I,2)),abs(signedFocusedW4Ps4I(ii,:)),2); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
              end
                precompute.twinGeneAxes4I_origUnits(isnan(precompute.twinGeneAxes4I_origUnits)) = 0; %only needed for very high NaN rate (>~90%)
                precompute.twinGeneAxes4I_sdUnits(isnan(precompute.twinGeneAxes4I_sdUnits)) = 0; %only needed for very high NaN rate (>~90%)
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.twinGeneAxes4I_origUnits = sign(precompute.twinGeneAxes4I_origUnits).*dampenOutliers(abs(precompute.twinGeneAxes4I_origUnits), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G(min(end,sL)), 1); 
                  precompute.twinGeneAxes4I_sdUnits = sign(precompute.twinGeneAxes4I_sdUnits).*dampenOutliers(abs(precompute.twinGeneAxes4I_sdUnits), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G(min(end,sL)), 1); 
                end
              %Compute vector norms:
                precompute.norm4twinGeneAxes4I_origUnits = BM.euclidW(precompute.twinGeneAxes4I_origUnits, 1, 1);          
                precompute.norm4twinGeneAxes4I_sdUnits = BM.euclidW(precompute.twinGeneAxes4I_sdUnits, 1, 1);          
            end
          %Twin sample axes for .geneAxes4J:
            if(true)
              geneAxes4J_origUnits = iif(bRestrictToPerformanceSubspace, precompute.geneAxes4J_origUnits(BsPerformanceSubspace4Correlations{1},:), precompute.geneAxes4J_origUnits);
              geneAxes4J_sdUnits = iif(bRestrictToPerformanceSubspace, precompute.geneAxes4J_sdUnits(BsPerformanceSubspace4Correlations{1},:), precompute.geneAxes4J_sdUnits);
              signedFocusedW4Gs4J = iif(bRestrictToPerformanceSubspace, precompute.signedFocusedW4Gs4J(BsPerformanceSubspace4Correlations{1},:), precompute.signedFocusedW4Gs4J);
              precompute.twinSampleAxes4J_sdUnits = zeros(size(signedFocusedW4Gs4J,2),nP);
              precompute.twinSampleAxes4J_origUnits = zeros(size(signedFocusedW4Gs4J,2),nP);
              for jj=1:size(signedFocusedW4Gs4J,2)
                precompute.twinSampleAxes4J_origUnits(jj,BsPerformanceSubspace4Correlations{2}) = projectW(...
                  reduced_currentData.L2Rs ... vectorsInSourceSpace...
                 ,geneAxes4J_origUnits(:,jj) ... axesInSourceSpace
                 ,abs(signedFocusedW4Gs4J(:,jj)) ... weights4sourceDims
                 ,1 ... sourceSpaceMatrixDim
                ,BM.meanW, BM.euclidW) / BM.euclidW(ones(size(signedFocusedW4Gs4J,1),1),abs(signedFocusedW4Gs4J(:,jj)),1); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
                precompute.twinSampleAxes4J_sdUnits(jj,BsPerformanceSubspace4Correlations{2}) = projectW(...
                  reduced_currentData.sdL2Rs ... vectorsInSourceSpace...
                 ,geneAxes4J_sdUnits(:,jj) ... axesInSourceSpace
                 ,abs(signedFocusedW4Gs4J(:,jj)) ... weights4sourceDims
                 ,1 ... sourceSpaceMatrixDim
                ,BM.meanW, BM.euclidW) / BM.euclidW(ones(size(signedFocusedW4Gs4J,1),1),abs(signedFocusedW4Gs4J(:,jj)),1); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
              end
                precompute.twinSampleAxes4J_origUnits(isnan(precompute.twinSampleAxes4J_origUnits)) = 0; %only needed for very high NaN rate (>~90%)
                precompute.twinSampleAxes4J_sdUnits(isnan(precompute.twinSampleAxes4J_sdUnits)) = 0; %only needed for very high NaN rate (>~90%)
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.twinSampleAxes4J_origUnits = sign(precompute.twinSampleAxes4J_origUnits).*dampenOutliers(abs(precompute.twinSampleAxes4J_origUnits), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                  precompute.twinSampleAxes4J_sdUnits = sign(precompute.twinSampleAxes4J_sdUnits).*dampenOutliers(abs(precompute.twinSampleAxes4J_sdUnits), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                end
              %Compute vector norms:
                precompute.norm4twinSampleAxes4J_origUnits = BM.euclidW(precompute.twinSampleAxes4J_origUnits, 1, 2)';
                precompute.norm4twinSampleAxes4J_sdUnits = BM.euclidW(precompute.twinSampleAxes4J_sdUnits, 1, 2)';
            end
        %% Initial signature focus for twin axes:
          precompute.signedFocusedW4Gs4I = calcSignatureFocus(precompute.twinGeneAxes4I_sdUnits, lc_currentData.noiseEstimation.PSignalStrength4G, [], [], 1, inInfo.searchStrategy.signatureFocus, sL);
          precompute.signedFocusedW4Ps4J = calcSignatureFocus(precompute.twinSampleAxes4J_sdUnits,   lc_currentData.noiseEstimation.PSignalStrength4P, [], [], 2, inInfo.searchStrategy.signatureFocus, sL);
          %Add .full2focusWeightRatio weights from the current and the twin space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the twin space:
            precompute.signedFocusedW4Gs4I_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4I, nan(0,nP), 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
            precompute.signedFocusedW4Ps4J_withPerpendicularSpace = addFlatWeights(nan(nG,0), precompute.signedFocusedW4Ps4J, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
          %Performance/no need to compute _sdUnits vectors hereafter (only needed for the initial signature focus); assert that they are not used hereafter by setting them to NaN:
            if(true)
              precompute.sampleAxes4I_sdUnits(:) = NaN;
              precompute.geneAxes4J_sdUnits(:) = NaN;
              precompute.twinGeneAxes4I_sdUnits(:) = NaN;
              precompute.twinSampleAxes4J_sdUnits(:) = NaN;
              precompute.norm4sampleAxes4I_sdUnits(:) = NaN;
              precompute.norm4geneAxes4J_sdUnits(:) = NaN;
              precompute.norm4twinSampleAxes4J_sdUnits(:) = NaN;
              precompute.norm4twinGeneAxes4I_sdUnits(:) = NaN;
            end

        %% Correlate the signal with primary and twin axes of signature candidates, using the initial signature focus:
          %% Correlate twin axes with all other genes respectively samples:
            if(bStatusOutput) SDCM_printStatus(5,'      - correlation with twin axes...\n'); end
            %Correlate .twinGeneAxes4I with all other samples:
              if(true)
                [precompute.R4Ps4Itwin, precompute.sampleSizes4R4Ps4Itwin] = uncenteredWeightedCorrelation(...
                   precompute.twinGeneAxes4I_origUnits, lc_currentData.current.(fn4L2Rs_origUnits), 1 ...
                  ,abs(precompute.signedFocusedW4Gs4I_withPerpendicularSpace) ...
                  ,BM.meanW ...
                  ,false ...bJustDiagonal ...
                  ,iif(bOnlyOneSidedPerfSubspace4Corr... 
                    ,{BsPerformanceSubspace4Correlations{1},true(1,nP)} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                    ,BsPerformanceSubspace4Correlations...
                   )...
                );
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.R4Ps4Itwin = sign(precompute.R4Ps4Itwin).*dampenOutliers(abs(precompute.R4Ps4Itwin), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                end
              end
            %Correlate .twinSampleAxes4J with all other genes:
              if(true)
                [precompute.R4Gs4Jtwin, precompute.sampleSizes4R4Gs4Jtwin] = uncenteredWeightedCorrelation(...
                   precompute.twinSampleAxes4J_origUnits, lc_currentData.current.(fn4L2Rs_origUnits), 2 ...
                  ,abs(precompute.signedFocusedW4Ps4J_withPerpendicularSpace) ...
                  ,BM.meanW ...
                  ,false ...bJustDiagonal ...
                  ,iif(bOnlyOneSidedPerfSubspace4Corr... 
                    ,{true(nG,1),BsPerformanceSubspace4Correlations{2}} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                    ,BsPerformanceSubspace4Correlations... 
                   )...
                );
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.R4Gs4Jtwin = sign(precompute.R4Gs4Jtwin).*dampenOutliers(abs(precompute.R4Gs4Jtwin), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G(min(end,sL)), 1); 
                end
              end
          %% Correlate primary (or double-indirect) axes with all other genes respectively samples:
            if(bStatusOutput) SDCM_printStatus(5,'      - correlation with primary/double-indirect signatures...\n'); end
            %Correlate .sampleAxes4I with all other genes:
              if(true)
                [precompute.R4Gs4I, precompute.sampleSizes4R4Gs4I] = uncenteredWeightedCorrelation(...
                   precompute.sampleAxes4I_origUnits, lc_currentData.current.(fn4L2Rs_origUnits), 2 ...
                  ,abs(precompute.signedFocusedW4Ps4I_withPerpendicularSpace) ...
                  ,BM.meanW ...
                  ,false ...bJustDiagonal ...
                  ,iif(bOnlyOneSidedPerfSubspace4Corr... 
                    ,{true(nG,1),BsPerformanceSubspace4Correlations{2}} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                    ,BsPerformanceSubspace4Correlations... 
                   )...
                );
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.R4Gs4I = sign(precompute.R4Gs4I).*dampenOutliers(abs(precompute.R4Gs4I), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G(min(end,sL)), 1); 
                end
              end
            %Correlate .geneAxes4J with all other samples:
              if(true)
                [precompute.R4Ps4J, precompute.sampleSizes4R4Ps4J] = uncenteredWeightedCorrelation(...
                   precompute.geneAxes4J_origUnits, lc_currentData.current.(fn4L2Rs_origUnits), 1 ...
                  ,abs(precompute.signedFocusedW4Gs4J_withPerpendicularSpace) ...
                  ,BM.meanW ...
                  ,false ...bJustDiagonal ...
                  ,iif(bOnlyOneSidedPerfSubspace4Corr... 
                    ,{BsPerformanceSubspace4Correlations{1},true(1,nP)} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                    ,BsPerformanceSubspace4Correlations...
                   )...
                );
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.R4Ps4J = sign(precompute.R4Ps4J).*dampenOutliers(abs(precompute.R4Ps4J), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                end
              end
          %% Correlation significance: estimate p values for correlations (used in the signature focus):
            if(bStatusOutput) SDCM_printStatus(5,'      - estimating correlation significance...\n'); end
            precompute.P4R4Gs4I = pValues4Correlations(...
              precompute.R4Gs4I, precompute.sampleSizes4R4Gs4I...
            );
            precompute.P4R4Ps4I = pValues4Correlations(...
              precompute.R4Ps4Itwin, precompute.sampleSizes4R4Ps4Itwin...
            );

            precompute.P4R4Gs4J = pValues4Correlations(...
              precompute.R4Gs4Jtwin, precompute.sampleSizes4R4Gs4Jtwin...
            );
            precompute.P4R4Ps4J = pValues4Correlations(...
              precompute.R4Ps4J, precompute.sampleSizes4R4Ps4J...
            );

            %Preformance/not needed here; will be done below in the focusing step:
              precompute.log10_p4Correlations4I = nan(1,length(precompute.I),numericTargetPrecision);
              precompute.log10_p4Correlations4J = nan(1,length(precompute.J),numericTargetPrecision);
          %% Extended signature focus for each gene and sample (needed e.g. for signature size estiamtions):
            if(bStatusOutput) SDCM_printStatus(5,'      - computing extended signature focus from correlations...\n'); end
            precompute.signedExtendedW4Ps4I = calcSignatureFocus([],[], precompute.R4Ps4Itwin, precompute.P4R4Ps4I, 2, inInfo.searchStrategy.extendedFocus, sL);
            precompute.signedExtendedW4Gs4J = calcSignatureFocus([],[], precompute.R4Gs4Jtwin, precompute.P4R4Gs4J, 1, inInfo.searchStrategy.extendedFocus, sL);
            precompute.signedExtendedW4Gs4Itwin = calcSignatureFocus([],[], precompute.R4Gs4I, precompute.P4R4Gs4I, 1, inInfo.searchStrategy.extendedFocus, sL);
            precompute.signedExtendedW4Ps4Jtwin = calcSignatureFocus([],[], precompute.R4Ps4J, precompute.P4R4Ps4J, 2, inInfo.searchStrategy.extendedFocus, sL);
          %% Signature size estimation by correlation sums:
            bIntegrateFlag = true;
            precompute.signatureSizeByCorrSum4G4I = BM.meanW(...
               abs(precompute.signedExtendedW4Gs4Itwin)...
              ,inInfo.reference.rowW ...
            ,1, bIntegrateFlag);
            precompute.signatureSizeByCorrSum4P4I = BM.meanW(...
               abs(precompute.signedExtendedW4Ps4I)...
              ,inInfo.reference.colW ...
            ,2, bIntegrateFlag)';
            precompute.signatureSizeByCorrSum2D4I = precompute.signatureSizeByCorrSum4G4I.*precompute.signatureSizeByCorrSum4P4I; %sum(sum(bsxfun(@times,W4NextR4Ps4Itwin,W4NextR4Gs4I))) = sum(W4NextR4Ps4Itwin)*sum(W4NextR4Gs4I)

            precompute.signatureSizeByCorrSum4G4J = BM.meanW(...
               abs(precompute.signedExtendedW4Gs4J)...
              ,inInfo.reference.rowW ...
            ,1, bIntegrateFlag);
            precompute.signatureSizeByCorrSum4P4J = BM.meanW(...
               abs(precompute.signedExtendedW4Ps4Jtwin)...
              ,inInfo.reference.colW ...
            ,2, bIntegrateFlag)';
            precompute.signatureSizeByCorrSum2D4J = precompute.signatureSizeByCorrSum4G4J.*precompute.signatureSizeByCorrSum4P4J; %sum(sum(bsxfun(@times,W4NextR4Ps4Itwin,W4NextR4Gs4I))) = sum(W4NextR4Ps4Itwin)*sum(W4NextR4Gs4I)
      end

    %% Focusing step: To focus the feature and sample space onto the signature, use correlations to update weights (refined signautre focus), update the twin axes and correlate again in the refined signature focus:
      if(true)
        stepOne = precompute;
        %% Refined signature focus based on correlations for primary axes:
          precompute.signedFocusedW4Ps4I = calcSignatureFocus([],[], stepOne.R4Ps4Itwin, stepOne.P4R4Ps4I, 2, inInfo.searchStrategy.signatureFocus, sL);
          precompute.signedFocusedW4Gs4J = calcSignatureFocus([],[], stepOne.R4Gs4Jtwin, stepOne.P4R4Gs4J, 1, inInfo.searchStrategy.signatureFocus, sL);
          %Add .full2focusWeightRatio weights from the current and the twin space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the twin space:
            precompute.signedFocusedW4Ps4I_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4I, precompute.signedFocusedW4Ps4I, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
            precompute.signedFocusedW4Gs4J_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4J, precompute.signedFocusedW4Ps4J, 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
        %% Twin axes based on the refined signature focus:
          if(bStatusOutput) SDCM_printStatus(5,'      - focused aggregation of twin axes...\n'); end
          %Twin gene axes for .sampleAxes4I:
            if(true)
              sampleAxes4I_origUnits = iif(bRestrictToPerformanceSubspace, precompute.sampleAxes4I_origUnits(:,BsPerformanceSubspace4Correlations{2}), precompute.sampleAxes4I_origUnits);
              signedFocusedW4Ps4I = iif(bRestrictToPerformanceSubspace, precompute.signedFocusedW4Ps4I(:,BsPerformanceSubspace4Correlations{2}), precompute.signedFocusedW4Ps4I);
              precompute.twinGeneAxes4I_origUnits = zeros(nG,size(signedFocusedW4Ps4I,1));
              for ii=1:size(signedFocusedW4Ps4I,1)
                precompute.twinGeneAxes4I_origUnits(BsPerformanceSubspace4Correlations{1},ii) = projectW(...
                  reduced_currentData.L2Rs ... vectorsInSourceSpace...
                 ,sampleAxes4I_origUnits(ii,:) ... axesInSourceSpace
                 ,abs(signedFocusedW4Ps4I(ii,:)) ... weights4sourceDims
                 ,2 ... sourceSpaceMatrixDim
                ,BM.meanW, BM.euclidW) / BM.euclidW(ones(1,size(signedFocusedW4Ps4I,2)),abs(signedFocusedW4Ps4I(ii,:)),2); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
              end
                precompute.twinGeneAxes4I_origUnits(isnan(precompute.twinGeneAxes4I_origUnits)) = 0; %only needed for very high NaN rate (>~90%)
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.twinGeneAxes4I_origUnits = sign(precompute.twinGeneAxes4I_origUnits).*dampenOutliers(abs(precompute.twinGeneAxes4I_origUnits), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G(min(end,sL)), 1); 
                end
              %Compute vector norms:
                precompute.norm4twinGeneAxes4I_origUnits = BM.euclidW(precompute.twinGeneAxes4I_origUnits, 1, 1);          
            end
          %Twin sample axes for .geneAxes4J:
            if(true)
              geneAxes4J_origUnits = iif(bRestrictToPerformanceSubspace, precompute.geneAxes4J_origUnits(BsPerformanceSubspace4Correlations{1},:), precompute.geneAxes4J_origUnits);
              signedFocusedW4Gs4J = iif(bRestrictToPerformanceSubspace, precompute.signedFocusedW4Gs4J(BsPerformanceSubspace4Correlations{1},:), precompute.signedFocusedW4Gs4J);
              precompute.twinSampleAxes4J_origUnits = zeros(size(signedFocusedW4Gs4J,2),nP);
              for jj=1:size(signedFocusedW4Gs4J,2)
                precompute.twinSampleAxes4J_origUnits(jj,BsPerformanceSubspace4Correlations{2}) = projectW(...
                  reduced_currentData.L2Rs ... vectorsInSourceSpace...
                 ,geneAxes4J_origUnits(:,jj) ... axesInSourceSpace
                 ,abs(signedFocusedW4Gs4J(:,jj)) ... weights4sourceDims
                 ,1 ... sourceSpaceMatrixDim
                ,BM.meanW, BM.euclidW) / BM.euclidW(ones(size(signedFocusedW4Gs4J,1),1),abs(signedFocusedW4Gs4J(:,jj)),1); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
              end
                precompute.twinSampleAxes4J_origUnits(isnan(precompute.twinSampleAxes4J_origUnits)) = 0; %only needed for very high NaN rate (>~90%)
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.twinSampleAxes4J_origUnits = sign(precompute.twinSampleAxes4J_origUnits).*dampenOutliers(abs(precompute.twinSampleAxes4J_origUnits), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                end
              %Compute vector norms:
                precompute.norm4twinSampleAxes4J_origUnits = BM.euclidW(precompute.twinSampleAxes4J_origUnits, 1, 2)';
            end
        %% Refined signature focus based on correlations for twin axes:
            precompute.signedFocusedW4Gs4I = calcSignatureFocus([],[], stepOne.R4Gs4I, stepOne.P4R4Gs4I, 1, inInfo.searchStrategy.signatureFocus, sL);
            precompute.signedFocusedW4Ps4J = calcSignatureFocus([],[], stepOne.R4Ps4J, stepOne.P4R4Ps4J, 2, inInfo.searchStrategy.signatureFocus, sL);
          %Add .full2focusWeightRatio weights from the current and the twin space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the twin space:
            precompute.signedFocusedW4Gs4I_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4I, precompute.signedFocusedW4Ps4I, 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
            precompute.signedFocusedW4Ps4J_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4J, precompute.signedFocusedW4Ps4J, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));

        %% Correlate the signal with primary and twin axes of signature candidates, using the refined signature focus:
          %% Correlate twin axes with all other genes respectively samples:
            if(bStatusOutput) SDCM_printStatus(5,'      - focused correlation with twin axes...\n'); end
            %Correlate .twinGeneAxes4I with all other samples:
              if(true)
                [precompute.R4Ps4Itwin, precompute.sampleSizes4R4Ps4Itwin] = uncenteredWeightedCorrelation(...
                   precompute.twinGeneAxes4I_origUnits, lc_currentData.current.(fn4L2Rs_origUnits), 1 ...
                  ,abs(precompute.signedFocusedW4Gs4I_withPerpendicularSpace) ...
                  ,BM.meanW ...
                  ,false ...bJustDiagonal ...
                  ,iif(bOnlyOneSidedPerfSubspace4Corr... 
                    ,{BsPerformanceSubspace4Correlations{1},true(1,nP)} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                    ,BsPerformanceSubspace4Correlations...
                   )...
                );
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.R4Ps4Itwin = sign(precompute.R4Ps4Itwin).*dampenOutliers(abs(precompute.R4Ps4Itwin), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                end
              end
            %Correlate .twinSampleAxes4J with all other genes:
              if(true)
                [precompute.R4Gs4Jtwin, precompute.sampleSizes4R4Gs4Jtwin] = uncenteredWeightedCorrelation(...
                   precompute.twinSampleAxes4J_origUnits, lc_currentData.current.(fn4L2Rs_origUnits), 2 ...
                  ,abs(precompute.signedFocusedW4Ps4J_withPerpendicularSpace) ...
                  ,BM.meanW ...
                  ,false ...bJustDiagonal ...
                  ,iif(bOnlyOneSidedPerfSubspace4Corr... 
                    ,{true(nG,1),BsPerformanceSubspace4Correlations{2}} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                    ,BsPerformanceSubspace4Correlations... 
                   )...
                );
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.R4Gs4Jtwin = sign(precompute.R4Gs4Jtwin).*dampenOutliers(abs(precompute.R4Gs4Jtwin), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G(min(end,sL)), 1); 
                end
              end
          %% Correlate primary (or double-indirect) axes with all other genes respectively samples:
            if(bStatusOutput) SDCM_printStatus(5,'      - focused correlation with primary/double-indirect signatures...\n'); end
            %Correlate .sampleAxes4I with all other genes:
              if(true)
                [precompute.R4Gs4I, precompute.sampleSizes4R4Gs4I] = uncenteredWeightedCorrelation(...
                   precompute.sampleAxes4I_origUnits, lc_currentData.current.(fn4L2Rs_origUnits), 2 ...
                  ,abs(precompute.signedFocusedW4Ps4I_withPerpendicularSpace) ...
                  ,BM.meanW ...
                  ,false ...bJustDiagonal ...
                  ,iif(bOnlyOneSidedPerfSubspace4Corr... 
                    ,{true(nG,1),BsPerformanceSubspace4Correlations{2}} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                    ,BsPerformanceSubspace4Correlations... 
                   )...
                );
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.R4Gs4I = sign(precompute.R4Gs4I).*dampenOutliers(abs(precompute.R4Gs4I), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G(min(end,sL)), 1); 
                end
              end
            %Correlate .geneAxes4J with all other samples:
              if(true)
                [precompute.R4Ps4J, precompute.sampleSizes4R4Ps4J] = uncenteredWeightedCorrelation(...
                   precompute.geneAxes4J_origUnits, lc_currentData.current.(fn4L2Rs_origUnits), 1 ...
                  ,abs(precompute.signedFocusedW4Gs4J_withPerpendicularSpace) ...
                  ,BM.meanW ...
                  ,false ...bJustDiagonal ...
                  ,iif(bOnlyOneSidedPerfSubspace4Corr... 
                    ,{BsPerformanceSubspace4Correlations{1},true(1,nP)} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                    ,BsPerformanceSubspace4Correlations...
                   )...
                );
                if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                  precompute.R4Ps4J = sign(precompute.R4Ps4J).*dampenOutliers(abs(precompute.R4Ps4J), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                end
              end
          %% Correlation significance: estimate p values for correlations (used in the signature focus):
            if(bStatusOutput) SDCM_printStatus(5,'      - estimating focused correlation significance...\n'); end
            [precompute.P4R4Gs4I, log10_p4Correlations4G4I] = pValues4Correlations(...
              precompute.R4Gs4I, precompute.sampleSizes4R4Gs4I...
            ,1);
            [precompute.P4R4Ps4I, log10_p4Correlations4P4I] = pValues4Correlations(...
              precompute.R4Ps4Itwin, precompute.sampleSizes4R4Ps4Itwin...
            ,2); log10_p4Correlations4P4I=log10_p4Correlations4P4I';
              precompute.log10_p4Correlations4I = min(log10_p4Correlations4G4I, log10_p4Correlations4P4I); %these are two indirectly dependent p values (same L2Rs base); the min should be an upper bound of the combined p value; too conservative?!

            [precompute.P4R4Gs4J, log10_p4Correlations4G4J] = pValues4Correlations(...
              precompute.R4Gs4Jtwin, precompute.sampleSizes4R4Gs4Jtwin...
            ,1);
            [precompute.P4R4Ps4J, log10_p4Correlations4P4J] = pValues4Correlations(...
              precompute.R4Ps4J, precompute.sampleSizes4R4Ps4J...
            ,2); log10_p4Correlations4P4J=log10_p4Correlations4P4J';
              precompute.log10_p4Correlations4J = min(log10_p4Correlations4G4J, log10_p4Correlations4P4J); %these are two indirectly dependent p values (same L2Rs base); the min should be an upper bound of the combined p value; too conservative?!
          %% Signature focus:
            if(bStatusOutput) SDCM_printStatus(5,'      - updating signature foci with correlations...\n'); end
            if(true)
              precompute.signedFocusedW4Ps4I = calcSignatureFocus([],[], precompute.R4Ps4Itwin, precompute.P4R4Ps4I, 2, inInfo.searchStrategy.signatureFocus, sL);
              precompute.signedFocusedW4Gs4J = calcSignatureFocus([],[], precompute.R4Gs4Jtwin, precompute.P4R4Gs4J, 1, inInfo.searchStrategy.signatureFocus, sL);
              %Add .full2focusWeightRatio weights from the current and the twin space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the twin space:
                precompute.signedFocusedW4Ps4I_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4I, precompute.signedFocusedW4Ps4I, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
                precompute.signedFocusedW4Gs4J_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4J, precompute.signedFocusedW4Ps4J, 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
            end
            if(true)
              precompute.signedFocusedW4Gs4I = calcSignatureFocus([],[], precompute.R4Gs4I, precompute.P4R4Gs4I, 1, inInfo.searchStrategy.signatureFocus, sL);
              precompute.signedFocusedW4Ps4J = calcSignatureFocus([],[], precompute.R4Ps4J, precompute.P4R4Ps4J, 2, inInfo.searchStrategy.signatureFocus, sL);
              %Add .full2focusWeightRatio weights from the current and the twin space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the twin space:
                precompute.signedFocusedW4Gs4I_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4I, precompute.signedFocusedW4Ps4I, 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
                precompute.signedFocusedW4Ps4J_withPerpendicularSpace = addFlatWeights(precompute.signedFocusedW4Gs4J, precompute.signedFocusedW4Ps4J, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
            end
          %% Extended signature focus (needed e.g. for signature size estiamtions):
            if(bStatusOutput) SDCM_printStatus(5,'      - computing extended signature focus from correlations...\n'); end
            precompute.signedExtendedW4Ps4I = calcSignatureFocus([],[], precompute.R4Ps4Itwin, precompute.P4R4Ps4I, 2, inInfo.searchStrategy.extendedFocus, sL);
            precompute.signedExtendedW4Gs4J = calcSignatureFocus([],[], precompute.R4Gs4Jtwin, precompute.P4R4Gs4J, 1, inInfo.searchStrategy.extendedFocus, sL);
            precompute.signedExtendedW4Gs4Itwin = calcSignatureFocus([],[], precompute.R4Gs4I, precompute.P4R4Gs4I, 1, inInfo.searchStrategy.extendedFocus, sL);
            precompute.signedExtendedW4Ps4Jtwin = calcSignatureFocus([],[], precompute.R4Ps4J, precompute.P4R4Ps4J, 2, inInfo.searchStrategy.extendedFocus, sL);
          %% Signature size estimation by correlation sums:
            bIntegrateFlag = true;
            precompute.signatureSizeByCorrSum4G4I = BM.meanW(abs(precompute.signedExtendedW4Gs4Itwin), inInfo.reference.rowW, 1, bIntegrateFlag);
            precompute.signatureSizeByCorrSum4P4I = BM.meanW(abs(precompute.signedExtendedW4Ps4I), inInfo.reference.colW, 2, bIntegrateFlag)';
            precompute.signatureSizeByCorrSum2D4I = precompute.signatureSizeByCorrSum4G4I.*precompute.signatureSizeByCorrSum4P4I; %sum(sum(bsxfun(@times,W4NextR4Ps4Itwin,W4NextR4Gs4I))) = sum(W4NextR4Ps4Itwin)*sum(W4NextR4Gs4I)

            precompute.signatureSizeByCorrSum4G4J = BM.meanW(abs(precompute.signedExtendedW4Gs4J), inInfo.reference.rowW, 1, bIntegrateFlag);
            precompute.signatureSizeByCorrSum4P4J = BM.meanW(abs(precompute.signedExtendedW4Ps4Jtwin), inInfo.reference.colW, 2, bIntegrateFlag)';
            precompute.signatureSizeByCorrSum2D4J = precompute.signatureSizeByCorrSum4G4J.*precompute.signatureSizeByCorrSum4P4J; %sum(sum(bsxfun(@times,W4NextR4Ps4Itwin,W4NextR4Gs4I))) = sum(W4NextR4Ps4Itwin)*sum(W4NextR4Gs4I)
      end

    %% Rate signature candidates (compute average signal and average correlation strength needed during qualification in step 1): 
      if(bStatusOutput) SDCM_printStatus(5,'      - computing average signal strengths of the candidate signatures...\n'); drawnow; end
      %Process gene-based signature candidates:
        %Initialize outputs:
          precompute.signatureAbsMean2D4I = nan(1,length(precompute.I),numericTargetPrecision);
          precompute.sampleSize4signatureAbsMean2D4I = nan(1,length(precompute.I),numericTargetPrecision);
          precompute.signatureAbsSD2D4I = nan(1,length(precompute.I),numericTargetPrecision);
          precompute.signatureCorrInExtendedFocus4I = nan(1,length(precompute.I),numericTargetPrecision);
          precompute.log10_p4SignalStrength4I = nan(1,length(precompute.I),numericTargetPrecision);
          precompute.log10_p4I = nan(1,length(precompute.I),numericTargetPrecision);
        for ii=1:length(precompute.I)
          %Performance: exclude zero weights to save RAM:
            BNonZero4G = precompute.signedFocusedW4Gs4I(:,ii)~=0;
            BNonZero4P = precompute.signedFocusedW4Ps4I(ii,:)~=0;
          %Compute signature statistics:
            [ precompute.signatureAbsMean2D4I(ii), precompute.sampleSize4signatureAbsMean2D4I(ii), precompute.signatureAbsSD2D4I(ii), precompute.log10_p4SignalStrength4I(ii)...
             ,precompute.signatureCorrInExtendedFocus4I(ii), precompute.log10_p4I(ii) ...
            ] = computeSignatureStatistics(...
               lc_currentData.current.L2Rs(BNonZero4G,BNonZero4P) ...
              ,precompute.signedExtendedW4Gs4Itwin(BNonZero4G,ii), precompute.signedExtendedW4Ps4I(ii,BNonZero4P) ...
              ,precompute.R4Gs4I(BNonZero4G,ii), precompute.R4Ps4Itwin(ii,BNonZero4P) ...
              ,precompute.signatureSizeByCorrSum4G4I(ii), precompute.signatureSizeByCorrSum4P4I(ii)...
              ,precompute.log10_p4Correlations4I(ii) ...
              ,lc_currentData.noiseEstimation ...
              ,BM ...
            );
        end
      %Process sample based signature candidates:
        %Initialize outputs:
          precompute.signatureAbsMean2D4J = nan(1,length(precompute.J),numericTargetPrecision);
          precompute.sampleSize4signatureAbsMean2D4J = nan(1,length(precompute.J),numericTargetPrecision);
          precompute.signatureAbsSD2D4J = nan(1,length(precompute.J),numericTargetPrecision);
          precompute.signatureCorrInExtendedFocus4J = nan(1,length(precompute.J),numericTargetPrecision);
          precompute.log10_p4SignalStrength4J = nan(1,length(precompute.J),numericTargetPrecision);
          precompute.log10_p4J = nan(1,length(precompute.J),numericTargetPrecision);
        for jj=1:length(precompute.J)
          %Performance: exclude zero weights to save RAM:
            BNonZero4G = precompute.signedFocusedW4Gs4J(:,jj)~=0;
            BNonZero4P = precompute.signedFocusedW4Ps4J(jj,:)~=0;
          %Compute signature statistics:
            [ precompute.signatureAbsMean2D4J(jj), precompute.sampleSize4signatureAbsMean2D4J(jj), precompute.signatureAbsSD2D4J(jj), precompute.log10_p4SignalStrength4J(jj)...
             ,precompute.signatureCorrInExtendedFocus4J(jj), precompute.log10_p4J(jj) ...
            ] = computeSignatureStatistics(...
               lc_currentData.current.L2Rs(BNonZero4G,BNonZero4P) ...
              ,precompute.signedExtendedW4Gs4J(BNonZero4G,jj), precompute.signedExtendedW4Ps4Jtwin(jj,BNonZero4P) ...
              ,precompute.R4Gs4Jtwin(BNonZero4G,jj), precompute.R4Ps4J(jj,BNonZero4P) ...
              ,precompute.signatureSizeByCorrSum4G4J(jj), precompute.signatureSizeByCorrSum4P4J(jj)...
              ,precompute.log10_p4Correlations4J(jj) ...
              ,lc_currentData.noiseEstimation ...
              ,BM ...
            );
        end
  end

