%ABSTRACT
% Subfunction for SDCM step 1: Detect the next signature by searching for a 
% qualifying initial representative gene or sample with high signature functional.

  function [signatureInfo, processingOrder4Step1] = SDCM_step1_searchStrategy(inInfo)
    %% Initialize and precalculate:
      global eDState;
      %Abbreviations:
        nG = size(eDState.current.L2Rs,1);
        nP = size(eDState.current.L2Rs,2);            
        sL = min(eDState.current.multiPassesLevel, inInfo.searchStrategy.nDefinedPasses);
      %Get basic math function handles according to bDataContainsNaNs:
        BM = getBasicMathFunctions(inInfo.preprocessing.bDataContainsNaNs);
      %Remember bestScoresSeenSoFar (initialize all fields here for status output):
        if(true)
          bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D = eDState.signatureTemplate;
            bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.signatureSizeByCorrSum2D = 0; %first comparison value.
            bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.label = 'NA';
            bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.cursor = NaN;
            
          bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus = eDState.signatureTemplate;
            bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.signatureCorrInExtendedFocus = 0; %first comparison value.
            bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.label = 'NA';
            bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.cursor = NaN;

          bestScoresSeenSoFar.highest_signatureAbsMean2D = eDState.signatureTemplate;
            bestScoresSeenSoFar.highest_signatureAbsMean2D.signatureAbsMean2D = 0; %first comparison value.
            bestScoresSeenSoFar.highest_signatureAbsMean2D.label = 'NA';
            bestScoresSeenSoFar.highest_signatureAbsMean2D.cursor = NaN;
          
          bestScoresSeenSoFar.lowest_log10_p = eDState.signatureTemplate;
            bestScoresSeenSoFar.lowest_log10_p.log10_p = 1; %first comparison value.
            bestScoresSeenSoFar.lowest_log10_p.label = 'NA';
            bestScoresSeenSoFar.lowest_log10_p.cursor = NaN;
        end
      %Collect previously detected signature norms (to compute normRatiosToStrongest signature statistics in status output):
        alreadyDetected = struct();
        if(eDState.current.k>1)
          alreadyDetected.norm4geneAxis_origUnits = [eDState.signatures(1:eDState.current.k-1).step2_finalSignatureAxes]; 
          alreadyDetected.norm4sampleAxis_origUnits = [eDState.signatures(1:eDState.current.k-1).step2_finalSignatureAxes]; 
          alreadyDetected.norm4geneAxis_origUnits = [alreadyDetected.norm4geneAxis_origUnits.norm4geneAxis_origUnits]; %these are always the norm4geneAxis_origUnits; standardized norms are in norm4sdSampleAxis in the public eigenOrder interface.
          alreadyDetected.norm4sampleAxis_origUnits = [alreadyDetected.norm4sampleAxis_origUnits.norm4sampleAxis_origUnits]; %these are always the norm4sampleAxis_origUnits; standardized norms are in norm4sdGeneAxis in the public eigenOrder interface.
        else
          alreadyDetected.norm4geneAxis_origUnits = NaN; 
          alreadyDetected.norm4sampleAxis_origUnits = NaN; 
          alreadyDetected.norm4geneAxis_origUnits = NaN;
          alreadyDetected.norm4sampleAxis_origUnits = NaN;
        end
    %% Local Greedy scores: To efficiently find genes or samples that have the most and most highly correlated partners, we test them in the order of descending greedy scores:
      processingOrder4Step1 = struct();
      %gene row signature scores:
        processingOrder4Step1.genes.greedyScores = localGreedyScores4ProcessingOrder(1);
      %sample column signature scores:
        processingOrder4Step1.samples.greedyScores = localGreedyScores4ProcessingOrder(2);
      %Combine row and column scores in a single processing order:
        if(true)
          processingOrder4Step1.global.greedyScores = [];
          processingOrder4Step1.global.signatureI = [];
          if(ismember('sampleAxis',inInfo.searchStrategy.allowedSignatureTypes))
            processingOrder4Step1.global.greedyScores = [processingOrder4Step1.global.greedyScores; processingOrder4Step1.genes.greedyScores];
            processingOrder4Step1.global.signatureI = [processingOrder4Step1.global.signatureI, 1:nG];
          end
          if(ismember('geneAxis',inInfo.searchStrategy.allowedSignatureTypes))
            processingOrder4Step1.global.greedyScores = [processingOrder4Step1.global.greedyScores; processingOrder4Step1.samples.greedyScores'];
            processingOrder4Step1.global.signatureI = [processingOrder4Step1.global.signatureI, -(1:nP)];                  
          end
          %Performance subspace: exclude noise based on weak signal strengths from the search for high correlations:
            BsPerformanceSubspace4Correlations = selectDimensionsInPerformanceSubspace(nG,nP ... %which genes/samples should be included when computing correlation below?
             ,nan(nG,1), nan(1,nP) ... %no signature has been detected yet => nothing is included into the performance subspace because of high correlations.
             ,eDState.noiseEstimation.log10PNoise4G, eDState.noiseEstimation.log10PNoise4P ... %deselect those with too weak signal
             ,inInfo.searchStrategy ... %contains cutoff params like .performance.alpha4Signal4inclusionInComputation4G
            );
            %=>Set the greedy scores of noise genes/samples to zero to also exclude them from candidate processing
              processingOrder4Step1.global.greedyScores([~BsPerformanceSubspace4Correlations{1}; ~BsPerformanceSubspace4Correlations{2}']) = 0;
          %Sort greedy scores to determine the order of candidate processing:
            [~,processingOrder4Step1.global.SII] = sort(illa(processingOrder4Step1.global.greedyScores,@isnan,0),'descend');
            processingOrder4Step1.global.SImJ = processingOrder4Step1.global.signatureI(processingOrder4Step1.global.SII);
        end
      %Performance: precompute constant values for all candidate test iterations:
        if(true)
          %initialize/clear precompute cache:
            if(ilv(gcp('nocreate'),@(pool)iif(isempty(pool),0,@()pool.NumWorkers))==0) %matlabpool('size')==0)
              warning('Hint: The pool of Matlab workers is currently closed; consider >>parpool(4) or >>parpool(8) depending on your CPU to speed up processing by parallelization.');
            end
            precomputeCandidatesInParallel([], inInfo); 
            %Cache logic: Performance/Restore candidate initial representatives computed in the last detection iteration that were not affected/invalidated by the last effect:
              if(inInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives(min(end,sL)) && isfield(eDState.current,'precompute_previousIterations') && ~isempty(eDState.current.precompute_previousIterations))
                eDState.current.precompute = eDState.current.precompute_previousIterations;
                cachedImJFromPreviousIterations = eDState.current.precompute.ImJ;
                if(~isempty(cachedImJFromPreviousIterations))
                  SDCM_printStatus(2,[
                       '   - restored %d precomputed initial representative candidates from the cache based on previous detection iteration(s).\n' ...
                     ]...
                    ,length(cachedImJFromPreviousIterations) ...
                  );
                end
%<-dev.note: maybe we should remove the cache, as it is hard to predict/invalidate cache entries that were affected "too much" by the last dissection and need to be recalculated... Additionally, this performance feature is incompatible with resume.
              else
                cachedImJFromPreviousIterations = [];
              end
          %Backup cache state at start of step 1 (for dev/resume):
            eDState.current.precompute_ImJ_atStartOfLastStep1 = eDState.current.precompute.ImJ;
          rowLabels = inInfo.reference.rowIDs;
            rowLabels = cellfun(@(c)iif(isnumeric(c),@()num2str(c),c), rowLabels, 'UniformOutput',false);
            rowLabels = strrep(rowLabels, '\color[rgb]{0 0 0.67}',''); %hotfix: do not display tex color formatting on console.
          colLabels = inInfo.reference.colIDs;
            colLabels = cellfun(@(c)iif(isnumeric(c),@()num2str(c),c), colLabels, 'UniformOutput',false);
            colLabels = strrep(colLabels, '\color[rgb]{0 0 0.67}',''); %hotfix: do not display tex color formatting on console.
          maxLabelLength = max(max(cellfun(@length,rowLabels)),max(cellfun(@length,colLabels))); %for status message alignment.
          bMoreDetectionPassesAhead = ...
              ~eDState.current.performance.bNoSignatureDetectedSinceLastBDidNotQualifyReset ... %if we found another signature in the last sweep
            || eDState.current.multiPassesLevel < inInfo.searchStrategy.nDefinedPasses ... %or if there are more sensitive configurations left to test:
          ;
        end
    %% Loop over gene and sample candidates in descending score of interest:
      nextValidCursor = 1; nVisitedSignatureCandidates = 0;
      bEarlyGlobalBreak_tooLongNoQualification = false;
      cursor = 0;
      initialSignatureScores4localRanking = nan(nG+nP,1); %just for status output.
      SDCM_printStatus(1 ...
        ,'   - starting initial representative search; current qualification thresholds are: .minAbsCorrSum4G=%0.2f, .minAbsCorrSum4P=%0.2f, .minAbsCorrSumSum=%0.2f, .minCorrInExtendedFocus=%0.2f, log10(.alpha)=(%6.1fcorr,%6.1fsignal,%6.1fcombined)\n'... 
        ,inInfo.searchStrategy.qualification.minAbsCorrSum4G(min(end,sL))...
        ,inInfo.searchStrategy.qualification.minAbsCorrSum4P(min(end,sL))...
        ,inInfo.searchStrategy.qualification.minAbsCorrSumSum(min(end,sL))...
        ,inInfo.searchStrategy.qualification.minCorrInExtendedFocus(min(end,sL))...
        ,log10(inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCorrP(min(end,sL)),nG+nP,1)) ...
        ,log10(inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForSignalP(min(end,sL)),nG+nP,1)) ...
        ,log10(inInfo.searchStrategy.qualification.alpha4combined(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCombinedP(min(end,sL)),nG+nP,1)) ...
      );
      while(cursor<length(processingOrder4Step1.global.SImJ)); cursor=cursor+1;
        %Loop logic: get next initial representative candidate:
          %loop forward to previously identified stronger lookahead candidate, if any.
            if(cursor<nextValidCursor) continue; end 
          %Get next candidate index:
            imj = processingOrder4Step1.global.SImJ(cursor);
          %Performance/skip genes and samples that were deselected already on greedy score level (e.g. outside of the BsPerformanceSubspace4Correlations) as per performance configuration:
            if(processingOrder4Step1.global.greedyScores(processingOrder4Step1.global.SII(cursor))==0)
              if(inInfo.export.nStatusOutputLevel >= 7) %performance: only compute the status meassage if requested.
                SDCM_printStatus(7 ...
                  ,['   - skipping candidate    #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]               (localZ=%4.2f), since it was not member of the performance subspace (cf. inInfo.searchStrategy.performance.alpha4Signal4inclusionInComputation4G)\n'] ...
                  ...Candidate indices and name:
                      ,cursor...
                      ,iif(imj>0, 'gene', 'sample') ...
                      ,iif(imj>0, +imj, -imj)...
                      ,iif(imj>0, @()rowLabels{+imj}, @()colLabels{-imj}) ...
                  ...Scores:
                    ...Local scores:
                        ,iif(imj>0, @()processingOrder4Step1.genes.greedyScores(+imj), @()processingOrder4Step1.samples.greedyScores(-imj)) ...
                );
              end
              continue;
            end
          %Performance/Skip, if the gene did not qualify before and if we search for a initial representative:
            if( imj>0 && eDState.current.performance.BGeneDidNotQualify(+imj)...
             || imj<0 && eDState.current.performance.BSampleDidNotQualify(-imj)...
            )
              if(inInfo.export.nStatusOutputLevel >= 5) %performance: only compute the status meassage if requested.
                SDCM_printStatus(5 ...
                  ,['   - skipping candidate    #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]               (localZ=%4.2f), since it did not qualify in a previous pass of this sensitivity level\n'] ...
                  ...Candidate indices and name:
                      ,cursor...
                      ,iif(imj>0, 'gene', 'sample') ...
                      ,iif(imj>0, +imj, -imj)...
                      ,iif(imj>0, @()rowLabels{+imj}, @()colLabels{-imj}) ...
                  ...Scores:
                    ...Local scores:
                        ,iif(imj>0, @()processingOrder4Step1.genes.greedyScores(+imj), @()processingOrder4Step1.samples.greedyScores(-imj)) ...
                );
              end
              continue;
            end
          %Performance: parallel precalculation of signatureCandidate*L2Rs covariances for chunks of the next initial representative candidates:
            if(~ismember(imj, eDState.current.precompute.ImJ))
              %Define ImJ to precalc: %<-Note: keep the performance/skip config in sync with the outer loop.
                precalcChunkSize = inInfo.searchStrategy.nLookahead4HigherScore(min(end,sL)); %at least the required lookahead range.              
                nextNeededImJ = [];
                precalcCursorOffset = -1;
                while(length(nextNeededImJ)<precalcChunkSize && cursor+precalcCursorOffset<length(processingOrder4Step1.global.SImJ))
                  precalcCursorOffset = precalcCursorOffset + 1;
                  %Performance/skip genes and samples that were deselected already on greedy score level as per performance configuration:
                    if(processingOrder4Step1.global.greedyScores(processingOrder4Step1.global.SII(cursor+precalcCursorOffset))==0) 
                      continue; 
                    end
                  imj_precalcCandidate = processingOrder4Step1.global.SImJ(cursor+precalcCursorOffset);
                  %Performance/skip genes and samples that are already in the cache (from previous iterations):
                    if(ismember(imj_precalcCandidate,eDState.current.precompute.ImJ))
                      continue; 
                    end
                  %Skip if already previously visited in thie detection pass and the candidate did not qualify:
                    if( imj_precalcCandidate>0 && eDState.current.performance.BGeneDidNotQualify(imj_precalcCandidate) ...
                     || imj_precalcCandidate<0 && eDState.current.performance.BSampleDidNotQualify(-imj_precalcCandidate) ...
                    ) 
                      continue; %performance.
                    end 
                  nextNeededImJ(end+1) = imj_precalcCandidate;
                end
              %Precalc next initial representative candidates:
                SDCM_printStatus(2,'   - Precomputing correlations with the next %d initial representative candidates (%d:%d+%d)/%d in parallel (performance subspace = [%d/%d,%d/%d])...\n'...
                  ,length(nextNeededImJ)...
                  ,cursor...
                  ,cursor,precalcCursorOffset...
                  ,length(processingOrder4Step1.global.SImJ)...
                  ,sum(BsPerformanceSubspace4Correlations{1}),length(BsPerformanceSubspace4Correlations{1})...
                  ,sum(BsPerformanceSubspace4Correlations{2}),length(BsPerformanceSubspace4Correlations{2})...
                );
                assert(ismember(imj, nextNeededImJ), 'code validation: current initial representative index imj=%d not present in nextNeededImJ', imj);
                %precomputeCandidatesInParallel(nextNeededImJ, inInfo); %no need to save this for step2 as step2 does not use double indirect signatures.
                precomputeCandidatesInParallel([eDState.current.precompute.ImJ,nextNeededImJ], inInfo, BsPerformanceSubspace4Correlations); %save step1 cahce for step2
                drawnow; %keep GUI responsive.
            end
          %Get precomputed signature for imj (twin axis, axes correlations to all genes/samples and foci):
            signatureCand = getSignatureAxesAndScores(imj);
              bIsGeneAxisCand = signatureCand.imj>0;
              bIsSampleAxisCand = signatureCand.imj<0;
              signatureCand.greedyScore = processingOrder4Step1.global.greedyScores(processingOrder4Step1.global.SII(cursor));
            %Status output: initial representative candidate and its scores:
              if(inInfo.export.nStatusOutputLevel>=4) %performance: only compute the status meassage if requested.
                SDCM_printStatus(4 ...
                  ,['   - processing candidate  #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]               (corr=(% 7.1fg,% 7.1fs, %0.2fr); signal=(%0.2ffocus=%1.2fSNR); log10(p)=(%7.1fcorr,%6.1fsignal))%s\n'] ... localZ=%4.2f; c4g=%+0.2f, c4p=%+0.2f; 
                  ...Candidate indices and name:
                      ,cursor...
                      ,iif(bIsGeneAxisCand, 'gene', 'sample') ...
                      ,iif(bIsGeneAxisCand, signatureCand.imj, -signatureCand.imj)...
                      ,iif(bIsGeneAxisCand, @()rowLabels{signatureCand.i}, @()colLabels{signatureCand.j}) ...
                  ...Scores:
                    ...Signature size estimates and average absolute correlation in the focus:
                        ,signatureCand.signatureSizeByCorrSum4G, signatureCand.signatureSizeByCorrSum4P...
                        ,signatureCand.signatureCorrInExtendedFocus ...
                    ...Signal strengths relative to initial noise level:
                        ,signatureCand.signatureAbsMean2D, signatureCand.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ... 
                    ...Signature significance:
                        ,signatureCand.log10_p4Correlations, signatureCand.log10_p4SignalStrength ...
                    ...Cached flag:
                        ,iif(ismember(signatureCand.imj, cachedImJFromPreviousIterations), ' (cached)', '') ...
                );
              end
              if(mod(cursor,25)==0) drawnow; end

        %Qualification: Threshold checks, whether the current initial representative candidate is eligible at all based on the sensitivity configuration:
          nVisitedSignatureCandidates = nVisitedSignatureCandidates + 1;
          bSkipAndContinue = false;
          bEarlyGlobalBreak_tooLongNoQualification = ...
                bMoreDetectionPassesAhead && nVisitedSignatureCandidates > inInfo.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification(min(end,sL)) ...
            || ~bMoreDetectionPassesAhead && nVisitedSignatureCandidates > inInfo.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification4lastPass(min(end,sL)) ...
          ;
            if(bEarlyGlobalBreak_tooLongNoQualification)
              break; %if there are 1000 genes or samples in greedyScore order in a row without any qualifying as initial representative, early break and do not visit the rest.
            end
          [bQualified4InitialSignature, bPerformanceExcludeInThisPass] = testQualificationAsInitialRepresenative(signatureCand);
            if(~bQualified4InitialSignature)
              if(bPerformanceExcludeInThisPass)
                if(bIsGeneAxisCand)
                  eDState.current.performance.BGeneDidNotQualify(signatureCand.i) = true;
                end
                if(bIsSampleAxisCand)
                  eDState.current.performance.BSampleDidNotQualify(signatureCand.j) = true;
                end
              end
              bSkipAndContinue= true;
            end
          %Cache logic: if a qualified representative's score were computed based on the signal of a previous iteration, we need to recompute its correlations and recheck the thresholds (only after qualification for performance reasons):
            originaCursor = cursor;
            if(inInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives(min(end,sL)) && bQualified4InitialSignature && ismember(signatureCand.imj, cachedImJFromPreviousIterations))
              SDCM_printStatus(4, '     <- retesting this qualified representative candidate, since it is a cached result from a previous detection iteration...\n');
              cachedImJFromPreviousIterations(cachedImJFromPreviousIterations==signatureCand.imj) = [];
              precomputeCandidatesInParallel(signatureCand.imj, 'clearFromCache');
              precomputeCandidatesInParallel([eDState.current.precompute.ImJ,signatureCand.imj], inInfo, BsPerformanceSubspace4Correlations);
              cursor = cursor-1;
              bSkipAndContinue = true;
              signatureCand = getSignatureAxesAndScores(imj); %update for bestScoresSeenSoFar updates below.
            end
          %Logging: For the global break case that no more initial representative candidates qualify, keep track of the best scores seen as feedback to the user (Are my parameters too tight? How much more sensitivity do I need to get xyz as another qualified signature for an signature in a rerun?):
            if(signatureCand.signatureSizeByCorrSum2D > bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.signatureSizeByCorrSum2D)
              bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D = signatureCand;
              bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.cursor = originaCursor;
              bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.label = iif(bIsGeneAxisCand,@()rowLabels{signatureCand.i},@()colLabels{signatureCand.j});
            end
            if(signatureCand.signatureCorrInExtendedFocus > bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.signatureCorrInExtendedFocus)
              bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus = signatureCand;
              bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.cursor = originaCursor;
              bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.label = iif(bIsGeneAxisCand,@()rowLabels{signatureCand.i},@()colLabels{signatureCand.j});
            end
            if(signatureCand.signatureAbsMean2D > bestScoresSeenSoFar.highest_signatureAbsMean2D.signatureAbsMean2D)
              bestScoresSeenSoFar.highest_signatureAbsMean2D = signatureCand;
              bestScoresSeenSoFar.highest_signatureAbsMean2D.cursor = originaCursor;
              bestScoresSeenSoFar.highest_signatureAbsMean2D.label = iif(bIsGeneAxisCand,@()rowLabels{signatureCand.i},@()colLabels{signatureCand.j});
            end
            if(signatureCand.log10_p < bestScoresSeenSoFar.lowest_log10_p.log10_p)
              bestScoresSeenSoFar.lowest_log10_p = signatureCand;
              bestScoresSeenSoFar.lowest_log10_p.cursor = originaCursor;
              bestScoresSeenSoFar.lowest_log10_p.label = iif(bIsGeneAxisCand,@()rowLabels{signatureCand.i},@()colLabels{signatureCand.j});
            end
          %if the signature did not qualify or has to be retested as it was a cached item from a previous iteration, continue:
            if(bSkipAndContinue) continue; end
        %LOOKAHEAD: Enter the lookahead phase and check if there are stronger initial representative candidates in sight:
          assert(isempty(setdiff(fieldnames(signatureCand),fieldnames(eDState.signatureTemplate))), 'code validation: the signature template does not include all actually used signatureCand fields');
          score_candidate = signatureFunctional(signatureCand);
            initialSignatureScores4localRanking(iif(signatureCand.imj<0,nG-signatureCand.imj,signatureCand.imj)) = score_candidate;
          %Status output: display qualified signature:
            if(true)
              SDCM_printStatus(1 ...
                ,['   - candidate             #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ] HAS QUALIFIED (corr=(% 7.1fg,% 7.1fs, %0.2fr); signal=(%0.2ffocus=%1.2fSNR); log10(p)=(%7.1fcorr,%6.1fsignal); score=%0.2f, local rank=%4d) => entering lookahead phase of length %d...\n']... localZ=%4.2f; c4g=%+0.2f, c4p=%+0.2f; 
                ...Candidate indices and name:
                    ,cursor...
                    ,iif(bIsGeneAxisCand, 'gene', 'sample') ...
                    ,iif(bIsGeneAxisCand, signatureCand.imj, -signatureCand.imj)...
                    ,iif(bIsGeneAxisCand, @()rowLabels{signatureCand.i}, @()colLabels{signatureCand.j}) ...
                ...Scores:
                  ...Signature size estimates and average absolute correlation in the focus:
                      ,signatureCand.signatureSizeByCorrSum4G, signatureCand.signatureSizeByCorrSum4P...
                      ,signatureCand.signatureCorrInExtendedFocus ...
                  ...Signal strengths relative to initial noise level:
                      ,signatureCand.signatureAbsMean2D, signatureCand.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ....
                  ...Signature significance:
                      ,signatureCand.log10_p4Correlations, signatureCand.log10_p4SignalStrength ...
                  ...initial representative score and local rank (to give the user a feedback about how many detection iterations to wait before a particular qualified signature will be selected)
                      ,initialSignatureScores4localRanking(iif(signatureCand.imj<0,nG-signatureCand.imj,signatureCand.imj)), sum(initialSignatureScores4localRanking>=initialSignatureScores4localRanking(iif(signatureCand.imj<0,nG-signatureCand.imj,signatureCand.imj))) ...
               ,inInfo.searchStrategy.nLookahead4HigherScore(min(end,sL)) ...
              );
            end
          %Check the next genes (in the greedy order based on single gene statistics) for better correlation properties than the current candidate:
            bLookaheadFoundStrongerSignature = false;
            lookaheadCursor=1; %count signatureCand as one (also needed performance-wise to sync with the precalcChunks)
            actualLookaheads=0; nActualLookaheadsSinceLastAcceptedAccu=0;
            while(true)
              %End of lookahead logic:
                if(true)
                  %commmon end conditions:
                    bNoMoreCandidatesAvailable = cursor+lookaheadCursor >= length(processingOrder4Step1.global.SImJ);
                  %initial representative search mode / The search for a stronger initial representative ends if no stronger candidate was found either in the .nLookahead4HigherScore width:
                    %Check if we reached the configured lookahead width:
                      bConfiguredLookaheadWidthReached = ...
                           actualLookaheads >= inInfo.searchStrategy.nLookahead4HigherScore(min(end,sL)) ...
                      ;
                    %If the end conditions are met:
                      if(bNoMoreCandidatesAvailable || bConfiguredLookaheadWidthReached)
                        %If the lookahead ends before its configured width, explain why:
                          if(~bConfiguredLookaheadWidthReached)
                            sMessage = sprintf(...
                               '  <- ending lookahead phase after only %d/%d actual lookaheads since no more candidates qualified'...
                              ,actualLookaheads, inInfo.searchStrategy.nLookahead4HigherScore(min(end,sL)) ...
                            );
                            nDidNotQualify = iif(bIsGeneAxisCand ...
                              ,@()sum(eDState.current.performance.BGeneDidNotQualify(ilsub(processingOrder4Step1.global.SImJ(cursor:end),@(X)X>0)))...
                              ,@()sum(eDState.current.performance.BSampleDidNotQualify(abs(ilsub(processingOrder4Step1.global.SImJ(cursor:end),@(X)X<0)))) ...
                            );
                            sMessage = [sMessage, sprintf(...
                              ' (of the %d possible lookahead candidates, %d did not qualify previously on the same sensitivity level and the rest was performance-skipped due to zeroed greedy scores)\n'...
                              ,lookaheadCursor ...
                              ,nDidNotQualify ...
                            )];
                            SDCM_printStatus(2, sMessage);
                          end
                        break;
                      end
                end

              %Loop logic: get next lookahead initial representative candidate:
                lookaheadCursor = lookaheadCursor+1;
                %Performance/skip genes and samples that were deselected already on greedy score level as per performance configuration:
                  if(processingOrder4Step1.global.greedyScores(processingOrder4Step1.global.SII(cursor+lookaheadCursor))==0)
                    continue; %keep single line for profiling.
                  end
                %Get next lookahead candidate index:
                  imj_lkhd = processingOrder4Step1.global.SImJ(cursor+lookaheadCursor);
                %Performance/Skip, if the gene did not qualify before and if we search for a initial representative (in the focused lookahead, always consider all genes as they are already presented in correlation-multiplied order and we can use them in the accumulated signature, even if theyselves are no valid initial representatives):
                  if( imj_lkhd>0 && eDState.current.performance.BGeneDidNotQualify(+imj_lkhd) ...
                   || imj_lkhd<0 && eDState.current.performance.BSampleDidNotQualify(-imj_lkhd) ... 
                  ) 
                    continue; %keep single line for profiling. 
                  end 
                %Performance: parallel precalculation of signatureCandidate*L2Rs covariances for chunks of the next initial representative candidates:
                  if(~ismember(imj_lkhd, eDState.current.precompute.ImJ))
                    %Define ImJ to precompute: %<-Note: keep the performance/skip config in sync with the lookahead loop.
                      nextNeededImJ = [];
                      precalcChunkSize = inInfo.searchStrategy.nLookahead4HigherScore(min(end,sL))-actualLookaheads; %at least the remaining required lookahead range.
                      precalcCursorOffset = -1;
                      while(length(nextNeededImJ)<max(precalcChunkSize,inInfo.searchStrategy.nLookahead4HigherScore(min(end,sL))-actualLookaheads) && cursor+lookaheadCursor+precalcCursorOffset<length(processingOrder4Step1.global.SImJ)) %we need at least inInfo.searchStrategy.nLookahead4HigherScore(min(end,sL))-actualLookaheads for the lookahead to complete; use this information when precomputing (prevent a precompute chunk to end just before the required lookahead size as this required another inefficient precompute chunk; still use max with precalcChunkSize as we might as well find another stronger lookahead in the next chunk)
                        precalcCursorOffset = precalcCursorOffset + 1;
                        %Performance/skip genes and samples that were deselected already on greedy score level as per performance configuration:
                          if(processingOrder4Step1.global.greedyScores(processingOrder4Step1.global.SII(cursor+lookaheadCursor+precalcCursorOffset))==0) 
                            continue; 
                          end
                        imj_precalcCandidate = processingOrder4Step1.global.SImJ(cursor+lookaheadCursor+precalcCursorOffset);
                        %Performance/skip genes and samples that are already in the cache (from previous iterations):
                          if(ismember(imj_precalcCandidate,eDState.current.precompute.ImJ))
                            continue; 
                          end
                        if(imj_precalcCandidate>0)
                          %Performance/skip genes or samples that did not qualify previously in this detection pass:
                            if(eDState.current.performance.BGeneDidNotQualify(imj_precalcCandidate)) continue; end %performance. do not skip, if we need all possible candidates, even those with <minCorrGene correlating genes.
                        else
                          %Performance/skip genes or samples that did not qualify previously in this detection pass:
                            if(eDState.current.performance.BSampleDidNotQualify(-imj_precalcCandidate)) continue; end %performance. do not skip, if we need all possible candidates, even those with <minCorrGene correlating genes.
                        end
                        nextNeededImJ(end+1) = imj_precalcCandidate;
                      end
                    %Precalc next initial representative candidates:
                      SDCM_printStatus(2,'   - Precomputing correlations with the next %d initial representative candidates (%d:%d+%d)/%d in parallel (performance subspace = [%d/%d,%d/%d])...\n'...
                        ,length(nextNeededImJ)...
                        ,cursor+lookaheadCursor...
                        ,cursor+lookaheadCursor,precalcCursorOffset...
                        ,length(processingOrder4Step1.global.SImJ)...
                        ,sum(BsPerformanceSubspace4Correlations{1}),length(BsPerformanceSubspace4Correlations{1})...
                        ,sum(BsPerformanceSubspace4Correlations{2}),length(BsPerformanceSubspace4Correlations{2})...
                      );
                      assert(ismember(imj_lkhd, nextNeededImJ), 'code validation: current lookahead candidate index imj=%d not present in nextNeededImJ', imj_lkhd);
                      precomputeCandidatesInParallel([eDState.current.precompute.ImJ,nextNeededImJ], inInfo, BsPerformanceSubspace4Correlations); %save step1 cahce for step2
                      drawnow; %keep GUI responsive.
                  end

              %Get precomputed signature for imj (twin axis, axes correlations to all genes/samples and foci):
                signatureLkhd = getSignatureAxesAndScores(imj_lkhd);
                  actualLookaheads = actualLookaheads+1; %only count those as members of the lookahead interval whose scores were actually computed.
                  signatureLkhd.greedyScore = processingOrder4Step1.global.greedyScores(processingOrder4Step1.global.SII(cursor+lookaheadCursor));
                  bIsGeneAxisCand_lkhd = signatureLkhd.imj>0;
                  bIsSampleAxisCand_lkhd = signatureLkhd.imj<0;
                nActualLookaheadsSinceLastAcceptedAccu = nActualLookaheadsSinceLastAcceptedAccu+1; %for the check bConfiguredLookaheadWidthReached in signature generalization mode.
                %Status output: lookahead candidate and its scores
                  if(inInfo.export.nStatusOutputLevel>=4) %performance: only compute the status meassage if requested.
                    SDCM_printStatus(4 ...
                      ,['      - processing lkhd. candidate #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]               (corr=(% 7.1fg,% 7.1fs, %0.2fr); signal=(%0.2ffocus=%1.2fSNR); log10(p)=(%7.1fcorr,%6.1fsignal))%s\n'] ... localZ=%4.2f; c4g=%+0.2f, c4p=%+0.2f; 
                      ...Candidate indices and name:
                          ,cursor+lookaheadCursor...
                          ,iif(bIsGeneAxisCand_lkhd, 'gene', 'sample') ...
                          ,iif(bIsGeneAxisCand_lkhd, signatureLkhd.imj, -signatureLkhd.imj)...
                          ,iif(bIsGeneAxisCand_lkhd, @()rowLabels{signatureLkhd.i}, @()colLabels{signatureLkhd.j}) ...
                      ...Scores:
                        ...Signature size estimates and average absolute correlation in the focus:
                            ,signatureLkhd.signatureSizeByCorrSum4G, signatureLkhd.signatureSizeByCorrSum4P...
                            ,signatureLkhd.signatureCorrInExtendedFocus ...
                        ...Signal strengths relative to initial noise level:
                            ,signatureLkhd.signatureAbsMean2D, signatureLkhd.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ....
                        ...Signature significance:
                            ,signatureLkhd.log10_p4Correlations, signatureLkhd.log10_p4SignalStrength ...
                        ...Cached flag:
                            ,iif(ismember(signatureLkhd.imj, cachedImJFromPreviousIterations), ' (cached)', '') ...
                    );
                  end
                  if(mod(cursor+lookaheadCursor,25)==0) drawnow; end
                  if(actualLookaheads==1)
                    assert(isempty(setdiff(fieldnames(signatureLkhd),fieldnames(eDState.signatureTemplate))), 'code validation: the signature template does not include all actually used signatureLkhd fields');
                  end

              %Qualification: Threshold checks, whether the current initial representative candidate is eligible at all based on the sensitivity configuration:
                [bQualified4InitialSignature_lkhd, bPerformanceExcludeInThisPass_lkhd] = testQualificationAsInitialRepresenative(signatureLkhd);
                  if(~bQualified4InitialSignature_lkhd)
                    if(bPerformanceExcludeInThisPass_lkhd)
                      if(bIsGeneAxisCand_lkhd)
                        eDState.current.performance.BGeneDidNotQualify(signatureLkhd.i) = true;
                      end
                      if(bIsSampleAxisCand_lkhd)
                        eDState.current.performance.BSampleDidNotQualify(signatureLkhd.j) = true;
                      end
                    end
                    continue;
                  end

              %Is signatureLkhd better than signatureCand?: Evaluate lookahead candidate and skip to it if the lookahead found a better starting signature for step 2:
                score_lkhd = signatureFunctional(signatureLkhd);
                  initialSignatureScores4localRanking(iif(signatureLkhd.imj<0,nG-signatureLkhd.imj,signatureLkhd.imj)) = score_lkhd;
                if(score_lkhd > score_candidate * inInfo.searchStrategy.performance.minScoreIncreaseByLookaheadCandidate(min(end,sL))) ... %ignore tiny score increases for performance reasons.
                  %Cache logic: if the qualified signatureLkhd's score was computed based on the signal of a previous iteration, we need to recompute its correlations and recheck the thresholds (only after qualification for performance reasons):
                    if(inInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives(min(end,sL)) && ismember(signatureLkhd.imj, cachedImJFromPreviousIterations))
                      SDCM_printStatus(4, '     <- retesting this stronger lookahead candidate, since it is a cached result from a previous detection iteration...\n');
                      cachedImJFromPreviousIterations(cachedImJFromPreviousIterations==signatureLkhd.imj) = [];
                      precomputeCandidatesInParallel(signatureLkhd.imj, 'clearFromCache');
                      precomputeCandidatesInParallel([eDState.current.precompute.ImJ,signatureLkhd.imj], inInfo, BsPerformanceSubspace4Correlations);
                      lookaheadCursor = lookaheadCursor-1;
                      continue;
                    end
                  bLookaheadFoundStrongerSignature = true;
                  break; %if a better initial representative candidate was found during lookahead, skip to it in the outer loop and restart the lookahead phase of this loop.
                else
                  if(inInfo.export.nStatusOutputLevel>=5) %performance: only compute the status meassage if requested.
                    SDCM_printStatus(5 ...
                      ,['        <- current candidate #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ] is stronger; lookahead candidate just has score=%0.2f (relative rank=%4d).\n'] ...
                      ...Candidate indices and name:
                          ,cursor...
                          ,iif(bIsGeneAxisCand, 'gene', 'sample') ...
                          ,iif(bIsGeneAxisCand, signatureCand.imj, -signatureCand.imj)...
                          ,iif(bIsGeneAxisCand, @()rowLabels{signatureCand.i}, @()colLabels{signatureCand.j}) ...
                      ...initial representative score and local rank (to give the user a feedback about how many detection iterations to wait before a particular qualified signature will be selected)
                          ,initialSignatureScores4localRanking(iif(signatureLkhd.imj<0,nG-signatureLkhd.imj,signatureLkhd.imj)), sum(initialSignatureScores4localRanking>=initialSignatureScores4localRanking(iif(signatureLkhd.imj<0,nG-signatureLkhd.imj,signatureLkhd.imj))) ...
                    );
                  elseif(inInfo.export.nStatusOutputLevel<4) %at least log the qualified candidates in the lookahead phase (to gie the user a rough indication of how many signatures are still ahead)
                    SDCM_printStatus(2 ... %this is an important info to see how many candidates there still are from the log => use status level 2.
                      ,['        <- is stronger than       #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]               (corr=(% 7.1fg,% 7.1fs, %0.2fr); signal=(%0.2ffocus=%1.2fSNR); log10(p)=(%7.1fcorr,%6.1fsignal); score=% 2.2f, local rank=%4d)%s that also qualified.\n'] ... %localZ=%4.2f; c4g=%+0.2f, c4p=%+0.2f; 
                      ...Candidate indices and name:
                          ,cursor+lookaheadCursor...
                          ,iif(bIsGeneAxisCand_lkhd, 'gene', 'sample') ...
                          ,iif(bIsGeneAxisCand_lkhd, signatureLkhd.imj, -signatureLkhd.imj)...
                          ,iif(bIsGeneAxisCand_lkhd, @()rowLabels{signatureLkhd.i}, @()colLabels{signatureLkhd.j}) ...
                      ...Scores:
                        ...Signature size estimates and average absolute correlation in the focus:
                            ,signatureLkhd.signatureSizeByCorrSum4G, signatureLkhd.signatureSizeByCorrSum4P...
                            ,signatureLkhd.signatureCorrInExtendedFocus ...
                        ...Signal strengths relative to initial noise level:
                            ,signatureLkhd.signatureAbsMean2D, signatureLkhd.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ....
                        ...Signature significance:
                            ,signatureLkhd.log10_p4Correlations, signatureLkhd.log10_p4SignalStrength ...
                        ...initial representative score and local rank (to give the user a feedback about how many detection iterations to wait before a particular qualified signature will be selected)
                            ,initialSignatureScores4localRanking(iif(signatureLkhd.imj<0,nG-signatureLkhd.imj,signatureLkhd.imj)), sum(initialSignatureScores4localRanking>=initialSignatureScores4localRanking(iif(signatureLkhd.imj<0,nG-signatureLkhd.imj,signatureLkhd.imj))) ...
                        ...Cached flag:
                            ,iif(ismember(signatureLkhd.imj, cachedImJFromPreviousIterations), ' (cached)', '') ...
                    );
                  end
                  continue; %the outer loop signatureCand is stronger than the current lookahead candidate => continue with the next lookahead candidate.
                end
            end
          
        %Evaluate results from lookahead phase and continue searching or return:
          %if bNoMoreCandidatesAvailable and this is not the last detection pass, continue to the next detection pass (thereby revisiting eDState.current.performance.BGeneDidNotQualify and eDState.current.performance.BSampleDidNotQualify in order to try getting a full lookahead interval)
            if(bNoMoreCandidatesAvailable && ~bConfiguredLookaheadWidthReached && bMoreDetectionPassesAhead)
              SDCM_printStatus(2 ...
                ,'    <-  since the configured lookahead length was not reached, but another detection pass is yet to perform, skipping to it now...\n'...
                ,nextValidCursor...
              );
              break;
            end
          %inital signature search mode: if the lookahead found a better gene, skip to it and continue searching; else RETURN the qualified inital signature now:
            if(bLookaheadFoundStrongerSignature) %skip to the next interesting gene in the signature search loop:
              if(inInfo.export.nStatusOutputLevel>=3) %performance: only compute the status meassage if requested.
                SDCM_printStatus(3 ...
                  ,['        <- lookahead candidate        #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ] is stronger than the current initial representative candidate => now skipping to it.\n'] ...
                  ...Candidate indices and name:
                      ,cursor+lookaheadCursor...
                      ,iif(bIsGeneAxisCand_lkhd, 'gene', 'sample') ...
                      ,iif(bIsGeneAxisCand_lkhd, signatureLkhd.imj, -signatureLkhd.imj)...
                      ,iif(bIsGeneAxisCand_lkhd, @()rowLabels{signatureLkhd.i}, @()colLabels{signatureLkhd.j}) ...
                );
              end
              nextValidCursor = cursor+lookaheadCursor; %-1;
              SDCM_printStatus(3 ...
                ,'     -> lookahead found stronger initial representative candidate at rank #%03d; now skipping to it.\n'...
                ,nextValidCursor...
              );
            else %return the found qualified initial representative
              %Status output:
                signatureCand.sInitialRepresenative = strtrim(strrep(strrep(SDCM_printStatus(0 ...
                  ,[ '  => INITIAL REPRES.:      #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]              '...
                     ' (corr=(% 7.1fg,% 7.1fs, %0.2fr)'...
                     '; signal=(%0.2ffocus=%1.2fSNR)'...
                     '; normRatiosToStrongest=(%0.2fg,%0.2fs)'...
                     '; log10(p)=(%7.1fcorr,%6.1fsignal)'...
                     '; config.level=%d/%d'......
                  ')\n'] ...
                  ...Candidate indices and name:
                      ,cursor...
                      ,iif(bIsGeneAxisCand, 'gene', 'sample') ...
                      ,iif(bIsGeneAxisCand, signatureCand.imj, -signatureCand.imj)...
                      ,iif(bIsGeneAxisCand, @()rowLabels{signatureCand.i}, @()colLabels{signatureCand.j}) ...
                  ...Scores:
                    ...Signature size estimates and average absolute correlation in the focus:
                        ,signatureCand.signatureSizeByCorrSum4G, signatureCand.signatureSizeByCorrSum4P...
                        ,signatureCand.signatureCorrInExtendedFocus ...
                    ...Signal strengths relative to initial noise level:
                        ,signatureCand.signatureAbsMean2D, signatureCand.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ...
                    ...Relative signature norm to strongest so far detected signature:
                        ,signatureCand.norm4sampleAxis_origUnits/max(alreadyDetected.norm4sampleAxis_origUnits), signatureCand.norm4geneAxis_origUnits/max(alreadyDetected.norm4geneAxis_origUnits) ...
                    ...Signature significance:
                        ,signatureCand.log10_p4Correlations, signatureCand.log10_p4SignalStrength ...
                    ...Detection sensitivity level:
                        ,sL,inInfo.searchStrategy.nDefinedPasses ...
                ),'<-',''),']               (','] ('));
              %For interactive results control and searching by file, create a file name with the core statistics of the qualified initial representative:
                if(inInfo.export.infoFile4initialSignatureStats.bEnabled)
                  sInitialRepresenativeStatistics = sprintf('%03d, initial representative = corr=(%2.1fg,%2.1fs,%0.2fr); p=(%1.0ecorr,%1.0esignal), sL=%dof%d'...
                     ,eDState.current.k ...
                     ...Signature size estimates and average absolute correlation in the focus:
                         ,signatureCand.signatureSizeByCorrSum4G, signatureCand.signatureSizeByCorrSum4P...
                         ,signatureCand.signatureCorrInExtendedFocus ...
                     ...Signature significance:
                         ,10^signatureCand.log10_p4Correlations, 10^signatureCand.log10_p4SignalStrength ...
                     ...Detection sensitivity level:
                         ,sL,inInfo.searchStrategy.nDefinedPasses ...
                  );
                  fileName = sprintf('%s.info', sInitialRepresenativeStatistics);
                  targetDir = inInfo.export.rootDir;
                    if(ispc() && ~(length(targetDir)>=3 && strcmp(targetDir(1:3),'\\?'))) %add \\?\ prefix to support long paths.
                      %If we already have a network UNC path line \\Server\Share\..., we need to replace \\ by UNC\ to get long path support in Windows:
                        bIsAlreadyANetworkPath = length(targetDir)>=2 && strcmp(targetDir(1:2),'\\');
                        if(bIsAlreadyANetworkPath)
                          targetDir = ['UNC',targetDir(2:end)];
                        end
                      targetDir = ['\\?\',targetDir];
                    end
                  fid = fopen(fullfile(targetDir,fileName),'w');
                  fprintf(fid,'%s',sInitialRepresenativeStatistics);
                  fclose(fid);
                end
              %memory performance: remove .BGeneDidNotQualify genes and .BSampleDidNotQualify from the precompute cache as they are skipped before getSignatureAxesAndScores(imj) anyway and never needed again until pass reset (in which case the precompute cache is cleared, too)
                if(inInfo.internal.bEnablePostProductionCodeOptimizations)
                  precomputeCandidatesInParallel(intersect(eDState.current.precompute.ImJ...
                    ,[find(eDState.current.performance.BGeneDidNotQualify)',-find(eDState.current.performance.BSampleDidNotQualify)]...
                  ),'clearFromCache');
                end
              %Performance/save candidate initial representatives computed in this detection iteration for later iterations:
                if(inInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives(min(end,sL)))
                  eDState.current.precompute_previousIterations = eDState.current.precompute;
                end
              %Output:
                signatureInfo = signatureCand;
                return;
            end
      end
    %% Global break case: no gene or sample candidate qualified as initial representative any more:
      if(true)
        if(bEarlyGlobalBreak_tooLongNoQualification)
          SDCM_printStatus(1 ...
            ,[ '  <- no initial representative candidate qualified for inInfo.searchStrategy.globalBreakConditions.%s=%0.0f candidates.\n'...
              ,'     <- qualification thresholds were: .minAbsCorrSum4G=%0.2f, .minAbsCorrSum4P=%0.2f, .minAbsCorrSumSum=%0.2f, .minCorrInExtendedFocus=%0.2f, log10(.alpha)=(%6.1fcorr,%6.1fsignal,%6.1fcombined)\n'... .minSignatureSNR=%0.1f, 
             ]...
            ,iif(bMoreDetectionPassesAhead,'maxLookaheadWithoutAnyQualification','maxLookaheadWithoutAnyQualification4lastPass') ...
            ,iif(bMoreDetectionPassesAhead,inInfo.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification(min(end,sL)),inInfo.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification4lastPass(min(end,sL))) ...
            ,inInfo.searchStrategy.qualification.minAbsCorrSum4G(min(end,sL))...
            ,inInfo.searchStrategy.qualification.minAbsCorrSum4P(min(end,sL))...
            ,inInfo.searchStrategy.qualification.minAbsCorrSumSum(min(end,sL))...
            ,inInfo.searchStrategy.qualification.minCorrInExtendedFocus(min(end,sL))...
            ,log10(inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCorrP(min(end,sL)),nG+nP,1)) ...
            ,log10(inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForSignalP(min(end,sL)),nG+nP,1)) ...
            ,log10(inInfo.searchStrategy.qualification.alpha4combined(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCombinedP(min(end,sL)),nG+nP,1)) ...
          );
        else
          SDCM_printStatus(1 ...
            ,[ '  <- no possible initial representative candidates qualified any more wrt. the inInfo.searchStrategy.qualification configuration.\n'...
              ,'     <- qualification thresholds were: .minAbsCorrSum4G=%0.2f, .minAbsCorrSum4P=%0.2f, .minAbsCorrSumSum=%0.2f, .minCorrInExtendedFocus=%0.2f, log10(.alpha)=(%6.1fcorr,%6.1fsignal,%6.1fcombined)\n'... , .minSignatureSNR=%0.1f
             ]...
            ,inInfo.searchStrategy.qualification.minAbsCorrSum4G(min(end,sL))...
            ,inInfo.searchStrategy.qualification.minAbsCorrSum4P(min(end,sL))...
            ,inInfo.searchStrategy.qualification.minAbsCorrSumSum(min(end,sL))...
            ,inInfo.searchStrategy.qualification.minCorrInExtendedFocus(min(end,sL))...
            ,log10(inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCorrP(min(end,sL)),nG+nP,1)) ...
            ,log10(inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForSignalP(min(end,sL)),nG+nP,1)) ...
            ,log10(inInfo.searchStrategy.qualification.alpha4combined(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCombinedP(min(end,sL)),nG+nP,1)) ...
          );
        end
        SDCM_printStatus(1 ...
          ,['  <- Best candidates for initial representatives were:\n'...
           ,'     - with maximal signal strength:       corr=(% 7.1fg,% 7.1fs, %0.2fr), signal=(%0.2ffocus=%1.2fSNR), normRatiosToStrongest=(%0.2fg,%0.2fs), log10(p)=(%7.1fcorr,%6.1fsignal,%6.1fcombined) for #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]\n'...
           ,'     - with max. signatureSizeByCorrSum2D: corr=(% 7.1fg,% 7.1fs, %0.2fr), signal=(%0.2ffocus=%1.2fSNR), normRatiosToStrongest=(%0.2fg,%0.2fs), log10(p)=(%7.1fcorr,%6.1fsignal,%6.1fcombined) for #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]\n'...
           ,'     - with maximal averag corr:           corr=(% 7.1fg,% 7.1fs, %0.2fr), signal=(%0.2ffocus=%1.2fSNR), normRatiosToStrongest=(%0.2fg,%0.2fs), log10(p)=(%7.1fcorr,%6.1fsignal,%6.1fcombined) for #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]\n'...
           ,'     - with minimal log10_p:               corr=(% 7.1fg,% 7.1fs, %0.2fr), signal=(%0.2ffocus=%1.2fSNR), normRatiosToStrongest=(%0.2fg,%0.2fs), log10(p)=(%7.1fcorr,%6.1fsignal,%6.1fcombined) for #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ]\n'...
           ]...
          ...,inInfo.searchStrategy.extendedFocus.absCorrSum(min(end,sL)), ils(inInfo.searchStrategy.extendedFocus.absCorrSum(min(end,sL)),1,'st',2,'nd',3,'rd','otherwise','th') ...           
            ,bestScoresSeenSoFar.highest_signatureAbsMean2D.signatureSizeByCorrSum4G, bestScoresSeenSoFar.highest_signatureAbsMean2D.signatureSizeByCorrSum4P, bestScoresSeenSoFar.highest_signatureAbsMean2D.signatureCorrInExtendedFocus...
            ,bestScoresSeenSoFar.highest_signatureAbsMean2D.signatureAbsMean2D, bestScoresSeenSoFar.highest_signatureAbsMean2D.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ....signaturePower/eDState.initialSignal.signalStrengthStatistics.signalPower ... 
            ,[bestScoresSeenSoFar.highest_signatureAbsMean2D.norm4sampleAxis_origUnits/max(alreadyDetected.norm4sampleAxis_origUnits), bestScoresSeenSoFar.highest_signatureAbsMean2D.norm4geneAxis_origUnits/max(alreadyDetected.norm4geneAxis_origUnits)] ...
            ,bestScoresSeenSoFar.highest_signatureAbsMean2D.log10_p4Correlations, bestScoresSeenSoFar.highest_signatureAbsMean2D.log10_p4SignalStrength, bestScoresSeenSoFar.highest_signatureAbsMean2D.log10_p...
            ,bestScoresSeenSoFar.highest_signatureAbsMean2D.cursor, iif(bestScoresSeenSoFar.highest_signatureAbsMean2D.imj>0, 'gene', 'sample'), abs(bestScoresSeenSoFar.highest_signatureAbsMean2D.imj), bestScoresSeenSoFar.highest_signatureAbsMean2D.label ...
          ...,inInfo.searchStrategy.extendedFocus.absCorrSum(min(end,sL)), ils(inInfo.searchStrategy.extendedFocus.absCorrSum(min(end,sL)),1,'st',2,'nd',3,'rd','otherwise','th') ...           
            ,bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.signatureSizeByCorrSum4G, bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.signatureSizeByCorrSum4P, bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.signatureCorrInExtendedFocus...
            ,bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.signatureAbsMean2D, bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ... .signaturePower/eDState.initialSignal.signalStrengthStatistics.signalPower  ... 
            ,[bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.norm4sampleAxis_origUnits/max(alreadyDetected.norm4sampleAxis_origUnits), bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.norm4geneAxis_origUnits/max(alreadyDetected.norm4geneAxis_origUnits)] ...
            ,bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.log10_p4Correlations, bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.log10_p4SignalStrength, bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.log10_p ...
            ,bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.cursor, iif(bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.imj>0, 'gene', 'sample'), abs(bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.imj), bestScoresSeenSoFar.highest_signatureSizeByCorrSum2D.label ...
          ...,inInfo.searchStrategy.extendedFocus.absCorrSum(min(end,sL)), ils(inInfo.searchStrategy.extendedFocus.absCorrSum(min(end,sL)),1,'st',2,'nd',3,'rd','otherwise','th') ...           
            ,bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.signatureSizeByCorrSum4G, bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.signatureSizeByCorrSum4P, bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.signatureCorrInExtendedFocus...
            ,bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.signatureAbsMean2D, bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ... .signaturePower/eDState.initialSignal.signalStrengthStatistics.signalPower  ... 
            ,[bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.norm4sampleAxis_origUnits/max(alreadyDetected.norm4sampleAxis_origUnits), bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.norm4geneAxis_origUnits/max(alreadyDetected.norm4geneAxis_origUnits)] ...
            ,bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.log10_p4Correlations, bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.log10_p4SignalStrength, bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.log10_p ...
            ,bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.cursor, iif(bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.imj>0, 'gene', 'sample'), abs(bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.imj), bestScoresSeenSoFar.highest_signatureCorrInExtendedFocus.label ...
          ...,inInfo.searchStrategy.extendedFocus.absCorrSum(min(end,sL)), ils(inInfo.searchStrategy.extendedFocus.absCorrSum(min(end,sL)),1,'st',2,'nd',3,'rd','otherwise','th') ...           
            ,bestScoresSeenSoFar.lowest_log10_p.signatureSizeByCorrSum4G, bestScoresSeenSoFar.lowest_log10_p.signatureSizeByCorrSum4P, bestScoresSeenSoFar.lowest_log10_p.signatureCorrInExtendedFocus ...
            ,bestScoresSeenSoFar.lowest_log10_p.signatureAbsMean2D, bestScoresSeenSoFar.lowest_log10_p.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ....signaturePower/eDState.initialSignal.signalStrengthStatistics.signalPower ... 
            ,[bestScoresSeenSoFar.lowest_log10_p.norm4sampleAxis_origUnits/max(alreadyDetected.norm4sampleAxis_origUnits), bestScoresSeenSoFar.lowest_log10_p.norm4geneAxis_origUnits/max(alreadyDetected.norm4geneAxis_origUnits)] ...
            ,bestScoresSeenSoFar.lowest_log10_p.log10_p4Correlations, bestScoresSeenSoFar.lowest_log10_p.log10_p4SignalStrength, bestScoresSeenSoFar.lowest_log10_p.log10_p...
            ,bestScoresSeenSoFar.lowest_log10_p.cursor, iif(bestScoresSeenSoFar.lowest_log10_p.imj>0, 'gene', 'sample'), abs(bestScoresSeenSoFar.lowest_log10_p.imj), bestScoresSeenSoFar.lowest_log10_p.label ...
        );
      end
      %Cache logic: Performance/save candidate initial representatives computed in this detection iteration for later iterations:
        if(eDState.current.multiPassesLevel < inInfo.searchStrategy.nDefinedPasses)
          eDState.current.precompute_previousIterations = []; %never cache across detection passes (with possibly changing sensitivity/score configurations.
        end
        %<-but do allow caching when rescanning one last time with same config (no emergent correlations assumption).
      %Output:
        signatureInfo = [];
        if(inInfo.internal.bDevEnableInteractiveBreaks)
          interactiveCheckPoint();
        end
        return;

    %% Score and qualification functions:
      function localScores = localGreedyScores4ProcessingOrder(dim)
        fcnZScore = @(X)(X-BM.meanW(X,1))/(nanstd(X)+eps); %+eps to define 0/0=0.
        %gene row signature scores:
          if(dim==1)
            ZprecenteredVarianceScores = fcnZScore(sqrt(BM.uncenteredVar(eDState.current.L2Rs,2)));
            ZpurityScores = fcnZScore(nanmax(abs(eDState.current.sdL2Rs),[],2)); %signatures with their norm concentrated on fewer pure directions/genes are less likely noise or highly overlapped than those with their norm spread out over all dimensions. Because often signatures extend to some genes for which they are the only signature (no superposition/overlapping), these genes are good candidates for initial representatives. Still, false positives and single-dimension outliers will also have high maximum norm of the standardized singal. Therefore qualification and scoring is needed to further evaluate the candidates.
            localScores = max(ZprecenteredVarianceScores, ZpurityScores);
          end
        %sample column signature scores:
          if(dim==2)
            ZprecenteredVarianceScores = fcnZScore(sqrt(BM.uncenteredVar(eDState.current.L2Rs,1)));
            ZpurityScores = fcnZScore(nanmax(abs(eDState.current.sdL2Rs),[],1)); %signatures with their norm concentrated on fewer pure directions/genes are less likely noise or highly overlapped than those with their norm spread out over all dimensions. Because often signatures extend to some genes where tey are the only signature (no overlapping), they are good candidates for dissection. Still, false positives and single-dimension outliers will also have high maximum norm of the standardized singal. Therefore qualification and scoring is needed to further evaluate the candidates.
            localScores = max(ZprecenteredVarianceScores, ZpurityScores);
          end
      end
      function [bQualifiedAsInitialRepresenativeCandidate, bPerformanceExcludeInThisDetectionPass] = testQualificationAsInitialRepresenative(signature)
        %Check signature size and correlation strength:
          bTooSmall = ...
              signature.signatureSizeByCorrSum4G + signature.signatureSizeByCorrSum4P < inInfo.searchStrategy.qualification.minAbsCorrSumSum(min(end,sL)) - 10*eps(inInfo.searchStrategy.qualification.minAbsCorrSumSum(min(end,sL))) ... %-eps for synthethic nonNoise-data where numeric effects may cause "1<1"
           || signature.signatureSizeByCorrSum4G                                      < inInfo.searchStrategy.qualification.minAbsCorrSum4G(min(end,sL))  - 10*eps(inInfo.searchStrategy.qualification.minAbsCorrSum4G(min(end,sL))) ...  %-eps for synthethic nonNoise-data where numeric effects may cause "1<1"
           || signature.signatureSizeByCorrSum4P                                      < inInfo.searchStrategy.qualification.minAbsCorrSum4P(min(end,sL))  - 10*eps(inInfo.searchStrategy.qualification.minAbsCorrSum4P(min(end,sL))) ...  %-eps for synthethic nonNoise-data where numeric effects may cause "1<1"
          ;
          %If signatures are narrow with respect to genes xor samples (stripe-like signautre signals), require a higher ratio of member samples xor genes as minimum signature size:
            if(~bTooSmall) 
              minSizeRatioIfDualIsSmall = 0.075; %TODO: include this as parameter in inInfo.
              bTooSmall = ...
                  signature.signatureSizeByCorrSum4G < nG * max(0, minSizeRatioIfDualIsSmall - signature.signatureSizeByCorrSum4P/nP) ...
               || signature.signatureSizeByCorrSum4P < nP * max(0, minSizeRatioIfDualIsSmall - signature.signatureSizeByCorrSum4G/nG) ...
              ;
            end
        %Compute significance threshold flags:
          if(true)
            bCorrSimpleQuali       = isnan(inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))) || ...
              10^signature.log10_p4Correlations   <= inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))...
            ;
            bCorrBonferroniQuali   = isnan(inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))) || ...
              10^signature.log10_p4Correlations   <= inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))/(nG+nP)...
            ;
            
            bSignalSimpleQuali     = isnan(signature.log10_p4SignalStrength) || isnan(inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL))) || ... %signal strength p value may be NaN for the first detection iteration as we have no reliable noise distribution, yet.
              10^signature.log10_p4SignalStrength <= inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL)) ...
            ;
            bSignalBonferroniQuali = isnan(signature.log10_p4SignalStrength) || isnan(inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL))) || ... %signal strength p value may be NaN for the first detection iteration as we have no reliable noise distribution, yet.
              10^signature.log10_p4SignalStrength <= inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL))/(nG+nP)...
            ;
            
            bCombinedSimpleQuali     = isnan(inInfo.searchStrategy.qualification.alpha4combined(min(end,sL))) || ...
              10^signature.log10_p <= inInfo.searchStrategy.qualification.alpha4combined(min(end,sL))...
            ;
            bCombinedBonferroniQuali = isnan(inInfo.searchStrategy.qualification.alpha4combined(min(end,sL))) || ...
              10^signature.log10_p <= inInfo.searchStrategy.qualification.alpha4combined(min(end,sL))/(nG+nP)...
            ;
          end
        %Check minimum average correlation strength:
          bSufficientCorrelationStrength = signature.signatureCorrInExtendedFocus >= inInfo.searchStrategy.qualification.minCorrInExtendedFocus(min(end,sL));
        %Check minimum average signal strength:
          bSufficientsiAbsoluteSignalStrength = signature.signatureAbsMean2D >= inInfo.searchStrategy.qualification.minAbsMean2DInExtendedFocus(min(end,sL));

        %Overall qualification flag:
          bQualifiedAsInitialRepresenativeCandidate = ~bTooSmall && bSufficientCorrelationStrength && bSufficientsiAbsoluteSignalStrength && bCorrSimpleQuali && bSignalSimpleQuali && bCombinedSimpleQuali;
          if(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCorrP(min(end,sL)))
            bQualifiedAsInitialRepresenativeCandidate = bQualifiedAsInitialRepresenativeCandidate && bCorrBonferroniQuali;
          end
          if(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForSignalP(min(end,sL)))
            bQualifiedAsInitialRepresenativeCandidate = bQualifiedAsInitialRepresenativeCandidate && bSignalBonferroniQuali;
          end
          if(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCombinedP(min(end,sL)))
            bQualifiedAsInitialRepresenativeCandidate = bQualifiedAsInitialRepresenativeCandidate && bCombinedBonferroniQuali;
          end

        %Performance: do not revisit the genes/samples for later k in this detection pass, if it is clearly disqualified:
          disqualificationFactorForPerformanceExclusion = 1/2;
          bPerformanceExcludeInThisDetectionPass = ...
              signature.signatureSizeByCorrSum4G + signature.signatureSizeByCorrSum4P    < disqualificationFactorForPerformanceExclusion*inInfo.searchStrategy.qualification.minAbsCorrSumSum(min(end,sL)) ...
           || signature.signatureSizeByCorrSum4G                                         < disqualificationFactorForPerformanceExclusion*inInfo.searchStrategy.qualification.minAbsCorrSum4G(min(end,sL)) ... %are there enough correlated pixels supporting the signature?
           || signature.signatureSizeByCorrSum4P                                         < disqualificationFactorForPerformanceExclusion*inInfo.searchStrategy.qualification.minAbsCorrSum4P(min(end,sL)) ... %are there enough correlated pixels supporting the signature?               
           || signature.signatureCorrInExtendedFocus                                     < disqualificationFactorForPerformanceExclusion*inInfo.searchStrategy.qualification.minCorrInExtendedFocus(min(end,sL)) ... %is the correlation strong enoug?
           || signature.log10_p4SignalStrength                                           >=disqualificationFactorForPerformanceExclusion*log10( inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL)) * iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCorrP(min(end,sL)),1/(nG+nP),1) ) ...
           || signature.log10_p4Correlations                                             >=disqualificationFactorForPerformanceExclusion*log10( inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL)) * iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForSignalP(min(end,sL)),1/(nG+nP),1) ) ...
           || signature.log10_p                                                          >=disqualificationFactorForPerformanceExclusion*log10( inInfo.searchStrategy.qualification.alpha4combined(min(end,sL)) * iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCombinedP(min(end,sL)),1/(nG+nP),1) ) ...
          ;
      end
  end

