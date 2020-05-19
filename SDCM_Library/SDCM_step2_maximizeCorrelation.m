%ABSTRACT
% Subfunction for SDCM step 2: Generalize from the initial representative 
% by maximizing the signature functional until convergence of signature axes.

  function generalizedSignature = SDCM_step2_maximizeCorrelation(initialSignature, inInfo)
    %% Initialize:
      global eDState;
      %Shortcuts:
        nG = size(eDState.current.L2Rs,1);
        nP = size(eDState.current.L2Rs,2);            
        sL = min(eDState.current.multiPassesLevel, inInfo.searchStrategy.nDefinedPasses);
      %Get basic math function handles:
        BM = getBasicMathFunctions(inInfo.preprocessing.bDataContainsNaNs);
      %Performance/precalculate constants:
        if(true)
          rowLabels = inInfo.reference.rowIDs;
            rowLabels = cellfun(@(c)iif(isnumeric(c),num2str(c),c), rowLabels, 'UniformOutput',false);
            rowLabels = cellfun(@(s)strrep(s,'\color[rgb]{0 0 0.67}',''), rowLabels, 'UniformOutput',false); %hotfix: do not display tex color formatting on console.
          colLabels = inInfo.reference.colIDs;
            colLabels = cellfun(@(c)iif(isnumeric(c),num2str(c),c), colLabels, 'UniformOutput',false);
            colLabels = cellfun(@(s)strrep(s,'\color[rgb]{0 0 0.67}',''), colLabels, 'UniformOutput',false); %hotfix: do not display tex color formatting on console.
          maxLabelLength = max(max(cellfun(@length,rowLabels)),max(cellfun(@length,colLabels))); %for status message alignment.
        end
      %Initialize generalizedSignature by starting with the initial representative from step 1 (single gene or single sample)
        if(true)
          %Initialize with the initial representative and its correlations:
            generalizedSignature = eDState.signatureTemplate;
            generalizedSignature.sampleAxis_origUnits = initialSignature.sampleAxis_origUnits;
            generalizedSignature.geneAxis_origUnits = initialSignature.geneAxis_origUnits;
            generalizedSignature.R4G = initialSignature.R4G;
            generalizedSignature.R4P = initialSignature.R4P;
            generalizedSignature.signedFocusedW4P = initialSignature.signedFocusedW4P;
            generalizedSignature.signedFocusedW4G = initialSignature.signedFocusedW4G;
            generalizedSignature.signedFocusedW4P_withPerpendicularSpace = initialSignature.signedFocusedW4P_withPerpendicularSpace;
            generalizedSignature.signedFocusedW4G_withPerpendicularSpace = initialSignature.signedFocusedW4G_withPerpendicularSpace;
            generalizedSignature.sampleSizes4R4G = initialSignature.sampleSizes4R4G;
            generalizedSignature.sampleSizes4R4P = initialSignature.sampleSizes4R4P;
          %Add the initial representative from step 1 as first member of the generalizedSignature:
            %precomputeCandidatesInParallel([], inInfo); %reset the precompute cache.
            %precomputeCandidatesInParallel(initialSignature.imj, inInfo);
            %since we cache all step 1 precomputation results, the initial representative should already be in the cache:
              assert(ismember(initialSignature.imj, eDState.current.precompute.ImJ), 'code validation: initial representative not present in eDState.current.precompute.ImJ');
          %Add initial representative as first member and accumulate all members: 
            nextMember = getSignatureAxesAndScores(initialSignature.imj);
            generalizedSignatureb4Accu = generalizedSignature;
            generalizedSignature = addMemberAndAccumulate(...
               generalizedSignature, initialSignature.imj...
              ,nextMember.sampleAxis_origUnits, nextMember.geneAxis_origUnits...
            );
          %Define the performance subspace:
            BsPerformanceSubspace4Correlations = selectDimensionsInPerformanceSubspace(nG,nP ...
             ,generalizedSignature.R4G, generalizedSignature.R4P ...
             ,eDState.noiseEstimation.log10PNoise4G,eDState.noiseEstimation.log10PNoise4P ...
             ,inInfo.searchStrategy ...
            );
          %Correlate the new generalized signatures with all genes respectively samples:
            generalizedSignature = correlationConvergenceStep(generalizedSignature); 
          %Compute combined score for selecting next accumulation members:
            initialCombinedScore = signatureFunctional(generalizedSignature);
            lastAcceptedCombinedScores = initialCombinedScore;
        end
      %Initialize .plots.focusConvergence:
        global sConvergenceFigure;
        if(inInfo.plots.focusConvergence.bEnabled)
          sConvergenceFigure.f = figure(...
             'Position',[25 295 1897 685]... 
            ,'Name', sprintf(...
               '%03d, signature focus convergence, initial corr=(%2.1fg,%2.1fs,%0.2fr); p=(%1.0ecorr,%1.0esignal)'...
               ,eDState.current.k ...
               ...Signature size estimates and average absolute correlation in the focus:
                   ,initialSignature.signatureSizeByCorrSum4G, initialSignature.signatureSizeByCorrSum4P...
                   ,initialSignature.signatureCorrInExtendedFocus ...
               ...Signature significance:
                   ,10^initialSignature.log10_p4Correlations, 10^initialSignature.log10_p4SignalStrength ...
             )...
            ,'Color','w'...
            ,'Visible',iif(inInfo.plots.bVisible,'on','off')...
          );
            sConvergenceFigure.a4CorrelationConvergence = subplot(3,3,1); set(sConvergenceFigure.a4CorrelationConvergence, 'YScale', 'log');
            sConvergenceFigure.a4SignatureConvergence = subplot(3,3,4); set(sConvergenceFigure.a4SignatureConvergence, 'YScale', 'log');
            sConvergenceFigure.a4R4P = subplot(3,3,2);
            sConvergenceFigure.a4R4G = subplot(3,3,3);
            sConvergenceFigure.a4sampleAxis = subplot(3,3,5);
            sConvergenceFigure.a4geneAxis = subplot(3,3,6);
            box(sConvergenceFigure.a4CorrelationConvergence, 'on');
            box(sConvergenceFigure.a4SignatureConvergence, 'on');
            box(sConvergenceFigure.a4R4P, 'on');
            box(sConvergenceFigure.a4R4G, 'on');
            box(sConvergenceFigure.a4sampleAxis, 'on');
            box(sConvergenceFigure.a4geneAxis, 'on');
            hold(sConvergenceFigure.a4CorrelationConvergence, 'on');
            hold(sConvergenceFigure.a4SignatureConvergence, 'on');
            hold(sConvergenceFigure.a4R4P, 'on');
            hold(sConvergenceFigure.a4R4G, 'on');
            hold(sConvergenceFigure.a4sampleAxis, 'on');
            hold(sConvergenceFigure.a4geneAxis, 'on');
            xlim(sConvergenceFigure.a4R4P,[1,nP]);
            xlim(sConvergenceFigure.a4R4G,[1,nG]);
            xlim(sConvergenceFigure.a4sampleAxis,[1,nP]);
            xlim(sConvergenceFigure.a4geneAxis,[1,nG]);
            grid(sConvergenceFigure.a4CorrelationConvergence, 'on');
            grid(sConvergenceFigure.a4SignatureConvergence, 'on');
            grid(sConvergenceFigure.a4R4P, 'on');
            grid(sConvergenceFigure.a4R4G, 'on');
            grid(sConvergenceFigure.a4sampleAxis, 'on');
            grid(sConvergenceFigure.a4geneAxis, 'on');
            title(sConvergenceFigure.a4CorrelationConvergence, sprintf('correlation detlas to the previous convegence iteration'));
            title(sConvergenceFigure.a4SignatureConvergence, sprintf('signature deltas to the previous convegence iteration'));
            title(sConvergenceFigure.a4R4P, sprintf('r_P(samples,a_G%s) (also defining the sample presentation order)',iif(false,'^S','')));
            title(sConvergenceFigure.a4R4G, sprintf('r_G(genes,a_P%s) (also defining the gene presentation order)',iif(false,'^S','')));
            title(sConvergenceFigure.a4sampleAxis,   sprintf('a_P%s (= sample signal strengths wrt. the signature, i.e. the sample axis)',iif(false,'^S','')));
            title(sConvergenceFigure.a4geneAxis, sprintf('a_G%s (= genes signal strengths wrt. the signature, i.e. the gene axis)',iif(false,'^S','')));
            xlabel(sConvergenceFigure.a4CorrelationConvergence, 'signature member #');
            xlabel(sConvergenceFigure.a4SignatureConvergence, 'signature member #');
            xlabel(sConvergenceFigure.a4R4P, 'current');
            xlabel(sConvergenceFigure.a4R4G, 'current');
            xlabel(sConvergenceFigure.a4sampleAxis, 'current');
            xlabel(sConvergenceFigure.a4geneAxis, 'current');
            ylabel(sConvergenceFigure.a4CorrelationConvergence, sprintf('1 - r_{focused}(previous, current)\nfor gene and sample correlations'));
            ylabel(sConvergenceFigure.a4SignatureConvergence, sprintf('1 - r_{focused}(previous, current)\nfor gene and gene axes'));
            ylabel(sConvergenceFigure.a4R4P, 'previous');
            ylabel(sConvergenceFigure.a4R4G, 'previous');
            ylabel(sConvergenceFigure.a4sampleAxis, 'previous');
            ylabel(sConvergenceFigure.a4geneAxis, 'previous');
          sConvergenceFigure.caAccepted = {};
          sConvergenceFigure.hInitialSignatureText = [];
          sConvergenceFigure.hAcceptedText = [];
          sConvergenceFigure.HmembersOfgeneAxis = [];
          sConvergenceFigure.HmembersOfsampleAxis = [];
          %Plot configured convergence thresholds into .a4SignatureConvergence:
            line([0,1000]...
              ,inInfo.correlationMaximization.convergenceEpsilon4deltaCorr4Rs(min(end,sL))*[1,1]...
              ,'LineStyle','-.','Color',[0 1 0],'LineWidth',1.5 ...
              ,'Parent',sConvergenceFigure.a4CorrelationConvergence...
            );
            line([0,1000]...
              ,inInfo.correlationMaximization.convergenceEpsilon4deltaCorr4AxesInNoiseSDs(min(end,sL))*eDState.noiseEstimation.sd4epsilon*[1,1]...
              ...
              ,'LineStyle','--','Color',[0 1 0],'LineWidth',1.5 ...
              ,'Parent',sConvergenceFigure.a4SignatureConvergence...
            );
        end
        if(inInfo.plots.focusConvergence.bEnabled) %plot initial representative.
          %Set plots from previous focusing iterations to gray:
            set(findobj(sConvergenceFigure.a4R4P,'Type','line'),'Color',[.5 .5 .5]);%,'LineStyle','-');
            set(findobj(sConvergenceFigure.a4R4G,'Type','line'),'Color',[.5 .5 .5]);%,'LineStyle','-');
            set(findobj(sConvergenceFigure.a4sampleAxis,'Type','line'),'Color',[.5 .5 .5]); %,'LineStyle','-');
            set(findobj(sConvergenceFigure.a4geneAxis,'Type','line'),'Color',[.5 .5 .5]); %,'LineStyle','-');
          plot(sConvergenceFigure.a4R4G, generalizedSignature.R4G, generalizedSignatureb4Accu.R4G, '.','Color',[0 0 0.85]); 
          plot(sConvergenceFigure.a4R4P, generalizedSignature.R4P, generalizedSignatureb4Accu.R4P, 'r.'); 
            %<-Note: dass previous==0 bei erster Aktualisierung liegt daran, dass das die in Step1 vie BsPerformanceSubspace ausgeschlossenen Punkte sind, die erst in Step2 eine Korrelation zu den aktuellen Effektachsen erhalten.
          plot(sConvergenceFigure.a4geneAxis, generalizedSignature.geneAxis_origUnits, generalizedSignatureb4Accu.geneAxis_origUnits, '.','Color',[0 0 0.85]);
          plot(sConvergenceFigure.a4sampleAxis, generalizedSignature.sampleAxis_origUnits, generalizedSignatureb4Accu.sampleAxis_origUnits, 'r.');
          %tighten axes:
            ylim(sConvergenceFigure.a4R4G,1.03*[-1,1]);
            ylim(sConvergenceFigure.a4R4P,1.03*[-1,1]);
            ylim(sConvergenceFigure.a4geneAxis, 1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.geneAxis_origUnits)));
            ylim(sConvergenceFigure.a4sampleAxis,   1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.sampleAxis_origUnits)));
            xlim(sConvergenceFigure.a4R4G,1.03*[-1,1]);
            xlim(sConvergenceFigure.a4R4P,1.03*[-1,1]);
            xlim(sConvergenceFigure.a4geneAxis, 1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.geneAxis_origUnits)));
            xlim(sConvergenceFigure.a4sampleAxis,   1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.sampleAxis_origUnits)));
          %Write initial representative as title (every iteration, since the ylim might change):
            if(ishandle(sConvergenceFigure.hInitialSignatureText))
              delete(sConvergenceFigure.hInitialSignatureText);
            end
            sConvergenceFigure.hInitialSignatureText = text(...
               min(xlim(sConvergenceFigure.a4R4P))-abs(diff(xlim(sConvergenceFigure.a4R4P)))*(1/4+1/4+1+1/4) ...
              ,max(ylim(sConvergenceFigure.a4R4P))+abs(diff(ylim(sConvergenceFigure.a4R4P)))/2.5 ...
              ,initialSignature.sInitialRepresenative ...
              ,'HorizontalAlignment','left', 'VerticalAlignment','top' ....
              ,'Color',[0 0 1]...
              ,'FontName',iif(ispc(),'Consolas','Courier') ...
              ,'FontWeight','normal' ...,'FontSize',8.5  ...
              ,'FontUnits','pixels', 'FontSize',13 ...
              ,'Parent',sConvergenceFigure.a4R4P...
            );
          drawnow;
        end      
    %% Successively maximize correlation and generalize the signature via incorporating genes/samples that optimally increase the signature functional:
      nMembers4PlotUpdateInterval = 5; nMembersInLastUpdate = -Inf;
      while(true)
        %Look ahead for the next best gene/sample to add to the generalized signature:
          %Initialize:
            lookaheadCursor = 0;
            %Determine processing order of gene and sample candidates for the lookahead loop based on their correlations with the current generalizedSignature and exclude already added members:
              if(true)
                %Configure scores for presenting genes or samples during signature generalization based on their R4G resp. R4P and combine row and column scores in a single processing order:
                  processingOrder.scores = [];
                  processingOrder.signatureI = [];
                  if(ismember('sampleAxis',inInfo.searchStrategy.allowedSignatureTypes))
                    %Presorting scores: require high correlation to the current generalized sampleAxis and high SD signal in the current generalized sampleAxis to become the top next gene member candidate .
                      focusedGreedyScores4G = abs(generalizedSignature.R4G);
                    processingOrder.scores = [processingOrder.scores; focusedGreedyScores4G];
                    processingOrder.signatureI = [processingOrder.signatureI, 1:nG];
                  end
                  if(ismember('geneAxis',inInfo.searchStrategy.allowedSignatureTypes))
                    %Presorting scores: require high correlation to the current generalized geneAxis and high SD signal in the current generalized geneAxis from the top next sample member candidate.
                      focusedGreedyScores4P = abs(generalizedSignature.R4P');
                    processingOrder.scores = [processingOrder.scores; focusedGreedyScores4P];
                    processingOrder.signatureI = [processingOrder.signatureI, -(1:nP)];                  
                  end
                %Remove those with too low relative correlation (we never want to accumulate these):
                  BLin = [BsPerformanceSubspace4Correlations{1}; BsPerformanceSubspace4Correlations{2}']; %the performance subspace and accu candidate qualification are tied now.
                  processingOrder.scores = processingOrder.scores(BLin);
                  processingOrder.signatureI = processingOrder.signatureI(BLin);
                %Remove genes/samples that are already part of the generalized signature:
                  BDelete = ismember(processingOrder.signatureI, generalizedSignature.corrMaximizingMembers.ImJ);
                  processingOrder.scores(BDelete) = [];
                  processingOrder.signatureI(BDelete) = [];
                %Performance: Boost the scores of already precomputed candidates a little (for large signatures there are often genes/samples with nearly identical scores; in this case it is sufficient to not use the best candidates for accumulation, but those with nearly equal scores that have already been precomputed):
                  processingOrder.scores(ismember(processingOrder.signatureI, eDState.current.precompute.ImJ)) = ...
                      inInfo.correlationMaximization.performance.scoreBoostFactore4alreadyPrecomputedMemberCandidates(min(end,sL)) ...
                    * processingOrder.scores(ismember(processingOrder.signatureI, eDState.current.precompute.ImJ)) ...
                  ;
                %Define lookahead processing order for signature generalization by descending relative correlations:
                  [~,processingOrder.SII] = sort(ilv(processingOrder.scores, @(X)iif(isnan(X),0,X)),'descend');
                  processingOrder.SImJ = processingOrder.signatureI(processingOrder.SII);
              end
            %Reset candidate accumulations for the previous generalizedSignature:
              signatureTemplate = eDState.signatureTemplate; signatureTemplate.corrMaximizingMembers = struct();
              eDState.current.precompute.candidateAccus = repmat(signatureTemplate, length(eDState.current.precompute.ImJ), 1);
              eDState.current.precompute.BcandidateAccusUpToDate = false(size(eDState.current.precompute.candidateAccus));
            %Reset lookahead knowledge:
              if(true)
                lookahead4bestScoreIncrease.ImJ = [];
                lookahead4bestScoreIncrease.signatureSizeByCorrSum4GRatios = [];
                lookahead4bestScoreIncrease.signatureSizeByCorrSum4PRatios = [];
                lookahead4bestScoreIncrease.signatureSizeByCorrSum4G = [];
                lookahead4bestScoreIncrease.signatureSizeByCorrSum4P = [];
                lookahead4bestScoreIncrease.signatureSizeByCorrSum2Ds = [];
                lookahead4bestScoreIncrease.signatureSharpnesses = []; 
                lookahead4bestScoreIncrease.signatureCorrInExtendedFocuss = [];
                lookahead4bestScoreIncrease.log10_ps = [];
                lookahead4bestScoreIncrease.SD4sampleAxis_origUnits = [];
                lookahead4bestScoreIncrease.SD4geneAxis_origUnits = [];
                lookahead4bestScoreIncrease.r2currentGeneralizedGeneAxis = [];
                lookahead4bestScoreIncrease.r2currentGeneralizedSampleAxis = [];
                lookahead4bestScoreIncrease.signatureAbsMean2Ds4accumulatedSignatures = [];
                lookahead4bestScoreIncrease.sampleAxisMeanAbsRatio = [];
                lookahead4bestScoreIncrease.geneAxisMeanAbsRatio = [];

                lookahead4bestScoreIncrease.combinedScore = [];
              end
            %Performance: Clear entries from the precompute cache that are unrelated to the current focus (Note: the cache entries for the next STEP1 are saved in eDState.current.precompute_previousIterations already):
              precomputeCandidatesInParallel(setdiff(eDState.current.precompute.ImJ, processingOrder.SImJ), 'clearFromCache');
          bAlreadyFoundASufficientScore = false;
          while( lookaheadCursor < length(processingOrder.SImJ) && (...
                 lookaheadCursor < inInfo.correlationMaximization.maxLookahead4NextBestAccuCandidate(min(end,sL)) ... %only relevant, if no score increase can be found
              || isempty(lookahead4bestScoreIncrease.ImJ) && lookaheadCursor < inInfo.correlationMaximization.maxLookahead4NextBestAccuCandidateFactorInCaseOfZeroCandidates(min(end,sL))*inInfo.correlationMaximization.maxLookahead4NextBestAccuCandidate(min(end,sL)) ... %if we did not find any qualifying candidate, yet, allow the doubled .maxLookahead4NextBestAccuCandidate temporarily.
          ))
            lookaheadCursor = lookaheadCursor + 1;
            imj = processingOrder.SImJ(lookaheadCursor);
            %Performance: parallel precomputation of the next required candidate members and their candidate-accumulation to the current generalizedSignature:
              if(~ismember(imj, eDState.current.precompute.ImJ)) %if not already in the cache:
                %Define ImJ to precalc: %<-Note: keep the performance/skip config in sync with the lookahead loop.
                  precalcChunkSize = min(inInfo.correlationMaximization.maxLookahead4NextBestAccuCandidate(min(end,sL))-lookaheadCursor+1, ceil(sum(eDState.current.precompute.BcandidateAccusUpToDate)/2)); %increase cahce by the amount still needed for the lookahead, but never more then half of its size.
                  nWorkers = max(1,inInfo.correlationMaximization.performance.nMinPrecomputeSignaturesPerWorker(min(end,sL))*ilv(gcp('nocreate'),@(pool)iif(isempty(pool),0,pool.NumWorkers))); %max(1,x) since we have at least our own instance as worker; 2* to not have too large overhad.
                  precalcChunkSize = max(nWorkers, precalcChunkSize); %in anticipation of more members iterations and similar presort order, we compute a bit more than needed to fill the cache, if the request would otherwise be too few (i.e. much overhead percentage)
                  nextNeededImJ = processingOrder.SImJ(lookaheadCursor:end);
                  nextNeededImJ(ismember(nextNeededImJ,eDState.current.precompute.ImJ)) = [];
                  nextNeededImJ = nextNeededImJ(1:min(precalcChunkSize,end));
                  assert(ismember(imj,nextNeededImJ), 'code validation: current imj is not in nextNeededImJ');
                %Precalc next initial representative candidates, keeping already precomputed ones:
                  SDCM_printStatus(2,'   - Precomputing correlations with the next %d member candidates in parallel (performance subspace = [%d/%d,%d/%d])...\n'...
                    ,length(nextNeededImJ)...
                    ,sum(BsPerformanceSubspace4Correlations{1}),length(BsPerformanceSubspace4Correlations{1})...
                    ,sum(BsPerformanceSubspace4Correlations{2}),length(BsPerformanceSubspace4Correlations{2})...
                  );
                  nextI = length(eDState.current.precompute.ImJ)+(1:length(nextNeededImJ));
                  precomputeCandidatesInParallel([eDState.current.precompute.ImJ,nextNeededImJ], inInfo, BsPerformanceSubspace4Correlations);
                %Get the precomputed gene and gene axes that have to be accumulated with the present generalizedSignature members for the next candidateAccus:
                  nextMembers.sampleAxes_origUnits = nan(length(nextNeededImJ),nP);
                  nextMembers.geneAxes_origUnits = nan(nG,length(nextNeededImJ));
                  nextMemberCursor = 0;
                  for k=1:length(nextNeededImJ)
                    %Get candidate:
                      nextMember = getSignatureAxesAndScores(nextNeededImJ(k));
                    %Store next member candidates:
                      nextMemberCursor = nextMemberCursor+1;
                      nextMembers.sampleAxes_origUnits(nextMemberCursor,:) = nextMember.sampleAxis_origUnits;
                      nextMembers.geneAxes_origUnits(:,nextMemberCursor) = nextMember.geneAxis_origUnits;
                  end
                  nextMembers.sampleAxes_origUnits = nextMembers.sampleAxes_origUnits(1:nextMemberCursor,:);
                  nextMembers.geneAxes_origUnits = nextMembers.geneAxes_origUnits(:,1:nextMemberCursor);
                %Accumulate all members of the next candidateAccus: 
                  eDState.current.precompute.candidateAccus(nextI) = addMemberAndAccumulate(...
                    generalizedSignature ...
                   ,nextNeededImJ ...
                   ,nextMembers.sampleAxes_origUnits, nextMembers.geneAxes_origUnits ...
                  );
                %Correlate the new generalized signatures with all genes respectively samples:
                  eDState.current.precompute.candidateAccus(nextI) = correlationConvergenceStep(...
                    eDState.current.precompute.candidateAccus(nextI) ...
                  ); 
                  eDState.current.precompute.BcandidateAccusUpToDate(nextI) = true;
              end
            %Performance: parallel update of the next required candidate accumulations and convergence steps for cached signatures from previous iterations:
              if(~eDState.current.precompute.BcandidateAccusUpToDate(eDState.current.precompute.ImJ==imj))
                %Define ImJ to precalc: %<-Note: keep the performance/skip config in sync with the lookahead loop.
                  nWorkers = max(1,inInfo.correlationMaximization.performance.nMinPrecomputeSignaturesPerWorker(min(end,sL))*ilv(gcp('nocreate'),@(pool)iif(isempty(pool),0,pool.NumWorkers))); %max(1, since we have at least our own instance as worker; 4* to not have too large overhad.
                  precalcChunkSize = max(nWorkers, ceil(sum(eDState.current.precompute.BcandidateAccusUpToDate)/2)); %precompute (maximally) half of the cache, but minimally nWorkers
                  precalcChunkSize = max(precalcChunkSize, ceil(inInfo.correlationMaximization.maxLookahead4NextBestAccuCandidate(min(end,sL))/2)); %pre-accumulation is relatively fast and should be done at maximum two times per member loop.
                  if(lookaheadCursor<inInfo.correlationMaximization.maxLookahead4NextBestAccuCandidate(min(end,sL)))
                    precalcChunkSize = min(precalcChunkSize, inInfo.correlationMaximization.maxLookahead4NextBestAccuCandidate(min(end,sL))-lookaheadCursor+1); %never increase cache by more than needed.
                  end
                  if(ilv(gcp('nocreate'),@(pool)iif(isempty(pool),0,pool.NumWorkers))>0)
                    precalcChunkSize = ceil(precalcChunkSize/ilv(gcp('nocreate'),@(pool)iif(isempty(pool),0,pool.NumWorkers)))*ilv(gcp('nocreate'),@(pool)iif(isempty(pool),0,pool.NumWorkers)); %use all workers.
                  end
                  %precalcChunkSize = 1; %performance: parallelization disabled here due to overhead.
                    %<-for large data the overhead without parallelization is much higher!
                  ImJ2update = processingOrder.SImJ(lookaheadCursor:end);
                  ImJ2update(ismember(ImJ2update,eDState.current.precompute.ImJ(eDState.current.precompute.BcandidateAccusUpToDate))) = [];
                  ImJ2update(~ismember(ImJ2update,eDState.current.precompute.ImJ)) = []; %must precompute/cache the signature before it can be used as candidate member for accumulation.
                  ImJ2update = ImJ2update(1:min(precalcChunkSize,end));
                  assert(ismember(imj,ImJ2update), 'code validation: current imj is not in ImJ2update');
                %Get the precomputed gene and gene axes that have to be accumulated with the present generalizedSignature members for the next candidateAccus:
                  membersCandidates.sampleAxes_origUnits = nan(length(ImJ2update),nP);
                  membersCandidates.geneAxes_origUnits = nan(nG,length(ImJ2update));
                  II2update = nan(size(ImJ2update));
                  nextMemberCursor = 0;
                  for k=1:length(ImJ2update)
                    %Get candidate:
                      nextMember = getSignatureAxesAndScores(ImJ2update(k));
                    %Store next member candidates:
                      nextMemberCursor = nextMemberCursor+1;
                      membersCandidates.sampleAxes_origUnits(nextMemberCursor,:) = nextMember.sampleAxis_origUnits;
                      membersCandidates.geneAxes_origUnits(:,nextMemberCursor) = nextMember.geneAxis_origUnits;
                    %Index in candidateAccus:
                      II2update(k) = find(eDState.current.precompute.ImJ==ImJ2update(k));
                  end
                  membersCandidates.sampleAxes_origUnits = membersCandidates.sampleAxes_origUnits(1:nextMemberCursor,:);
                  membersCandidates.geneAxes_origUnits = membersCandidates.geneAxes_origUnits(:,1:nextMemberCursor);
                %Accumulate all members of the next candidateAccus: 
                  eDState.current.precompute.candidateAccus(II2update) = addMemberAndAccumulate(...
                    generalizedSignature ...
                   ,ImJ2update ...
                   ,membersCandidates.sampleAxes_origUnits, membersCandidates.geneAxes_origUnits ...
                  );
                %Correlate the new generalized signatures with all genes respectively samples:
                  eDState.current.precompute.candidateAccus(II2update) = correlationConvergenceStep(...
                    eDState.current.precompute.candidateAccus(II2update) ...
                  ); 
                  eDState.current.precompute.BcandidateAccusUpToDate(II2update) = true;
              end
            %Get the precomputed candidateAccu:
              candidateAccu = eDState.current.precompute.candidateAccus(eDState.current.precompute.ImJ==imj);
            %Qualification:
              if(candidateAccu.signatureSizeByCorrSum2D==0) continue; end %special case: zero focus candidate.
              %Note: no further qualification here, as correlation maximization is the only criterium.
            %Record score increase/decrease:
              if(true)
                lookahead4bestScoreIncrease.ImJ(end+1) = imj;
                  assert(length(unique(lookahead4bestScoreIncrease.ImJ))==length(lookahead4bestScoreIncrease.ImJ), 'code validation: duplicate imj in lookahead4bestScoreIncrease.ImJ');
                lookahead4bestScoreIncrease.signatureSizeByCorrSum4G(end+1) = candidateAccu.signatureSizeByCorrSum4G;
                lookahead4bestScoreIncrease.signatureSizeByCorrSum4P(end+1) = candidateAccu.signatureSizeByCorrSum4P;
                lookahead4bestScoreIncrease.signatureSizeByCorrSum4GRatios(end+1) = candidateAccu.signatureSizeByCorrSum4G / generalizedSignature.signatureSizeByCorrSum4G;
                lookahead4bestScoreIncrease.signatureSizeByCorrSum4PRatios(end+1) = candidateAccu.signatureSizeByCorrSum4P / generalizedSignature.signatureSizeByCorrSum4P;
                lookahead4bestScoreIncrease.signatureSizeByCorrSum2Ds(end+1) = candidateAccu.signatureSizeByCorrSum2D;
                lookahead4bestScoreIncrease.signatureCorrInExtendedFocuss(end+1) = candidateAccu.signatureCorrInExtendedFocus;
                lookahead4bestScoreIncrease.log10_ps(end+1) = candidateAccu.log10_p;
                lookahead4bestScoreIncrease.SD4sampleAxis_origUnits(end+1) = nanstd(candidateAccu.sampleAxis_origUnits);
                lookahead4bestScoreIncrease.SD4geneAxis_origUnits(end+1) = nanstd(candidateAccu.geneAxis_origUnits);
                lookahead4bestScoreIncrease.r2currentGeneralizedGeneAxis(end+1) = candidateAccu.corrMaximizingMembers.R4sampleAxisMembers(end);
                lookahead4bestScoreIncrease.r2currentGeneralizedSampleAxis(end+1) = candidateAccu.corrMaximizingMembers.R4geneAxisMembers(end);
                lookahead4bestScoreIncrease.signatureAbsMean2Ds4accumulatedSignatures(end+1) = candidateAccu.signatureAbsMean2D;
                lookahead4bestScoreIncrease.sampleAxisMeanAbsRatio(end+1) = BM.meanW(abs(candidateAccu.sampleAxis_origUnits))/BM.meanW(abs(generalizedSignature.sampleAxis_origUnits));
                lookahead4bestScoreIncrease.geneAxisMeanAbsRatio(end+1) = BM.meanW(abs(candidateAccu.geneAxis_origUnits))/BM.meanW(abs(generalizedSignature.geneAxis_origUnits));
                
                lookahead4bestScoreIncrease.combinedScore(end+1) = signatureFunctional(candidateAccu);
              end
              %Status output: accumulation candidate and its scores:
                if(inInfo.export.nStatusOutputLevel>=4) %performance: only compute the status meassage if requested.
                  SDCM_printStatus(4 ...
                   ,['   - processing candidate            #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ] with (r2g=%+0.2f, r2p=%+0.2f)  => refocused corr=(% 7.1fg,% 7.1fs, %0.2fr)->(% 7.1fg,% 7.1fs, %0.2fr)\n'...not computed for performance reasons/; signature/noise=%0.1f->%0.1f, log10(p)=%5.1f->%5.1f\n'...
                    ] ...
                   ...New member indices and name:
                       ,lookaheadCursor ...
                       ,iif(imj>0, 'gene', 'sample') ...
                       ,iif(imj>0, imj, -imj) ...
                       ,iif(imj>0, @()rowLabels{imj}, @()colLabels{-imj}) ...
                   ...Scores of the new member:
                       ,lookahead4bestScoreIncrease.r2currentGeneralizedGeneAxis(end)...
                       ,lookahead4bestScoreIncrease.r2currentGeneralizedSampleAxis(end)...
                   ...Score changed for the new focus:
                       ,generalizedSignature.signatureSizeByCorrSum4G, generalizedSignature.signatureSizeByCorrSum4P ...
                       ,generalizedSignature.signatureCorrInExtendedFocus...
                       ,candidateAccu.signatureSizeByCorrSum4G, candidateAccu.signatureSizeByCorrSum4P ...
                       ,candidateAccu.signatureCorrInExtendedFocus...
                  );
                end
            %Performance: If we have a score increase and reached the minimum lookahead width, break to directly accept it:
              bAlreadyFoundASufficientScore = bAlreadyFoundASufficientScore || (~isempty(lastAcceptedCombinedScores) ...
                && max(lookahead4bestScoreIncrease.combinedScore) > inInfo.correlationMaximization.earlyBreakLookahead.fcnMinScoreRatio2QualifyForEarlyBreak4CurrentLookahead(...
                    lookaheadCursor, length(lookahead4bestScoreIncrease.ImJ) ...
                  )*max(lastAcceptedCombinedScores)...
                && lookaheadCursor >= inInfo.correlationMaximization.earlyBreakLookahead.minLookahead4MemberCount(min(end,length(lookahead4bestScoreIncrease.ImJ))) ...
              );
              if(bAlreadyFoundASufficientScore)
                if(inInfo.export.nStatusOutputLevel>=3)
                  SDCM_printStatus(3 ...
                    ,'      <- breaking the lookahead loop for member %d early after %d iterations, since we already found a sufficient score increase as per .correlationMaximization.earlyBreakLookahead.fcnMinScoreRatio2QualifyForEarlyBreak4CurrentLookahead\n' ...
                    ,length(generalizedSignature.corrMaximizingMembers.ImJ) ...
                    ,lookaheadCursor ...
                  );
                end
                break; %configured lookahead width for this accu pass already reached and score increase was found.
              end
          end

        %Choose the next member for the generalizedSignature from the computed lookahead phase:
          if(isempty(lookahead4bestScoreIncrease.ImJ))
            sWarning = sprintf('WARNING: No more signatures qualified as member for accumulation => breaking generalization before convergence.\n(This is normal for signatures based on single or very few genes, but these signatures should be validated carefully.)');
              warning(sWarning);
            if(inInfo.plots.focusConvergence.bEnabled)
              sConvergenceFigure.caAccepted{end+1} = ''; %insert a blank line.
              sConvergenceFigure.caAccepted{end+1} = sWarning;
              %Write warning below plot:
                if(ishandle(sConvergenceFigure.hAcceptedText))
                  delete(sConvergenceFigure.hAcceptedText);
                end
                sConvergenceFigure.hAcceptedText = text(...
                   min(xlim(sConvergenceFigure.a4sampleAxis))-abs(diff(xlim(sConvergenceFigure.a4sampleAxis)))*(1/4+1/4+1+1/4) ...
                  ,min(ylim(sConvergenceFigure.a4sampleAxis))-abs(diff(ylim(sConvergenceFigure.a4sampleAxis)))/3 ...
                  ,sConvergenceFigure.caAccepted(max(1,end-13):end) ...
                  ,'HorizontalAlignment','left', 'VerticalAlignment','top' ....
                  ,'Color',[0 0 0]...
                  ,'FontName',iif(ispc(),'Consolas','Courier') ...
                  ,'FontWeight','normal' ...,'FontSize',8.5  ...
                  ,'FontUnits','pixels', 'FontSize',13 ...
                  ,'Parent',sConvergenceFigure.a4sampleAxis...
                );
                drawnow;
            end
            break;
          end
          if(true)
            combinedScores = lookahead4bestScoreIncrease.combinedScore;
            [ScombinedScores,SK] = sort(combinedScores, 'descend');
            %Performance/Batch-accept logic:
              if(true)
                minBestScoreRatio = 0.99;
                if(ScombinedScores(1)>0)
                  BSelectedNewMembers = ScombinedScores/ScombinedScores(1) >= minBestScoreRatio;
                else
                  BSelectedNewMembers = false(size(ScombinedScores));
                  BSelectedNewMembers(1) = true;
                  warning('best score was <=0; just taking the first one');
                end
                ImJ2accumulate = lookahead4bestScoreIncrease.ImJ(SK(BSelectedNewMembers));
                lastAcceptedCombinedScores = combinedScores(SK(BSelectedNewMembers));
              end
          end
          assert(~isempty(ImJ2accumulate), 'code validation: did not select a next candidate from a non-empty lookahead4bestScoreIncrease');
        %Add the selected member to the generalizedSignature and update correlations and scores:
          generalizedSignatureb4Accu = generalizedSignature;
          if(true)
            %Add selected members to generalizedSignature.corrMaximizingMembers:
              nextMemberCursor = length(generalizedSignature.corrMaximizingMembers.ImJ);
              for k=1:length(ImJ2accumulate)
                %Get candidate:
                  nextMember = getSignatureAxesAndScores(ImJ2accumulate(k));
                  nextMemberCursor = nextMemberCursor+1;
                %Store next member candidates:
                  generalizedSignature.corrMaximizingMembers.ImJ(end+1) = nextMember.imj;
                  generalizedSignature.corrMaximizingMembers.sampleAxes_origUnits(nextMemberCursor,:) = nextMember.sampleAxis_origUnits;
                  generalizedSignature.corrMaximizingMembers.geneAxes_origUnits(:,nextMemberCursor) = nextMember.geneAxis_origUnits;
              end
              generalizedSignature = accumulateMembers(generalizedSignature); %reaccumulate with new members.
            %Correlate the new generalized signatures with all genes respectively samples:
              bOnlyOneSidedPerfSubspaceToFindCorrelatedButWeak = true; %<-dev.note: for every accepted member, compute the corrs for all genes/samples (and not just in the corr subspace) to allow member selection from highly correlated but signal-weak genes/samples as per .minAbsR4inclusionInCorrComputation config, too.
              generalizedSignature = correlationConvergenceStep(generalizedSignature, bOnlyOneSidedPerfSubspaceToFindCorrelatedButWeak);
            %Update the performance subspace:
              BsPerformanceSubspace4Correlations = selectDimensionsInPerformanceSubspace(nG,nP ...
               ,generalizedSignature.R4G, generalizedSignature.R4P ...
               ,eDState.noiseEstimation.log10PNoise4G,eDState.noiseEstimation.log10PNoise4P ...
               ,inInfo.searchStrategy ...
              );
          end
          
        %Convergence check: Break accumulation, if the correlations R4G and R4P and axes of the generalizedSignature have converged and if the signature has been sufficiently generalized (i.e. is based on sufficiently many members):
          if(true)
            %Check for axes convergence:
              deltaCorr4R4G =        1-uncenteredWeightedCorrelation(double(generalizedSignatureb4Accu.R4G),                  double(generalizedSignature.R4G),                  1, double(abs(generalizedSignature.signedFocusedW4G)));
              deltaCorr4R4P =        1-uncenteredWeightedCorrelation(double(generalizedSignatureb4Accu.R4P),                  double(generalizedSignature.R4P),                  2, double(abs(generalizedSignature.signedFocusedW4P)));
              deltaCorr4SampleAxis = 1-uncenteredWeightedCorrelation(double(generalizedSignatureb4Accu.geneAxis_origUnits),   double(generalizedSignature.geneAxis_origUnits),   1, double(abs(generalizedSignature.signedFocusedW4G)));
              deltaCorr4GeneAxis =   1-uncenteredWeightedCorrelation(double(generalizedSignatureb4Accu.sampleAxis_origUnits), double(generalizedSignature.sampleAxis_origUnits), 2, double(abs(generalizedSignature.signedFocusedW4P)));
              %In case of multi accept, distribute the deltas over all selected candidates:
                if(length(ImJ2accumulate)>1)
                  deltaCorr4R4G = deltaCorr4R4G/length(ImJ2accumulate);
                  deltaCorr4R4P = deltaCorr4R4P/length(ImJ2accumulate);
                  deltaCorr4SampleAxis = deltaCorr4SampleAxis/length(ImJ2accumulate);
                  deltaCorr4GeneAxis = deltaCorr4GeneAxis/length(ImJ2accumulate);
                end
              bCorrelationsConverged = ...
                   (deltaCorr4R4G+deltaCorr4R4P)/2             < inInfo.correlationMaximization.convergenceEpsilon4deltaCorr4Rs(min(end,sL)) ...
                && (deltaCorr4SampleAxis+deltaCorr4GeneAxis)/2 < inInfo.correlationMaximization.convergenceEpsilon4deltaCorr4AxesInNoiseSDs(min(end,sL))*eDState.noiseEstimation.sd4epsilon ...  %<-TODO/dev.note: consider removing *eDState.noiseEstimation.sd4epsilon from the break condition and allow spevifying an absolute epsilon to make this decision independent from the noise estimation.
                ...&& (deltaCorr4SampleAxis+deltaCorr4GeneAxis)/2 < inInfo.correlationMaximization.convergenceEpsilon4deltaCorr4Axes(min(end,sL)) ...
              ;
            %Check for sufficient representatives for generalization:
              nDistinctGenesSoFar = sum(inInfo.reference.rowW(+generalizedSignature.corrMaximizingMembers.ImJ(generalizedSignature.corrMaximizingMembers.ImJ>0))); %respect external weights (e.g. in case of multiple probesets for the same gene)
              nDistinctSamplesSoFar = sum(inInfo.reference.colW(-generalizedSignature.corrMaximizingMembers.ImJ(generalizedSignature.corrMaximizingMembers.ImJ<0))); %respect external weights (e.g. in case of multiple probesets for the same gene)
              %Get estimated signature size:
                nSignatureSize4G = generalizedSignature.signatureSizeByCorrSum4G;
                nSignatureSize4P = generalizedSignature.signatureSizeByCorrSum4P;
              minDistinctMembers = min([
                inInfo.correlationMaximization.nSufficient4Generalization(min(end,sL));
                inInfo.correlationMaximization.minRatioOfSignatureSize4Generalization(min(end,sL))*(nSignatureSize4G+nSignatureSize4P);
              ]);
              bGeneralized = ...
                  nDistinctGenesSoFar >= nSignatureSize4G*inInfo.correlationMaximization.minRatioOfSignatureSizeInSingleDim4Generalization(min(end,sL)) ... 
               || nDistinctSamplesSoFar >= nSignatureSize4P*inInfo.correlationMaximization.minRatioOfSignatureSizeInSingleDim4Generalization(min(end,sL)) ... 
               || nDistinctGenesSoFar+nDistinctSamplesSoFar >= minDistinctMembers ...
              ;
            %Emergency break in case of no convergence even after incorporating many members:
              if(~bCorrelationsConverged && nDistinctGenesSoFar+nDistinctSamplesSoFar>=10*inInfo.correlationMaximization.nSufficient4Generalization(min(end,sL)))
                warning('emergency break of signature axes convergence: no convergence relative to configured epsilons reached with 10*.correlationMaximization.nSufficient4Generalization=%0.1f members; no accepting AS IS; inspect or use >>dbcont to continue', 10*inInfo.correlationMaximization.nSufficient4Generalization(min(end,sL)));
                bCorrelationsConverged = true;
                bGeneralized = true;
%keyboard; %never happened so far => check interactively, once/if this ever happens.
              end

            %Status output: display accepted signature members:
              if(inInfo.export.nStatusOutputLevel>=2 || inInfo.plots.focusConvergence.bEnabled)
                for k=1:length(ImJ2accumulate)
                  sAccepted = SDCM_printStatus(2 ...
                   ,['   ->accepted as representative %03d: #%04d=[% 7s#%05d % ',num2str(maxLabelLength),'s ] with (r2g=%+0.2f,r2p=%+0.2f) => corr=(% 7.1fg,% 7.1fs, %0.2fr)->(% 7.1fg,% 7.1fs, %0.2fr)\n'...not computed for performance reasons/; signature/noise=%0.1f->%0.1f, log10(p)=%5.1f->%5.1f\n'...
                    ] ...
                   ...New member indices and name:
                       ,length(generalizedSignature.corrMaximizingMembers.ImJ)-length(ImJ2accumulate)+k ...
                       ...,iif(ImJ2accumulate(1)>0,'g','p'), ilv(processingOrder.SImJ(1:find(processingOrder.SImJ==ImJ2accumulate(1),1)),@(ImJprocessedBefore)iif(ImJ2accumulate(1)>0,sum(ImJprocessedBefore>0),sum(ImJprocessedBefore<0))) ...
                       ,find(processingOrder.SImJ==ImJ2accumulate(k),1) ...
                       ,iif(ImJ2accumulate(k)>0, 'gene', 'sample') ...
                       ,iif(ImJ2accumulate(k)>0, ImJ2accumulate(k), -ImJ2accumulate(k)) ...
                       ,iif(ImJ2accumulate(k)>0, @()rowLabels{ImJ2accumulate(k)}, @()colLabels{-ImJ2accumulate(k)}) ...
                   ...Scores of the new member:
                       ,lookahead4bestScoreIncrease.r2currentGeneralizedGeneAxis(lookahead4bestScoreIncrease.ImJ==ImJ2accumulate(k))...
                       ,lookahead4bestScoreIncrease.r2currentGeneralizedSampleAxis(lookahead4bestScoreIncrease.ImJ==ImJ2accumulate(k))...
                   ...Score changed for the new focus:
                      ...Signature size estimates and average absolute correlation in the focus:
                           ,generalizedSignatureb4Accu.signatureSizeByCorrSum4G, generalizedSignatureb4Accu.signatureSizeByCorrSum4P ...
                           ,generalizedSignatureb4Accu.signatureCorrInExtendedFocus...
                           ,generalizedSignature.signatureSizeByCorrSum4G, generalizedSignature.signatureSizeByCorrSum4P ...
                           ,generalizedSignature.signatureCorrInExtendedFocus...
                      ...Signature significance:
                           ...not computed for performance reasons/,generalizedSignatureb4Accu.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k), generalizedSignature.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ...
                           ...not computed for performance reasons/,generalizedSignatureb4Accu.log10_p, generalizedSignature.log10_p ...
                  );
                  if(inInfo.plots.focusConvergence.bEnabled)
                    sConvergenceFigure.caAccepted{end+1} = strtrim(strrep(sAccepted,'<-',''));
                  end
                end
                if(inInfo.internal.bDevEnableInteractiveBreaks)
                  interactiveCheckPoint();
                end
              end
            %Update plot for correlation maximization/convergence:
              if(inInfo.plots.focusConvergence.bEnabled)
                %Write convergence deltas in the plots:
                  plot(sConvergenceFigure.a4CorrelationConvergence, length(generalizedSignature.corrMaximizingMembers.ImJ), deltaCorr4R4P+eps, 'r.');
                  plot(sConvergenceFigure.a4CorrelationConvergence, length(generalizedSignature.corrMaximizingMembers.ImJ), deltaCorr4R4G+eps, '.','Color',[0 0 0.85]);
                  plot(sConvergenceFigure.a4CorrelationConvergence, length(generalizedSignature.corrMaximizingMembers.ImJ), (deltaCorr4R4G+deltaCorr4R4P+eps)/2, 'k.');
                  plot(sConvergenceFigure.a4SignatureConvergence, length(generalizedSignature.corrMaximizingMembers.ImJ), deltaCorr4GeneAxis+eps, 'r.');
                  plot(sConvergenceFigure.a4SignatureConvergence, length(generalizedSignature.corrMaximizingMembers.ImJ), deltaCorr4SampleAxis+eps, '.','Color',[0 0 0.85]);
                  plot(sConvergenceFigure.a4SignatureConvergence, length(generalizedSignature.corrMaximizingMembers.ImJ), (deltaCorr4GeneAxis+deltaCorr4SampleAxis+eps)/2, 'k.');
                  %Adapt axes:
                    if(inInfo.plots.bVisible && length(generalizedSignature.corrMaximizingMembers.ImJ)-nMembersInLastUpdate>=nMembers4PlotUpdateInterval) %performance (only every tenth addition and if visible)
                      ylim(sConvergenceFigure.a4geneAxis, 1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.geneAxis_origUnits)));
                      ylim(sConvergenceFigure.a4sampleAxis,   1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.sampleAxis_origUnits)));
                      xlim(sConvergenceFigure.a4geneAxis, 1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.geneAxis_origUnits)));
                      xlim(sConvergenceFigure.a4sampleAxis,   1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.sampleAxis_origUnits)));
                      xlim(sConvergenceFigure.a4CorrelationConvergence, [1,length(generalizedSignature.corrMaximizingMembers.ImJ)+nMembers4PlotUpdateInterval-1]+1/2);
                      xlim(sConvergenceFigure.a4SignatureConvergence, [1,length(generalizedSignature.corrMaximizingMembers.ImJ)+nMembers4PlotUpdateInterval-1]+1/2);
                    end
                %Set plots from previous focusing iterations to gray:
                  set(findobj(sConvergenceFigure.a4R4P,'Type','line'),'Color',[.5 .5 .5]);%,'LineStyle','-');
                  set(findobj(sConvergenceFigure.a4R4G,'Type','line'),'Color',[.5 .5 .5]);%,'LineStyle','-');
                  set(findobj(sConvergenceFigure.a4sampleAxis,'Type','line'),'Color',[.5 .5 .5]); %,'LineStyle','-');
                  set(findobj(sConvergenceFigure.a4geneAxis,'Type','line'),'Color',[.5 .5 .5]); %,'LineStyle','-');
                  plot(sConvergenceFigure.a4R4G, generalizedSignature.R4G, generalizedSignatureb4Accu.R4G, '.','Color',[0 0 0.85]); 
                  plot(sConvergenceFigure.a4R4P, generalizedSignature.R4P, generalizedSignatureb4Accu.R4P, 'r.'); 
                    %<-Note: dass previous==0 bei erster Aktualisierung liegt daran, dass das die in Step1 vie BsPerformanceSubspace ausgeschlossenen Punkte sind, die erst in Step2 eine Korrelation zu den aktuellen Effektachsen erhalten.
                  plot(sConvergenceFigure.a4geneAxis, generalizedSignature.geneAxis_origUnits, generalizedSignatureb4Accu.geneAxis_origUnits, '.','Color',[0 0 0.85]);
                  plot(sConvergenceFigure.a4sampleAxis, generalizedSignature.sampleAxis_origUnits, generalizedSignatureb4Accu.sampleAxis_origUnits, 'r.');
                %Write accepted genes/samples below the plot:
                  if(inInfo.plots.bVisible && length(generalizedSignature.corrMaximizingMembers.ImJ)-nMembersInLastUpdate>=nMembers4PlotUpdateInterval) %performance (only every tenth addition and if visible)
                    if(ishandle(sConvergenceFigure.hAcceptedText))
                      delete(sConvergenceFigure.hAcceptedText);
                    end
                    sConvergenceFigure.hAcceptedText = text(...
                       min(xlim(sConvergenceFigure.a4sampleAxis))-abs(diff(xlim(sConvergenceFigure.a4sampleAxis)))*(1/4+1/4+1+1/4) ...
                      ,min(ylim(sConvergenceFigure.a4sampleAxis))-abs(diff(ylim(sConvergenceFigure.a4sampleAxis)))/3 ...
                      ,sConvergenceFigure.caAccepted(max(1,end-13):end) ...
                      ,'HorizontalAlignment','left', 'VerticalAlignment','top' ....
                      ,'Color',[0 0 0]...
                      ,'FontName',iif(ispc(),'Consolas','Courier') ...
                      ,'FontWeight','normal' ...,'FontSize',8.5  ...
                      ,'FontUnits','pixels', 'FontSize',13 ...
                      ,'Parent',sConvergenceFigure.a4sampleAxis...
                    );
                  end
              end
            %keep GUI responsive:
              if(inInfo.plots.bVisible && length(generalizedSignature.corrMaximizingMembers.ImJ)-nMembersInLastUpdate>=nMembers4PlotUpdateInterval) %performance (only every tenth addition and if visible)
                nMembersInLastUpdate = length(generalizedSignature.corrMaximizingMembers.ImJ);
                drawnow;
              end
            %Status output: display convergence status:
              if(bCorrelationsConverged && bGeneralized)
                sConverged = SDCM_printStatus(1 ...
                  ,'  <- MAXIMIZED and GENERALIZED: breaking the accumulation phase for correlation maximization and signature generalization after %d representatives, as they are sufficient for generalization as per .correlationMaximization.nSufficient4Generalization and .minRatioOfSignatureSize4Generalization configuration.\n' ...
                  ,length(generalizedSignature.corrMaximizingMembers.ImJ) ...
                );
                if(inInfo.plots.focusConvergence.bEnabled)
                  sConvergenceFigure.caAccepted{end+1} = strtrim(sConverged);
                  %Write below plot:
                    if(ishandle(sConvergenceFigure.hAcceptedText))
                      delete(sConvergenceFigure.hAcceptedText);
                    end
                    sConvergenceFigure.hAcceptedText = text(...
                       min(xlim(sConvergenceFigure.a4sampleAxis))-abs(diff(xlim(sConvergenceFigure.a4sampleAxis)))*(1/4+1/4+1+1/4) ...
                      ,min(ylim(sConvergenceFigure.a4sampleAxis))-abs(diff(ylim(sConvergenceFigure.a4sampleAxis)))/3 ...
                      ,sConvergenceFigure.caAccepted(max(1,end-13):end) ...
                      ,'HorizontalAlignment','left', 'VerticalAlignment','top' ....
                      ,'Color',[0 0 0]...
                      ,'FontName',iif(ispc(),'Consolas','Courier') ...
                      ,'FontWeight','normal' ...,'FontSize',8.5  ...
                      ,'FontUnits','pixels', 'FontSize',13 ...
                      ,'Parent',sConvergenceFigure.a4sampleAxis...
                    );
                    drawnow;
                end
                break; %accept accumulation result as generalized signature 
              else
                SDCM_printStatus(2 ...
                 ,[   '     <- convergence: deltas=(%0.1e dR4G, %0.1e dSampleAxis, %0.1e dR4P, %0.1e dGeneAxis) = (%0.1fcorr, %0.1faxes)*epsilon'...
                   ,'\n                     independent members mass for generalization = (%0.1fg+%0.1fs)/%0.1f=%0.1f%% => not converged or not generalized, yet.\n' ...
                  ]...
                 ...Convergence:
                     ,deltaCorr4R4G, deltaCorr4SampleAxis,  deltaCorr4R4P, deltaCorr4GeneAxis...
                     ,(deltaCorr4R4G+deltaCorr4R4P)/2 / inInfo.correlationMaximization.convergenceEpsilon4deltaCorr4Rs(min(end,sL)) ...
                     ,(deltaCorr4SampleAxis+deltaCorr4GeneAxis)/2 / (inInfo.correlationMaximization.convergenceEpsilon4deltaCorr4AxesInNoiseSDs(min(end,sL))*eDState.noiseEstimation.sd4epsilon) ...
                     ...
                 ...Generalization:
                    ,nDistinctGenesSoFar, nDistinctSamplesSoFar, minDistinctMembers, 100*(nDistinctGenesSoFar+nDistinctSamplesSoFar)/minDistinctMembers ...
                );
              end
          end
      end
      %Plot correlation maximization/convergence, if enabled:
        if(inInfo.plots.focusConvergence.bEnabled)
          %Finalize axes view:
            %Adapt axes for correlation convergence overview:
              axis(sConvergenceFigure.a4CorrelationConvergence, 'tight');
              axis(sConvergenceFigure.a4SignatureConvergence, 'tight');
              ylim(sConvergenceFigure.a4CorrelationConvergence, ylim(sConvergenceFigure.a4CorrelationConvergence).*[1/1.21,1.21]);
              ylim(sConvergenceFigure.a4SignatureConvergence, ylim(sConvergenceFigure.a4SignatureConvergence).*[1/1.21,1.21]);
              xlim(sConvergenceFigure.a4CorrelationConvergence, [1,length(generalizedSignature.corrMaximizingMembers.ImJ)+eps]+1/2); %+eps in case only one signature qualifies for accumulation.
              xlim(sConvergenceFigure.a4SignatureConvergence, [1,length(generalizedSignature.corrMaximizingMembers.ImJ)+eps]+1/2);
              set(sConvergenceFigure.a4CorrelationConvergence,'YTick',10.^(max(-100,floor(log10(min(ylim(sConvergenceFigure.a4CorrelationConvergence))))):1:0));
              set(sConvergenceFigure.a4SignatureConvergence,'YTick',10.^(max(-100,floor(log10(min(ylim(sConvergenceFigure.a4SignatureConvergence))))):1:0));
            %tighten axes for signatures:
              ylim(sConvergenceFigure.a4R4G,1.03*[-1,1]);
              ylim(sConvergenceFigure.a4R4P,1.03*[-1,1]);
              ylim(sConvergenceFigure.a4geneAxis, 1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.geneAxis_origUnits)));
              ylim(sConvergenceFigure.a4sampleAxis,   1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.sampleAxis_origUnits)));
              xlim(sConvergenceFigure.a4R4G,1.03*[-1,1]);
              xlim(sConvergenceFigure.a4R4P,1.03*[-1,1]);
              xlim(sConvergenceFigure.a4geneAxis, 1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.geneAxis_origUnits)));
              xlim(sConvergenceFigure.a4sampleAxis,   1.03*[-1,1]*nanmax(eps+abs(generalizedSignature.sampleAxis_origUnits)));
            %Update text (position has changed due to changed xlim):
              %Write initial representative as title (every iteration, since the ylim might change):
                if(ishandle(sConvergenceFigure.hInitialSignatureText)) 
                  delete(sConvergenceFigure.hInitialSignatureText);
                end
                sConvergenceFigure.hInitialSignatureText = text(...
                   min(xlim(sConvergenceFigure.a4R4P))-abs(diff(xlim(sConvergenceFigure.a4R4P)))*(1/4+1/4+1+1/4) ...
                  ,max(ylim(sConvergenceFigure.a4R4P))+abs(diff(ylim(sConvergenceFigure.a4R4P)))/2.5 ...
                  ,initialSignature.sInitialRepresenative ...
                  ,'HorizontalAlignment','left', 'VerticalAlignment','top' ....
                  ,'Color',[0 0 1]...
                  ,'FontName',iif(ispc(),'Consolas','Courier') ...
                  ,'FontWeight','normal'...,'FontSize',8.5  ...
                  ,'FontUnits','pixels','FontSize',13 ...
                  ,'Parent',sConvergenceFigure.a4R4P...
                );
              %Write accepted genes/samples below the plot:
                if(ishandle(sConvergenceFigure.hAcceptedText))
                  delete(sConvergenceFigure.hAcceptedText);
                end
                sConvergenceFigure.hAcceptedText = text(...
                   min(xlim(sConvergenceFigure.a4sampleAxis))-abs(diff(xlim(sConvergenceFigure.a4sampleAxis)))*(1/4+1/4+1+1/4) ...
                  ,min(ylim(sConvergenceFigure.a4sampleAxis))-abs(diff(ylim(sConvergenceFigure.a4sampleAxis)))/3 ...
                  ,sConvergenceFigure.caAccepted(max(1,end-13):end) ...
                  ,'HorizontalAlignment','left', 'VerticalAlignment','top' ....
                  ,'Color',[0 0 0]...
                  ,'FontName',iif(ispc(),'Consolas','Courier') ...
                  ,'FontWeight','normal' ...,'FontSize',8.5  ...
                  ,'FontUnits','pixels', 'FontSize',13 ...
                  ,'Parent',sConvergenceFigure.a4sampleAxis...
                );
          %Export, if requested:
            if(inInfo.export.plots.bEnabled)
              try
                %Hotfix for truncated fonts in PNG export:
                  set(sConvergenceFigure.hInitialSignatureText,'FontUnits','pixels','FontSize',10);
                  set(sConvergenceFigure.hAcceptedText,'FontUnits','pixels','FontSize',10);
                    figPos = get(sConvergenceFigure.f,'Position');
                    exportPos = get(0,'ScreenSize');
                      exportPos(4) = figPos(4)*1.09;
                    saveFig(sConvergenceFigure.f,inInfo.export.rootDir,sprintf(...
                       '%03d, step 2, signature focus convergence'...
                       ,eDState.current.k ...
                    ),{'png'},[],[],[],exportPos); %use png because the overlay plotting ansatz of all iterations produces huge eps files...
                  set(sConvergenceFigure.hInitialSignatureText,'FontUnits','pixels','FontSize',13);
                  set(sConvergenceFigure.hAcceptedText,'FontUnits','pixels','FontSize',13);
                if(inInfo.export.plots.bCloseAfterExport) close(sConvergenceFigure.f); else set(sConvergenceFigure.f,'Visible','on'); end          
              catch ex
                warning('Error during export of the focus convergence figure; details: %s\nInteractive break point for manual retry/inspection; use >>dbcont to ignore this and continue.', ex.message);
                keyboard;
              end
            end
        end
    %% Finalize the signature for regression:
      %Performance: combined p values for Rs were omitted in each correlationConvergenceStep for performance reasons; now compute them just once for the converged axes/corrs:
        [~, log10_p4Correlations4G] = pValues4Correlations(...
          generalizedSignature.R4G ... 
         ,generalizedSignature.sampleSizes4R4G ... 
        ,1);
        [~, log10_p4Correlations4P] = pValues4Correlations(...
          generalizedSignature.R4P ... 
         ,generalizedSignature.sampleSizes4R4P ... 
        ,2);
        generalizedSignature.log10_p4Correlations = min(log10_p4Correlations4G, log10_p4Correlations4P); %these are two indirectly dependent p values (same L2Rs base); the min should be an upper bound of the combined p value; too conservative?!
      %Compute signature strengths (i.e. projections on axes) for all points (these signature strenghts determine the signautre eigenOrder for regression in the default):
        generalizedSignature.signatureStrengths4G = projectW(...
          eDState.current.L2Rs ... vectorsInSourceSpace...
         ,generalizedSignature.sampleAxis_origUnits ... axesInSourceSpace
         ,abs(generalizedSignature.signedFocusedW4P) ... weights4sourceDims
         ,2 ... sourceSpaceMatrixDim
        ,BM.meanW, BM.euclidW) / BM.euclidW(ones(1,nP),abs(generalizedSignature.signedFocusedW4P),2); %normalize such that every compoenent is a representative expression independent of the size of the aggregated space.
        generalizedSignature.signatureStrengths4P = projectW(...
          eDState.current.L2Rs ... vectorsInSourceSpace...
         ,generalizedSignature.geneAxis_origUnits ... axesInSourceSpace
         ,abs(generalizedSignature.signedFocusedW4G) ... weights4sourceDims
         ,1 ... sourceSpaceMatrixDim
        ,BM.meanW, BM.euclidW) / BM.euclidW(ones(nG,1),abs(generalizedSignature.signedFocusedW4G),1); %normalize such that every compoenent is a representative expression independent of the size of the aggregated space.
      %Status output: print the signature-generalizing genes and sample names and the top of the signature as per inInfo.signatureDefinition.generalization.minRelAbsCorr*:
        if(true) %final scores for the signature.
          SDCM_printStatus(2, [...
            '  <- After correlation maximization, signature statistics increased by corr=(% 7.1fg,% 7.1fs, %0.2fr)->(% 7.1fg,% 7.1fs, %0.2fr), signal=(%0.2ffocus=%1.2fSNR)->(%0.2ffocus=%1.2fSNR), log10(p)=(%7.1fcorr,%6.1fsignal)->(%6.1fcorr,%6.1fsignal)\n' ...
            ]...
           ...Signature size estimates and average absolute correlation in the focus:
               ,initialSignature.signatureSizeByCorrSum4G, initialSignature.signatureSizeByCorrSum4P ...
               ,initialSignature.signatureCorrInExtendedFocus...
               ,generalizedSignature.signatureSizeByCorrSum4G, generalizedSignature.signatureSizeByCorrSum4P ...
               ,generalizedSignature.signatureCorrInExtendedFocus...
           ...Signal strengths relative to initial noise level:
               ,initialSignature.signatureAbsMean2D, initialSignature.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ... 
               ,generalizedSignature.signatureAbsMean2D, generalizedSignature.signatureAbsMean2D/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) ... 
           ...Signature significance:
               ,initialSignature.log10_p4Correlations, initialSignature.log10_p4SignalStrength ...
               ,generalizedSignature.log10_p4Correlations, generalizedSignature.log10_p4SignalStrength ...
          );
          if(  10^generalizedSignature.log10_p4Correlations   > inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))... 
            || 10^generalizedSignature.log10_p4SignalStrength > inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL))...
          )
            warning('<- After focusing the p value increased above the configured significance threshold; maybe the signal is to weak for stable focusing. Depending on the context, this signature should be validated in a separate dataset.');
          end
        end
        if(true) %signature generalizing/spanning genes/samples:
          SDCM_printStatus(2, ['   - Generalization finished using the following genes and samples as signature representatives:\n']);
          I = generalizedSignature.corrMaximizingMembers.ImJ(generalizedSignature.corrMaximizingMembers.ImJ>0);
            I = I(ilnth(2,@sort,abs(generalizedSignature.R4G(I)),'descend'));
            B = ismember((1:nG)',I) & generalizedSignature.R4G>0;
              SDCM_printStatus(2 ...
                ,'        <- % 3d correlated genes used for generalization:       %s\n' ...
                ,sum(B) ...
                ,cellstring2separatedList(rowLabels(B),' | ') ...
              );
            B = ismember((1:nG)',I) & generalizedSignature.R4G<0;
              SDCM_printStatus(2 ...
                ,'        <- % 3d anticorrelated genes used for generalization:   %s\n' ...
                ,sum(B) ...
                ,cellstring2separatedList(rowLabels(B),' | ') ...
              );

          J = -generalizedSignature.corrMaximizingMembers.ImJ(generalizedSignature.corrMaximizingMembers.ImJ<0);
            J = J(ilnth(2,@sort,abs(generalizedSignature.R4P(J)),'descend'));
            B = ismember(1:nP,J) & generalizedSignature.R4P>0;
              SDCM_printStatus(2 ...
                ,'        <- % 3d correlated samples used for generalization:     %s\n' ...
                ,sum(B) ...
                ,cellstring2separatedList(colLabels(B),' | ') ...
              );
            B = ismember(1:nP,J) & generalizedSignature.R4P<0;
              SDCM_printStatus(2 ...
                ,'        <- % 3d anticorrelated samples used for generalization: %s\n' ...
                ,sum(B) ...
                ,cellstring2separatedList(colLabels(B),' | ') ...
              );
        end
        if(true) %Information about the top correlated genes/samples (cutoff via inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio*)
          SDCM_printStatus(1, ['   ->The following genes and samples are top correlated to the signature (wrt. inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds):\n']);

          B = generalizedSignature.R4G/max(abs(generalizedSignature.R4G)) >= +inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G;
          I = ilsub(find(B),ilnth(2,@sort,abs(generalizedSignature.R4G(B)),'descend')); %present in descending order of their absolute correlations.
            labelsWithCorrs = cellfun(@(s,r)sprintf('%s (r=%0.2f)',s,r), rowLabels(I), num2cell(generalizedSignature.R4G(I)), 'UniformOutput', false);
            sStatusMessage = SDCM_printStatus(1 ...
              ,'        <- % 5d genes   with relative correlation >=%+0.2f: %s\n'...
              ,length(I)...
              ,+inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G ...
              ,cellstring2separatedList(labelsWithCorrs,' | ') ...
            );
              if(inInfo.export.infoFile4topMembers.bEnabled4topCorrGenes && any(B)) %For interactive results control and searching by file, create a file name with the core statistics of the qualified initial representative:
                nMaxLengh = 100;
                fileName = sprintf('%03d, top corr. genes = %s.info', eDState.current.k, string2fieldname(ilv(cellstring2separatedList(labelsWithCorrs,'; '),@(s)s(1:min(end,nMaxLengh))),true,true,true));
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
                fprintf(fid,'%s',sStatusMessage);
                fclose(fid);
              end
          
          B = generalizedSignature.R4G/max(abs(generalizedSignature.R4G)) <= -inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G;
          I = ilsub(find(B),ilnth(2,@sort,abs(generalizedSignature.R4G(B)),'descend')); %present in descending order of their absolute correlations.
            labelsWithCorrs = cellfun(@(s,r)sprintf('%s (r=%0.2f)',s,r), rowLabels(I), num2cell(generalizedSignature.R4G(I)), 'UniformOutput', false);
            sStatusMessage = SDCM_printStatus(1 ...
              ,'        <- % 5d genes   with relative correlation <=%+0.2f: %s\n'...
              ,length(I)...
              ,-inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G ...
              ,cellstring2separatedList(labelsWithCorrs,' | ') ...
            );
            if(inInfo.export.infoFile4topMembers.bEnabled4topAntiCorrGenes && any(B)) %For interactive results control and searching by file, create a file name with the core statistics of the qualified initial representative:
              nMaxLengh = 100;
              fileName = sprintf('%03d, top anti-corr. genes = %s.info', eDState.current.k, string2fieldname(ilv(cellstring2separatedList(labelsWithCorrs,'; '),@(s)s(1:min(end,nMaxLengh))),true,true,true));
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
              fprintf(fid,'%s',sStatusMessage);
              fclose(fid);
            end

          B = generalizedSignature.R4P/max(abs(generalizedSignature.R4P)) >= +inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P;
          J = ilsub(find(B),ilnth(2,@sort,abs(generalizedSignature.R4P(B)),'descend')); %present in descending order of their absolute correlations.
            labelsWithCorrs = cellfun(@(s,r)sprintf('%s (r=%0.2f)',s,r), colLabels(J), num2cell(generalizedSignature.R4P(J)), 'UniformOutput', false);
            sStatusMessage = SDCM_printStatus(1 ...
              ,'        <- % 5d samples with relative correlation >=%+0.2f: %s\n'...
              ,length(J)...
              ,+inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P ...
              ,cellstring2separatedList(labelsWithCorrs,' | ') ...
            );
            if(inInfo.export.infoFile4topMembers.bEnabled4topCorrSamples && any(B)) %For interactive results control and searching by file, create a file name with the core statistics of the qualified initial representative:
              nMaxLengh = 100;
              fileName = sprintf('%03d, top corr. samples = %s.info', eDState.current.k, string2fieldname(ilv(cellstring2separatedList(labelsWithCorrs,'; '),@(s)s(1:min(end,nMaxLengh))),true,true,true));
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
              fprintf(fid,'%s',sStatusMessage);
              fclose(fid);
            end
          B = generalizedSignature.R4P/max(abs(generalizedSignature.R4P)) <= -inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P;
          J = ilsub(find(B),ilnth(2,@sort,abs(generalizedSignature.R4P(B)),'descend')); %present in descending order of their absolute correlations.
            labelsWithCorrs = cellfun(@(s,r)sprintf('%s (r=%0.2f)',s,r), colLabels(J), num2cell(generalizedSignature.R4P(J)), 'UniformOutput', false);
            sStatusMessage = SDCM_printStatus(1 ...
              ,'        <- % 5d samples with relative correlation <=%+0.2f: %s\n'...
              ,length(J)...
              ,-inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P ...
              ,cellstring2separatedList(labelsWithCorrs,' | ') ...
            );
            if(inInfo.export.infoFile4topMembers.bEnabled4topAntiCorrSamples && any(B)) %For interactive results control and searching by file, create a file name with the core statistics of the qualified initial representative:
              nMaxLengh = 100;
              fileName = sprintf('%03d, top anti-corr. samples = %s.info', eDState.current.k, string2fieldname(ilv(cellstring2separatedList(labelsWithCorrs,'; '),@(s)s(1:min(end,nMaxLengh))),true,true,true));
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
              fprintf(fid,'%s',sStatusMessage);
              fclose(fid);
            end
        end
      %2020-02/TODO: if signature statistics decreased (no stable local maximum, wandering axes in high-dim), add an option to just accept the (in this case better) estimation based on the initial representative from the search strategy.
        %TODO.
      if(inInfo.internal.bDevEnableInteractiveBreaks)
        interactiveCheckPoint();
      end
      return; %return accumulation result as generalized signature 

    %% Accumulation of members and correlation update/convergence step:
      function candidateAccus = addMemberAndAccumulate(generalizedSignature, ImJ, nextGeneAxes_origUnits, nextSampleAxes_origUnits)
        %Initialize:
          nCurrentMembers = length(generalizedSignature.corrMaximizingMembers.ImJ);
          nCandidates4Accumulation = length(ImJ);
          assert(size(nextGeneAxes_origUnits,1)==nCandidates4Accumulation,'code validation: length(ImJ)~=size(nextGeneAxes,1)');
          assert(size(nextSampleAxes_origUnits,2)==nCandidates4Accumulation,'code validation: length(ImJ)~=size(nextSampleAxes,2)');
        %Accumulation:
          %Get all current and candidate members and compute their standard deviations:
            allGeneAxes_origUnits = [generalizedSignature.corrMaximizingMembers.sampleAxes_origUnits; nextGeneAxes_origUnits];
            allSampleAxes_origUnits = [generalizedSignature.corrMaximizingMembers.geneAxes_origUnits, nextSampleAxes_origUnits];
          %Compute member correlations with the signature of the current generalized focus (for accumulation weights below):
            if(true)
              R4sampleAxisMembers = uncenteredWeightedCorrelation(... %!must use the correct signs for constructive accumulation!
                 generalizedSignature.sampleAxis_origUnits...
                ,allGeneAxes_origUnits ...
                ,2 ...
                ,abs(generalizedSignature.signedFocusedW4P_withPerpendicularSpace) ... 
                ,BM.meanW ...
              ); %<-Note.depthTest: must only use sign instead of Rs for depth dissection!
              R4geneAxisMembers = uncenteredWeightedCorrelation(... %!must use the correct signs for constructive accumulation!
                 generalizedSignature.geneAxis_origUnits...
                ,allSampleAxes_origUnits ...
                ,1 ...
                ,abs(generalizedSignature.signedFocusedW4G_withPerpendicularSpace) ...
                ,BM.meanW ...
              ); %<-Note.depthTest: must only use sign instead of Rs for depth dissection!
            end
          %Build/accumulate candidates:
            candidateAccus = repmat(generalizedSignature,nCandidates4Accumulation,1);
            for k=1:nCandidates4Accumulation
              %Add member:
                candidateAccus(k).corrMaximizingMembers.ImJ = [generalizedSignature.corrMaximizingMembers.ImJ, ImJ(k)];
                candidateAccus(k).corrMaximizingMembers.sampleAxes_origUnits = [generalizedSignature.corrMaximizingMembers.sampleAxes_origUnits; nextGeneAxes_origUnits(k,:)];
                candidateAccus(k).corrMaximizingMembers.geneAxes_origUnits = [generalizedSignature.corrMaximizingMembers.geneAxes_origUnits, nextSampleAxes_origUnits(:,k)];
              %Accumulate sample axis:
                candidateAccus(k).corrMaximizingMembers.R4sampleAxisMembers = R4sampleAxisMembers([1:nCurrentMembers,nCurrentMembers+k]);
                candidateAccus(k).sampleAxis_origUnits = BM.meanW(...
                   candidateAccus(k).corrMaximizingMembers.sampleAxes_origUnits ...
                  ,candidateAccus(k).corrMaximizingMembers.R4sampleAxisMembers ... %signs for constructive average.
                ,1);
                %It is important to keep r4g and geneAxis positively correlated for later dissection monotonicity direction:
                  if(sign(uncenteredWeightedCorrelation(candidateAccus(k).R4P,candidateAccus(k).sampleAxis_origUnits,2,1,BM.meanW))<0)
                    warning('code validation: candidateAccus(k).sampleAxis is not positively correlated with candidateAccus(k).R4P; flipping .sampleAxis');
                    candidateAccus(k).sampleAxis_origUnits = -candidateAccus(k).sampleAxis_origUnits;
                  end
              %Accumulate gene axis:
                candidateAccus(k).corrMaximizingMembers.R4geneAxisMembers = R4geneAxisMembers([1:nCurrentMembers,nCurrentMembers+k]);
                candidateAccus(k).geneAxis_origUnits = BM.meanW(... 
                   candidateAccus(k).corrMaximizingMembers.geneAxes_origUnits ...
                  ,candidateAccus(k).corrMaximizingMembers.R4geneAxisMembers ... %signs for constructive average
                ,2);
                %It is important to keep r4g and geneAxis positively correlated for later dissection monotonicity direction:
                  if(sign(uncenteredWeightedCorrelation(candidateAccus(k).R4G,candidateAccus(k).geneAxis_origUnits,1,1,BM.meanW))<0)
                    warning('code validation: candidateAccus(k).geneAxis is not positively correlated with candidateAccus(k).R4G; flipping .geneAxis');
                    candidateAccus(k).geneAxis_origUnits = -candidateAccus(k).geneAxis_origUnits;
                  end
              %Compute signature vector norms: 
                candidateAccus(k).norm4geneAxis_origUnits = BM.euclidW(candidateAccus(k).geneAxis_origUnits, 1, 1);          
                candidateAccus(k).norm4sampleAxis_origUnits = BM.euclidW(candidateAccus(k).sampleAxis_origUnits, 1, 2);
            end
      end
        function generalizedSignature = accumulateMembers(generalizedSignature)
          %Reuse addMemberAndAccumulate:
            lastmember.imj = generalizedSignature.corrMaximizingMembers.ImJ(end);
            lastmember.sampleAxis_origUnits = generalizedSignature.corrMaximizingMembers.sampleAxes_origUnits(end,:);
            lastmember.geneAxis_origUnits = generalizedSignature.corrMaximizingMembers.geneAxes_origUnits(:,end);
            generalizedSignature_allButLast = generalizedSignature;
              generalizedSignature_allButLast.corrMaximizingMembers.ImJ = generalizedSignature_allButLast.corrMaximizingMembers.ImJ(1:end-1);
              generalizedSignature_allButLast.corrMaximizingMembers.sampleAxes_origUnits = generalizedSignature_allButLast.corrMaximizingMembers.sampleAxes_origUnits(1:end-1,:);
              generalizedSignature_allButLast.corrMaximizingMembers.geneAxes_origUnits = generalizedSignature_allButLast.corrMaximizingMembers.geneAxes_origUnits(:,1:end-1);
            generalizedSignature = addMemberAndAccumulate(generalizedSignature_allButLast, lastmember.imj, lastmember.sampleAxis_origUnits, lastmember.geneAxis_origUnits);
        end
      function signatures = correlationConvergenceStep(signatures, bOnlyOneSidedPerfSubspaceToFindCorrelatedButWeak)
        %Initialize:
          fn4L2Rs_origUnits = 'L2Rs';
          if(nargin<2) bOnlyOneSidedPerfSubspaceToFindCorrelatedButWeak = false; end
        %Correlate the sample axis with all genes while focusing on the actually affected samples with high abs(generalizedSignature.sampleAxis) and on the already learned current signature via abs(generalizedSignature.R4P):
          if(true)
            [newR4Gs, newsampleSizes4R4Gs] = uncenteredWeightedCorrelation(... %nanmean for the first step where generalizedSignature.R4G ist still nan(nG,1)
               cat(1,signatures.sampleAxis_origUnits)...
              ,eDState.current.(fn4L2Rs_origUnits) ...
              ,2 ...
              ,abs(cat(1,signatures.signedFocusedW4P_withPerpendicularSpace)) ...+100*eps(inInfo.preprocessing.numericTargetPrecision) ... %+100*eps against degenerate all=1 or all=-1 correlations (due to a 1:0 binary signature focus and synthethic noise-free data); in this case we need a bit information from outside the signature focus to get a proper correlation eigenOrder (not all exactly equal to 1 or -1).
              ,BM.meanW ...
              ,false ...bJustDiagonal ...
              ,iif(bOnlyOneSidedPerfSubspaceToFindCorrelatedButWeak...
                ,{true(nG,1),BsPerformanceSubspace4Correlations{2}} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                  ...%<-dev.note: this is important to let the performance subspace grow to significantly correlated genes/samples as per .performance.minAbsR4inclusionInCorrComputation in the next iteration; otherwise, the perf subspace will never grow and might be a too narrow view on the signal.
                ,BsPerformanceSubspace4Correlations... %retry with full perf subspace in step2 (especially means that no corr but low sig genes can be found as members! still those found in corrMoments of step1 should all be available.
               )...
            );
            [signatures.R4G] = ilC2A(num2cell(...
               mean(cat(3, cat(2,signatures.R4G), newR4Gs),3)... %converge smoothly without oscillations.
            ,1)); 
              if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                [signatures.R4G] = ilC2A(num2cell(...
                  sign(cat(2,signatures.R4G)).*dampenOutliers(abs(cat(2,signatures.R4G)), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G(min(end,sL)), 1) ... %<-dampen abs(X) (for less total dampening interval)
                ,1));
              end

            [signatures.sampleSizes4R4G] = ilC2A(num2cell(...
               mean(cat(3, cat(2,signatures.sampleSizes4R4G), newsampleSizes4R4Gs),3) ... %estimated DOFs for the averaged Rs
            ,1)); 
          end
        %Correlate the gene axis with all samples while focusing on the actually affected genes with high abs(generalizedSignature.geneAxis) and on the already learned current signature via abs(generalizedSignature.R4G):
          if(true)
            [newR4Ps, newsampleSizes4R4Ps] = uncenteredWeightedCorrelation(... %nanmean for the first step where generalizedSignature.R4P ist still nan(1,nP)
               cat(2,signatures.geneAxis_origUnits) ...
              ,eDState.current.(fn4L2Rs_origUnits) ...
              ,1 ...
              ,abs(cat(2,signatures.signedFocusedW4G_withPerpendicularSpace)) ...+100*eps(inInfo.preprocessing.numericTargetPrecision) ... %+100*eps against degenerate all=1 or all=-1 correlations (due to a 1:0 binary signature focus and synthethic noise-free data); in this case we need a bit information from outside the signature focus to get a proper correlation eigenOrder (not all exactly equal to 1 or -1).
              ,BM.meanW ...
              ,false ...bJustDiagonal ...
              ,iif(bOnlyOneSidedPerfSubspaceToFindCorrelatedButWeak...
                ,{BsPerformanceSubspace4Correlations{1},true(1,nP)} ... %use subspace in corrDim, but compute R4P for all samples including noisy ones; 
                  ...%<-dev.note: this is important to let the performance subspace grow to significantly correlated genes/samples as per .performance.minAbsR4inclusionInCorrComputation in the next iteration; otherwise, the perf subspace will never grow and might be a too narrow view on the signal.
                ,BsPerformanceSubspace4Correlations... %retry with full perf subspace in step2 (especially means that no corr but low sig genes can be found as members! still those found in corrMoments of step1 should all be available.
               )...
            );
            [signatures.R4P] = ilC2A(num2cell(...
              mean(cat(3, cat(1,signatures.R4P), newR4Ps),3)... %converge smoothly without oscillations.
            ,2)); 
              if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                [signatures.R4P] = ilC2A(num2cell(...
                  sign(cat(1,signatures.R4P)).*dampenOutliers(abs(cat(1,signatures.R4P)), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2) ... %<-dampen abs(X) (for less total dampening interval)
                ,2));
              end
              
            [signatures.sampleSizes4R4P] = ilC2A(num2cell(...
               mean(cat(3, cat(1,signatures.sampleSizes4R4P), newsampleSizes4R4Ps),3) ... %estimated DOFs for the averaged Rs
            ,2)); 
          end
        %Significance of the weighted correlations: 
          [signatures.P4R4G] = ilC2A(num2cell(pValues4Correlations(...
            cat(2,signatures.R4G), cat(2,signatures.sampleSizes4R4G)...
          ),1));
          [signatures.P4R4P] = ilC2A(num2cell(pValues4Correlations(...
            cat(1,signatures.R4P), cat(1,signatures.sampleSizes4R4P)...
          ),2));
        %Update signature focus based on current correlations:
          newSignedW4G = calcSignatureFocus([],[], cat(2,signatures.R4G), cat(2,signatures.P4R4G), 1, inInfo.searchStrategy.signatureFocus, sL);
          newSignedW4P = calcSignatureFocus([],[], cat(1,signatures.R4P), cat(1,signatures.P4R4P), 2, inInfo.searchStrategy.signatureFocus, sL);
          %Add .full2focusWeightRatio weights from the current and the twin space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the twin space:
            newSignedW4P_withPerpendicularSpace = addFlatWeights(newSignedW4G, newSignedW4P, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
            newSignedW4G_withPerpendicularSpace = addFlatWeights(newSignedW4G, newSignedW4P, 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
          [signatures.signedFocusedW4G] = ilC2A(num2cell(newSignedW4G,1));
          [signatures.signedFocusedW4P] = ilC2A(num2cell(newSignedW4P,2));
          [signatures.signedFocusedW4G_withPerpendicularSpace] = ilC2A(num2cell(newSignedW4G_withPerpendicularSpace,1));
          [signatures.signedFocusedW4P_withPerpendicularSpace] = ilC2A(num2cell(newSignedW4P_withPerpendicularSpace,2));
        %Update extended signature focus:
          [signatures.signedExtendedW4G] = ilC2A(num2cell(...
            calcSignatureFocus([],[], cat(2,signatures.R4G), cat(2, signatures.P4R4G), 1, inInfo.searchStrategy.extendedFocus, sL) ...
          ,1));
          [signatures.signedExtendedW4P] = ilC2A(num2cell(...
            calcSignatureFocus([],[], cat(1,signatures.R4P), cat(1, signatures.P4R4P), 2, inInfo.searchStrategy.extendedFocus, sL) ...
          ,2));
        %Update signature size estimations:
          bIntegrateFlag = true;
          [signatures.signatureSizeByCorrSum4G] = ilC2A(num2cell(BM.meanW(...
             abs(cat(2,signatures.signedExtendedW4G))...
            ,inInfo.reference.rowW...
          ,1, bIntegrateFlag), 1));
          [signatures.signatureSizeByCorrSum4P] = ilC2A(num2cell(BM.meanW(...
             abs(cat(1,signatures.signedExtendedW4P))...
            ,inInfo.reference.colW...
          ,2, bIntegrateFlag), 2));
          [signatures.signatureSizeByCorrSum2D] = ilC2A(num2cell([signatures.signatureSizeByCorrSum4G].*[signatures.signatureSizeByCorrSum4P]));
        %Update signature abs correlations in their new extended foci:
          for k=1:length(signatures)
            signatures(k).signatureCorrInExtendedFocus = (...
               signatures(k).signatureSizeByCorrSum4P*BM.meanW(abs(signatures(k).R4G), abs(signatures(k).signedExtendedW4G), 1) ...
             + signatures(k).signatureSizeByCorrSum4G*BM.meanW(abs(signatures(k).R4P), abs(signatures(k).signedExtendedW4P), 2) ...
            )/(signatures(k).signatureSizeByCorrSum4G+signatures(k).signatureSizeByCorrSum4P);
          end
      end
  end

