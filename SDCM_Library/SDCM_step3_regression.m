%ABSTRACT
% Subfunction for SDCM step 3: Regress the 2D bimonotonic 
% baseline signal for the detected signature.

  function [regressionResults, forPlots] = SDCM_step3_regression(signatureDefinition, inInfo)
    %% Initialize:
      global eDState;
      nG = signatureDefinition.reference.nG;
      nP = signatureDefinition.reference.nP;
      sL = min(eDState.current.multiPassesLevel, inInfo.searchStrategy.nDefinedPasses);
      %Get basic math function handles:
        BM = getBasicMathFunctions(inInfo.preprocessing.bDataContainsNaNs);
      %Define the performance subspace:
        BsPerformanceSubspace4Correlations = selectDimensionsInPerformanceSubspace(nG,nP ...
         ,signatureDefinition.step2_finalSignatureAxes.R4G, signatureDefinition.step2_finalSignatureAxes.R4P ...
         ,eDState.noiseEstimation.log10PNoise4G,eDState.noiseEstimation.log10PNoise4P ...
         ,inInfo.searchStrategy ...
        );
      %Initialize outputs:
        signatureDefinition.step3_regression = struct(); %contains .eigenSI, .eigenSJ based on configured order metrics.
        if(~isfield(signatureDefinition,'forPlots')) 
          signatureDefinition.forPlots = struct(); 
        end
        %Copy needed fields form step 2:
          if(true)
            signatureDefinition.step3_regression.R4G = signatureDefinition.step2_finalSignatureAxes.R4G;
            signatureDefinition.step3_regression.P4R4G = signatureDefinition.step2_finalSignatureAxes.P4R4G;
            signatureDefinition.step3_regression.signedFocusedW4G = signatureDefinition.step2_finalSignatureAxes.signedFocusedW4G;
            signatureDefinition.step3_regression.signedFocusedW4G_withPerpendicularSpace = signatureDefinition.step2_finalSignatureAxes.signedFocusedW4G_withPerpendicularSpace;
            signatureDefinition.step3_regression.signedExtendedW4G = signatureDefinition.step2_finalSignatureAxes.signedExtendedW4G;
            signatureDefinition.step3_regression.signatureSizeByCorrSum4G = signatureDefinition.step2_finalSignatureAxes.signatureSizeByCorrSum4G;

            signatureDefinition.step3_regression.R4P = signatureDefinition.step2_finalSignatureAxes.R4P;
            signatureDefinition.step3_regression.P4R4P = signatureDefinition.step2_finalSignatureAxes.P4R4P;
            signatureDefinition.step3_regression.signedFocusedW4P = signatureDefinition.step2_finalSignatureAxes.signedFocusedW4P;
            signatureDefinition.step3_regression.signedFocusedW4P_withPerpendicularSpace = signatureDefinition.step2_finalSignatureAxes.signedFocusedW4P_withPerpendicularSpace;
            signatureDefinition.step3_regression.signedExtendedW4P = signatureDefinition.step2_finalSignatureAxes.signedExtendedW4P;
            signatureDefinition.step3_regression.signatureSizeByCorrSum4P = signatureDefinition.step2_finalSignatureAxes.signatureSizeByCorrSum4P;

            signatureDefinition.step3_regression.log10_p4Correlations = signatureDefinition.step2_finalSignatureAxes.log10_p4Correlations;
            signatureDefinition.step3_regression.signatureSizeByCorrSum2D = signatureDefinition.step2_finalSignatureAxes.signatureSizeByCorrSum2D;
          end
        %Copy order metrics to step3_regression, if not already present (i.e. if not the default signatureStrengths4* are configured):
          if(isfield(signatureDefinition.step2_finalSignatureAxes, signatureDefinition.reference.metric4geneOrder))
            signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder) = signatureDefinition.step2_finalSignatureAxes.(signatureDefinition.reference.metric4geneOrder);
          end
          if(isfield(signatureDefinition.step2_finalSignatureAxes, signatureDefinition.reference.metric4sampleOrder))
            signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder) = signatureDefinition.step2_finalSignatureAxes.(signatureDefinition.reference.metric4sampleOrder);
          end

    %% Calculate the signature eigensignal via isotonic regression and Fourier2D-based smoothening in the space spanned by configured signature order metrics (i.e. by projections on signature axes in the default):
      SDCM_printStatus(2 ...
        ,'   - Use the configured signature order metrics .%s and .%s to span a 2D space with an equidistant grid of the configured resolution %d*%d, then downscale the current signal to it for fast processing.\n'...
        ,signatureDefinition.reference.metric4geneOrder, signatureDefinition.reference.metric4sampleOrder...
        ,inInfo.dissection.signatureSpace.resolution4G, inInfo.dissection.signatureSpace.resolution4P...
      );
      bRegressionAndSignatureSelectionConverged = false; nRegressionIteration = 0; 
      while(~bRegressionAndSignatureSelectionConverged) nRegressionIteration = nRegressionIteration + 1;
        %Get order metrics, the current raw signal in the signature's eigenOrder and determine signature weights for regression:
          %Refine signature strenghts by projecting on the previously regressed gene/sample curves (instead of onto the linear gene/sample axes):
            if(nRegressionIteration > 1)
              %Get the latest regressed signature eigensignal in original signal units:
                signatureEigensignalOSOrigUnits = signatureEigensignalOS;
              signatureDefinition.step3_regression.signatureStrengths4G = projectW(...
                eDState.current.L2Rs ... vectorsInSourceSpace...
               ,bsxfun(@times, sign(signatureDefinition.step3_regression.R4G), signatureEigensignalOSOrigUnits) ... axesInSourceSpace
               ,abs(signatureDefinition.step3_regression.signedFocusedW4P) ... weights4sourceDims
               ,2 ... sourceSpaceMatrixDim
              ,BM.meanW, BM.euclidW) / BM.euclidW(ones(1,nP),abs(signatureDefinition.step3_regression.signedFocusedW4P),2); %normalize such that every compoenent is a representative expression independent of the size of the aggregated space.
              signatureDefinition.step3_regression.signatureStrengths4P = projectW(...
                eDState.current.L2Rs ... vectorsInSourceSpace...
               ,bsxfun(@times, sign(signatureDefinition.step3_regression.R4P), signatureEigensignalOSOrigUnits) ... axesInSourceSpace
               ,abs(signatureDefinition.step3_regression.signedFocusedW4G) ... weights4sourceDims
               ...,W2DOR ...
               ,1 ... sourceSpaceMatrixDim
              ,BM.meanW, BM.euclidW) / BM.euclidW(ones(nG,1),abs(signatureDefinition.step3_regression.signedFocusedW4G),1); %normalize such that every compoenent is a representative expression independent of the size of the aggregated space.
            end
          %Update correlations and the signature focus:
            if(inInfo.dissection.bimonotonicRegression.bUpdateSignatureFocusInEachSignatureCurveRegressionIteration) %TODO: this is unneeded; just use step 2 results and deprecate this for speedup and ease of algorithm. 
              %Update correlations:
                if(true)
                  signatureDefinition.step3_regression.R4G = uncenteredWeightedCorrelation(... %nanmean for the first step where generalizedSignature.R4G ist still nan(nG,1)
                     signatureDefinition.step2_finalSignatureAxes.sampleAxis_origUnits ...
                    ,eDState.current.L2Rs ...
                    ,2 ...
                    ,abs(signatureDefinition.step3_regression.signedFocusedW4P_withPerpendicularSpace) ...
                    ,BM.meanW ...
                    ,false ...bJustDiagonal ...
                    ,{true(nG,1),BsPerformanceSubspace4Correlations{2}} ... %use subspace in corrDIm, but compute R4P for all samples including noisy ones.
                  );

                  signatureDefinition.step3_regression.R4P = uncenteredWeightedCorrelation(... %nanmean for the first step where generalizedSignature.R4P ist still nan(1,nP)
                     signatureDefinition.step2_finalSignatureAxes.geneAxis_origUnits ...
                    ,eDState.current.L2Rs ...
                    ,1 ...
                    ,abs(signatureDefinition.step3_regression.signedFocusedW4G_withPerpendicularSpace) ...
                    ,BM.meanW ...
                    ,false ...bJustDiagonal ...
                    ,{BsPerformanceSubspace4Correlations{1},true(1,nP)} ... %use subspace in corrDIm, but compute R4P for all samples including noisy ones.
                  );
                end
                %<-TODO/performance: we use signature axes as correlation partners here => initial Rs are identical to the ones already computed => no need to recompute. Besides, correlations are only used for initial slope signs in step 3...
              %Comptue P4Rs:
                if(true)
                    [signatureDefinition.step3_regression.P4R4G, log10_p4Correlations4G] = pValues4Correlations(...
                      signatureDefinition.step3_regression.R4G ... 
                     ,signatureDefinition.step2_finalSignatureAxes.sampleSizes4R4G ... 
                    ,1);

                    [signatureDefinition.step3_regression.P4R4P, log10_p4Correlations4P] = pValues4Correlations(...
                      signatureDefinition.step3_regression.R4P ... 
                     ,signatureDefinition.step2_finalSignatureAxes.sampleSizes4R4P ... 
                    ,2);
                      signatureDefinition.step3_regression.log10_p4Correlations = min(log10_p4Correlations4G, log10_p4Correlations4P); %these are two indirectly dependent p values (same L2Rs base); the min should be an upper bound of the combined p value; too conservative?!
                end
              %Update the signature focus (.signedFocusedW4*):
                if(true)
                  %Update extended signature focus (also needed for signature size estimation):
                    if(true)
%                     signatureDefinition.step3_regression.signedExtendedW4G = corrMoments(...
%                        10.^eDState.noiseEstimation.log10PNoise4G ...
%                       ,signatureDefinition.step3_regression.R4G, signatureDefinition.step3_regression.P4R4G ...
%                     ,1, inInfo.searchStrategy.correlationMoments, sL);
                    signatureDefinition.step3_regression.signedExtendedW4G = calcSignatureFocus(...
                       [],[]...
                      ,signatureDefinition.step3_regression.R4G, signatureDefinition.step3_regression.P4R4G ...
                    ,1, inInfo.searchStrategy.extendedFocus, sL);

%                     signatureDefinition.step3_regression.signedExtendedW4P = corrMoments(...
%                        10.^eDState.noiseEstimation.log10PNoise4P ...
%                       ,signatureDefinition.step3_regression.R4P, signatureDefinition.step3_regression.P4R4P ...
%                     ,2, inInfo.searchStrategy.correlationMoments, sL);
                    signatureDefinition.step3_regression.signedExtendedW4P = calcSignatureFocus(...
                       [],[] ...
                      ,signatureDefinition.step3_regression.R4P, signatureDefinition.step3_regression.P4R4P ...
                    ,2, inInfo.searchStrategy.extendedFocus, sL);
                    end
                    %Update Correlation moments of the effects (i.e. effect size and its correlation strength) via simple integration:
                      signatureDefinition.step3_regression.signatureSizeByCorrSum4G = BM.meanW(abs(signatureDefinition.step3_regression.signedExtendedW4G), inInfo.reference.rowW, 1, true);
                      signatureDefinition.step3_regression.signatureSizeByCorrSum4P = BM.meanW(abs(signatureDefinition.step3_regression.signedExtendedW4P), inInfo.reference.colW, 2, true);
                      signatureDefinition.step3_regression.signatureSizeByCorrSum2D = signatureDefinition.step3_regression.signatureSizeByCorrSum4G * signatureDefinition.step3_regression.signatureSizeByCorrSum4P;
                  %Update the signature focus:
                    if(true)
                        signatureDefinition.step3_regression.signedFocusedW4G = calcSignatureFocus(...
                           [], [] ... signatureDefinition.step3_regression.signatureStrengths4G, 10.^eDState.noiseEstimation.log10PNoise4G...
                          ,signatureDefinition.step3_regression.R4G, signatureDefinition.step3_regression.P4R4G...
                          ...,signatureDefinition.step2_finalSignatureAxes.corrMoment4G ...
                        ,1, inInfo.searchStrategy.signatureFocus, sL);

                        signatureDefinition.step3_regression.signedFocusedW4P = calcSignatureFocus(...
                           [],[] ...signatureDefinition.step3_regression.signatureStrengths4P, 10.^eDState.noiseEstimation.log10PNoise4P...
                          ,signatureDefinition.step3_regression.R4P, signatureDefinition.step3_regression.P4R4P...
                          ...,signatureDefinition.step2_finalSignatureAxes.corrMoment4P ...
                        ,2, inInfo.searchStrategy.signatureFocus, sL);
                    end
                    assert(sum(abs(signatureDefinition.step3_regression.signedFocusedW4G))>0, 'code validation: all .signedFocusedW4G==0');
                    assert(sum(abs(signatureDefinition.step3_regression.signedFocusedW4P))>0, 'code validation: all .signedFocusedW4P==0');
                  %Add .full2focusWeightRatio weights from the current and the dual space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the dual space:
                    if(true)
                      signatureDefinition.step3_regression.signedFocusedW4P_withPerpendicularSpace = addFlatWeights(signatureDefinition.step3_regression.signedFocusedW4G, signatureDefinition.step3_regression.signedFocusedW4P, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
                      signatureDefinition.step3_regression.signedFocusedW4G_withPerpendicularSpace = addFlatWeights(signatureDefinition.step3_regression.signedFocusedW4G, signatureDefinition.step3_regression.signedFocusedW4P, 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
                    end
                end
            end
          %Compute/update the gene and sample eigenOrder (updates in nRegressionIteration>1 are needed in case signature strengths updated above are used as eigenOrder metric):
            if(true)
              signatureDefinition.reference.metric4geneOrder = inInfo.dissection.metric4geneOrder;
              signatureDefinition.reference.metric4sampleOrder = inInfo.dissection.metric4sampleOrder;
              %Map any NaNs in the eigenOrder metrics to zero with warning (e.g. signatureStrengths for allNaNRows):
                if(any(isnan(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder))))
                  warning('%d genes had NaN as eigenOrder metric %s; now mapping to zero.', sum(isnan(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder))), signatureDefinition.reference.metric4geneOrder);
                  signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder)(isnan(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder))) = 0;
                end
                if(any(isnan(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder))))
                  warning('%d samples had NaN as eigenOrder metric %s; now mapping to zero.', sum(isnan(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder))), signatureDefinition.reference.metric4sampleOrder);
                  signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder)(isnan(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder))) = 0;
                end

              signatureDefinition.step3_regression.eigenSI = ilv(...
                signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder) ...
               ,@(M)ilnth(2,@sort,M)...
              );
              signatureDefinition.step3_regression.eigenSJ = ilv(...
                signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder) ...
               ,@(M)ilnth(2,@sort,M)...
              );
            end

        %!Bimonotonically regress the signature signal in signature strength order:
          if(true)
            %Transform source t values, if configured:
              trafo4G = inInfo.dissection.signatureSpace.trafo4G;
              trafo4P = inInfo.dissection.signatureSpace.trafo4P;
            %axes and meshgrid:
              orderMetric4GES = trafo4G(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder)(signatureDefinition.step3_regression.eigenSI));
              orderMetric4PES = trafo4P(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder)(signatureDefinition.step3_regression.eigenSJ));
              assert(all(diff(orderMetric4GES)>=0), 'code validation: inInfo.dissection.signatureSpace.trafo4G(signatureDefinition.reference.metric4geneOrder) is not monotonic in signatureDefinition.step3_regression.eigenSI');
              assert(all(diff(orderMetric4PES)>=0), 'code validation: inInfo.dissection.signatureSpace.trafo4P(signatureDefinition.reference.metric4sampleOrder) is not monotonic in signatureDefinition.step3_regression.eigenSJ');
            %Sorted prepared signal:
              L2Rs4dissection_preparedOS = eDState.current.L2Rs;
              L2Rs4dissection_preparedES = L2Rs4dissection_preparedOS(signatureDefinition.step3_regression.eigenSI, signatureDefinition.step3_regression.eigenSJ);
              %Save the previous regression result for convergence checking:
                if(nRegressionIteration==1)
                  signatureEigensignalOS_previous = L2Rs4dissection_preparedOS;
                else
                  signatureEigensignalOS_previous = signatureEigensignalOS;
                end
            %enforce monotonicity via isotonic regression to prevent information loss about superposed signatures:
              assert(sum(signatureDefinition.step3_regression.signedFocusedW4G~=0)>=2, 'code validation: less than two genes have nonzero weight => cannot regress');
              assert(sum(signatureDefinition.step3_regression.signedFocusedW4P~=0)>=2, 'code validation: less than two samples have nonzero weight => cannot regress');
              bMexVersionAvailable = ~isempty(which('GPAV_static_mex'));
              SDCM_printStatus(2-double(~bMexVersionAvailable),'   - Starting bimonotonic regression iteration %d...\n', nRegressionIteration);
              %Performance: define subsampling grid for isotonic regression adaptive to signature strengths (no need to regress in the full signal space, if the signal does not change much over the signature order anyway...):
                if(true)
                  %Distribute nG4subSampling many subII4GES via an equidistant grid over abs(.signatureStrengths4G) in eigenOrder excluding BsPerformanceSubspace4Correlations{1}:
                    [signatureSignalStrengths4GES,SI] = sort(abs(signatureDefinition.step3_regression.signatureStrengths4G)); 
                      if(inInfo.internal.bEnablePostProductionCodeOptimizations)
                        [signatureSignalStrengths4GES,SI] = sort(abs(signatureDefinition.step3_regression.signatureStrengths4G).*(signatureDefinition.step3_regression.signedFocusedW4G~=0)); 
                        %<-TODO/toCheck: can genes despite having zero weight in the signature be useful support points for regression down to / near zero?
                      end
                      signatureSignalStrengths4GES(isnan(signatureSignalStrengths4GES)) = 0;
                    absSignalDeltas4G = abs(diff(signatureSignalStrengths4GES));
                      absSignalDeltas4G(isnan(absSignalDeltas4G)) = 0;
                    adaptiveIntervalMasses4G = cumsum([0;absSignalDeltas4G]);
                      adaptiveIntervalMasses4G = adaptiveIntervalMasses4G/(eps+adaptiveIntervalMasses4G(end))*nG;
                    nG4subSampling = max(...
                      inInfo.dissection.bimonotonicRegression.maxResolutionOnSignalDeltaGrid4G ...
                     ,signatureDefinition.step3_regression.signatureSizeByCorrSum4G ... %for very broad signatures with signatureSizeByCorrSum>>resolution it may be unprecise to just use the configured resolution => adapt to signature size.
                    ); 
                    BDistinct = [true;diff(adaptiveIntervalMasses4G)>eps(adaptiveIntervalMasses4G(end))];
                    if(sum(BDistinct&BsPerformanceSubspace4Correlations{1}(SI))>2)
                      subII4G = interp1(...
                         adaptiveIntervalMasses4G(BDistinct&BsPerformanceSubspace4Correlations{1}(SI))...
                        ,ilsub(1:nG,BDistinct&BsPerformanceSubspace4Correlations{1}(SI))...
                        ,linspace(min(adaptiveIntervalMasses4G(BDistinct&BsPerformanceSubspace4Correlations{1}(SI))), max(adaptiveIntervalMasses4G(BDistinct&BsPerformanceSubspace4Correlations{1}(SI))), nG4subSampling) ...
                        ,'linear' ...
                      );
                    else
                      warning('BDistinct&BsPerformanceSubspace4Correlations{1}(SI) has <=2 nonzero elems; trying just BDistinct');
                      subII4G = interp1(...
                         adaptiveIntervalMasses4G(BDistinct)...
                        ,ilsub(1:nG,BDistinct)...
                        ,linspace(min(adaptiveIntervalMasses4G(BDistinct)), max(adaptiveIntervalMasses4G(BDistinct)), nG4subSampling)...
                        ,'linear' ...
                      );
                    end
                      subII4G = unique(round(subII4G)); %Remove duplicates.
                      subB4GOS = ismember((1:nG)',SI(subII4G));
                      subII4GES = find(subB4GOS(signatureDefinition.step3_regression.eigenSI));

                  %Distribute nP4subSampling many subII4P via an equidistant grid over abs(.signatureStrengths4P) in eigenOrder excluding BsPerformanceSubspace4Correlations{2}:
                    [signatureSignalStrengths4PES,SJ] = sort(abs(signatureDefinition.step3_regression.signatureStrengths4P));
                      if(inInfo.internal.bEnablePostProductionCodeOptimizations)
                        [signatureSignalStrengths4PES,SJ] = sort(abs(signatureDefinition.step3_regression.signatureStrengths4P).*(signatureDefinition.step3_regression.signedFocusedW4P~=0));
                        %<-TODO/toCheck: can genes despite having zero weight in the signature be useful support points for regression down to / near zero?
                      end
                      signatureSignalStrengths4PES(isnan(signatureSignalStrengths4PES)) = 0;
                    absSignalDeltas4P = abs(diff(signatureSignalStrengths4PES));
                      absSignalDeltas4P(isnan(absSignalDeltas4P)) = 0;
                    adaptiveIntervalMasses4P = cumsum([0,absSignalDeltas4P]);
                      adaptiveIntervalMasses4P = adaptiveIntervalMasses4P/adaptiveIntervalMasses4P(end)*nP;
                    nP4subSampling = max(...
                      inInfo.dissection.bimonotonicRegression.maxResolutionOnSignalDeltaGrid4P ...
                     ,signatureDefinition.step3_regression.signatureSizeByCorrSum4P ... %for very broad signatures with signatureSizeByCorrSum>>resolution it may be unprecise to just use the configured resolution => adapt to signature size.
                    ); 
                    BDistinct = [true,diff(adaptiveIntervalMasses4P)>eps(absSignalDeltas4P(end))];
                    if(sum(BDistinct&BsPerformanceSubspace4Correlations{2}(SJ))>2)
                      subII4P = interp1(...
                         adaptiveIntervalMasses4P(BDistinct&BsPerformanceSubspace4Correlations{2}(SJ))...
                        ,ilsub(1:nP,BDistinct&BsPerformanceSubspace4Correlations{2}(SJ))...
                        ,linspace(min(adaptiveIntervalMasses4P(BDistinct&BsPerformanceSubspace4Correlations{2}(SJ))),max(adaptiveIntervalMasses4P(BDistinct&BsPerformanceSubspace4Correlations{2}(SJ))),nP4subSampling)...
                        ,'linear' ...
                      );
                    else
                      warning('BDistinct&BsPerformanceSubspace4Correlations{1}(SI) has <=2 nonzero elems; trying just BDistinct');
                      subII4P = interp1(...
                         adaptiveIntervalMasses4P(BDistinct)...
                        ,ilsub(1:nP,BDistinct)...
                        ,linspace(min(adaptiveIntervalMasses4P(BDistinct)),max(adaptiveIntervalMasses4P(BDistinct)),nP4subSampling)...
                        ,'linear' ...
                      );
                    end
                      subII4P = unique(round(subII4P)); %Remove duplicates.
                      subB4POS = ismember(1:nP,SJ(subII4P));
                      subII4PES = find(subB4POS(signatureDefinition.step3_regression.eigenSJ));

                  %Exclude all-NaN rows/cols, if any:
                    subII4GES = setdiff(subII4GES,find(all(isnan(L2Rs4dissection_preparedES),2)));
                    subII4PES = setdiff(subII4PES,find(all(isnan(L2Rs4dissection_preparedES),1)));
                else
                  warning('dev/performance: reenable performance subspace for monotonic regression');
                  subII4GES = 1:nG;
                  subII4PES = 1:nP;
                  subII4GES = setdiff(subII4GES,find(all(isnan(L2Rs4dissection_preparedES),2)));
                  subII4PES = setdiff(subII4PES,find(all(isnan(L2Rs4dissection_preparedES),1)));
                end
              %!Compute monotonic regression in the subsampled signal strength space:
                if(true)
                  %Prescribe directions of monotonicity (increasing/deacreasing) based on correlations, if sufficiently certain:
                    slopeSigns4GES = signatureDefinition.step3_regression.R4G(signatureDefinition.step3_regression.eigenSI);
                    slopeSigns4PES = signatureDefinition.step3_regression.R4P(signatureDefinition.step3_regression.eigenSJ);
                    minAbsRForSlopePrescription = 0.5;
                      slopeSigns4GES(abs(slopeSigns4GES)<minAbsRForSlopePrescription) = NaN; %let the regression find the best for unsure corrs (=>takes twice as long).
                      slopeSigns4PES(abs(slopeSigns4PES)<minAbsRForSlopePrescription) = NaN; %let the regression find the best for unsure corrs (=>takes twice as long).
                    slopeSigns4GES = sign(slopeSigns4GES);
                    slopeSigns4PES = sign(slopeSigns4PES);
                  %Get 2D regression weights:
                    if(true)
                      weights2DES = bsxfun(@times...
                        ,abs(signatureDefinition.step3_regression.signedFocusedW4G(signatureDefinition.step3_regression.eigenSI)) .* BsPerformanceSubspace4Correlations{1}(signatureDefinition.step3_regression.eigenSI) ...
                        ,abs(signatureDefinition.step3_regression.signedFocusedW4P(signatureDefinition.step3_regression.eigenSJ)) .* BsPerformanceSubspace4Correlations{2}(signatureDefinition.step3_regression.eigenSJ) ...
                      );
                      if(inInfo.internal.bEnablePostProductionCodeOptimizations)
                        weights2DES = bsxfun(@times...
                          ,    abs(signatureDefinition.step3_regression.signedFocusedW4G(signatureDefinition.step3_regression.eigenSI))...
                            .* (BsPerformanceSubspace4Correlations{1}(signatureDefinition.step3_regression.eigenSI) | ismember((1:nG)',subII4GES)) ... %do not exclude elements of the subsampling grid for regression.
                          ,    abs(signatureDefinition.step3_regression.signedFocusedW4P(signatureDefinition.step3_regression.eigenSJ))...
                            .* (BsPerformanceSubspace4Correlations{2}(signatureDefinition.step3_regression.eigenSJ) | ismember(1:nP, subII4PES)) ... %do not exclude elements of the subsampling grid for regression.
                        );
                      end
                    end
                  [signatureEigensignalES, signatureDefinition] = bimonotonicRegression(...
                       L2Rs4dissection_preparedES ...
                   ...prescribed directions of monotonicity:
                      ,slopeSigns4GES, slopeSigns4PES ... 
                   ...%order metrics (used for smoothening and linear interpolation at NaNs)
                      ,orderMetric4GES, orderMetric4PES ... 
                   ...%regression weights:
                      ,weights2DES ...
                   ...%performance subgrid:
                      ,subII4GES, subII4PES ... 
                   ,signatureDefinition, inInfo ...
                  ); %iterative monotonic regression of each gene and sample in parallel and mixing results until convergence.
                end
              %Unsort back to the global reference order:
                signatureEigensignalOS(...
                  signatureDefinition.step3_regression.eigenSI, signatureDefinition.step3_regression.eigenSJ...
                ) = signatureEigensignalES;                  
                clear signatureEigensignalES; %always use the eigensignal in reference order (to not mix up eigenOrder iterations).
          end
          %Check regression convergence:
            if(nRegressionIteration>1) %use all points before computing convergence break option.
              %Get the extended product focus in 2D, as when computing signature statistics:
                W2D = sqrt(bsxfun(@times, abs(signatureDefinition.step3_regression.signedExtendedW4G), abs(signatureDefinition.step3_regression.signedExtendedW4P)));
              %Performance/only test the top half of the signature for convergence (no need to wait for a dangling curve near zero where the eigenOrder may be noisy/changing with every iteration...)
                BTopSignatureMembers2D = W2D > nanmax(W2D(:))/2; 
              r4TopOfSignatureDissection = uncenteredWeightedCorrelation(... %nanmean for the first step where generalizedSignature.R4P ist still nan(1,nP)
                 signatureEigensignalOS_previous(BTopSignatureMembers2D) ...(BTopSignatureMembers) ...
                ,signatureEigensignalOS(BTopSignatureMembers2D) ...(BTopSignatureMembers) ...
                ,1 ...
                ,W2D(BTopSignatureMembers2D) ...
                ,@weightedMean_supportsNaNs ... %some points may be nan due to inInfo.dissection.ignoreSignalFromNonSignatureMembers.bEnabled of the current or previous regression iteration.
              );
              SDCM_printStatus(1,'  <- Bimonotonic regression iteration %d finished. Correlation to previous iteration =%0.4f (configured threshold for convergence =%0.4f).\n', nRegressionIteration, r4TopOfSignatureDissection, inInfo.dissection.bimonotonicRegression.convergenceThreshold4r4TopOfSignatureDissection); 
            else
              SDCM_printStatus(1,'  <- Bimonotonic regression iteration %d finished.\n', nRegressionIteration); 
              r4TopOfSignatureDissection = NaN;
            end
            bRegressionAndSignatureSelectionConverged = r4TopOfSignatureDissection > inInfo.dissection.bimonotonicRegression.convergenceThreshold4r4TopOfSignatureDissection;
            %in application mode, we apply the externally provided signature just once (no axes or focus propagation is allowed) => a single regression is enough already:
              bRegressionAndSignatureSelectionConverged = bRegressionAndSignatureSelectionConverged || inInfo.applicationMode.bEnabled; 
            %respect .bimonotonicRegression.maxRegressionIterations config:
              if(nRegressionIteration>inInfo.dissection.bimonotonicRegression.maxRegressionIterations && ~bRegressionAndSignatureSelectionConverged)
                warning('Even after %d regression refinement iterations the convergence criterion was not met; now accepting AS IS (correlation to last r=%0.4f). Visually control the signature membership diagram.', inInfo.dissection.bimonotonicRegression.maxRegressionIterations, r4TopOfSignatureDissection); 
                bRegressionAndSignatureSelectionConverged = true;
              end

        %Compute a final Fourier2D-smoothening to polish the regressed signature curve (this local smoothening cannot change monotonicity):
          if(bRegressionAndSignatureSelectionConverged)
            signatureEigensignalES_smoothed = smoothenInSignatureStrengthsSpace(...
              orderMetric4GES, orderMetric4PES ...
             ,signatureEigensignalOS(signatureDefinition.step3_regression.eigenSI, signatureDefinition.step3_regression.eigenSJ)...
             ,inInfo, false ... %inInfo, bSupressStatusmessages
             ,signatureEigensignalOS(signatureDefinition.step3_regression.eigenSI, signatureDefinition.step3_regression.eigenSJ), signatureDefinition... 
            );
            signatureEigensignalOS_smoothed(...
              signatureDefinition.step3_regression.eigenSI, signatureDefinition.step3_regression.eigenSJ...
            ) = signatureEigensignalES_smoothed;
              clear signatureEigensignalES; %always use the eigensignal in reference order (to not mix up eigenOrder iterations).
            signatureEigensignalOS = signatureEigensignalOS_smoothed;
          end
        %Dissection strengths (i.e. the removal mask) for the regressed signature eigensignal:
          if(bRegressionAndSignatureSelectionConverged)
            %Signature dissection strengths (product signature focus):
              signatureDissectionStrengthsOS = sqrt(abs(bsxfun(@times ...
                ,signatureDefinition.step3_regression.signedFocusedW4G ...
                ,signatureDefinition.step3_regression.signedFocusedW4P ...
              )));
            %Apply dissection strengths:
              signatureEigensignalTimesMaskOS = signatureEigensignalOS .* signatureDissectionStrengthsOS;
            %Optionally save some plotting info in the signatureDefinition:
              if(inInfo.export.matFile4eachSignatureDefinition.bIncludeDefPlotPanels || inInfo.plots.signatureDefinitions.bEnabled)
                signatureDefinition.forPlots.overview.current.signatureEigensignalb4applyingSignatureDissectionStrengths = signatureEigensignalOS;
              end
              if(isa(inInfo.export.matFile4eachSignatureDefinition.bIncludeRegressedSignalEstimationSteps,'function_handle'))
                current_bSaveIntermediateResultsForDevPlots = inInfo.export.matFile4eachSignatureDefinition.bIncludeRegressedSignalEstimationSteps(eDState.current.k);
              else
                current_bSaveIntermediateResultsForDevPlots = inInfo.export.matFile4eachSignatureDefinition.bIncludeRegressedSignalEstimationSteps;
              end
                if(current_bSaveIntermediateResultsForDevPlots) %<-save additional data maybe needed by some plots later:
                  signatureDefinition.forPlots.dissectionPipelines.nonEqdInput.trafoG = trafo4G;
                  signatureDefinition.forPlots.dissectionPipelines.nonEqdInput.trafoP = trafo4P;
                  signatureDefinition.forPlots.dissectionPipelines.nonEqdInput.T4P = orderMetric4PES;
                  signatureDefinition.forPlots.dissectionPipelines.nonEqdInput.T4G = orderMetric4GES;
                  signatureDefinition.forPlots.dissectionPipelines.nonEqdInput.L2Rs = L2Rs4dissection;
                  signatureDefinition.forPlots.dissectionPipelines.eqdInterp.T4P = orderMetric4P_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.eqdInterp.T4G = orderMetric4G_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.eqdInterp.L2Rs = signal_eqd;
                  padded_T4P = [orderMetric4P_eqd(end/2+1:end), orderMetric4P_eqd, orderMetric4P_eqd(1:end/2)];
                  padded_T4G = [orderMetric4G_eqd(end/2+1:end); orderMetric4G_eqd; orderMetric4G_eqd(1:end/2)];
                  signatureDefinition.forPlots.dissectionPipelines.paddedEqd.T4P = padded_T4P;
                  signatureDefinition.forPlots.dissectionPipelines.paddedEqd.T4G = padded_T4G;
                  signatureDefinition.forPlots.dissectionPipelines.paddedEqd.L2Rs = padded_convolution_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.paddedEqdConvoluted.T4P = padded_T4P;
                  signatureDefinition.forPlots.dissectionPipelines.paddedEqdConvoluted.T4G = padded_T4G;
                  signatureDefinition.forPlots.dissectionPipelines.paddedEqdConvoluted.L2Rs = padded_convolution_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.unpaddedEqdConvoluted.T4P = orderMetric4P_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.unpaddedEqdConvoluted.T4G = orderMetric4G_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.unpaddedEqdConvoluted.L2Rs = convolution_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.unpaddedEqdConvolutedIsotonic.T4P = orderMetric4P_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.unpaddedEqdConvolutedIsotonic.T4G = orderMetric4G_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.unpaddedEqdConvolutedIsotonic.L2Rs = bimonotonic_convolution_eqd;
                  signatureDefinition.forPlots.dissectionPipelines.basline4origT.T4P = orderMetric4PES;
                  signatureDefinition.forPlots.dissectionPipelines.basline4origT.T4G = orderMetric4GES;
                  signatureDefinition.forPlots.dissectionPipelines.basline4origT.L2Rs = signatureEigensignalOS;
                  signatureDefinition.forPlots.dissectionPipelines.remainingSignal.T4P = orderMetric4PES;
                  signatureDefinition.forPlots.dissectionPipelines.remainingSignal.T4G = orderMetric4GES;
                  signatureDefinition.forPlots.dissectionPipelines.remainingSignal.L2Rs = L2Rs4dissection-signatureEigensignalOS;
                end
          end
      end
    %% Finish and return the signature's eigensignal:
      %Save the signature''s eigensignal in the signatureDefinition (i.e. the dissection/shift results):
        signatureDefinition.step3_regression.signatureEigensignal = signatureEigensignalTimesMaskOS;
        signatureDefinition.step3_regression.signatureDissectionStrengthsOS = signatureDissectionStrengthsOS; %save for noise estiamtion.
          assert(~any(isnan(signatureDefinition.step3_regression.signatureEigensignal(:))), 'data validation: some values in signatureDefinition.step3_regression.signatureEigensignal were NaN.');
        %Performance: save memory by compressing/removing zero rows/cols:
          if(true)
            signatureDefinition.step3_regression.BNonzeroRows4signatureEigensignal = ~all(signatureDefinition.step3_regression.signatureEigensignal==0,2);
            signatureDefinition.step3_regression.BNonzeroCols4signatureEigensignal = ~all(signatureDefinition.step3_regression.signatureEigensignal==0,1);
            signatureDefinition.step3_regression.signatureEigensignal(~signatureDefinition.step3_regression.BNonzeroRows4signatureEigensignal,:) = [];
            signatureDefinition.step3_regression.signatureEigensignal(:,~signatureDefinition.step3_regression.BNonzeroCols4signatureEigensignal) = [];
            signatureDefinition.step3_regression.signatureDissectionStrengthsOS(~signatureDefinition.step3_regression.BNonzeroRows4signatureEigensignal,:) = [];
            signatureDefinition.step3_regression.signatureDissectionStrengthsOS(:,~signatureDefinition.step3_regression.BNonzeroCols4signatureEigensignal) = [];
          end
      %Use a precise/narrow output interface (instead of returning the whole signature definition structure):
        regressionResults = signatureDefinition.step3_regression;
        forPlots = signatureDefinition.forPlots;
      SDCM_printStatus(2,'  <- finished calculating the signature''s bi-monotonic eigensignal.\n');
        if(inInfo.internal.bDevEnableInteractiveBreaks)
          interactiveCheckPoint();
        end
  end

