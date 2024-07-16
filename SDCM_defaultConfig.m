%ABSTRACT
% Library function for SDCM. Declares the developed general purpose
% default parameter set. Optimizes with the documented test
% signals to minimize false negatives as primary and to minimize
% false positives as secondary priority. Some parameters are
% adaptive to the signal size [nG,nP].
% TODO: QC of explanation texts...
function defaultInfo = SDCM_defaultConfig(nG,nP)
  %Initialize:
    explanations = struct();
  %% Reference info:
    defaultInfo.reference = struct(); %.reference at field position one.
    %Gene and sample IDs/names (for status output and plots):
      explanations.reference.rowIDs = 'n*1 cell array of row identifiers (typically gene or probeset names). Default: gene index numbers starting with an i)';
      explanations.reference.rowGroups = 'n*1 numeric class array determining row groups. If for example eight probesets measure the same gene and a gene level interpretation of parameters like .minAbsCorrSum4G is desired, set the numeric class of these rows to the same value. This will cause them to be weighted by 1/8 each when determining row counts. Default: (1:nG)'', i.e. every row gets the same weight.';
      explanations.reference.colIDs = '1*m cell array of column identifiers (typically sample IDs or sample names). Default: sample index numbers starting with a j';
      explanations.reference.colGroups = '1*m numeric class array determining column groups. If for example four arrays measure the same sample and a sample level interpretation of parameters like .minAbsCorrSum4P is desired, set the numeric class of these columns to the same value. This will cause them to be weighted one fourth each when determining column counts. Default: -(1:nP), i.e. every column gets the same weight.';
      if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
        defaultInfo.reference.rowIDs = cellfun(@(n)sprintf('i%d',n),num2cell((1:nG)'),'UniformOutput',false);
        defaultInfo.reference.rowGroups = (1:nG)';
        defaultInfo.reference.colIDs = cellfun(@(n)sprintf('j%d',n),num2cell((1:nP)),'UniformOutput',false);
        defaultInfo.reference.colGroups = -(1:nP); %use negative numbers to discern cols from rows as the default naming convention.
      end
    %Version info:
      defaultInfo.reference.version = 3.61;
      defaultInfo.reference.versionText = '(c) Michael Grau, April 2016';
  %% Preprocessing options:
    explanations.preprocessing.numericTargetPrecision = 'The numeric precision for the signal and all computations. Must be either single or double. Default: ''double''; use single in order to be more memory-efficient for very large data sets.';
      defaultInfo.preprocessing.numericTargetPrecision = 'double'; 

    explanations.preprocessing.duplicateRows.bErrorIfFound = 'Rows with exactly the same expressions do not carry any information and are typically errors from mapping measurement values to their annotations. Additionally they cause signatures to be found although the acutal information base is weak. Set this flag to prevent accepting such artefacts in the input signal. Default: true.';
      defaultInfo.preprocessing.duplicateRows.bErrorIfFound = true;
    explanations.preprocessing.duplicateRows.bCheckAndRemove = 'Rows with exactly the same expressions do not carry any information and are typically errors from mapping measurement values to their annotations. They can be removed with this flag, but it is recommended to handle this on caller level. Default: false.';
      defaultInfo.preprocessing.duplicateRows.bCheckAndRemove = false;

    explanations.preprocessing.dampenOutliers4G.bEnabled = 'Many statistics like standard deviations or correlations are not robust against far outliers. Therefore it is recommended to check for them and dampen them before signature dissection. This flag enables outlier dampening for each gene row of the input signal. Default: true.';
      defaultInfo.preprocessing.dampenOutliers4G.bEnabled = true;
      explanations.preprocessing.dampenOutliers4G.nTopRanksToCheck = 'Number of top and bottom samples with respect to the sort order of the respective gene row that are checked for being an outlier. Default: min(4,floor(nG/20))';
        if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
          defaultInfo.preprocessing.dampenOutliers4G.nTopRanksToCheck = min(4,floor(nG/20)); %Note: must not be too large to prevent smoothing out true narrow and sharp signatures.
        end
      explanations.preprocessing.dampenOutliers4G.maxAllowedNeighboursRatio = 'Maximum allowed signal ratio relative to the next sample in the sort order of the respective gene row. If it is larger, the sample that is nearer to the sort border is set to the maximum (respectively minimum) expression allowed by this maximum ratio parameter. Default: 2^(1/.dampenOutliers4G.nTopRanksToCheck), i.e. allow a log2(ratio) doubling at the sort border only over the interval of .dampenOutliers4G.nTopRanksToCheck samples.';
        if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
          defaultInfo.preprocessing.dampenOutliers4G.maxAllowedNeighboursRatio = 2^(1/defaultInfo.preprocessing.dampenOutliers4G.nTopRanksToCheck);
        end

    explanations.preprocessing.dampenOutliers4P.bEnabled = 'Many statistics like standard deviations or correlations are not robust against far outliers. Therefore it is recommended to check for them and dampen them before signature dissection. This flag enables outlier dampening for each sample column of the input signal. Default: true.';
      defaultInfo.preprocessing.dampenOutliers4P.bEnabled = true;
      explanations.preprocessing.dampenOutliers4P.nTopRanksToCheck = 'Number of top and bottom genes with respect to the sort order of the respective sample column that are checked for being an outlier. Default: min(4,floor(nP/20))';
        if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
          defaultInfo.preprocessing.dampenOutliers4P.nTopRanksToCheck = min(4,floor(nP/20)); %Note: must not be too large to prevent smoothing out true narrow and sharp signatures.
        end
      explanations.preprocessing.dampenOutliers4P.maxAllowedNeighboursRatio = 'Maximum allowed signal ratio relative to the next gene in the sort order of the respective sample column. If it is larger, the gene that is nearer to the sort border is set to the maximum (respectively minimum) expression allowed by this maximum ratio parameter. Default: 2^(1/.dampenOutliers4P.nTopRanksToCheck), i.e. allow a log2(ratio) doubling at the sort border only over the interval of .dampenOutliers4P.nTopRanksToCheck genes.';
        if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
          defaultInfo.preprocessing.dampenOutliers4P.maxAllowedNeighboursRatio = 2^(1/defaultInfo.preprocessing.dampenOutliers4P.nTopRanksToCheck); %maximal ein L2Rs-Anstieg von 2 innerhalb des Randbereichs.
        end

    explanations.preprocessing.bQuantilenormalizeSamples = 'Quantile normalization (cf. doc quantilenorm) of the input signal columns. For microarray experiments that were not yet normalized to each other, it is recommended to enable this option; otherwise global brightness differences between arrays will be found and need to be identified as lab signatures. Note that this option can be left enabled even if the inputs were already quantile-normalized, since quantile normalization is a fixpoint operator. However, this operation generally breaks the transpose symmetry, i.e. a re-run with flipped roles between genes and samples will return different results (quantilenorm(X)~=quantilenorm(X'')''); therefore it is disabled by default. Default: false.';
      defaultInfo.preprocessing.bQuantilenormalizeSamples = false;

    explanations.preprocessing.bConvertIntoCohortRelativeL2RsAtStart = 'If the inputs are log2(intensities), convert them into log2(ratio)s in order to focus on signatures that differentiate samples. If this flag is enabled, log2(ratio)s are defined relative to cohort mean intensities. The cohort mean baseline is only well-defined for large cohorts; for very small cohorts an external baseline on caller level should be used to directly pass log2(ratio)s; in this case. Note that this preprocessing step generally breaks the transpose symmetry, i.e. a re-run with flipped roles between genes and samples will return different results; therefore it is disabled by default. It is recommended to enable it for raw microarray data (in log2(intensities) units). Otherwise the first detected signature will most likely represent and dissect platform-specific technical offsets. Default: false.';
      defaultInfo.preprocessing.bConvertIntoCohortRelativeL2RsAtStart = false;

    explanations.preprocessing.bSetNaNsToZeroL2R = 'Use this option to set missing values to zero log2(ratios). This allows for faster computation but has lower imputation quality than to detect signatures with NaN-robust aggregation functions and then to impute missing values by nearest 2D-neighbours in the detected signature orders, which is the default for missing values. Default: false.';
      defaultInfo.preprocessing.bSetNaNsToZeroL2R = false;
  %% Configuration of STEP 1) Search strategy: find the gene or sample that is the most promosing first representative of a new signature:
    %Declare defined sensitivity levels in this paramter set:
      explanations.searchStrategy.nDefinedPasses = 'Declares the number of defined sensistivity levels. Every detection parameter can be of size [1,nDefinedPasses] and SDCM uses .parameterName(sensitivityLevel). If a parameter specified less than nDefinedPasses values, .parameterName(min(end,sensitivityLevel)) is used.  Default: 1, i.e. only one common parameter set (one sensitivity level) is defined for all detection passes; this is usually sufficient.';
        defaultInfo.searchStrategy.nDefinedPasses = 1;
    %Define allowed signature dimensions:
      explanations.searchStrategy.allowedSignatureTypes = 'Cell string with allowed elements: ''sampleAxis'' and ''geneAxis''. Declares the allowed types of initial representative for a signature. Allow both if you want a signal dissection that treats genes and samples symmetrically. Default: {''sampleAxis'', ''geneAxis''}, i.e. symmetric.';
        defaultInfo.searchStrategy.allowedSignatureTypes = {'sampleAxis', 'geneAxis'};
    %Search for the next initial representative:
      %Lookahead length for a better score:
        explanations.searchStrategy.nLookahead4HigherScore = 'The number of genes and samples in presort order that are still tested for a being a better initial representative after a qualifying initial representative has been found. Set this number to >nG+nP to essentially force searching for the global optimum in every detection iteration k (very slow and not recommended as the detection order of signatures is usually not important). Default: 200';
          defaultInfo.searchStrategy.nLookahead4HigherScore = 200;
        explanations.searchStrategy.performance.minScoreIncreaseByLookaheadCandidate = 'To not unnecessarily prolong the lookahead phase by candidates with only a slightly better initial representative score, set this to a value of >1. Default: 1.015-';
          defaultInfo.searchStrategy.performance.minScoreIncreaseByLookaheadCandidate = 1.015;
      %Global break conditions (additional to "no gene or sample qualified as initial representative" in step 1):
        explanations.searchStrategy.globalBreakConditions.maxSignaturesToDetect = 'Maximum number of signatures to detect and return. In the default, no hard limit is set; instead qualification thresholds and the signal''s complexity deterine the number of identified signatures. WHen using a hard limit, detection may stop, before the full signal is dissected, i.e. the remaining signal may still contain significant correlations patterns. Default: Inf, i.e. no hard limit';
          defaultInfo.searchStrategy.globalBreakConditions.maxSignaturesToDetect = Inf;
        explanations.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification = 'The number of candidates in presort order tested for qualification as initial representative of a signature before breaking the current pass. If no qualifying candidate is found within this looakahead window, SDCM advances to the next pass (and, if configured, the next sensitivity level; see .nDefinedPasses). Default: 2*.nLookahead4HigherScore';
          defaultInfo.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification = 2*defaultInfo.searchStrategy.nLookahead4HigherScore; %if there are 1000 genes or samples in greedyScore order in a row without any qualifying as initial representative, early break and do not visit the rest.
        explanations.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification4lastPass = 'The number of candidates in presort order tested for qualification in the last pass. If no qualifying candidate is found within this looakahead window, SDCM terminates (global break condition). Set this to Inf to test every gene and sample in the last pass (slow as it tests all noise genes). Default: 0.2*(nG+nP), i.e. 20% of all gene and sample candidates suffice.';
          if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
            defaultInfo.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification4lastPass = min(floor(0.2*(nG+nP)), 2.5*defaultInfo.searchStrategy.globalBreakConditions.maxLookaheadWithoutAnyQualification-1); %2.5 times the lookahead or 20% of the whole dataset, whichever comes first; -1 to not precomute a full chunk for the last needed initial representative in the last iteration.
          end
      %Qualification thresholds for signatures based on a candidate initial representative:
        if(true)
          explanations.searchStrategy.qualification.minAbsCorrSum4G = 'Specifies the minimum amount of genes participating in a signature for the signature to be considered "interesting". Specified in units of the abs(correlations) sum, i.e. as equivalent amount of fully correlated (|r|=1) genes. A correspondingly higher amount of less strongly correlated (|r|<1) genes can also qualify. We define the minimum size of interest adaptive to the signal size, but never require more than 20% of all genes to be correlated stronger than |r|>0.67. Default: min(ceil(0.2*nG)*0.67, 0.5*log2(nP))';
            if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
              defaultInfo.searchStrategy.qualification.minAbsCorrSum4G = min(ceil(0.2*nG)*0.67, 0.5*log2(nP));
            end
          explanations.searchStrategy.qualification.minAbsCorrSum4P = 'Specifies the minimum amount of samples participating in a signature for the signature to be considered "interesting". Specified in units of the abs(correlations) sum, i.e. as equivalent amount of fully correlated (|r|=1) samples. A correspondingly higher amount of less strongly correlated (|r|<1) samples can also qualify. We define the minimum size of interest adaptive to the signal size, but never require more than 20% of all samples to be correlated stronger than |r|>0.67. Default: min(ceil(0.2*nP)*0.67, 0.5*log2(nG))';
            if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
              defaultInfo.searchStrategy.qualification.minAbsCorrSum4P = min(ceil(0.2*nP)*0.67, 0.5*log2(nG));
            end
          explanations.searchStrategy.qualification.minAbsCorrSumSum = 'Specified the minimum amount of genes plus samples participating in a signature. Usually signatures of minimum size in both order dimensions are not of interest. In this case, this parameter should be configured higher than just .minAbsCorrSum4G + .minAbsCorrSum4P. Default: 2.5*(.minAbsCorrSum4G + .minAbsCorrSum4P)';
            if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
              defaultInfo.searchStrategy.qualification.minAbsCorrSumSum = 2.5*(defaultInfo.searchStrategy.qualification.minAbsCorrSum4G + defaultInfo.searchStrategy.qualification.minAbsCorrSum4P);
            end

          explanations.searchStrategy.qualification.alpha4correlations = 'The maximum p-value acceptable for all correlations to a signature. All gene correlations and all sample correlations are each combined (by Kolmogorov Smirnov tests against t CDFs for expected noise correlations). Hence, this p-value is much lower the correlation p-value for a single gene or a single sample already for only moderately sized signatures. Therefore, this alpha can be configured very conservatively (and should be configured so to prevent false positive in form of ery small signatures or signatures with many participants, but very weak correlations). Default: 1e-10';
            defaultInfo.searchStrategy.qualification.alpha4correlations = 1e-10; 
          explanations.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCorrP = 'If true, .alpha4correlations is effectively divided by nG+nP, i.e. by the maximum amount of hypothesis (i.e. initial representative candidates) tested. (Bonferroni correction for multiple hypothesis tests.) Default: true';
            defaultInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCorrP = true;

          explanations.searchStrategy.qualification.alpha4signalStrength = 'The maximum p-value acceptable for the signature''s signal strength. This does not control the p-values of a single gene or sample, but for all (gene, sample) values in the signature focus. Tested by a weighted t-test comparing (gene, sample) values in the signature focus with the estimated noise distribution (only tested for k>1, as no noise distribution is available before the first dissection). Default: 1e-5';
            defaultInfo.searchStrategy.qualification.alpha4signalStrength = 1e-3; %1e-5; 
          explanations.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForSignalP = 'If true, .alpha4signalStrength is effectively divided by nG+nP, i.e. by the maximum amount of hypothesis (i.e. initial representative candidates) tested. (Bonferroni correction for multiple hypothesis tests.) Default: true';
            defaultInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForSignalP = true;

          explanations.searchStrategy.qualification.alpha4combined = 'The maximum p-value acceptable for the signature p-value that combined p-values based on correlations and signal strengths via Fisher''s method. Default: 1 (individual thresholds .alpha4correlations and .alpha4signalStrength are used instead)';
            defaultInfo.searchStrategy.qualification.alpha4combined = 1;
          explanations.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCombinedP = 'If true, .alpha4combined is effectively divided by nG+nP, i.e. by the maximum amount of hypothesis (i.e. initial representative candidates) tested. (Bonferroni correction for multiple hypothesis tests.) Default: true';
            defaultInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCombinedP = true;

          explanations.searchStrategy.qualification.minCorrInExtendedFocus = 'Minimum demanded mean absolute correlation of all participating genes and samples in the signature''s extended focus. Default: 0.4';
            defaultInfo.searchStrategy.qualification.minCorrInExtendedFocus = 0.4;

          %202301: in some input contexts (e.g. quantile axes units), we may want to specify an absolute minSignalStrengthInSignatureFocus (the same as is printed in status outpus as %0.2ffocus)
            explanations.searchStrategy.qualification.minAbsMean2DInExtendedFocus = 'Minimum demanded mean absolute signal in the signature focus (in your input units). Default: 0, i.e. disabled.';
              defaultInfo.searchStrategy.qualification.minAbsMean2DInExtendedFocus = 0; %test condition: signatureCand.signatureAbsMean2D>=.minAbsMean2DInExtendedFocus.
        end
    %Signature computation:
      %Signature focus (initial and refined weights):
        if(true)
          %Initial weights (based on standardized signal strengths):
            explanations.searchStrategy.signatureFocus.relSignal.threshold = 'Threshold paramter in the initial weights formula; see calcSignatureFocus. Default: 0, i.e. no threshold'; 
              defaultInfo.searchStrategy.signatureFocus.relSignal.threshold = 0; 
            explanations.searchStrategy.signatureFocus.relSignal.sufficient4fullWeight = 'Relative signal mapped to full initial weight; see calcSignatureFocus. Default: 0.8, i.e. 80% of the maximum already suffices for full weight'; 
              defaultInfo.searchStrategy.signatureFocus.relSignal.sufficient4fullWeight = 0.8;
            explanations.searchStrategy.signatureFocus.relSignal.exponent = 'Exponent in the initial weights formula; see calcSignatureFocus. Default: 1, i.e. a simple linear ramp.'; 
              defaultInfo.searchStrategy.signatureFocus.relSignal.exponent  = 1;
          %Refined weights (correlation based):
            explanations.searchStrategy.signatureFocus.relCorr.bEnabled = 'Controls whether relative or absolute correlations are used in the (refined) weights formula; see calcSignatureFocus. Default: true, i.e. rescale correlations relative to max(abs(Rs))'; 
              defaultInfo.searchStrategy.signatureFocus.relCorr.bEnabled = true;
            explanations.searchStrategy.signatureFocus.relCorr.threshold = 'Threshold paramter in the refined weights formula; see calcSignatureFocus. Default: 0, i.e. no threshold'; 
              defaultInfo.searchStrategy.signatureFocus.relCorr.threshold = 0; 
            explanations.searchStrategy.signatureFocus.relCorr.sufficient4fullWeight = 'Relative absolute correlations (times (1-p)^2) mapped to full weight; see calcSignatureFocus. Default: 0.5, i.e. 50% of the maximum already suffices for full weight'; 
              defaultInfo.searchStrategy.signatureFocus.relCorr.sufficient4fullWeight = 0.5;
            explanations.searchStrategy.signatureFocus.relCorr.exponent = 'Exponent in the refined weights formula; see calcSignatureFocus. Default: 1, i.e. a simple linear ramp.'; 
              defaultInfo.searchStrategy.signatureFocus.relCorr.exponent = 1;
            explanations.searchStrategy.signatureFocus.relCorr.cutNoiseBelowQuantileAxisRatio = 'To exclude any influence of (potentially very many) noise gene (repsectively samples), we set all weights below the cut point of the ascendingly sorted weights curve with the linear quantile ramp times this .cutNoiseBelowQuantileAxisRatio paramter to zero exactly. (Opposed to a hard cutoff, this cut point is adaptive to the actual signature size, as determined by the distribution of all computed weights.) Default: 0.67'; 
              defaultInfo.searchStrategy.signatureFocus.relCorr.cutNoiseBelowQuantileAxisRatio = 0.67;
            explanations.searchStrategy.signatureFocus.bWarnOnZeroWeights = 'Development flag only. When adjusting signature focus paramters, it might be useful to warn if >5% of all genes (or samples) get all-zero weight vectors (indicative for too tight focusing). Default: false';
              defaultInfo.searchStrategy.signatureFocus.bWarnOnZeroWeights = false;
            explanations.searchStrategy.correlation.full2focusWeightRatio = 'To stabilize correlation computation in very low-dimensional scenarios, add weights to dimensions that would otherwise have zero weight. (If we only have three genes, do not exclude one or two just because a signature mainly extends along the remaining third basis gene, as the two other genes still carry useful information opposed to computing correlations over only one or two points. Signature focusing should mainly occur via the much higher-dimensional twin/sample space in these scenarios.) The .full2focusWeightRatio defines the ratio of sum(abs(weights)) that is distributed over all dimensions; see addFlatWeights for details. Default: 0.1';
              defaultInfo.searchStrategy.correlation.full2focusWeightRatio = 0.1;
        end
      %Extended signature focus for signature size estimation:
        if(true)
          explanations.searchStrategy.extendedFocus.relCorr.bEnabled = 'Controls whether relative or absolute correlations are used in the (refined) weights formula; see calcSignatureFocus. Default: false, i.e. do not compute correlations relative to max(abs(Rs))'; 
            defaultInfo.searchStrategy.extendedFocus.relCorr.bEnabled = false;
          explanations.searchStrategy.extendedFocus.relCorr.threshold = 'Threshold paramter in the extended (refined) weights formula; see calcSignatureFocus. Default: 0, i.e. no threshold'; 
            defaultInfo.searchStrategy.extendedFocus.relCorr.threshold = 0; 
          explanations.searchStrategy.extendedFocus.relCorr.sufficient4fullWeight = 'Relative absolute correlations (times (1-p)^2) mapped to full extended weight; see calcSignatureFocus. Default: 1, i.e. only r==1 with p->0 is mapped to full extended weight'; 
            defaultInfo.searchStrategy.extendedFocus.relCorr.sufficient4fullWeight = 1;
          explanations.searchStrategy.extendedFocus.relCorr.exponent = 'Exponent in the extended (refined) weights formula; see calcSignatureFocus. Default: 1, i.e. a simple linear ramp.'; 
            defaultInfo.searchStrategy.extendedFocus.relCorr.exponent = 1;
          explanations.searchStrategy.extendedFocus.relCorr.cutNoiseBelowQuantileAxisRatio = 'To exclude any influence of (potentially very many) noise gene (repsectively samples), we set all weights below the cut point of the ascendingly sorted weights curve with the linear quantile ramp times this .cutNoiseBelowQuantileAxisRatio paramter to zero exactly. (Opposed to a hard cutoff, this cut point is adaptive to the actual signature size, as determined by the distribution of all computed weights.) Default: 0.4'; 
            defaultInfo.searchStrategy.extendedFocus.relCorr.cutNoiseBelowQuantileAxisRatio = 0.4;
          explanations.searchStrategy.extendedFocus.bWarnOnZeroWeights = 'Development flag only. When adjusting signature focus paramters, it might be useful to warn if >5% of all genes (or samples) get all-zero weight vectors (indicative for too tight focusing). Default: false';
            defaultInfo.searchStrategy.extendedFocus.bWarnOnZeroWeights = false;
        end
      %Optionally dampen outliers (to become robust against their influence on statistics like signature size estiamtions):
        if(true)
          explanations.searchStrategy.dampenOutliers.bEnabled = 'Many statistics like standard deviations or correlations are not robust against far outliers. Therefore it is recommended to dampen them. This flag enables outlier dampening for signature axes and correlation vectors of signature candidates. Default: true (If .preprocessing.dampenOutliers4? is enabled, this can be disabled to speed up processing.)';
            defaultInfo.searchStrategy.dampenOutliers.bEnabled = true;
          explanations.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio = 'Maximum allowed ratio relative to the next gene or sample in the sort order by the respective axis. If the ratio between neighbours is larger, the neighbour that is nearer to the sort border is set to the maximum (respectively minimum) value allowed by this maximum ratio parameter. Default: 1.03';
            defaultInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio = 1.03;
          if(exist('nG','var') && exist('nP','var'))
            explanations.searchStrategy.dampenOutliers.nTopRanksToCheck4G = 'Number of top and bottom genes with respect to the sort order of the respective sample column that are checked for being an outlier. Default: min(4,floor(nG/20))';
              defaultInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4G = min(4,floor(nG/20)); %Note: must not be too large to prevent smoothing out true narrow and sharp signatures.
            explanations.searchStrategy.dampenOutliers.nTopRanksToCheck4P = 'Number of top and bottom genes with respect to the sort order of the respective sample column that are checked for being an outlier. Default: min(4,floor(nP/20))';
              defaultInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P = min(4,floor(nP/20)); %Note: must not be too large to prevent smoothing out true narrow and sharp signatures.
          end
        end
    %Performance enhancements:
      if(true)
        explanations.searchStrategy.performance.bSkipGenesThatDidNotQualifyInPreviousRuns = 'If a gene that was visited during a lookahead phase in a previous detection iteration and did not qualify wrt. .searchStrategy.qualification.*, it will be skipped. This prevents recalculating signatures for genes (as initial representatives) having high greedy score (for example due to strong folding) but few or none correlated partners. (Since dissections may uncover previously outshined correlations, these genes are testes again irrespective of this flag, once a full pass of all genes and samples has been completed with no qualifying candidates.) Default: true';
          defaultInfo.searchStrategy.performance.bSkipGenesThatDidNotQualifyInPreviousRuns = true;
        explanations.searchStrategy.performance.bSkipSamplesThatDidNotQualifyInPreviousRuns = 'If a sample that was visited during a lookahead phase in a previous detection iteration and did not qualify wrt. .searchStrategy.qualification.*, it will be skipped. This prevents recalculating signatures samples (as initial representatives) having high greedy score (for example due to strong folding) but few or none correlated partners. (Since dissections may uncover previously outshined correlations, these samples are tested agin irrespective of this flag, once a full pass of all genes and samples has been completed with no qualifying candidates.) Default: true';
          defaultInfo.searchStrategy.performance.bSkipSamplesThatDidNotQualifyInPreviousRuns = true;

        explanations.searchStrategy.performance.bCachePotentialInitialRepresentatives = 'Cache initial representatives from step 1 that were not affected by dissections of the detected signatures(s). Use this for large datasets to speedup computation. Default: true if nG>20000 && nP>=200.';
          if(exist('nG','var') && exist('nP','var')) %data dependent defaults require the signal size
            defaultInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives = nG>20000 && nP>=200;
          end
        explanations.searchStrategy.performance.maxAllowedDissectionStrengthToKeepInCache = 'If .bCachePotentialInitialRepresentatives, this defiend tha maximum allowed dissection strengths for an initial representative gene or sample to remain in the step 1 cache. If it was affected with a higher dissection strengths, it will be recomputed. Default: 0.05.';
          defaultInfo.searchStrategy.performance.maxAllowedDissectionStrengthToKeepInCache = 0.05;

        explanations.searchStrategy.performance.alpha4Signal4inclusionInComputation4G = 'When determining over which genes the correlation of samples with sample axes needs to be computed for a precise result, genes with p_Signal<.alpha4Signal4inclusionInComputation4G are always included; others may be excluded for performance reasons. (See selectDimensionsInPerformanceSubspace. for details.) Useful to speed up correlation computation in large spaces with small signatures. Default: 0.1';
          defaultInfo.searchStrategy.performance.alpha4Signal4inclusionInComputation4G = 0.1;
        explanations.searchStrategy.performance.alpha4Signal4inclusionInComputation4P = 'When determining over which samples the correlation of genes with gene axes needs to be computed for a precise result, samples with p_Signal<.alpha4Signal4inclusionInComputation4P are always included; others may be excluded for performance reasons. (See selectDimensionsInPerformanceSubspace. for details.) Useful to speed up correlation computation in large spaces with small signatures. Default: 0.1';
          defaultInfo.searchStrategy.performance.alpha4Signal4inclusionInComputation4P = 0.1;
        explanations.searchStrategy.performance.minAbsR4inclusionInCorrComputation = 'When determining over which genes/samples the correlation of samples/genes with sample/gene axes needs to be computed for a precise result, samples with |r| >.minAbsR4inclusionInCorrComputation are always included; others may be excluded for performance reasons. (See selectDimensionsInPerformanceSubspace. for details.) Useful to speed up correlation computation in large spaces with small signatures. Default: 0.8';
          defaultInfo.searchStrategy.performance.minAbsR4inclusionInCorrComputation = 0.8;
        explanations.searchStrategy.performance.minPerformanceSpaceSize4G = 'Never reduce the number of genes over which sample correlations are computed below this threshold. (Set this to >nG to always include all genes of the signal when computing correlations; not recommended as this is unnecessarily slow for large datasets.) Default: 1000';
          defaultInfo.searchStrategy.performance.minPerformanceSpaceSize4G = 1000;
        explanations.searchStrategy.performance.minPerformanceSpaceSize4P = 'Never reduce the number of samples over which gene correlations are computed below this threshold. (Set this to >nP to always include all samples of the signal when computing correlations; not recommended as this is unnecessarily slow for large datasets.) Default: 1000';
          defaultInfo.searchStrategy.performance.minPerformanceSpaceSize4P = 1000;
      end
  %% Configuration of STEP 2) Correlation maximization for signature generalization and definition (climb/orient axes to the center of the correlation-defined pattern in order to get independent of the initial representative)
    %Convergence criteria for the correlation-maximized and generalized signature:
      explanations.correlationMaximization.convergenceEpsilon4deltaCorr4AxesInNoiseSDs = 'Difference from one for the correlation (of signature axes to before the incorporation of the last representative into the signature) that is considered sufficiently small for axes convergence. Specified in units of the noise SD. Default: 5e-3';
        defaultInfo.correlationMaximization.convergenceEpsilon4deltaCorr4AxesInNoiseSDs = 5e-3;
      %TODO: allow absolute specification of explanations.correlationMaximization.convergenceEpsilon4deltaCorr4Axes = 'Difference from one for the correlation (of signature axes to before the incorporation of the last representative into the signature) that is considered sufficiently small for axes convergence. Default: 5e-3';
      %  defaultInfo.correlationMaximization.convergenceEpsilon4deltaCorr4Axes = 5e-3;
      explanations.correlationMaximization.convergenceEpsilon4deltaCorr4Rs = 'Difference from one for the correlation (of signature correlations to before the incorporation of the last representative into the signature) that is considered sufficiently small for convergence. Default: 1e-4';
        defaultInfo.correlationMaximization.convergenceEpsilon4deltaCorr4Rs = 1e-4;
    %Generalization criteria: Check for sufficiently many representatives for a generalized signature definition (that is independent of characteristics from individual representatives):
      explanations.correlationMaximization.nSufficient4Generalization = 'If more genes and samples have been selected as representatives for the current signature than .nSufficient4Generalization it is considered sufficiently generalized (i.e. sufficiently independent of features from individual representatives). Default: 15';
        defaultInfo.correlationMaximization.nSufficient4Generalization = 15; 
      explanations.correlationMaximization.minRatioOfSignatureSize4Generalization = 'If more genes and samples have been selected as representatives for the current signature than .minRatioOfSignatureSize4Generalization times the current signature size (correlation sums of genes plus samples), it is considered sufficiently generalized (i.e. sufficiently independent of features from individual representatives). Default: 0.2';
        defaultInfo.correlationMaximization.minRatioOfSignatureSize4Generalization = 0.2;
      explanations.correlationMaximization.minRatioOfSignatureSizeInSingleDim4Generalization = 'If more genes repsectively samples have been selected as representatives for the current signature than .minRatioOfSignatureSizeInSingleDim4Generalization times the current signature size in the respective order dimension (correlation sum of genes respectively samples), it is considered sufficiently generalized (i.e. sufficiently independent of features from individual representatives). Default: 0.8';
        defaultInfo.correlationMaximization.minRatioOfSignatureSizeInSingleDim4Generalization = 0.8;
    %Lookahead for next signature member towards optimal correlation/consistency score:
      explanations.correlationMaximization.maxLookahead4NextBestAccuCandidate = 'Minimum number of candidate representatives (genes or samples) to be tested from the presort order by correlation to signature axes for the optimal step towards maximal signature functional. Default: 20'; 
        defaultInfo.correlationMaximization.maxLookahead4NextBestAccuCandidate = 20; 
      explanations.correlationMaximization.maxLookahead4NextBestAccuCandidateFactorInCaseOfZeroCandidates = 'If no signature functional maximizing candidate representative has been found within the .maxLookahead4NextBestAccuCandidate, but the signature still requires additional representatives for sufficient generalization, temporarily increase the lookahead by this factor. Default: 2.5';
        defaultInfo.correlationMaximization.maxLookahead4NextBestAccuCandidateFactorInCaseOfZeroCandidates = 2.5;
    %Performance enhancements:
      explanations.correlationMaximization.performance.nMinPrecomputeSignaturesPerWorker = 'Configures the minimum workload per Matlab worker to minimize the paralellization overhead percentage. Default: 4';
        defaultInfo.correlationMaximization.performance.nMinPrecomputeSignaturesPerWorker = 4;
      explanations.correlationMaximization.performance.scoreBoostFactore4alreadyPrecomputedMemberCandidates = 'When selecting additional members towards signature functional maximization, this factor is applied to slightly increase scores for already precomputed genes/samples in order to select them instead of nearly score-equal alternatives without cache entry. Default: 1.15';
        defaultInfo.correlationMaximization.performance.scoreBoostFactore4alreadyPrecomputedMemberCandidates = 1.15;
      %Allow early breaks of lookahead iterations, if we already found a good candidate for the next representative:
        explanations.correlationMaximization.earlyBreakLookahead.fcnMinScoreRatio2QualifyForEarlyBreak4CurrentLookahead = 'We allow to break the search for the next representative early (before testing all .maxLookahead4NextBestAccuCandidate candidates), if we find a candidate in presort order by descending correlations that leads to a signature functional score increase >= this ratio. Default: 1.01-0.02*nRepresentativesSoFar/nSufficient4Generalization, i.e. require score increases for all representatives up tp 50%*nSufficient4Generalization; after that we allow early-breaks even for candidates that keep the score constant or are only slightly score-decreasing (their addition still needs to maintain >=99% of the score).';
          lc_nSufficient4Generalization = defaultInfo.correlationMaximization.nSufficient4Generalization;
          defaultInfo.correlationMaximization.earlyBreakLookahead.fcnMinScoreRatio2QualifyForEarlyBreak4CurrentLookahead = ...
             @(nCurrentLookaheadRange,nRepresentativesSoFar) 1.01 - 0.02*nRepresentativesSoFar/lc_nSufficient4Generalization ... %in the beginning, require score increases, later during generalization 99% of the previous score is fine for perfomance.
          ;
        explanations.correlationMaximization.earlyBreakLookahead.minLookahead4MemberCount = 'We allow to break the search for the next representative early (before testing all .maxLookahead4NextBestAccuCandidate candidates), but we require at least .minLookahead4MemberCount candidates to be tested dependent on the number of so far selected representatives. Default: 2.^(4:-1:0), i.e. at the beginning we do not early break after the first score increase, but lookahead for a few candidates in presort order by descending correlation to signature axes to find the optimal increase of the signature functional';
          defaultInfo.correlationMaximization.earlyBreakLookahead.minLookahead4MemberCount = 2.^(4:-1:0); 
  %% Application mode) Replacement for STEP 1 and STEP 2 in the signature application phase (i.e. when signature axes are alreadig known and should just be applied for bimonotonic regression to, e.g., a new sample cohort):
    explanations.applicationMode.bEnabled = 'If true, the search and correlation maximization strategy for detection (steps 1 and 2) is replaced by externally provided signatures via SDCM''s third input parameter previouslyDetectedSignaturesForApplicationMode. Default: false';
      defaultInfo.applicationMode.bEnabled = false;
    explanations.applicationMode.externalInitialShifts = 'Optional nG*nP matrix of log2(ratio)s that are treated as "already explained signal parts" (they will be subtracted from the original signal and be visible in panels b of signature plots). Default: 0';
      defaultInfo.applicationMode.externalInitialShifts = 0;
    explanations.applicationMode.bApplyEveryExternalSignatureToInitialSignal = 'If true, every signature in previouslyDetectedSignaturesForApplicationMode is applied to the initial signal. This has the benefit of validation results for applied signatures being independent of the application order. If false, they are applied and dissected sequentially in the externally provided order. For validation purposes this flag is recommended (even in case of several superposed signatures, correlations with provided provided signature axes should be significant, if the signature also exists in the validation input dataset). Default: true';
      defaultInfo.applicationMode.bApplyEveryExternalSignatureToInitialSignal = true;
    explanations.applicationMode.kOffset4filenames = 'The k of the first processed externally provided signature in previouslyDetectedSignaturesForApplicationMode is 1 plus this offset. (Useful for controlling exported file names when passing single signatures as previouslyDetectedSignaturesForApplicationMode interactively). Default: 0';
      defaultInfo.applicationMode.kOffset4filenames = 0;
    explanations.applicationMode.maxk_externallyValidatedSignatures = 'The total number of external signatures being applied (used in plot titles). Default: length(previouslyDetectedSignaturesForApplicationMode)';
      defaultInfo.applicationMode.maxk_externallyValidatedSignatures = 'DEFAULT';
  %% Configuration of STEP 3) Bimonotonic regression of all genes/arrays in the signature's eigenOrder of genes and samples to determine the signal parts contributed by the current signature (i.e. by its underlying interaction) to the signal sum:
    %Ordering metrics: define the order of genes and samples for the signature:
      explanations.dissection.metric4geneOrder = 'The signature fieldname used to determine the gene order of the signature for regression. Default: signatureStrengths4G, i.e. weighted projections of all genes on the signature''s gene axis';
        defaultInfo.dissection.metric4geneOrder = 'signatureStrengths4G'; %'R4G'; 
      explanations.dissection.metric4sampleOrder = 'The signature fieldname used to determine the sample order of the signature for regression. Default: signatureStrengths4P, i.e. weighted projections of all samples on the signature''s sample axis';
        defaultInfo.dissection.metric4sampleOrder = 'signatureStrengths4P'; %'R4P'; 
    %Bimonotonic regression to make sure no signal of superposed signatures are included in (or partly explained by) the current signature's eigensignal:
      explanations.dissection.bimonotonicRegression.convergenceThreshold4r4TopOfSignatureDissection = 'The minimum correlation between the current and previous bi-monotonically regressed signature eigensignal in the top signature focus required for bimonotonic regression to be considered converged in the (outer) iteration in SDCM_step3_regression. Default: 0.99';
        defaultInfo.dissection.bimonotonicRegression.convergenceThreshold4r4TopOfSignatureDissection = 0.99; 
      explanations.dissection.bimonotonicRegression.maxRegressionIterations = 'The maximum number of (outer) bi-monotonic regressions in SDCM_step3_regression performed to obtain the signature eigensignal. Use Inf for no limitation. Default: 7 (Usually only 2-3 are needed. This is just a failsave parameter; if this limit should ever be reached, a warning is issued and SDCM continues detecting the next signature.)';
        defaultInfo.dissection.bimonotonicRegression.maxRegressionIterations = 7;
      explanations.dissection.bimonotonicRegression.epsilon4convergenceInSDs = 'Maximum allowed remaining standard deviations in the signature focus between current and previous inner iteration in bimonotonicRegression. Specified in units of the estimated noise standard deviation of the signal. Once falling below this threshold, the matrix is numerically considered bi-monotonic. Default: 1e-3';
        defaultInfo.dissection.bimonotonicRegression.epsilon4convergenceInSDs = 1e-3;
      explanations.dissection.bimonotonicRegression.bUpdateSignatureFocusInEachSignatureCurveRegressionIteration = 'If enabled, correlations .R4G and .R4P are refined based on the current signature focus and then the signature focus is refined based on these refined correlations in each outer bimonotonic regression iteration until convergence. Default: true';
        defaultInfo.dissection.bimonotonicRegression.bUpdateSignatureFocusInEachSignatureCurveRegressionIteration = true; %false;
      %Performance enhancements:
        explanations.dissection.bimonotonicRegression.maxResolutionOnSignalDeltaGrid4G = 'Instead of isotonically regressing all genes in the signature focus, we limit regresion to .maxResolutionOnSignalDeltaGrid4G genes that are equidistantly spaced with respect to signatureStrengths4G. Others are interpolated. For large datasets, this may be considerably faster than regressing every single gene. Default: 1000';
          defaultInfo.dissection.bimonotonicRegression.maxResolutionOnSignalDeltaGrid4G = 1000;
        explanations.dissection.bimonotonicRegression.maxResolutionOnSignalDeltaGrid4P = 'Instead of isotonically regressing all samples in the signature focus, we limit regresion to .maxResolutionOnSignalDeltaGrid4P samples that are equidistantly spaced with respect to signatureStrengths4P. For large datasets, this may be considerably faster than regressing every single sample. Default: 1000';
          defaultInfo.dissection.bimonotonicRegression.maxResolutionOnSignalDeltaGrid4P = 1000;
    %Adaptive Smoothening: For (gene,sample) points outside of the signature focus, the signature signal is defined based on smoothening of the signal in signature order using a constant smoothing window in the space of equidistant signature strengths (corresponding to an adaptive window in the original index space).
      %Configure resolution for the rescaled&downscaled signature strengths space for smoothing:
        explanations.dissection.signatureSpace.resolution4G = 'Number of grid cells in gene direction for downscaling and rescaling into the the signature strengths space. (From minimum to maximum of .metric4geneOrder, .resolution4G equidistantly spaced grid cells are used.) Default: 512 (this is sufficient to resolve monotonic non-linearities)';
          defaultInfo.dissection.signatureSpace.resolution4G = 512;
        explanations.dissection.signatureSpace.resolution4P = 'Number of grid cells in sample direction for downscaling and rescaling into the the signature strengths space. (From minimum to maximum of .metric4sampleOrder, .resolution4P equidistantly spaced grid cells are used.) Default: 512 (this is sufficient to resolve monotonic non-linearities)';
          defaultInfo.dissection.signatureSpace.resolution4P = 512;
        explanations.dissection.signatureSpace.trafo4G = 'Optional transformation function of .trafo4G(.dissection.metric4geneOrder) values before rescaling/interpolation on equidistant axes. (Use e.g. @(X)X.^3 to assign more relative resolution to values far from zero.) Default: identity';
          defaultInfo.dissection.signatureSpace.trafo4G = @(X)X;
        explanations.dissection.signatureSpace.trafo4P = 'Optional transformation function of .trafo4P(.dissection.metric4sampleOrder) values before rescaling/interpolation on equidistant axes. (Use e.g. @(X)X.^3 to assign more relative resolution to values far from zero.) Default: identity';
          defaultInfo.dissection.signatureSpace.trafo4P = @(X)X;
      %Smoothing kernel used to equalize neighbouring expressions (with effectively adaptive kernel size in the original index space):
        explanations.dissection.gaussianKernel.sigma4G = 'Determines the width of the Gaussian smoothing kernel for the .trafo4G(.dissection.metric4geneOrder) values. Default: NaN; .eqdResolutionInSigma4G is used isntead to become independent of signal units';
          defaultInfo.dissection.gaussianKernel.sigma4G = NaN;
        explanations.dissection.gaussianKernel.sigma4P = 'Determines the width of the Gaussian smoothing kernel for the .trafo4P(.dissection.metric4sampleOrder) values. Default: NaN; .eqdResolutionInSigma4G is used isntead to become independent of signal units';
          defaultInfo.dissection.gaussianKernel.sigma4P = NaN;
        %Alternatively specify sigmas in terms of the neighbouring pixels in .signatureSpace that it should comprise, i.e. a specification relative to .resolution4?:
          explanations.dissection.gaussianKernel.eqdResolutionInSigma4G = 'Smoothening kernel sigma/width in gene direction in numbers of pixels in the .signatureSpace.resolution4G grid. Default: 8';
            defaultInfo.dissection.gaussianKernel.eqdResolutionInSigma4G = 8;
          explanations.dissection.gaussianKernel.eqdResolutionInSigma4P = 'Smoothening kernel sigma/width in sample direction in numbers of pixels in the .signatureSpace.resolution4P grid. Default: 8';
            defaultInfo.dissection.gaussianKernel.eqdResolutionInSigma4P = 8;
  %% Plot configuration:
    explanations.plots.bVisible = 'Controls if generated figures are visible. Default: ~inInfo.export.plots.bEnabled';
      defaultInfo.plots.bVisible = 'DEFAULT';
    explanations.plots.signatureDefinitions.bEnabled = 'Enables or disables the main definition plot of a detected signature (showing original expressions in signature order + already explained signal parts = current signal giving rise to the detected signature = bimonotonically regressed signal + remaining signal after dissection). Note that the plot code is provided AS IS and may not work on your local platform as intended. You are free to adapt the plot code to your needs, but this is not supported. In case of plotting errors, SDCM shows a warning, but continues with detection. Default: true';
      defaultInfo.plots.signatureDefinitions.bEnabled = true;
      explanations.plots.signatureDefinitions.nInitialSignaturesNotToPlot = 'Controls how many signatures in detection order k are skipped when plotting the already explained signal in panel b. (Sometimes e.g. large superposed lab signature are so uninteresting that we do not want to see them in the "already explained signal" panel of definition plots.) Default: 0, i.e. plot all previously detected signatures';
        defaultInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot = 0;
      %<-Note: many more plot options are available and defined in plotDetectedSignature.m
    explanations.plots.focusConvergence.bEnabled = 'Enables or disables the plot illustrating selected representatives and the convergence towards signature functional maximization in step 2. Default: false'; 
      defaultInfo.plots.focusConvergence.bEnabled = false;
    explanations.plots.postprocessingComparisons.bEnabled = 'Enables or disables the plot visualizing correlations between all detected signature axes. Allows studying partial correlations between detected signatures. Default: false';
      defaultInfo.plots.postprocessingComparisons.bEnabled = false;
    %Colormap used by heatmap plots:
      if(true)
        colorResolution=252;
        ramp = ((1:colorResolution/2)/(colorResolution/2))';  
        CM = flipud([[flipud(ramp) zeros(size(ramp)) zeros(size(ramp))]; [0 0 0]; [zeros(size(ramp)) zeros(size(ramp)), ramp]]); 
        CM = bsxfun(@times,sum(CM,2).^(1/3),(1-bsxfun(@times,1-sum(CM,2)/3,exp(-3*CM))).^1); %increase dynamic color range a bit.
        CM = CM.^0.67;
        defaultInfo.plots.colormap = CM;
        defaultInfo.plots.colormap_withNaNColor = [[.5 .5 .5];defaultInfo.plots.colormap]; %use gray as NaN color.
      end
  %% Postprocessing options:
    %cutoffs for signature statistics and plots:
      explanations.postprocessing.cutoffs4statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot = 'Minimum dissection strength anywhere in a gene row or sample column to be considered as influenced by the detected signature. This will determine the signatureDefinition.statistics.signatureMembership.* counts. This cutoff may also determine genes shown in the default plot (if .plots.signatureDefinitions.selection.bExcludeIfBelowMinRequireddDissectionStrength). Default: 0.05';
        defaultInfo.postprocessing.cutoffs4statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot = 0.05;

      explanations.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G = 'Minimum ratio relative to the maximum absolute correlation to consider a gene being part of the "top" of the signature. This will determine the signatureDefinition.statistics.relCorrGtTopThreshold.* counts like .nTopCoregGenes. This cutoff also determines the gene lists generated by .export.infoFile4topMembers.*. Default: 0.67';
        defaultInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G = 0.67;
      explanations.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P = 'Minimum ratio relative to the maximum absolute correlation to consider a sample being part of the "top" of the signature. This will determine the signatureDefinition.statistics.relCorrGtTopThreshold.* counts like .nTopAffectedSamples. This cutoff also determines the gene lists generated by .export.infoFile4topMembers.*. Default: 0.67';
        defaultInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P = 0.67;

      explanations.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4G = 'Minimum ratio relative to the maximum absolute correlation to consider a gene being part of the "outreach" of the signature. This will determine the signatureDefinition.statistics.relCorrGtNoiseThreshold.* counts like .nRepresentativeCoregGenes. This cutoff may also determine genes shown in the default plot (if .plots.signatureDefinitions.selection.bExcludeIfBelowRelCorrNoiseThresholds). Default: 0.5';
        defaultInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4G = 0.5;
      explanations.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4P = 'Minimum ratio relative to the maximum absolute correlation to consider a sample being part of the "outreach" of the signature. This will determine the signatureDefinition.statistics.relCorrGtNoiseThreshold.* counts like .nRepresentativeCoregSamples. This cutoff may also determine samples shown in the default plot (if .plots.signatureDefinitions.selection.bExcludeIfBelowRelCorrNoiseThresholds). Default: 0.5';
        defaultInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4P = 0.5;

    %alternative p-value based cutoffs for defining "top" members of a signature:
      defaultInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G = 0.001;
        explanations.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G = 'Maximum p value for the correlation of a gene with the gene axis to count it as top correlated / bottom anti-correlated gene of the signature. Default: 1e-3';
      defaultInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P = 0.001;
        explanations.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P = 'Maximum p value for the correlation of a sample with the samples axis to count it as top correlated / bottom anti-correlated sample of the signature. Default: 1e-3';
      defaultInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4G = 0.2;
        explanations.postprocessing.cutoffs4statistics.noiseThresholds.alpha4G = 'Maximum p value for the correlation of a gene with the gene axis to not count it as noise gene with respect to the signature. Note that the absolue correlation strength should be considered too, especially in large datasets where only moderate correlations may already be significant due to the amount of supporting points. Default: 0.2';
      defaultInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P = 0.2;
        explanations.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P = 'Maximum p value for the correlation of a sample with the samples axis to not count it as noise sample with respect to the signature. Note that the absolue correlation strength should be considered too, especially in large datasets where only moderate correlations may already be significant due to the amount of supporting points. Default: 0.2';
  %% Export configuration:
    explanations.export.rootDir = 'The export root directory wherein all files and subfolders for detected signatures are created. Default: none, i.e. a required user parameter (if none is specified, any enabled exports are disabled)';
      defaultInfo.export.rootDir = 'REQUIRED';
    %Status output: Detail level of status messages exported/printed to the console:
      explanations.export.nStatusOutputLevel = 'Configures the detail level of status messages printed to the console during detection. Default: 2 (set to inf to see all detail messages)';
        defaultInfo.export.nStatusOutputLevel = 2;
    %Reproducability:
      %Logging:
        explanations.export.bCreateDetectionLogFile = 'If this flag is set, the console output of SDCM will be saved in .rootDir\detection.log. Warning: If a detection.log file already exists in the configured export directory, it is overwritten. Default: true';
          defaultInfo.export.bCreateDetectionLogFile = true;
      %Configuration and input data:
        explanations.export.matFile4configurationAndInputData.bEnabled = 'If enabled, the used detection configuration and input data will be exported to disk after initialization for reference and reproducability. Default: true';
          defaultInfo.export.matFile4configurationAndInputData.bEnabled = true;
      %SDCM copy for reference:
        explanations.export.copyCurrentAlgorithm.bEnabled = 'If enabled, a SDCM is copied to the export dir for reference and reproducability. Default: false';
          defaultInfo.export.copyCurrentAlgorithm.bEnabled = false;
    %Signature definitions as Matlab structure:
      explanations.export.matFile4eachSignatureDefinition.bEnabled = 'If enabled, each detected signature will be exported to disk at the end of the each detection iteration as Matlab structure. Default: true';
        defaultInfo.export.matFile4eachSignatureDefinition.bEnabled = true; 
      %Optional additional information to be saved in MAT files for each detected signature:
        explanations.export.matFile4eachSignatureDefinition.bIncludeDefPlotPanels = 'If enabled, MAT-exported signature definitions and SDCM''s outInfo structure will contain .forPlots.overview.* fields with the results for each panel in the signature definition plot. Disable for speed and lower file size. Default: true';
          defaultInfo.export.matFile4eachSignatureDefinition.bIncludeDefPlotPanels = true;
        explanations.export.matFile4eachSignatureDefinition.bIncludeRegressedSignalEstimationSteps = 'If enabled, MAT-exported signature definitions and SDCM''s outInfo structure will contain .forPlots.dissectionPipelines.* fields with intermediate results in the equidistant signature space. Warning: consumes lots of memory for large datasets. Default: false';
          defaultInfo.export.matFile4eachSignatureDefinition.bIncludeRegressedSignalEstimationSteps = false;
        explanations.export.matFile4eachSignatureDefinition.bIncludeBimonotonicConvergenceSteps = 'If enabled, MAT-exported signature definitions and SDCM''s outInfo structure will contain .forPlots.bimonotonicConvergence.* fields with intermediate results for the convergence iterations towards bi-monotonicity. Warning: consumes lots of memory for large datasets. Default: false';
          defaultInfo.export.matFile4eachSignatureDefinition.bIncludeBimonotonicConvergenceSteps = false;
      %Flags for disabling standard exports:
        explanations.export.matFile4eachSignatureDefinition.bDisableNoiseDistributionFileExport = 'If enabled, SDCM does not export "*noise distribution.mat" files. Even if not exported, this information is still removed from RAM (if .bRemoveNoiseEstimationResumeDataAfterExport is set) and a rerun is necessary to get these information again. No resume is available if this option is enabled. Use this for runs when these files and the resume feature are not needed. Default: false';
          defaultInfo.export.matFile4eachSignatureDefinition.bDisableNoiseDistributionFileExport = false;
        explanations.export.matFile4eachSignatureDefinition.bDisablePipelineStateForPlotsFileExport = 'If enabled, SDCM does not export "*pipeline state for plots.mat" files. Even if not exported, this information is still removed from RAM (if .bRemoveNoncoreInfoFromRAMandOutputAfterExport is set) and a rerun is necessary to get these information again. Use this when these files or post-detection plots are not of interest. Default: false';
          defaultInfo.export.matFile4eachSignatureDefinition.bDisablePipelineStateForPlotsFileExport = false;
      %Performance enhancements / RAM savers:
        explanations.export.matFile4eachSignatureDefinition.bRemoveNoiseEstimationResumeDataAfterExport = 'If enabled, saves RAM by removing information that is not crucial for dissection to continue. Resume data must/can be loaded from disk afterwards, if needed (the .forResume.noiseEstimation field points to the correct exported file, if ~.matFile4eachSignatureDefinition.bDisableNoiseDistributionFileExport). Default: true';
          defaultInfo.export.matFile4eachSignatureDefinition.bRemoveNoiseEstimationResumeDataAfterExport = true;
        explanations.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport = 'If enabled, saves RAM by removing information that is not crucial for dissection. Plotting info for individual signatures will not be present in SDCM''s outInfo structure and will have to be loaded from disk (assuming .export.matFile4eachSignatureDefinition.bEnabled is enabled, otherwise these additional information are lost). Default: true';
          defaultInfo.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport = true;
    %Optional text files with simple signature stats (for searching in the file system):
      explanations.export.infoFile4initialSignatureStats = 'Creates a file with the core stats of a detected signature''s initial representative in the export dir (for manual control and searching by file name). Default: true';
        defaultInfo.export.infoFile4initialSignatureStats.bEnabled = true;
      explanations.export.infoFile4focusedSignatureStats = 'Creates a file with the core stats of a focused signature after correlation maximization in the export dir (for manual control and searching by file name). Default: true';
        defaultInfo.export.infoFile4focusedSignatureStats.bEnabled = true;
      explanations.export.infoFile4topMembers.bEnabled4topCorrGenes = 'Creates a file listing the top correlated genes (as defined in .reference.rowIDs) of a detected signature in the export dir (as determined by .postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G). Useful for searching for particular gene names with text find tools (e.g. findstr in Windows). Default: true';
        defaultInfo.export.infoFile4topMembers.bEnabled4topCorrGenes = true;
      explanations.export.infoFile4topMembers.bEnabled4topAntiCorrGenes = 'Creates a file listing the top anti-correlated genes (as defined in .reference.rowIDs) of a detected signature in the export dir (as determined by .postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G). Useful for searching for particular gene names with text find tools (e.g. findstr in Windows). Default: true';
        defaultInfo.export.infoFile4topMembers.bEnabled4topAntiCorrGenes = true;
      explanations.export.infoFile4topMembers.bEnabled4topCorrSamples = 'Creates a file listing the top correlated sample IDs (as defined in .reference.colIDs) of a detected signature in the export dir (as determined by .postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P). Useful for searching for particular sample names with text find tools (e.g. findstr in Windows). Default: true';
        defaultInfo.export.infoFile4topMembers.bEnabled4topCorrSamples = true;
      explanations.export.infoFile4topMembers.bEnabled4topAntiCorrSamples = 'Creates a file listing the top anti-correlated sample IDs (as defined in .reference.colIDs) of a detected signature in the export dir (as determined by .postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P). Useful for searching for particular sample names with text find tools (e.g. findstr in Windows). Default: true';
        defaultInfo.export.infoFile4topMembers.bEnabled4topAntiCorrSamples = true;
    %Signature definitions as tables:
      explanations.export.table4eachSignatureDefinition.bEnabled = 'If enabled, a table for each detected signature with the defining signature axes, signature focus, correlations and strengths is exported to disk at the end of the each detection iteration. Default: true';
        defaultInfo.export.table4eachSignatureDefinition.bEnabled = true;
    %Signature definitions as figures and other figure exports:
      explanations.export.plots.bEnabled = 'If enabled, the signature definition plot for the detected signature will be exported to disk at the end of the each detection iteration. Default: true';
        defaultInfo.export.plots.bEnabled = true;
      explanations.export.plots.exportFormats = 'Specifies the picture formats used for figure exports as cellstring of file extensions. Default: {''eps''}';
        defaultInfo.export.plots.exportFormats = {'eps'};
      explanations.export.plots.bCloseAfterExport = 'If .export.plots.bEnabled and this flag is set, generated figures will be invisible during computation and will be closed after export. Default: true, except if inInfo.plots.bVisible==true';
        defaultInfo.export.plots.bCloseAfterExport = 'DEFAULT';
  %% Internal flags used during algorithm development:
    explanations.internal.bResume = '(for development only) If the global eDState variable already contains detected signatures of the current L2Is (by a previous run or loaded from disk), this flag saves times by reusing them. Default: false';
      defaultInfo.internal.bResume = false;
    explanations.internal.bResumeFromDisk = '(for development only) Enables loading previously detected and exported signatures from the configured export directory. Replaces global eDState that can subsequently be used by .bResume. Default: false';
      defaultInfo.internal.bResumeFromDisk = false;
    explanations.internal.bReuseCacheOnEndResume = '(for development only) If global eDState.current.precompute_previous contains the current/valid cache for the last/next signature, this option reuses it instead of clearing it. Default: false';
      defaultInfo.internal.bReuseCacheOnEndResume = false;
    explanations.internal.bDevEnableInteractiveBreaks = '(for development only) If enabled, MATLAB''s interactive mode is entered after the first detection iteration to allow inspection of variables (via the >>keyboard command). Continue execution via >>dbcont. A small figure window will be opened. After subsequent detection iterations the interactive mode will only be re-entered if this figure is closed during detection. Default: false';
      defaultInfo.internal.bDevEnableInteractiveBreaks = false;
    explanations.export.nStatusOutputLevel4checkpoints = '(for development only) Execution interrupts at all status messages <= this level, if .bDevEnableInteractiveBreaks and the trigger figure has been closed. Default: 2';
      defaultInfo.export.nStatusOutputLevel4checkpoints = Inf;
    explanations.internal.bDevPlots = '(for development only) Several smaller quick and dirty visualizations of intermediate results are enabled by this flag (only useful during debugging). Default: false';
      defaultInfo.internal.bDevPlots = false;
    explanations.internal.bConsistencyCheckAfterEachDetectionIteration = '(for development only) Numerically assert after each detection iteration that the invariant sum(explained signals) + remaining signal = input signal. Default: false';
      defaultInfo.internal.bConsistencyCheckAfterEachDetectionIteration = false; %true;
    explanations.internal.bEnablePostProductionCodeOptimizations = 'Some code optimizations were implemented after production. For backwards compatibility, they can be centrally disabled with this flag. Default: true';
      defaultInfo.internal.bEnablePostProductionCodeOptimizations = true;
    explanations.internal.bAllowDissectionOfTinyDatasets = 'Originally, numerics of SDCM were developed for high-dimensional datasets. SDCM''s concepts extends mathematically to low-dimensional scenarios as well. However, application of SDCM to low-dimensional data has not been sufficiently tested yet and, therefore, is still considered experimental. To override this (e.g. for generating 3D concept plots or for algorithm development) set this flag. Default: false, i.e. only allow sufficently large input datasets.';
      defaultInfo.internal.bAllowDissectionOfTinyDatasets = false;
  %Output: add field explanations to the defaultInfo config structure:
    defaultInfo.explanations = explanations;
end

