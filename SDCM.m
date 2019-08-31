%ABSTRACT
%  Takes a 2D matrix as input signal (e.g. a log2-transformed gene expression dataset 
%  with all measured genes as matrix rows and all measured samples as columns). 
%  Dissects the signal into signatures as defined in [1]. Each signature has a 
%  bi-monotonic signal for its gene and sample order in a weighted subset of 
%  genes and samples. For each detected signature, the regressed bi-monotonic signal
%  is saved in the output and dissected from the signal. Signatures are detected
%  iteratively until statistically only noise seems to be left, as per configured 
%  significance thresholds.
%  Signatures are comprised of correlated genes and samples. Therefore, these signatures 
%  are usually biologically specific and meaningful. All detected signatures 
%  and the residual signal are available in the output for further analyses (for 
%  example, for a cross-cohort validation or for pipelining into statistical tests 
%  like gene set enrichment analyses for signature interpretation).
%SYNTAX
%  outInfo = SDCM(L2Is, inInfo)
%  outInfo = SDCM(L2Is, inInfo, previouslyDetectedSignatures)
%INPUT
%  L2Is: nG*nP numeric matrix (may be single precision to save RAM). For example, L2Is 
%    may be log2(gene expressions) of nG measured genes (or probesets) for nP samples.
%  inInfo: Configuration structure:
%  - The default parameter set is returned by:
%      >> defaultConfig = SDCM_defaultConfig();
%  - All available parameters like defaultConfig.x.y are explained by:
%      >> defaultConfig.explanations.x.y
%  - Configuration categories are:
%      .reference: Declares row and column IDs and their external weights, if any.
%      .preprocessing: Configures data preprocessing steps before dissection.
%      .searchStrategy: Detection parameters for the search strategy that looks for 
%         initial signature representatives.
%      .correlationMaximization: Configures correlation maximization phase.
%      .dissection: Parameters for regressing the bi-monotonic signal of a 
%         signature and its dissection from the signal.
%      .applicationMode: Parameters for the application phase, for which the input 
%         parameter previouslyDetectedSignatures must be provided 
%         (see below). This replaces the search strategy by provided previously learnt
%         signature axes. Useful, e.g., to check for signatures in validation datasets.
%      .plots: Configures plots of detected signatures.
%      .export: Sets if, what and where to export results to disk.
%  previouslyDetectedSignatures: Optional array of signature axes in the structure of
%    outInfo.signatures.step2_finalSignatureAxes from previous SDCM detection runs. 
%    Used, e.g., for validation purposes; see .applicationMode above.
%OUTPUT
%  outInfo: 
%  - All available output parameters like outInfo.x.y are explained by:
%      >> outInfo.explanations.x.y
%  - In .signatures all detected signatures are returned. Each one contains in particular:
%     .step2_finalSignatureAxes: signature axes and correlations with them
%     .step3_regression: the signature's contributions to the overall signal sum.
%DEPENDENCIES
%  Requires Matlab's Statistics Toolbox, Bioinformatics Toolbox and Matlab's Parallel 
%  Computing Toolbox for parallelization. (Additional required tools are provided in
%  the SDCM_Library folder for standalone deployment and execution.)
%HINTS
%  The source code is structured like a document where comments are headers for the
%  following code lines. The document's hierarchy is realized by indentation and code 
%  folding. When working with the source code, it is recommended to configure Matlab's 
%  code folding to also include if-structures, as I often use them to semantically 
%  group lines of code.
%REFERENCES
%  For a complete documentation of signal dissection, see the following article:
%  [1] M. Grau, G. Lenz and P. Lenz, "SDCM: Inference of interactions from high-
%      dimensional datasets and survival prediction in diffuse large B-cell lymphoma", 
%      Nature Communications (submitted).
%  The GPAV algorithm is used for 1D monotonic regressions; for details see:
%  [2] O. Burdakov, O. Sysoev, A. Grimvall and M. Hussian. "An O(n2) algorithm for 
%      isotonic regression." In: G. Di Pillo and M. Roma (Eds) Large-Scale Nonlinear 
%      Optimization. Series: Nonconvex Optimization and Its Applications, Springer-
%      Verlag, (2006) 83, pp. 25-33.
%SUPPORT AND UPDATES
%  Please contact me via sciene dot Michael dot Grau at gmail. For support requests, 
%  use "SDCM support: <problem title>" as your eMail subject. Include or upload a 
%  ZIP with your L2Is and inInfo configuration saved as .mat file and a script that 
%  reproduces or demonstrates the problem. I will try to debug any reproducible problems 
%  on the basis of available time resources. (Questions that are answered by [1],
%  its SIs, the source code or by MATLAB's documentation have to be ignored; sorry.)
%LICENSE (colloquial)
%  Free to use and modify in science, as long as [1] is properly cited. Redistribution 
%  of the code or offerring indirect access to it (e.g. via a web service) requires  
%  written consent. Free to test for anyone. Free for personal use. A license is required 
%  as soon as you directly or indirectly generate revenue with it; in this case contact me  
%  for a commercial license. The code is provided AS IS with no warranties whatsoever.
%  (For a juridical formulation, see license.txt.)
%AUTHOR
%  (c) Michael Grau, 2012-2016; version 3.55, April 2016.

function outInfo = SDCM(L2Is, inInfo, previouslyDetectedSignatures)
  %%.Initialization and parameters logic:
  
    %Warn on missing dependencies:
      if(true)
        installedToolboxesInfo = ver;
        csRequiredToolboxes = {
          {'Statistics Toolbox','Statistics and Machine Learning Toolbox'} %'Statistics Toolbox' is the older name of this toolbox
          {'Bioinformatics Toolbox'}
          {'Parallel Computing Toolbox'}
        };
        for i=1:length(csRequiredToolboxes)
          if(~any(ismember(csRequiredToolboxes{i}, {installedToolboxesInfo.Name})))
            if(isscalar(csRequiredToolboxes{i}))
              sMissingToolbox = csRequiredToolboxes{i}{1};
            else
              sMissingToolbox = ['with a name in {',strtrim(evalc('disp(csRequiredToolboxes{i}(:))')),'}'];
            end
            warning('Dependencies check: MISSING TOOLBOX %s\n(errors will occur, once functions of this toolbox are invoked by SDCM)', sMissingToolbox);
          end
        end
      end
    %Check access to library functions:
      if(exist('SDCM_printStatus','file')~=2)
        sLibrarySubdir = 'SDCM_Library';
        sPathToSelf = fileparts(mfilename('fullpath'));
        if(exist(fullfile(sPathToSelf,sLibrarySubdir,'tools','SDCM_printStatus.m'),'file')>0)
          warning('Now adding SDCM library functions in [%s] to the Matlab path.', fullfile(sPathToSelf,sLibrarySubdir));
          addpath(genpath(fullfile(sPathToSelf,sLibrarySubdir)));
          assert(exist('SDCM_printStatus','file')==2, 'Cannot find SDCM library functions such as SDCM_printStatus.m even after adding %s to the Matlab path. Please check your installation and add the library functions manually to your Matlab path.', fullfile(sPathToSelf,sLibrarySubdir));
        else
          error('Cannot find SDCM library functions subfolder at %s. Please add the library functions manually to your Matlab path.', sPathToSelf);
        end
      end
    %Get configuration and fill in defaults for parameters not provided by the user:
      nG = size(L2Is,1);
      nP = size(L2Is,2);
      if(nargin('SDCM')<2) inInfo = struct(); end
      assert(isstruct(inInfo) && isscalar(inInfo), 'inInfo must be a scalar structure');
      defConfig = SDCM_defaultConfig(nG,nP);
        if(~isempty(setdiff(fieldnames(inInfo),fieldnames(defConfig))))
          warning('### The following unknown options in inInfo will be IGNORED: %s', strtrim(strrep(evalc('setdiff(fieldnames(inInfo),fieldnames(defConfig))'),'ans =','')));
        end
      inInfo = mergeStructs(inInfo, defConfig);
    %Initialize status output:
      if(inInfo.internal.bDevEnableInteractiveBreaks)
        interactiveCheckPoint(true); %initialize.
      end
      SDCM_printStatus('INTERNAL_configureStatusOutputLevelForThisSession', inInfo.export.nStatusOutputLevel);
      SDCM_printStatus('INTERNAL_configureCheckpointLevelForThisSession', iif(inInfo.internal.bDevEnableInteractiveBreaks, inInfo.export.nStatusOutputLevel4checkpoints, -Inf));
      %Also initialize in workers, if any:
        if(ilv(gcp('nocreate'),@(pool)iif(isempty(pool),0,@()pool.NumWorkers))>0) %matlabpool('size')==0)
          inInfo_internal_nStatusOutputLevel = inInfo.export.nStatusOutputLevel;
          spmd
            SDCM_printStatus('INTERNAL_configureStatusOutputLevelForThisSession', inInfo_internal_nStatusOutputLevel);
          end
        end
      SDCM_printStatus(0,'################################################################################################################################################################\n');
      SDCM_printStatus(0,'### INITILIZATION) Checking configuration...\n'); drawnow;
    %Perform some basic configuration checks:
      function checkParams()
        %Note: more sophisticated parameter validations will be added here to protect against possible user errors in future versions...

        %Input signal:
          %Check for real values:
            if(~isreal(L2Is))
              error('L2Is must be real-valued.');
            end
        
        %Metadata:
          assert(all(size(inInfo.reference.rowIDs)==[nG,1]), 'inInfo.reference.rowIDs must be of size nG*1');
          assert(all(size(inInfo.reference.colIDs)==[1,nP]), 'inInfo.reference.colIDs must be of size 1*nP');

        %Preprocessing:
          if(inInfo.preprocessing.dampenOutliers4G.bEnabled && inInfo.preprocessing.dampenOutliers4G.nTopRanksToCheck==0)
            if(nG>=5) %only warn, if this is not the default.
              warning('.preprocessing.dampenOutliers4G is enabled, but the .dampenOutliers4G.nTopRanksToCheck was zero which essentially disables this preprocessing step. This can happen in the default config for less than 5 input genes (currenty nG=%d), because in this case there is not enough information to dampen outliers along the genes dimension without possibly destroying valuable information about true signatures.', nG);
            end
          end
          if(inInfo.preprocessing.dampenOutliers4P.bEnabled && inInfo.preprocessing.dampenOutliers4P.nTopRanksToCheck==0)
            if(nP>=5) %only warn, if this is not the default.
              warning('.preprocessing.dampenOutliers4P is enabled, but the .dampenOutliers4P.nTopRanksToCheck was zero which essentially disabled this preprocessing step. This can happen in the default config for less than 5 input samples (currenty nP=%d), because in this case there is not enough information to dampen outliers along the samples dimension without possibly destroying valuable information about true signatures.', nP);
            end
          end

        %for STEP 3:
          if(~ismember(inInfo.preprocessing.numericTargetPrecision,{'single','double'})) 
            error('inInfo.preprocessing.numericTargetPrecision must be either ''single'' or ''double''.'); 
          end

        %application mode checks:
          if(exist('previouslyDetectedSignatures','var') && ~inInfo.applicationMode.bEnabled)
            error('The third input argument previouslyDetectedSignatures is only allowed if inInfo.applicationMode.bEnabled');
          end
          if(inInfo.applicationMode.bEnabled)
            assert(nargin('SDCM')>=3 && isstruct(previouslyDetectedSignatures), 'If inInfo.applicationMode.bEnabled, the third input argument previouslyDetectedSignatures must be a array of signature axes structures to apply.');
          end

        %Plot settings:
          if(ischar(inInfo.plots.bVisible) && strcmp(inInfo.plots.bVisible, 'DEFAULT'))
            inInfo.plots.bVisible = ~inInfo.export.plots.bEnabled';
          end
          if(~inInfo.export.plots.bEnabled && ~inInfo.plots.bVisible)
            if(inInfo.plots.focusConvergence.bEnabled)
              warning('neither .export.plots.bEnabled export nor .plots.bVisible => disabling inInfo.plots.focusConvergence.bEnabled');
              inInfo.plots.focusConvergence.bEnabled = false;
            end
            if(inInfo.plots.signatureDefinitions.bEnabled)
              warning('neither .export.plots.bEnabled export nor .plots.bVisible => disabling inInfo.plots.signatureDefinitions.bEnabled');
              inInfo.plots.signatureDefinitions.bEnabled = false;
            end
          end
        %Export settings:
          if(ischar(inInfo.export.plots.bCloseAfterExport) && strcmp(inInfo.export.plots.bCloseAfterExport, 'DEFAULT'))
            if(inInfo.plots.bVisible)
              inInfo.export.plots.bCloseAfterExport = false;
            elseif(~inInfo.plots.bVisible)
              inInfo.export.plots.bCloseAfterExport = true;
            end
          end
          if(~islogical(inInfo.export.matFile4eachSignatureDefinition.bEnabled) || ~isscalar(inInfo.export.matFile4eachSignatureDefinition.bEnabled))
            error('inInfo.export.matFile4eachSignatureDefinition.bEnabled must be a scalar logical.')
          end
          inInfo.export.bAnyFileExportEnabled = ...
               inInfo.export.matFile4configurationAndInputData.bEnabled ...
            || inInfo.export.matFile4eachSignatureDefinition.bEnabled ...
            || inInfo.export.copyCurrentAlgorithm.bEnabled ...
            || inInfo.export.bCreateDetectionLogFile ...
            || inInfo.export.plots.bEnabled ...
            || inInfo.export.postprocessingTables.bEnabled ...
            || inInfo.export.infoFile4initialSignatureStats.bEnabled ...
            || inInfo.export.infoFile4focusedSignatureStats.bEnabled ...
            || inInfo.export.infoFile4topMembers.bEnabled4topCorrGenes ...
            || inInfo.export.infoFile4topMembers.bEnabled4topAntiCorrGenes ...
            || inInfo.export.infoFile4topMembers.bEnabled4topCorrSamples ...
            || inInfo.export.infoFile4topMembers.bEnabled4topAntiCorrSamples ...
          ;
          if(inInfo.export.bAnyFileExportEnabled)
            if(~ischar(inInfo.export.rootDir) || strcmp(inInfo.export.rootDir,'REQUIRED'))
              warning('inInfo.export.rootDir must be a configured for any file export option (e.g. .export.matFile4eachSignatureDefinition.bEnabled || .export.copyCurrentAlgorithm.bEnabled || .export.bCreateDetectionLogFile || .export.plots.bEnabled || .export.postprocessingTables.bEnabled); now DISABLING all configured file export options.');
              %Disable all file export options:
                inInfo.export.matFile4configurationAndInputData.bEnabled = false;
                inInfo.export.matFile4eachSignatureDefinition.bEnabled = false;
                inInfo.export.copyCurrentAlgorithm.bEnabled = false;
                inInfo.export.bCreateDetectionLogFile = false;
                inInfo.export.plots.bEnabled = false;
                inInfo.export.table4eachSignatureDefinition.bEnabled = false;
                inInfo.export.postprocessingTables.bEnabled = false;
                inInfo.export.infoFile4initialSignatureStats.bEnabled = false;
                inInfo.export.infoFile4focusedSignatureStats.bEnabled = false;
                inInfo.export.infoFile4topMembers.bEnabled4topCorrGenes = false;
                inInfo.export.infoFile4topMembers.bEnabled4topAntiCorrGenes = false;
                inInfo.export.infoFile4topMembers.bEnabled4topCorrSamples = false;
                inInfo.export.infoFile4topMembers.bEnabled4topAntiCorrSamples = false;
                inInfo.export.bAnyFileExportEnabled = false;
              %make figures visible (no point in plotting invisible figures and closing them without export):
                inInfo.plots.bVisible = true;
                inInfo.export.plots.bCloseAfterExport = false; 
            elseif(~exist(inInfo.export.rootDir,'dir'))
              mkdir(inInfo.export.rootDir);
            end
          end
          if(inInfo.export.plots.bEnabled && (~iscellstr(inInfo.export.plots.exportFormats) || isempty(inInfo.export.plots.exportFormats)))
            error('inInfo.export.plots.exportFormats must be a cell string of file extensions for figure exports, if plot exports are enabled.')
          end          

        %Configruation hints:
          %Small sample cohort / general warning:
            if((nG<50 || nP<50) && ~inInfo.internal.bAllowDissectionOfTinyDatasets)
              warning('DATASET SMALL: SDCM was originally designed only for large datasets like 50000 probesets times 500 samples, but the current input size is rather small (<50 samples or genes). Signature detection is based on correlations. For small sample sets it is increasingly easier for noise to disturb the correlation footprint of true interactions or to spawn noise signatures. Therefore the signatures returned by this detection run should not be used without any proper validation in independent datasets. You also might consider increasing e.g. the default inInfo.searchStrategy.qualification.minAbsCorrSum4P to prevent false positives.');
            end
          %Tiny sample cohort / error:
            if((nG<10 || nP<10) && ~inInfo.internal.bAllowDissectionOfTinyDatasets)
              error(sprintf('DATASET TOO SMALL, dissection denied: SDCM was originally designed only for large datasets like 50000 probesets times 500 samples. The current input size is too small (<10 samples or <10 genes). This could lead to artefacts, because some numeric parts of the current implementation are still in experimental state for such small numbers of samples or genes.\n(For scientific purposes in the context of algorithm development, you may continue anyway by setting inInfo.internal.bAllowDissectionOfTinyDatasets=true.)'));
            end
      end
      checkParams(); 
    %Initialize detection log, if enabled:
      if(inInfo.export.bCreateDetectionLogFile)
        sLogFileName = iif(inInfo.applicationMode.bEnabled, sprintf('signatureApplication_kStart=%03d.log', inInfo.applicationMode.kOffset4filenames+1), 'detection.log');
        SDCM_printStatus(0, '   - initializing log file [%s]...\n', sLogFileName);
        diary off;
        if(exist(fullfile(inInfo.export.rootDir,sLogFileName),'file') && ~inInfo.internal.bResume) %reset log:
          delete(fullfile(inInfo.export.rootDir,sLogFileName)); 
        end
        diary(fullfile(inInfo.export.rootDir,sLogFileName));
      end
      SDCM_printStatus(1, '   - version: %0.2f %s, executed from: [%s].\n', inInfo.reference.version, inInfo.reference.versionText, mfilename('fullpath'));
    %Initialize output structure:
      outInfo = struct();  
        outInfo.explanations.reference.initialSignal.originalL2IsBeforePreprocessing = 'original input log2(intensities) before any changes by enabled preprocessing components';
                     outInfo.reference.initialSignal.originalL2IsBeforePreprocessing = L2Is;

  %% Preparations: preprocess the initial signal and prepare the state variable for detection and dissection:
    SDCM_printStatus(0,'################################################################################################################################################################\n');
    SDCM_printStatus(0,'### PREPROCESSING) Prepare raw input signal for signature detection and dissection:\n'); drawnow;
      outInfo.explanations.reference.timing = 'time stamps (in the format of Matlab''s now function) at the various steps of the algorithm';
        outInfo.reference.timing.a_startOfPreprocessing = now;
    %Save memory by conversion to single precision:
      if(~isa(L2Is,inInfo.preprocessing.numericTargetPrecision))
        switch(inInfo.preprocessing.numericTargetPrecision)
          case 'single'
            L2Is = single(L2Is);
          case 'double'
            L2Is = double(L2Is);
        end
      end
    %Exclude duplicate rows, if any:
      if(inInfo.preprocessing.duplicateRows.bCheckAndRemove)
        SDCM_printStatus(0, '   - checking for identical duplicate input rows...\n');        
        [M,SI] = sortrows(L2Is);
        BValueDuplicateRows = sum(diff(M,1),2)==0;
        if(any(BValueDuplicateRows))
          outInfo.explanations.preprocessing.IRemovedDuplicateValueRows = 'rows removed from L2Is since they were duplicates of other rows';
                       outInfo.preprocessing.IRemovedDuplicateValueRows = SI(BValueDuplicateRows);
          if(inInfo.preprocessing.duplicateRows.bErrorIfFound)
            error(' - preprocessing/length(outInfo.preprocessing.IRemovedDuplicateValueRows)=%d source data rows\n            are identical value duplicates of other source data rows! This is an error per inInfo.preprocessing.duplicateRows.bErrorIfFound configuration. Remove the duplicates on caller level.', length(outInfo.preprocessing.IRemovedDuplicateValueRows));
          end
          warning(' - preprocessing/REMOVING length(outInfo.preprocessing.IRemovedDuplicateValueRows)=%d source data rows\n            as they are identical value duplicates of other source data rows!', length(outInfo.preprocessing.IRemovedDuplicateValueRows));
          L2Is(outInfo.preprocessing.IRemovedDuplicateValueRows,:) = [];
        end
      end
    %Dampen outliers, if enabled:
      if(inInfo.preprocessing.dampenOutliers4G.bEnabled)
        SDCM_printStatus(0, '   - searching and dampening outliers in the sample eigenOrder of every gene...\n');        
        L2Is = dampenOutliers(L2Is...
          ,inInfo.preprocessing.dampenOutliers4G.maxAllowedNeighboursRatio...
          ,inInfo.preprocessing.dampenOutliers4G.nTopRanksToCheck...
        ,2);
      end
      if(inInfo.preprocessing.dampenOutliers4P.bEnabled)
        SDCM_printStatus(0, '   - searching and dampening outliers in the gene eigenOrder of every sample...\n');        
        L2Is = dampenOutliers(L2Is...
          ,inInfo.preprocessing.dampenOutliers4G.maxAllowedNeighboursRatio...
          ,inInfo.preprocessing.dampenOutliers4G.nTopRanksToCheck...
        ,1);
      end
    %Check and warn if data contains Infs; replace them by NaNs:
      nInfs = sum(isinf(L2Is(:)));
      if(nInfs>0)
        warning('Your signal contains %d infinity values (Inf ratio = %0.1f%%). Correlations cannot count Infs and thus they will be replaced by NaNs (missing values) now.'...
          ,nInfs, nInfs/numel(L2Is)*100 ...
        );
        L2Is(isinf(L2Is)) = NaN;
      end
    %Check if data contains NaNs:
      nNaNs = sum(isnan(L2Is(:)));
      if(nNaNs==numel(L2Is))
        warning('The input data is COMPRISED OF NaNs; nothing to detect!');
        inInfo.preprocessing.bQuantilenormalizeSamples = false; %Matlab's quantilenorm function crashes on allNaNs inputs; disable as it is irrelevant anyway in this case.
      end
      inInfo.preprocessing.bDataContainsNaNs = nNaNs>0;
      if(~inInfo.preprocessing.bSetNaNsToZeroL2R)
        SDCM_printStatus(0, '   - Your signal contains %d missing values (NaN ratio = %0.1f%%). Note: This requires the use of slower NaN-robust aggregation functions like nanmean instead of mean. You might want to consider using inInfo.preprocessing.bSetNaNsToZeroL2R for speed. If quality maximization is of primary interest, keep this configuration and let the algorithm impute missing values based on the eigenorders detected for the non-missing values.'...
          ,nNaNs, nNaNs/numel(L2Is)*100 ...
        );
      end
    %Get basic math function handles according to bDataContainsNaNs:
      BM = getBasicMathFunctions(inInfo.preprocessing.bDataContainsNaNs);
    %Quantile-normalize all samples, if enabled:
      if(inInfo.preprocessing.bQuantilenormalizeSamples)
        SDCM_printStatus(0, '   - quantile-normalize all sample arrays...\n');        
        if(isempty(which('quantilenorm')))
          error('  <- quantilenorm NOT AVAILABLE: inInfo.preprocessing.bQuantilenormalizeSamples was enabled, but the local Matlab installation does not know the quantilenorm function. Either intall the Bioinformatics Toolbox or normalize the L2Is input matrix on caller level and disable inInfo.preprocessing.bQuantilenormalizeSamples.');
        end
        try
          BAllNaNColumns = all(isnan(L2Is),1); %Workaround as quantilenorm crashes for allNaN rows/cols.
          BAllNaNRows = all(isnan(L2Is),2);
          L2Is(~BAllNaNRows,~BAllNaNColumns) = quantilenorm(L2Is(~BAllNaNRows,~BAllNaNColumns));
          if(~inInfo.preprocessing.bDataContainsNaNs && ~all(~isnan(L2Is(:))))
            warning('<- quantilenorm produced %d NaNs; setting them to zero', sum(isnan(L2Is(:)))); %this happened once with NaN free L2Is (only 100 of 55000*100, one NaN in every column, but in different rows; cause unknown)
            L2Is(isnan(L2Is)) = 0;
          end
        catch ex  %Note: quantilenorm cannot handle very high NaN rates near 100%.
          rethrow(addCause(ex,MException('MATLAB:quantilenorm','Error detected in Matlab''s quantilenorm function; consider disabling inInfo.preprocessing.bQuantilenormalizeSamples')));
        end
      end
    %Precompute inInfo.reference.rowW and inInfo.reference.colW to correctly count genes in case of multiple rows (probesets) per gene:
      if(true)
        SDCM_printStatus(0, '   - precomputing row and column weights based on inInfo.reference.rowGroups and inInfo.reference.colGroups respectively (so that multiple probesets for the same gene count as one and not as many genes)\n');
        uniqueRowGroups = unique(inInfo.reference.rowGroups);
        if(length(uniqueRowGroups)<nG)
          groupSizes = joinAndFill(inInfo.reference.rowGroups,(1:nG)',uniqueRowGroups,@(I)length(I));
          inInfo.reference.rowW = 1./joinAndFill(uniqueRowGroups,groupSizes,inInfo.reference.rowGroups);
            inInfo.reference.rowW(isnan(inInfo.reference.rowW)) = 1; %NaNs are treated as separate groups.
        else
          inInfo.reference.rowW = ones(nG,1);
        end
        uniqueColGroups = unique(inInfo.reference.colGroups);
        if(length(uniqueColGroups)<nP)
          groupSizes = joinAndFill(inInfo.reference.colGroups,1:nP,uniqueColGroups,@(I)length(I));
          inInfo.reference.colW = 1./joinAndFill(uniqueColGroups,groupSizes,inInfo.reference.colGroups);
          inInfo.reference.colW(isnan(inInfo.reference.colW)) = 1; %NaNs are treated as separate groups.
        else
          inInfo.reference.colW = ones(1,nP);
        end    
      end

    %Get the link to the global state structure (use a global var for the resume feature and also for easier development):
      global eDState; 
      if(isempty(eDState))
        eDState = struct();
      end
    %Store final initial signal for reference: 
      eDState.initialSignal.L2Is = L2Is; %start with (possibly preprocessed) log2(intensities)
    %Compute initial log2(ratio)s relative to the cohort average (subtract per-gene brightness/sensitivity offset signature):
      SDCM_printStatus(0, '   - define initial log2(ratios) for detection\n');
      if(inInfo.preprocessing.bConvertIntoCohortRelativeL2RsAtStart)
        SDCM_printStatus(0, '     <- subtracting gene offsets to transform log2(intensities) into log2(ratios)...\n');
        eDState.initialSignal.base4G = BM.mean(eDState.initialSignal.L2Is,2);
        eDState.initialSignal.base4P = zeros(1,nP);
      else
        SDCM_printStatus(0, '     <- leaving any possible brightness offsets per gene in the signal (will be subtracted mainly by the first detected signature)\n');
        eDState.initialSignal.base4G = zeros(nG,1);
        eDState.initialSignal.base4P = zeros(1,nP);
        %eDState.initialSignal.base4P = fcnMedian(eDState.initialSignal.L2Is,1); %common reference brightness after quantile normalization for all samples.
      end
      eDState.initialSignal.L2Rs = bsxfun(@plus,-eDState.initialSignal.base4G,bsxfun(@plus,-eDState.initialSignal.base4P,eDState.initialSignal.L2Is)); %needed for quality scores.
      %Zero NaNs in L2Rs, if this imputation is enabled:
        if(inInfo.preprocessing.bSetNaNsToZeroL2R)
          SDCM_printStatus(0, '   - replacing NaNs with a log2(ratio) of zero\n');
          eDState.initialSignal.L2Rs(isnan(eDState.initialSignal.L2Rs)) = 0;
          inInfo.preprocessing.bDataContainsNaNs =  false;
          %Update basic math functions (get faster functions without NaN support):
            BM = getBasicMathFunctions(inInfo.preprocessing.bDataContainsNaNs);
        end
    %Initialize the signal of the current/first detection iteration with the initial signal the thus-far explained signal with zeros:
      eDState.current.L2Rs = eDState.initialSignal.L2Rs;
      eDState.current.explainedSignal = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision);
      eDState.current.k = 0; %we are before the first detected signature.
      eDState.current.multiPassesLevel = 1; %support for different detection parameters per pass; start detection with configuration level 1.
      %Support for inInfo.applicationMode.externalInitialShifts (subtract externally explained or previously learnt signal parts):
        if(inInfo.applicationMode.bEnabled)
          assert(all(~isnan(inInfo.applicationMode.externalInitialShifts(:))), 'some inInfo.applicationMode.externalInitialShifts were NaN.');
          eDState.current.explainedSignal = eDState.current.explainedSignal + inInfo.applicationMode.externalInitialShifts;
          eDState.current.L2Rs = eDState.current.L2Rs + inInfo.applicationMode.externalInitialShifts;
        end
    %Initialize/reset the performance state:
      eDState.current.performance.bNoSignatureDetectedSinceLastBDidNotQualifyReset = true;
      eDState.current.performance.BGeneDidNotQualify = false(nG,1);
      eDState.current.performance.BSampleDidNotQualify = false(1,nP);
      if(~inInfo.internal.bResume)
        eDState.current.precompute_previousIterations = []; %clear cache.
      else
        %will be cleared below.
      end
    %Initialize/reset the noise estimation:
      eDState.noiseEstimation = struct();

    %RESUME feature: optionally restore previously detected signatures for the same input data and resume detection&dissection:
      %Optionally load previously detected signatures from disk:
        if(inInfo.internal.bResumeFromDisk)
          warning('inInfo.internal.bResumeFromDisk feature: about to REPLACE the current eDState.signatures by loading all ''###, definition.mat'' files in [%s]; confirm with >>dbcont', inInfo.export.rootDir);
          keyboard;
          eDState.signatures = [];
          k=0; while(true) k=k+1;
            fileName = fullfile(inInfo.export.rootDir, sprintf('%03d, definition.mat',k));
            if(~exist(fileName,'file')) break; end
            temp = load(fileName);
              signatureDef = temp.signatureDef;
              %mem savers:
                %Mem saver / remove redundant reference info from signature definition:
                  if(isfield(signatureDef.reference,'inInfo'))
                    signatureDef.reference = rmfield(signatureDef.reference,'inInfo'); %this is only included to make the signature saved to disc become self-contained.
                  end
                %Mem saver / remove plotting info from RAM, if enabled:
                  if(inInfo.export.matFile4eachSignatureDefinition.bEnabled && inInfo.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport)
                    if(isfield(signatureDef,'forPlots')) %make structs array compatible wrt. fields, if we resumed from different save plotsInfo config.
                      signatureDef = rmfield(signatureDef, 'forPlots');
                    end
                  end
              clear temp;
            if(k==1)
              eDState.signatures = signatureDef;
            else
              eDState.signatures(k) = signatureDef;
            end
            fprintf('imported previousy exported signature %03d from [%s]\n', k, inInfo.export.rootDir); drawnow;
          end
        end
      if(inInfo.internal.bResume)
        bCanResume = ...
             ~isempty(eDState) && isstruct(eDState) && isfield(eDState,'signatures') && ~isempty(eDState.signatures) ...
          && length(eDState.signatures(1).step3_regression.BNonzeroRows4signatureEigensignal)==nG ...
          && length(eDState.signatures(1).step3_regression.BNonzeroCols4signatureEigensignal)==nP ...
          && (abs( nansum(abs(eDState.current.L2Rs(:))) - eDState.signatures(1).reference.forResume.reconstrutionChecksums.L2Rs_before ) < 10^-2) ...
        ;
        if(bCanResume)
          %Select number of previously detected signatures to recover:
            if(inInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives) 
              warning('Resuming with inInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives enabled may lead to not perfectly reproducible results, because the cache state cannot be restored and detection will resume with an empty cache (the cache is not saved with the other resume information to disk, because it is generally too large and this is too inefficient).');
            end
            maxkInMem = find(cellfun(@(s)isfield(s.step3_regression,'signatureEigensignal')&&~isempty(s.step3_regression.signatureEigensignal),num2cell(eDState.signatures)),1,'last');
            eDState.current.k = maxkInMem;
              %<-this k plus one will become the next signature that will be searched for.
            %Let the user check k:
              SDCM_printStatus(0,'# - RESUME: will restore eDState.current.k = %d signatures from memory and then search for the next; reduce eDState.current.k, if this is not correct, otherwise use dbcont.\n', eDState.current.k);
              keyboard;
            %eDState.signatures(eDState.current.k+1:end) = []; %no need to delete; will be overwritten.
          %Cache resume logic (only available for end-resume)
            if(inInfo.internal.bReuseCacheOnEndResume && ~isempty(eDState.current.precompute_previousIterations) && eDState.current.k==maxkInMem)
              warning('REUSING eDState.current.precompute_previousIterations for this resume-from-last-detected');
            elseif(inInfo.internal.bReuseCacheOnEndResume)
              warning('inInfo.internal.bReuseCacheOnEndResume was true, but either no cache is present in memory or eDState.current.k<maxkInMem was chosen; the cache will now be cleared; use >>dbcont to confirm');
              keyboard;
              eDState.current.precompute_previousIterations = []; %clear cache.
            else
              warning('inInfo.internal.bReuseCacheOnEndResume is false => about to clear eDState.current.precompute_previousIterations; use >>dbcont to confirm');
              keyboard;
              eDState.current.precompute_previousIterations = []; %clear cache.
            end
          %Reset current state to the initial minus all selected so far detected signatures:
            if(eDState.current.k>0)
              %Define the current signal for the the detection start as the initial signal minus so far detected signatures:
                currentShifts = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision); %preallocate.
                for k=1:eDState.current.k
                  %%Check initial L2Rs:
                  %  if(k==2 || k==eDState.current.k) %performance/only check after first reconstruction and before last last:
                  %    assert(abs( nansum(abs(eDState.current.L2Rs(:))) - eDState.signatures(k).reference.forResume.reconstrutionChecksums.L2Rs_before ) < 10^-2, 'code validation: L2Is_before reconstruction checksum mismatch. Probably a different detection run has overwritten the global eDState store in the meantime since the last detection run on the current input data.');
                  %  end
                  %Initialize:
                    currentShifts(:) = 0;
                  %Add dissection shifts:
                    currentShifts(...
                      eDState.signatures(k).step3_regression.BNonzeroRows4signatureEigensignal...
                     ,eDState.signatures(k).step3_regression.BNonzeroCols4signatureEigensignal...
                    ) = -eDState.signatures(k).step3_regression.signatureEigensignal;
                  %Check reconstructed shifts:
                    if(k==eDState.current.k) %performance/only check for last:
                      assert(abs( nansum(abs(currentShifts(:))) - eDState.signatures(k).reference.forResume.reconstrutionChecksums.shifts ) < 10^-2, 'shifts reconstruction checksum mismatch');
                    end
                  %Apply shifts and check result:
                    if(~(inInfo.applicationMode.bEnabled && inInfo.applicationMode.bApplyEveryExternalSignatureToInitialSignal))
                      eDState.current.L2Rs = eDState.current.L2Rs + currentShifts;
                      if(  (k==1 || k==eDState.current.k) ... %performance/only check for first/last.
                        && isfield(eDState.signatures(k).reference.forResume.reconstrutionChecksums,'L2Is_after') ... %field not available, if loaded from disk.
                      )
                        assert(abs( nansum(abs(eDState.current.L2Rs(:))) - eDState.signatures(k).reference.forResume.reconstrutionChecksums.L2Rs_after ) < 10^-2, 'code validation: L2Is_after reconstruction checksum mismatch. Probably a different detection run has overwritten the global eDState store in the meantime since the last detection run on the current input data.');
                      end
                      eDState.current.explainedSignal = eDState.current.explainedSignal + currentShifts;
                    end
                  SDCM_printStatus(0, '      <- reconstructed shifts for already detected signature %d/%d\n',k,eDState.current.k);
                end
              %Start detection with the configuration level of the last resumed signature:
                eDState.current.multiPassesLevel = eDState.signatures(eDState.current.k).reference.configurationLevel;
              %Restore performance state:
                eDState.current.performance.BGeneDidNotQualify = eDState.signatures(eDState.current.k).reference.forResume.BGeneDidNotQualifyAfterDetectionOfThisSignature; %needed for resume feature.
                eDState.current.performance.BSampleDidNotQualify = eDState.signatures(eDState.current.k).reference.forResume.BSampleDidNotQualifyAfterDetectionOfThisSignature; %needed for resume feature.
                eDState.current.performance.bNoSignatureDetectedSinceLastBDidNotQualifyReset = eDState.signatures(eDState.current.k).reference.forResume.bNoSignatureDetectedSinceLastBDidNotQualifyReset; %needed for resume feature.
                if(isstruct(eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation))
                  eDState.noiseEstimation = eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation; 
                elseif(iscell(eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation) && ischar(eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation{1}) && strcmp(eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation{1},'loadFromFile')) %support for mem saver
                  try 
                    temp = load(eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation{2});
                    eDState.noiseEstimation = ilsub(temp,eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation{3});
                  catch ex
                    warning('ERROR loading noise estimation from [%s]; use >>dbcont to ignore and reinitialize noise estimation from scratch', eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation{2});
                    keyboard;
                    eDState.noiseEstimation = struct(); %let noise estiamtion be restarted/reinitialized below.
                  end
                  %eDState.signatures(eDState.current.k).reference.forResume.noiseEstimation = eDState.noiseEstimation;
                else
                  error('eDState.signatures(%d).reference.forResume.noiseEstimation is not of an expected format; cannot resume.', eDState.current.k);
                end
            else
              eDState.signatures = []; %reset to clean state.            
            end
        else
          warning('inInfo.internal.bResume is true, but no compatible previously detected signatures in the global eDState were found to resume => use >>dbcont to overwrite the global eDState and detect from scratch or >>dbquit to abort');
          keyboard;
          eDState.signatures = []; %reset to clean state.
        end
      else
        eDState.signatures = []; %reset to start detection with a clean state.
      end
      %<-Note: for exact resume reproductivity inInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives must be false since eDState.current.precompute_previousIterations is too large to save with the resume info of every detected eigenorder.

    %Declare the signature template (contains fields for all signature information in steps 1 and 2, except for the signature's eigenOrder and regressed signal from step 3):
      if(true)
        eDState.signatureTemplate = struct();
        %step 1:
          eDState.signatureTemplate.greedyScore = NaN;
          eDState.signatureTemplate.imj = NaN;
          eDState.signatureTemplate.i = NaN;
          eDState.signatureTemplate.j = NaN;
        %Gene space:
          eDState.signatureTemplate.sampleAxis_origUnits = nan(1,nP);
            eDState.signatureTemplate.norm4sampleAxis_origUnits = nan;
          eDState.signatureTemplate.R4G = nan(nG,1);
            eDState.signatureTemplate.sampleSizes4R4G = nan(nG,1);
          eDState.signatureTemplate.P4R4G = nan(nG,1);
          eDState.signatureTemplate.signedFocusedW4G = nan(nG,1);
            eDState.signatureTemplate.signedFocusedW4G_withPerpendicularSpace = nan(nG,1);
          eDState.signatureTemplate.signedExtendedW4G = nan(nG,1);
        %sample space:
          eDState.signatureTemplate.geneAxis_origUnits = nan(nG,1);
            eDState.signatureTemplate.norm4geneAxis_origUnits = nan;
          eDState.signatureTemplate.R4P = nan(1,nP);
            eDState.signatureTemplate.sampleSizes4R4P = nan(1,nP);
          eDState.signatureTemplate.P4R4P = nan(1,nP);
          eDState.signatureTemplate.signedFocusedW4P = nan(1,nP);
            eDState.signatureTemplate.signedFocusedW4P_withPerpendicularSpace = nan(1,nP);
          eDState.signatureTemplate.signedExtendedW4P = nan(1,nP);
        %signature statistics:
          %Size:
            eDState.signatureTemplate.signatureSizeByCorrSum2D = NaN;
              eDState.signatureTemplate.signatureSizeByCorrSum4G = NaN;
              eDState.signatureTemplate.signatureSizeByCorrSum4P = NaN;
          %Signal strength:
            eDState.signatureTemplate.signatureAbsMean2D = NaN;
              eDState.signatureTemplate.sampleSize4signatureAbsMean2D = NaN;
            eDState.signatureTemplate.signatureAbsSD2D = NaN;
          %Average correlation:
            eDState.signatureTemplate.signatureCorrInExtendedFocus = NaN;            
          %Significance:
            eDState.signatureTemplate.log10_p = NaN;
              eDState.signatureTemplate.log10_p4Correlations = NaN;
              eDState.signatureTemplate.log10_p4SignalStrength = NaN;
        %step 2:
          eDState.signatureTemplate.corrMaximizingMembers = struct();
            eDState.signatureTemplate.corrMaximizingMembers.ImJ = [];
            eDState.signatureTemplate.corrMaximizingMembers.sampleAxes_origUnits = nan(0,nP);
            eDState.signatureTemplate.corrMaximizingMembers.geneAxes_origUnits = nan(nG,0);
        %step 3:
          eDState.signatureTemplate.signatureStrengths4G = nan(nG,1);
          eDState.signatureTemplate.signatureStrengths4P = nan(1,nP);
        %<-Note: explanations of all signature fields are provided below in this file at the end of step 2 when assembling the output.
      end
    %Initialize fields to record the timing of every detection iteration:
      if(true)
        outInfo.reference.timing.b_startOfDetectionOrApplication = now;
        outInfo.reference.timing.b1_startOfStep1 = [];
        outInfo.reference.timing.b2_startOfStep2 = [];
        outInfo.reference.timing.b3_startOfStep3 = [];
        outInfo.reference.timing.b4_startOfExportsAndPlotting = [];
      end
    %Export reference information, if requested:
      if(true)
        %save the active inInfo configuration and the original.L2Rs used for detection:
          if(inInfo.export.matFile4configurationAndInputData.bEnabled)
            inputData = struct();
              inputData.original = eDState.initialSignal;
              inputData.inInfo = inInfo; %includes annotations4cols/rows, if provided.
            SDCM_printStatus(2,'  - reference export: saving configuration and the initial signal for reference to [%s]...\n', fullfile(inInfo.export.rootDir,'reference','configuration and input data.mat'));
            if(~exist(fullfile(inInfo.export.rootDir,'reference'),'dir')) mkdir(fullfile(inInfo.export.rootDir,'reference')); end;
            save(fullfile(inInfo.export.rootDir,'reference','configuration and input data.mat'), 'inputData');
            clear inputData;
          end        
        %copy the current version of this algorithm to the export dir for reference:
          if(inInfo.export.copyCurrentAlgorithm.bEnabled)
            targetRootDir = fullfile(inInfo.export.rootDir,'reference',[ilnth(2,@fileparts,mfilename),'_v',num2str(inInfo.reference.version,'%0.2f')]);
              if(~exist(targetRootDir,'dir')) mkdir(targetRootDir); end;
            SDCM_printStatus(2,'  - reference export: saving a copy of this version of the dissection algorithm to [%s]...\n', targetRootDir);
            %Copy main function:
              sourceRootDir = fileparts(mfilename('fullpath'));
              copyfile(fullfile(sourceRootDir,'SDCM*.m'), [targetRootDir,filesep]);            
            %Copy library functions:
              csLibraryPaths = strsplit(genpath(fullfile(sourceRootDir,'SDCM_Library')),pathsep)';
                csLibraryPaths(cellfun(@isempty,csLibraryPaths)) = [];
              for i4paths=1:length(csLibraryPaths)
                sSourceDir = csLibraryPaths{i4paths};
                sTargetDir = fullfile(targetRootDir,strrep(sSourceDir,sourceRootDir,''));
                if(~isempty(dir(fullfile(sSourceDir,'*.m'))))
                  copyfile(fullfile(sSourceDir,'*.m'), [sTargetDir,filesep]);
                end
                if(~isempty(dir(fullfile(sSourceDir,'*.mex*'))))
                  copyfile(fullfile(sSourceDir,'*.mex*'), [sTargetDir,filesep]);
                end
              end
          end
      end
    SDCM_printStatus(0, '  ->ready for detection.\n');

  %% SDCM: Iteratively detect signatures and remove their own regressed bimontonic signal, until statistically only noise seems to be left:
    while(true)
      eDState.current.k = eDState.current.k + 1;
      %% STEP 1) search strategy: find the gene or sample that is the most promosing first representative of a new signature:
        %Get config level: (support for different detection parameters per all-genes/all-samples pass)
          sL = min(eDState.current.multiPassesLevel, inInfo.searchStrategy.nDefinedPasses); 
          %<-note: for each per parameter, only one value is defined by the default parameter set that is used for all passes, i.e. .nDefinedPasses=1 and hence sL==1 always for the default paramter set.
        %Application mode break condition: once all previouslyDetectedSignatures have been applied and dissected, return (useful for application of previously detected signatures to new cohorts):
          if(inInfo.applicationMode.bEnabled && eDState.current.k == length(previouslyDetectedSignatures)+1)
            SDCM_printStatus(0,'#<-GLOBAL BREAK CONDITION: all %d provided externally provided signatures were applied and dissected.\n', length(previouslyDetectedSignatures));
            break;
          end
        %Check an optional global break condition (maximum configured number of detected signatures reached?):
          if(eDState.current.k > inInfo.searchStrategy.globalBreakConditions.maxSignaturesToDetect(min(end,sL)))
            warning('#<- GLOBAL BREAK CONDITION: inInfo.searchStrategy.globalBreakConditions.maxSignaturesToDetect(min(end,sL)) limit reached!'); %this is a warning rather than a usual log entry, because usually signature dissection decides the stop condition based on the remaining signal (no more correlations qualify as signature).
            break;
          end
        %Status Output:
          SDCM_printStatus(0,'################################################################################################################################################################\n');
          if(~inInfo.applicationMode.bEnabled)
            SDCM_printStatus(0,'### DETECTING signature %02d (at %s with dataset size [%d,%d]; detection sensitivity level %d/%d):\n',eDState.current.k, datestr(now), nG,nP, sL,inInfo.searchStrategy.nDefinedPasses);
            SDCM_printStatus(0,'### STEP 1) Search strategy: Find the next initial representative, i.e. ideally the gene or sample with maximal signature functional in the current signal:\n'); drawnow;
          end
          if(inInfo.applicationMode.bEnabled)
            SDCM_printStatus(0,'### APPLYING signature %02d (at %s with dataset size [%d,%d]; detection sensitivity level %d/%d):\n',eDState.current.k, datestr(now), nG,nP, sL,inInfo.searchStrategy.nDefinedPasses);
            SDCM_printStatus(0,'### STEP 1+2) Prepare externally provided previously learnt signatures for application:\n'); drawnow;              
          end
          outInfo.reference.timing.b1_startOfStep1(end+1) = now;
        %Initialize:
          %Reset .BGeneDidNotQualify and .BSampleDidNotQualify, if disabled:
            if(~inInfo.searchStrategy.performance.bSkipGenesThatDidNotQualifyInPreviousRuns(min(end,sL)))
              eDState.current.performance.BGeneDidNotQualify(:) = false;
            end
            if(~inInfo.searchStrategy.performance.bSkipSamplesThatDidNotQualifyInPreviousRuns(min(end,sL)))
              eDState.current.performance.BSampleDidNotQualify(:) = false;
            end
          %Noise estimation:
            if(true)
              bInitializeNoiseEstimation = eDState.current.k==1 || ~isfield(eDState,'noiseEstimation') || ~isstruct(eDState.noiseEstimation) || ~isfield(eDState.noiseEstimation,'preciseNoiseSD');
                if(bInitializeNoiseEstimation && eDState.current.k~=1) %may be okay when interactively using the resume feature.
                  warning('code validation: bInitializeNoiseEstimation, but eDState.current.k~=1; inspect manually and >>dbcont to confirm');
                  keyboard;
                end
              %get the epsilon used for convergence detection:
                if(bInitializeNoiseEstimation)
                  eDState.noiseEstimation.sd4epsilon = nanstd(eDState.current.L2Rs(:)); %not used for noise estimation, but only for precision epsilons for convergence detection (here it is uncritical, if noise is overestimated even by a factor 2 or so).
                else
                  eDState.noiseEstimation.sd4epsilon = eDState.noiseEstimation.preciseNoiseSD(eDState.current.k); %use the .preciseNoiseSD estimated based on the empirical noise distribution from the last iteration (this is NaN for k==1).
                end
              %Initialize noise distribution fields:
                if(bInitializeNoiseEstimation)
                  %Disable signal strength based p values for qualification in iteration k==1: (It is better to find some false positives in a noisy signal if there are significant correlations, than to not qualify anything in a superposition-heavy signal where the signal strength of the top signature is about (or below) the overestimated noise level...)
                    eDState.noiseEstimation.preciseNoiseSD(eDState.current.k) = NaN;
                    eDState.noiseEstimation.halfNormalAbsNoiseMedian(eDState.current.k) = NaN;
                    eDState.noiseEstimation.absMean2D(eDState.current.k) = NaN;
                    eDState.noiseEstimation.nsBehindAbsMean2D(eDState.current.k) = NaN;
                    eDState.noiseEstimation.SD4absMean2D(eDState.current.k) = NaN;
                    %performance for large datasets/exclude at least some noise genes based on .noiseEstimation.sd4epsilon=std(L2Rs(:)) in k==1 (the std estimator may overestimate the noise level in signals with many superposed signatures, but this estimation levels-in over k and there is no need to respect all the weak noise genes in large datasets when computing correlations)
                      if(true)
                        %Notes:
                        % - die Radialnormen sqrt(sum(x_i^2)) folgen einer chi-Verteilung mit nG bzw. nP Freiheitsgraden
                        % - die Radialvarianzen sum(x_i^2) folgen einer chi^2-Verteilung mit nG bzw. nP Freiheitsgraden
                        Z = eDState.current.L2Rs / (eDState.noiseEstimation.sd4epsilon/2); %factor 1/2 usually overcompensated noise overestimation to noise underestimation which is desired for k==1 to get a good first signature (underestimation causes larger selected performance subspace and makes it slower, but this is better than a too tight performance subspace that might influence results.
                        chiSquared4G = BM.meanW(Z.^2,1,2,true);
                        nu4G = BM.meanW(~isnan(Z),1,2,true); %make ready for missing values.
                          eDState.noiseEstimation.log10PNoise4G = exactUpperChi2CDF(chiSquared4G,nu4G)/log(10);
                        chiSquared4P = BM.meanW(Z.^2,1,1,true);
                        nu4P = BM.meanW(~isnan(Z),1,1,true); %make ready for missing values.
                          eDState.noiseEstimation.log10PNoise4P = exactUpperChi2CDF(chiSquared4P,nu4P)/log(10);
                          %figure;scatter(log10(eDState.noiseEstimation.log10PNoise4P),exactUpperChi2CDF(chiSquared4P,nG)/log(10));
                      end
                      %<-Note: this reduces the BsPerformanceSubspace4Correlations for step 1 to genes above the estiamted noise level according to inInfo.searchStrategy.performance.alpha4Signal4inclusionInComputation4G,
                      %        which can speed up computation dramatically in iteration k==1 for large datasets.
  %                   eDState.noiseEstimation.nullCDF.commonX4AbsL2RsAboveThreshold = linspace(...
  %                      0 ...
  %                     ,min(nanmax(abs(eDState.initialSignal.L2Rs(:))),nanmean(eDState.initialSignal.L2Rs(:))+5*nanstd(eDState.initialSignal.L2Rs(:))) ...
  %                     ,1000 ...
  %                   );
                else
                  %Note: .preciseNoiseSD for the current k has already been updated after dissection in the last iteration.
                  SDCM_printStatus(2,'   - info: current estimated pure noise SD = %0.2f\n', eDState.noiseEstimation.preciseNoiseSD(eDState.current.k));
                end
            end
          %Initialize a new signature definition structure:
            signatureDefinition = struct();
              signatureDefinition.explanations = struct();
              signatureDefinition.reference = struct();
        %Signal standardization (used, e.g., as part of the greedy scores for presorting candidates in step 1):
          %Standardize in 2D (both rows and cols) to get rid of tumor load and gene affinity signatures (as preparation for smooth dissection and for computation of standardized signatures like .sdGeneAxis that are used for comparison of signatures and their valiation):
            SDCM_printStatus(2,'   - standardize the signal in 2D to get rid of tumor load and gene affinity signatures before the 2D bimonotonic baseline estimation...\n'); drawnow;
            eDState.current.stds4L2Rs = std2D(eDState.current.L2Rs, BM.uncenteredVar, eDState.noiseEstimation.sd4epsilon/1000); %/10000); %use a small epsilon to avoid accumulating numeric differences that break the conceptual transpose symmetry dissectSignal(X)=dissectSignal(X')'.
            eDState.current.sdL2Rs = bsxfun(@times, eDState.current.L2Rs, 1./eDState.current.stds4L2Rs);
        %SEARCH for an initial representative for a new/thereby detected signature:
          if(~inInfo.applicationMode.bEnabled)
            %Find correlating genes/samples for the initial representative candidates and stop as soon as we find enough and if there is no stronger signature in the lookahead-sight:
              signatureDefinition.explanations.step1_initialRepresentative = 'Initial signature structure resulting from step 1. Based on an identified initial representative gene or sample. The subfield .imj identifies the particular gene (positions 1:nG) respectively sample (positions -(1:nP)). The structure contains the signature''s .sampleAxis and .geneAxis, correlations .R4G and .R4P of all genes and samples with these signature axes and the signature focus/weights. For detailed explanations see .step2_finalSignatureAxes.';
                signatureDefinition.step1_initialRepresentative = SDCM_step1_searchStrategy(inInfo);
            %Check the primary global break condition (is no gene or sample qualifying any longer?):
              if(isempty(signatureDefinition.step1_initialRepresentative))
                %Genes/Samples that were not qualifying as initial representatives earlier may now have become eligible due to data cleaning by dissection of already detected signatures:
                  %=>we need to reset BGeneDidNotQualify as long as signatures get detected and search again without any exclusions:
                  if(  ~eDState.current.performance.bNoSignatureDetectedSinceLastBDidNotQualifyReset ... %if we found another signature in the last sweep
                    || eDState.current.multiPassesLevel < inInfo.searchStrategy.nDefinedPasses ... %or if there are more sensitive configurations left to test:
                  )
                    SDCM_printStatus(0,'#<- No more initial representative candidates availble. RESETTING .performance.BGeneDidNotQualify, .performance.BSampleDidNotQualify and the precompute cache to check for emerging correlations due to previous dissections and to retest global break conditions over all genes and samples once more...\n'); drawnow;
                    eDState.current.performance.bNoSignatureDetectedSinceLastBDidNotQualifyReset = true;
                    eDState.current.k = eDState.current.k-1;
                    %Reset performance state, since performed dissections may have "enabled" previosly uncorrelated genes:
                      eDState.current.performance.BGeneDidNotQualify = false(nG,1);                    
                      eDState.current.performance.BSampleDidNotQualify = false(1,nP);
                    %Reset/Clear the precompute cache:
                      precomputeCandidatesInParallel(eDState.current.precompute.ImJ, 'clearFromCache');
                      eDState.current.precompute_previousIterations = [];
                      %<-Note: it is important to reset the cache: even if all cached items taking part in the neutralistion mask of detected signatures have already been invalidated, other genes/samples may be indirectly affted (changing correlation scores due to changes in the signal of their genes/samples) 
                    %Multi passes feature:
                      eDState.current.multiPassesLevel = eDState.current.multiPassesLevel + 1;
                      sL = min(eDState.current.multiPassesLevel, inInfo.searchStrategy.nDefinedPasses);
                      SDCM_printStatus(0 ...
                        ,[ '#-> Advancing to pass %d (configuration level %d/%d) with detection parameters:\n'...
                          ,'     - initial representative qualification: .minAbsCorrSum4G=%0.2f, .minAbsCorrSum4P=%0.2f, .minAbsCorrSumSum=%0.2f, .minCorrInExtendedFocus=%0.2f, log10(.alpha)=(%6.1fcorr,%6.1fsignal,%6.1fcombined)\n'... .minSignatureSNR=%0.1f, 
                          ...
                         ]...
                        ,eDState.current.multiPassesLevel, sL, inInfo.searchStrategy.nDefinedPasses ...
                        ,inInfo.searchStrategy.qualification.minAbsCorrSum4G(min(end,sL))...
                        ,inInfo.searchStrategy.qualification.minAbsCorrSum4P(min(end,sL))...
                        ,inInfo.searchStrategy.qualification.minAbsCorrSumSum(min(end,sL))...
                        ,inInfo.searchStrategy.qualification.minCorrInExtendedFocus(min(end,sL))...
                        ,log10(inInfo.searchStrategy.qualification.alpha4correlations(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCorrP(min(end,sL)),nG+nP,1)) ...
                        ,log10(inInfo.searchStrategy.qualification.alpha4signalStrength(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForSignalP(min(end,sL)),nG+nP,1)) ...
                        ,log10(inInfo.searchStrategy.qualification.alpha4combined(min(end,sL))/iif(inInfo.searchStrategy.qualification.bConservativeBonferroniMultipleHypothesisCorrectionForCombinedP(min(end,sL)),nG+nP,1)) ...
                      );
                    continue;
                  end
                SDCM_printStatus(0,'#<- GLOBAL BREAK CONDITION: No gene and no sample are eligible as signature any more.\n');
                break;
              end
              eDState.current.performance.bNoSignatureDetectedSinceLastBDidNotQualifyReset = false;
          end
      %% STEP 2) correlation maximization for signature generalization and definition (climb to the center of the correlation-defined pattern in order to get independent of the initial representative)
        %MAXIMIZE the correlations represented by the signature by selecting additional representives until signature axes convergence:
          if(~inInfo.applicationMode.bEnabled)
            SDCM_printStatus(0,'### STEP 2) Select additional representives via correlation maximization until signature axes convergence, thereby generalizing the signature:\n'); drawnow;
              outInfo.reference.timing.b2_startOfStep2(end+1) = now;
            signatureDefinition.step2_finalSignatureAxes = SDCM_step2_maximizeCorrelation(signatureDefinition.step1_initialRepresentative, inInfo);
              %<-Note: all fields are explained in detail for the output below.
          end
        %Application stage) replace steps 1&2 by previously learnt and externally provided signatures and just apply them to the current signal:
          if(inInfo.applicationMode.bEnabled)
            %Support for resume feature:
              if(eDState.current.k == length(previouslyDetectedSignatures)+1)
                eDState.current.k = eDState.current.k - 1;
                warning('All previouslyDetectedSignatures were already processed by the resumed data => the main detection/application loop will not make any iteration!');
                break;
              end
              outInfo.reference.timing.b2_startOfStep2(end+1) = now;
            %Build signature from external information (typically a signature axis detected in a training dataset or a combined axis from multiple training datasets):
              %Initialize:
                externallyProvidedSignature = previouslyDetectedSignatures(eDState.current.k);
                assert(isfield(externallyProvidedSignature,'step2_finalSignatureAxes'), 'externally provided signature %d does neither define a .geneAxis_origUnits (and .R4G and .sampleSizes4R4P fields) nor a .sampleAxis_origUnits (and .R4P and .sampleSizes4R4P fields)', eDState.current.k);
                externallyProvidedSignatureAxes = externallyProvidedSignature.step2_finalSignatureAxes;
              %For which space is an axis provided?:
                bGeneAxisProvided = isfield(externallyProvidedSignatureAxes, 'geneAxis_origUnits') && isfield(externallyProvidedSignatureAxes, 'R4G') && isfield(externallyProvidedSignatureAxes, 'sampleSizes4R4G'); 
                bSampleAxisProvided = isfield(externallyProvidedSignatureAxes, 'sampleAxis_origUnits') && isfield(externallyProvidedSignatureAxes, 'R4P') && isfield(externallyProvidedSignatureAxes, 'sampleSizes4R4P');
                assert(bGeneAxisProvided || bSampleAxisProvided, 'the externally provided signature %d does neither define a .geneAxis_origUnits (and .R4G and .sampleSizes4R4P fields) nor a .sampleAxis_origUnits (and .R4P and .sampleSizes4R4P fields) in substructure .step2_finalSignatureAxes', eDState.current.k);
              %Symmetrize the signature by computing the axis in the twin space of the current input signal, if only one axis is provided:
                fcnNumericTargetPrecision = iif(strcmp(inInfo.preprocessing.numericTargetPrecision,'single'),@single,@double);
                if(bGeneAxisProvided && bSampleAxisProvided)
                  %completely specified signature; nothing to compute here (e.g. for replaying detected signatures on the detection cohort with updated dissection or plotting settings)
                elseif(bGeneAxisProvided && ~bSampleAxisProvided) %gene axis provided => compute twin sample axis for the current input signal to complete/symmetrize the signature for regression:
                  %Assert size consistency between provided signature axis and the corresponding order dimension in the input signal:
                    assert(all(size(externallyProvidedSignatureAxes.geneAxis_origUnits)==[nG,1]), 'the externally provided signature %d defines a gene axis .geneAxis_origUnits of %d elements in .step2_finalSignatureAxes, but the current input dataset has %d genes', eDState.current.k, length(externallyProvidedSignatureAxes.geneAxis_origUnits), nG);
                    assert(all(size(externallyProvidedSignatureAxes.R4G)==[nG,1]), 'the externally provided signature %d defines gene correlations .R4G with %d elements in .step2_finalSignatureAxes, but the current input dataset has %d genes', eDState.current.k, length(externallyProvidedSignatureAxes.R4G), nG);
                  %Enforce numeric type compatibility for numeric computations:
                    externallyProvidedSignatureAxes.geneAxis_origUnits = fcnNumericTargetPrecision(externallyProvidedSignatureAxes.geneAxis_origUnits);
                      externallyProvidedSignatureAxes.norm4geneAxis_origUnits = BM.euclidW(externallyProvidedSignatureAxes.geneAxis_origUnits, 1, 1)';
                    externallyProvidedSignatureAxes.R4G = fcnNumericTargetPrecision(externallyProvidedSignatureAxes.R4G);
                  %Compute gene focus:
                    externallyProvidedSignatureAxes.P4R4G = pValues4Correlations(externallyProvidedSignatureAxes.R4G, externallyProvidedSignatureAxes.sampleSizes4R4G);
                    externallyProvidedSignatureAxes.signedFocusedW4G = calcSignatureFocus([],[], externallyProvidedSignatureAxes.R4G, externallyProvidedSignatureAxes.P4R4G, 1, inInfo.searchStrategy.signatureFocus, sL);
                    externallyProvidedSignatureAxes.signedFocusedW4G_withPerpendicularSpace = addFlatWeights(externallyProvidedSignatureAxes.signedFocusedW4G, nan(0,nP), 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
                  %Compute twin sample axis (analogous to precomputeCandidatesInParallel):
                    externallyProvidedSignatureAxes.sampleAxis_origUnits = projectW(...
                      eDState.current.L2Rs ... vectorsInSourceSpace...
                     ,externallyProvidedSignatureAxes.geneAxis_origUnits ... axesInSourceSpace
                     ,abs(externallyProvidedSignatureAxes.signedFocusedW4G) ... weights4sourceDims
                     ,1 ... sourceSpaceMatrixDim
                    ,BM.meanW, BM.euclidW) / BM.euclidW(ones(size(externallyProvidedSignatureAxes.signedFocusedW4G,1),1),abs(externallyProvidedSignatureAxes.signedFocusedW4G),1); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
                      externallyProvidedSignatureAxes.sampleAxis_origUnits(isnan(externallyProvidedSignatureAxes.sampleAxis_origUnits)) = 0;
                      if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                        externallyProvidedSignatureAxes.sampleAxis_origUnits = sign(externallyProvidedSignatureAxes.sampleAxis_origUnits).*dampenOutliers(abs(externallyProvidedSignatureAxes.sampleAxis_origUnits), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                      end
                    externallyProvidedSignatureAxes.norm4sampleAxis_origUnits = BM.euclidW(externallyProvidedSignatureAxes.sampleAxis_origUnits, 1, 2)';
                  %Compute correlations of samples in the current dataset (analogous to precomputeCandidatesInParallel):
                    [externallyProvidedSignatureAxes.R4P, externallyProvidedSignatureAxes.sampleSizes4R4P, externallyProvidedSignatureAxes.P4R4P] = uncenteredWeightedCorrelation(...
                       externallyProvidedSignatureAxes.geneAxis_origUnits, eDState.current.L2Rs, 1 ...
                      ,abs(externallyProvidedSignatureAxes.signedFocusedW4G_withPerpendicularSpace) ...
                      ,BM.meanW ...
                      ,false ...
                    );
                      if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                        externallyProvidedSignatureAxes.R4P = sign(externallyProvidedSignatureAxes.R4P).*dampenOutliers(abs(externallyProvidedSignatureAxes.R4P), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 2); 
                      end
                elseif(bSampleAxisProvided && ~bGeneAxisProvided) %sample axis provided => compute twin gene axis for the current input signal to complete/symmetrize the signature for regression:
                  %Assert size consistency between provided signature axis and the corresponding order dimension in the input signal:
                    assert(all(size(externallyProvidedSignatureAxes.sampleAxis_origUnits)==[1,nP]), 'the externally provided signature %d defines a sample axis .sampleAxis_origUnits of %d elements in .step2_finalSignatureAxes, but the current input dataset has %d samples', eDState.current.k, length(externallyProvidedSignatureAxes.sampleAxis_origUnits), nP);
                    assert(all(size(externallyProvidedSignatureAxes.R4P)==[1,nP]), 'the externally provided signature %d defines sample correlations .R4P with %d elements in .step2_finalSignatureAxes, but the current input dataset has %d samples', eDState.current.k, length(externallyProvidedSignatureAxes.R4P), nP);
                  %Enforce numeric type compatibility for numeric computations:
                    externallyProvidedSignatureAxes.sampleAxis_origUnits = fcnNumericTargetPrecision(externallyProvidedSignatureAxes.sampleAxis_origUnits);
                      externallyProvidedSignatureAxes.norm4sampleAxis_origUnits = BM.euclidW(externallyProvidedSignatureAxes.sampleAxis_origUnits, 1, 2)';
                    externallyProvidedSignatureAxes.R4P = fcnNumericTargetPrecision(externallyProvidedSignatureAxes.R4P);
                  %Compute sample focus:
                    externallyProvidedSignatureAxes.P4R4P = pValues4Correlations(externallyProvidedSignatureAxes.R4P, externallyProvidedSignatureAxes.sampleSizes4R4P);
                    externallyProvidedSignatureAxes.signedFocusedW4P = calcSignatureFocus([],[], externallyProvidedSignatureAxes.R4P, externallyProvidedSignatureAxes.P4R4P, 2, inInfo.searchStrategy.signatureFocus, sL);
                    externallyProvidedSignatureAxes.signedFocusedW4P_withPerpendicularSpace = addFlatWeights(nan(nG,0), externallyProvidedSignatureAxes.signedFocusedW4P, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
                  %Compute twin gene axis (analogous to precomputeCandidatesInParallel):
                    externallyProvidedSignatureAxes.geneAxis_origUnits = projectW(...
                      eDState.current.L2Rs ... vectorsInSourceSpace...
                     ,externallyProvidedSignatureAxes.sampleAxis_origUnits ... axesInSourceSpace
                     ,abs(externallyProvidedSignatureAxes.signedFocusedW4P) ... weights4sourceDims
                     ,2 ... sourceSpaceMatrixDim
                    ,BM.meanW, BM.euclidW) / BM.euclidW(ones(1,size(externallyProvidedSignatureAxes.signedFocusedW4P,2)),abs(externallyProvidedSignatureAxes.signedFocusedW4P),2); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
                      externallyProvidedSignatureAxes.geneAxis_origUnits(isnan(externallyProvidedSignatureAxes.geneAxis_origUnits)) = 0;
                      if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                        externallyProvidedSignatureAxes.geneAxis_origUnits = sign(externallyProvidedSignatureAxes.geneAxis_origUnits).*dampenOutliers(abs(externallyProvidedSignatureAxes.geneAxis_origUnits), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 1); 
                      end
                    externallyProvidedSignatureAxes.norm4geneAxis_origUnits = BM.euclidW(externallyProvidedSignatureAxes.geneAxis_origUnits, 1, 1)';
                  %Compute correlations of genes in the current dataset (analogous to precomputeCandidatesInParallel):
                    [externallyProvidedSignatureAxes.R4G, externallyProvidedSignatureAxes.sampleSizes4R4G, externallyProvidedSignatureAxes.P4R4G] = uncenteredWeightedCorrelation(...
                       externallyProvidedSignatureAxes.sampleAxis_origUnits, eDState.current.L2Rs, 2 ...
                      ,abs(externallyProvidedSignatureAxes.signedFocusedW4P_withPerpendicularSpace) ...
                      ,BM.meanW ...
                      ,false ...
                    );
                      if(inInfo.searchStrategy.dampenOutliers.bEnabled(min(end,sL)))
                        externallyProvidedSignatureAxes.R4G = sign(externallyProvidedSignatureAxes.R4G).*dampenOutliers(abs(externallyProvidedSignatureAxes.R4G), inInfo.searchStrategy.dampenOutliers.maxAllowedNeighboursRatio(min(end,sL)), inInfo.searchStrategy.dampenOutliers.nTopRanksToCheck4P(min(end,sL)), 1); 
                      end
                end
              %Complete the signature by updating its foci and computing its stats:
                if(true)
                  %Signature focus:
                    externallyProvidedSignatureAxes.signedFocusedW4G = calcSignatureFocus([],[], externallyProvidedSignatureAxes.R4G, externallyProvidedSignatureAxes.P4R4G, 1, inInfo.searchStrategy.signatureFocus, sL);
                    externallyProvidedSignatureAxes.signedFocusedW4P = calcSignatureFocus([],[], externallyProvidedSignatureAxes.R4P, externallyProvidedSignatureAxes.P4R4P, 2, inInfo.searchStrategy.signatureFocus, sL);
                    %Add .full2focusWeightRatio weights from the current and the twin space to respect the information from perpendicular dimensions and equalize in case of stark size differences between the current and the twin space:
                      externallyProvidedSignatureAxes.signedFocusedW4G_withPerpendicularSpace = addFlatWeights(externallyProvidedSignatureAxes.signedFocusedW4G, externallyProvidedSignatureAxes.signedFocusedW4P, 1, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
                      externallyProvidedSignatureAxes.signedFocusedW4P_withPerpendicularSpace = addFlatWeights(externallyProvidedSignatureAxes.signedFocusedW4G, externallyProvidedSignatureAxes.signedFocusedW4P, 2, inInfo.searchStrategy.correlation.full2focusWeightRatio(min(end,sL)));
                  %Extended signature focus:
                    externallyProvidedSignatureAxes.signedExtendedW4G = calcSignatureFocus([],[], externallyProvidedSignatureAxes.R4G, externallyProvidedSignatureAxes.P4R4G, 1, inInfo.searchStrategy.extendedFocus, sL);
                    externallyProvidedSignatureAxes.signedExtendedW4P = calcSignatureFocus([],[], externallyProvidedSignatureAxes.R4P, externallyProvidedSignatureAxes.P4R4P, 2, inInfo.searchStrategy.extendedFocus, sL);
                  %Signature size estimation by correlation sums:
                    externallyProvidedSignatureAxes.signatureSizeByCorrSum4G = BM.meanW(abs(externallyProvidedSignatureAxes.signedExtendedW4G), inInfo.reference.rowW, 1, true);
                    externallyProvidedSignatureAxes.signatureSizeByCorrSum4P = BM.meanW(abs(externallyProvidedSignatureAxes.signedExtendedW4P), inInfo.reference.colW, 2, true);
                    externallyProvidedSignatureAxes.signatureSizeByCorrSum2D = externallyProvidedSignatureAxes.signatureSizeByCorrSum4G * externallyProvidedSignatureAxes.signatureSizeByCorrSum4P;
                  %combined p values for Rs:
                    [~, log10_p4Correlations4G] = pValues4Correlations(...
                      externallyProvidedSignatureAxes.R4G, externallyProvidedSignatureAxes.sampleSizes4R4G ...
                    ,1);
                    [~, log10_p4Correlations4P] = pValues4Correlations(...
                      externallyProvidedSignatureAxes.R4P, externallyProvidedSignatureAxes.sampleSizes4R4P ...
                    ,2);
                    externallyProvidedSignatureAxes.log10_p4Correlations = min(log10_p4Correlations4G, log10_p4Correlations4P); %these are two indirectly dependent p values (same L2Rs base); the min should be an upper bound of the combined p value; too conservative?!
                  %Compute signature strengths (i.e. projections on axes) for all points (these signature strenghts determine the signautre eigenOrder for regression in the default):
                    externallyProvidedSignatureAxes.signatureStrengths4G = projectW(...
                      eDState.current.L2Rs ... vectorsInSourceSpace...
                     ,externallyProvidedSignatureAxes.sampleAxis_origUnits ... axesInSourceSpace
                     ,abs(externallyProvidedSignatureAxes.signedFocusedW4P) ... weights4sourceDims
                     ,2 ... sourceSpaceMatrixDim
                    ,BM.meanW, BM.euclidW) / BM.euclidW(ones(1,nP),abs(externallyProvidedSignatureAxes.signedFocusedW4P),2); %normalize such that every compoenent is a representative expression independent of the size of the aggregated space.
                    externallyProvidedSignatureAxes.signatureStrengths4P = projectW(...
                      eDState.current.L2Rs ... vectorsInSourceSpace...
                     ,externallyProvidedSignatureAxes.geneAxis_origUnits ... axesInSourceSpace
                     ,abs(externallyProvidedSignatureAxes.signedFocusedW4G) ... weights4sourceDims
                     ,1 ... sourceSpaceMatrixDim
                    ,BM.meanW, BM.euclidW) / BM.euclidW(ones(nG,1),abs(externallyProvidedSignatureAxes.signedFocusedW4G),1); %normalize such that every compoenent is a representative expression independent of the size of the aggregated space.
                end
            %Set final signature axes to be used for regression:
              signatureDefinition.explanations.step2_finalSignatureAxes = sprintf('Externally provided signature %d/%d, applied to the current input signal.', eDState.current.k, length(previouslyDetectedSignatures));
                signatureDefinition.step2_finalSignatureAxes = externallyProvidedSignatureAxes;
          end
        %Assemble and explain outputs for the signature definition:
          %Reproducibility/store reference information in the detected signature:
            if(true)
              signatureDefinition.explanations.reference.k = 'detection rank of this signature';
                           signatureDefinition.reference.k = eDState.current.k + inInfo.applicationMode.kOffset4filenames;
              signatureDefinition.explanations.reference.nG = 'number of genes, i.e. size(L2Is,1)';
                           signatureDefinition.reference.nG = nG;
              signatureDefinition.explanations.reference.nP = 'number of samples, i.e. size(L2Is,2)';
                           signatureDefinition.reference.nP = nP;
              signatureDefinition.explanations.reference.configurationLevel = 'The configuration level, the signature was detected on.';
                           signatureDefinition.reference.configurationLevel = sL;
              %Save configured order metrics used in step 3 to the output for reference (default is using projections .signatureStrengths4G respectively .signatureStrengths4P):
                signatureDefinition.explanations.reference.metric4geneOrder = 'The field used to define the gene order of the signature; typically the gene strengths .signatureStrengths4G.';
                             signatureDefinition.reference.metric4geneOrder = inInfo.dissection.metric4geneOrder;
                signatureDefinition.explanations.reference.metric4sampleOrder = 'The field used to define the sample order of the signature; typically the sample strengths .signatureStrengths4P.';
                             signatureDefinition.reference.metric4sampleOrder = inInfo.dissection.metric4sampleOrder;
              %Backup the complete config in the output for reference:
                signatureDefinition.explanations.reference.inInfo = 'Backup of the used detection configuration that led to this signature.';
                             signatureDefinition.reference.inInfo = inInfo;
            end
          %Explain all fields of the converged/focused final signature from step 2 in the output:
            signatureDefinition.explanations.step2_finalSignatureAxes = struct();
            for fn=fieldnames(signatureDefinition.step2_finalSignatureAxes)'; fn=fn{1};
              switch(fn)
                %step1/Initial representative index: (Note: unfold all case headers to see header comments)
                  case 'imj'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Only used for signatureDefinition.step1_initialRepresentative: There, .imj identifies the particular gene (positions 1:nG) respectively sample (positions -(1:nP)). For .step2_finalSignatureAxes, look at .step2_finalSignatureAxes.corrMaximizingMembers instead.';
                  case 'i'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Only used for signatureDefinition.step1_initialRepresentative: There, .i identifies the particular gene (positions 1:nG), if the initial representative is a gene. For .step2_finalSignatureAxes, look at .step2_finalSignatureAxes.corrMaximizingMembers instead.';
                  case 'j'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Only used for signatureDefinition.step1_initialRepresentative: There, .j identifies the particularsample (positions 1:nP), if the initial representative is a sample. For .step2_finalSignatureAxes, look at .step2_finalSignatureAxes.corrMaximizingMembers instead.';
                  case 'greedyScore'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Only used for signatureDefinition.step1_initialRepresentative: There, .greedyScore determines the rank of this particular gene or sample in the presort order of all candidates for initial representatives of some signature.';
                %step2/Correlation maximizing representatives:
                  case 'corrMaximizingMembers'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Structure lising all indices and step 1 axes of genes and samples used as representatives to define the converged signature axes. (Cf. explanation for the field .imj for index definition details and to fields like .sampleAxis_origUnits for details on axes.)';

                %Primary and twin axes:
                  case 'sampleAxis_origUnits'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The converged/focused sample axis in signal units for the signature from step 2.';
                  case 'geneAxis_origUnits'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The converged/focused gene axis in signal units for the signature from step 2.';
                  case 'norm4sampleAxis_origUnits'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Euclidean norm of the convergedi/focused sample axis in signal units for the signature from step 2.';
                  case 'norm4geneAxis_origUnits'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The converged/focused gene axis in signal units for the signature from step 2.';
                  case 'signatureStrengths4P'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The projections of samples onto .geneAxis, i.e. the Euclidean distance in R^nG of the input samples parallel to the gene axis. (Refined in step 3 by projecting onto the sample-specific point on the regressed gene curve instead of the gene axis.)';
                  case 'signatureStrengths4G'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The projection of genes onto .sampleAxis, i.e. the Euclidean distance in R^nP of the input genes parallel to the sample axis. (Refined in step 3 by projecting onto the gene-specific point on the regressed sample curve instead of the sample axis.)';

                %Correlations with the signature axes and their significance:
                  case 'R4G'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The correlations of all genes with the sample axis (pre-centered weighted correlations of the log2(ratio)s).';
                  case 'R4P'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Correlations of all samples with the gene axis (pre-centered weighted correlations of the log2(ratio)s).';
                  case 'sampleSizes4R4G'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The signatureive number of samples in the signature focus whose signals were used to compute .R4G.';
                  case 'sampleSizes4R4P'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The signatureive number of genes in the signature focus whose signals were used to compute .R4P.';
                  case 'P4R4G'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The estimated p values for the gene correlations .R4G with the sample axis (estimation based on permutation tests).';
                  case 'P4R4P'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Estimated p values for the sample correlations .R4P with the gene axis (estimation based on permutation tests).';

                %Signature focus (weights):
                  case 'signedFocusedW4G'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Converged signed gene weights defining the gene focus of the signature. Based on standardized signal strengths and correlations with the signature. The signs indicate co/anti-correlation with the sample axis of the signature.';
                  case 'signedFocusedW4P'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Converged signed sample weights defining the sample focus of the signature. Based on standardized signal strengths and correlations with the signature. The signs indicate co/anti-correlation with the gene axis of the signature.';
                  case 'signedFocusedW4G_withPerpendicularSpace'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Converged signed gene weights defining the gene focus of the signature. Based on standardized signal strengths and correlations with the signature. The signs indicate co/anti-correlation with the sample axis of the signature.';
                  case 'signedFocusedW4P_withPerpendicularSpace'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Converged signed sample weights defining the sample focus of the signature. Based on standardized signal strengths and correlations with the signature. The signs indicate co/anti-correlation with the gene axis of the signature.';
                %Extended signature focus and correaltion sums as signature size estimates:
                  case 'signedExtendedW4G'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Signed extended signature focus for every gene.';
                  case 'signedExtendedW4P'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Signed extended signature focus for every sample.';
                  case 'signatureSizeByCorrSum4G'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Estimation of the signature size in gene space. Computed as sum over abs(.signedExtendedW4G).';
                  case 'signatureSizeByCorrSum4P'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Estimation of the signature size in sample space. Computed as sum mover abs(.signedExtendedW4P).';
                  case 'signatureSizeByCorrSum2D'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Estimation of the signature size in units of independent (gene, sample) measurement values. Computed as product of .signatureSizeByCorrSum4G and .signatureSizeByCorrSum4P.';
                %p values of the signature as a whole:
                  case 'log10_p'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The overall log10(p value) for the signature. This is the combination of .log10_p4Correlations .log10_p4SignalStrength by Fisher''s method.';
                  case 'log10_p4Correlations'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'Combined p values for the correlations of the genes and sampels with the signature based on a Kolmogorov Smirnov test against the t(r) distribution expected for noise-based correlations.';
                  case 'log10_p4SignalStrength'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'log_10(p value) of the two-sample t-test of the absolute signal of the weighted (gene,sample) pixels in the signature focus versus the current estimated noise distribution of absolute signal values.';

                %Signature statistics/signal and correlation averages:
                  case 'signatureAbsMean2D'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The mean absolute current unexplained signal in the signature focus in signal units.';
                  case 'sampleSize4signatureAbsMean2D'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The mass of the signature focus (i.e. the integral over all absolute weights used for computing the mean absolute signal in the signature current/remaining signal in the signature focus).';
                  case 'signatureAbsSD2D'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The standard deviation of the current/remaining absolute signal in the signature focus in signal units.';
                  case 'signatureCorrInExtendedFocus'
                    signatureDefinition.explanations.step2_finalSignatureAxes.(fn) = 'The mean absolute correlation in the signature focus.';
                otherwise
                  error('code validation: unexplained field in signature structure signatureDefinition.step2_finalSignatureAxes: %s', fn);
              end
            end
      %% STEP 3) Bimonotonic regression of all genes/arrays in the signature's eigenOrder of genes and samples to determine the signal parts contributed by the current signature (i.e. by its underlying interaction) to the signal sum:
        SDCM_printStatus(0,'### STEP 3) Bimonotonic regression of the signal in the signature''s gene and sample order to determine those signal parts contributed by the current signature (caused by the underlying interaction represented by it):\n'); drawnow;
          outInfo.reference.timing.b3_startOfStep3(end+1) = now;
        %Bimonotonic regression and dissection strengths:
          bIsEmptySignature = ...
              sum(~all(isnan(eDState.current.L2Rs(signatureDefinition.step2_finalSignatureAxes.signedFocusedW4G~=0,:)),2)) <= 2 ...
           || sum(~all(isnan(eDState.current.L2Rs(:,signatureDefinition.step2_finalSignatureAxes.signedFocusedW4P~=0)),1)) <= 2 ...
          ;
          assert(~bIsEmptySignature, 'code validation: the current signature has an empty focus and should not have been detected; check .qualification.minAbsCorrSum4* configuration.');
          [signatureDefinition.step3_regression, signatureDefinition.forPlots] = SDCM_step3_regression(signatureDefinition, inInfo);
          %Explain all fields of the regression/from step 3 in the output:
            for fn=fieldnames(signatureDefinition.step3_regression)'; fn=fn{1};
              switch(fn)
                %Signature order (along which the signature signal has been bimonotonically regressed):
                  case 'eigenSI'
                    signatureDefinition.explanations.step3_regression.(fn) = 'The sort index vector of the signature igene order (denoted I_k in the paper).';
                  case 'eigenSJ'
                    signatureDefinition.explanations.step3_regression.(fn) = 'The sort index vector of the signature sample order (denoted J_k in the paper).';
                %Order metric fields:
                  case 'signatureStrengths4P'
                    signatureDefinition.explanations.step3_regression.(fn) = 'The projections of samples onto the regressed gene curve, i.e. the Euclidean distance in R^nG of the input samples parallel to the respective vector in the gene curve.';
                  case 'signatureStrengths4G'
                    signatureDefinition.explanations.step3_regression.(fn) = 'The projection of genes onto the regressed sample curve, i.e. the Euclidean distance in R^nP of the input genes parallel to the respective vector in the sample curve.';
                %Regression result:
                  case 'signatureEigensignal'
                    signatureDefinition.explanations.step3_regression.(fn) = 'The (nonzero parts of the) bimonotonically regressed signal contributions of this signature to the overall signal sum. (Execute "signatureSignal = zeros(nG,nP); signatureSignal(.step3_regression.BNonzeroRows4signatureEigensignal, .step3_regression.BNonzeroCols4signatureEigensignal) = .step3_regression.signatureEigensignal;" to obtain it as a full-sized (nG,nP) matrix.)';
                  case 'BNonzeroRows4signatureEigensignal'
                    signatureDefinition.explanations.step3_regression.(fn) = 'Boolean vector of size nG*1 determining the nonzero rows of this signature''s eigensignal; see .signatureEigensignal.';
                  case 'BNonzeroCols4signatureEigensignal'
                    signatureDefinition.explanations.step3_regression.(fn) = 'Boolean vector of size 1*nP determining the nonzero columns of this signature''s eigensignal; see .signatureEigensignal.';
                  case 'signatureDissectionStrengthsOS'
                    signatureDefinition.explanations.step3_regression.(fn) = 'The (nonzero parts of the) dissection strengths matrix. Can be used to estimate gene/sample memberships in this signature. Divide .signatureEigensignal by .signatureDissectionStrengthsOS to obtain the bimonotonically regressed signature signal before weakining parts of uncertain membership, i.e. to obtain the maximum amount of the signal sum that is consistent with and explainable by current signature axes (irrespective of weak correlation strength and significance).';
                %Copied fields from step 2:
                  case {'R4G','R4P','P4R4G','P4R4P','signedFocusedW4G','signedFocusedW4P','signedFocusedW4G_withPerpendicularSpace','signedFocusedW4P_withPerpendicularSpace','signedExtendedW4G','signedExtendedW4P','signatureSizeByCorrSum4G','signatureSizeByCorrSum4P','log10_p4Correlations','signatureSizeByCorrSum2D'}
                    signatureDefinition.explanations.step3_regression.(fn) = [signatureDefinition.explanations.step2_finalSignatureAxes.(fn), ' (Optionally updated for each regressed signature curve; otherwise idential to step 2.)'];
                otherwise
                  error('code validation: unexplained field in regression structure signatureDefinition.step3_regression: %s', fn);
              end
            end
        %Get the signatures's eigensignal:
          eigenSignal = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision);
          eigenSignal(...
            signatureDefinition.step3_regression.BNonzeroRows4signatureEigensignal...
           ,signatureDefinition.step3_regression.BNonzeroCols4signatureEigensignal...
          ) = signatureDefinition.step3_regression.signatureEigensignal;
          eDState.current.shifts = -eigenSignal; %code backwards compatibility: we "shift" the negative eigensignal of the signature to dissect the signature from the signal.
        %Estimate noise distribution based on signatureStrengths for the nonSignaturePoints in order to find the stop condition wrt. signature signal strengths:
          if(true)
            %Initialize:
              %for k==1, initialize the noise distribution storage fields:
                if(bInitializeNoiseEstimation) %reset/overwrite temporary overestimating distribution
                  eDState.noiseEstimation.nullL2Rs = nan(nG,nP);
                  eDState.noiseEstimation.noiseW2D = zeros(nG,nP);
                  eDState.noiseEstimation.sumDissectionStrength2DoverAllSignatures = zeros(nG,nP); %needed for correction (using pixels from which eigensignals were dissected in an uncorrected way is a bias to too low noise; every projection/dissection counts, even if it dissects a small zero-near offset, a dimension is lost)

                  %just needed for k==1:
                    eDState.noiseEstimation.absW4G4K = nan(nG,0);
                    eDState.noiseEstimation.absW4P4K = nan(0,nP);
                    eDState.noiseEstimation.pixelSDViaRadialNorms4G4K = nan(nG,0);
                    eDState.noiseEstimation.pixelSDViaRadialNorms4P4K = nan(0,nP);
                end
              %weights: reuse dissection strengths for noise estimation (i.e. for determining in-signature and out-of-signature points below):
                %W2D = abs(bsxfun(@times, signatureDefinition.step2_finalSignatureAxes.signedFocusedW4G, signatureDefinition.step2_finalSignatureAxes.signedFocusedW4P));
                dissectionStrengths2DOS = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision);
                  dissectionStrengths2DOS(...
                    signatureDefinition.step3_regression.BNonzeroRows4signatureEigensignal...
                   ,signatureDefinition.step3_regression.BNonzeroCols4signatureEigensignal...
                  ) = signatureDefinition.step3_regression.signatureDissectionStrengthsOS; %we shift the negative eigensignal of the signature to dissect the signature from the signal.
                W2D = dissectionStrengths2DOS; 
            %Estimate SD around the signature axis to exclude points having strong signals from perpendicular signatures:
              %Correction factors due to dimension loss by previous dissections (via E[chi(k)] = E[euclid(k dims)] and expectation ratio after loosing dims due to projections)
                if(true)
                  soFarProjectedGeneDims4P = max(eDState.noiseEstimation.sumDissectionStrength2DoverAllSignatures,[],1); %max since we are interested in the number of projections (a single projection may affect/dissect singals in many genes, but can only ever lose one space dimension)
                  soFarProjectedSampleDims4G = max(eDState.noiseEstimation.sumDissectionStrength2DoverAllSignatures,[],2); %max since we are interested in the number of projections (a single projection may affect/dissect singals in many genes, but can only ever lose one space dimension)
                    euclidNoiseNormUnderestimationFactorsDueToPreviousDissections4P = ...
                        (sqrt(2)*exp(gammaln((max(0,nG-soFarProjectedGeneDims4P)+1)/2)-gammaln(max(0,nG-soFarProjectedGeneDims4P)/2))) ...
                      / (sqrt(2)*exp(gammaln((nG+1)/2)-gammaln(nG/2))) ...
                    ;
                    euclidNoiseNormUnderestimationFactorsDueToPreviousDissections4G = ...
                        (sqrt(2)*exp(gammaln((max(0,nP-soFarProjectedSampleDims4G)+1)/2)-gammaln(max(0,nP-soFarProjectedSampleDims4G)/2))) ...
                      / (sqrt(2)*exp(gammaln((nP+1)/2)-gammaln(nP/2))) ...
                    ;
                    %<-A point should not be affected by so many projections/dissections that it has already lost >80% of its norm:
                      minRemainingNormRatioAfterProjections = 0.2;
                      euclidNoiseNormUnderestimationFactorsDueToPreviousDissections4G = max(minRemainingNormRatioAfterProjections, euclidNoiseNormUnderestimationFactorsDueToPreviousDissections4G);
                      euclidNoiseNormUnderestimationFactorsDueToPreviousDissections4P = max(minRemainingNormRatioAfterProjections, euclidNoiseNormUnderestimationFactorsDueToPreviousDissections4P);

                  eDState.noiseEstimation.sumDissectionStrength2DoverAllSignatures = eDState.noiseEstimation.sumDissectionStrength2DoverAllSignatures + dissectionStrengths2DOS;

                  soFarProjectedGeneDims4P = max(eDState.noiseEstimation.sumDissectionStrength2DoverAllSignatures,[],1); %max since we are interested in the number of projections (a single projection may affect/dissect singals in many genes, but can only ever lose one space dimension)
                  soFarProjectedSampleDims4G = max(eDState.noiseEstimation.sumDissectionStrength2DoverAllSignatures,[],2); %max since we are interested in the number of projections (a single projection may affect/dissect singals in many genes, but can only ever lose one space dimension)
                    underestimationFactorsDueToPreviousPlusCurrentDissections4P = ...
                        (sqrt(2)*exp(gammaln((max(0,nG-soFarProjectedGeneDims4P)+1)/2)-gammaln(max(0,nG-soFarProjectedGeneDims4P)/2))) ...
                      / (sqrt(2)*exp(gammaln((nG+1)/2)-gammaln(nG/2))) ...
                    ;
                    underestimationFactorsDueToPreviousPlusCurrentDissections4G = ...
                        (sqrt(2)*exp(gammaln((max(0,nP-soFarProjectedSampleDims4G)+1)/2)-gammaln(max(0,nP-soFarProjectedSampleDims4G)/2))) ...
                      / (sqrt(2)*exp(gammaln((nP+1)/2)-gammaln(nP/2))) ...
                    ;
                    %<-A point should not be affected by so many projections/dissections that it has already lost >80% of its norm:
                      BSufficientRemainingNormRatio4G = underestimationFactorsDueToPreviousPlusCurrentDissections4G > minRemainingNormRatioAfterProjections;
                      BSufficientRemainingNormRatio4P = underestimationFactorsDueToPreviousPlusCurrentDissections4P > minRemainingNormRatioAfterProjections;
                      
                  underestimationFactorsDueToPreviousPlusCurrentDissections4G = max(minRemainingNormRatioAfterProjections, underestimationFactorsDueToPreviousPlusCurrentDissections4G);
                  underestimationFactorsDueToPreviousPlusCurrentDissections4P = max(minRemainingNormRatioAfterProjections, underestimationFactorsDueToPreviousPlusCurrentDissections4P);
                end
              %estimate BCertainNonNoiseInPerpendicularDirections2D based on pixelSD:
                if(true)
                  %Comptue radial norms after dissection:
                    radialNorms4G = sqrt(BM.meanW((eDState.current.L2Rs + eDState.current.shifts).^2,1,2,true));
                    radialNorms4P = sqrt(BM.meanW((eDState.current.L2Rs + eDState.current.shifts).^2,1,1,true));
                    %Normalize per pixel: Berechne aus Euklidnormverteilung die mittlere SD pro Dimension:
                      if(true)
                        pixelSDViaRadialNorms4G = radialNorms4G / (sqrt(2)*exp(gammaln((nP+1)/2)-gammaln(nP/2)));
                        pixelSDViaRadialNorms4P = radialNorms4P / (sqrt(2)*exp(gammaln((nG+1)/2)-gammaln(nG/2)));
                      end
                    %Correct for dimension loss:
                      pixelSDViaRadialNorms4G = pixelSDViaRadialNorms4G ./ underestimationFactorsDueToPreviousPlusCurrentDissections4G;
                      pixelSDViaRadialNorms4P = pixelSDViaRadialNorms4P ./ underestimationFactorsDueToPreviousPlusCurrentDissections4P;
                  %Obtain pixelSDInSignaturePerpendicularDirections:
                    if(bInitializeNoiseEstimation) %in this case no .preciseNoiseSD is yet available, but we already need an estimate for the perpendicular extents of the signature around the signature curve; use the distribution of pixelSDViaRadialNorms4? to get a start:
                      %Select and save pixelSDsViaRadialNorms for current TOP signature members:
                        if(true)
                          %Selection criteria:
                            %The point should be a top member in the signature, because else it might be affected by perpendicular signatures and overestimate noise (implicit assumption: top members of signatures are not also superposed):
                              minSignatureWeight = 0.8; %0.9; %0.8; %0.95; %require certain signature membership in order to exclude points that are not pure members, but also in some perpendicular signatures. (needed for versa)
                              eDState.noiseEstimation.absW4G4K(:,eDState.current.k) = max(0, ...
                                ...double(BSufficientRemainingNormRatio4G).*abs(signatureDefinition.step2_finalSignatureAxes.signedFocusedW4G)-minSignatureWeight...
                                double(BSufficientRemainingNormRatio4G).*abs(signatureDefinition.step3_regression.signedFocusedW4G)-minSignatureWeight... %.step3_regression.signedFocusedW4G is identical to step 2 except if a backwards compatibility option is enabled.
                              )/(1-minSignatureWeight);
                              eDState.noiseEstimation.absW4P4K(eDState.current.k,:) = max(0, ...
                                ...double(BSufficientRemainingNormRatio4P).*abs(signatureDefinition.step2_finalSignatureAxes.signedFocusedW4P)-minSignatureWeight...
                                double(BSufficientRemainingNormRatio4P).*abs(signatureDefinition.step3_regression.signedFocusedW4P)-minSignatureWeight... %.step3_regression.signedFocusedW4P is identical to step 2 except if a backwards compatibility option is enabled.
                              )/(1-minSignatureWeight);
                          %Save selected:
                            eDState.noiseEstimation.pixelSDViaRadialNorms4G4K(:,eDState.current.k) = NaN;
                              eDState.noiseEstimation.pixelSDViaRadialNorms4G4K(eDState.noiseEstimation.absW4G4K(:,eDState.current.k)>0, eDState.current.k) = ...
                                pixelSDViaRadialNorms4G(eDState.noiseEstimation.absW4G4K(:,eDState.current.k)>0);
                            eDState.noiseEstimation.pixelSDViaRadialNorms4P4K(eDState.current.k,:) = NaN;
                              eDState.noiseEstimation.pixelSDViaRadialNorms4P4K(eDState.noiseEstimation.absW4P4K(eDState.current.k,:)>0) = ...
                                pixelSDViaRadialNorms4P(eDState.noiseEstimation.absW4P4K(eDState.current.k,:)>0);
                        end
                      %Estimate pixelSD:
                        if(true)
                          flatten = @(C)ilv(cellfun(@(c)c(:),C,'UniformOutput',false),@(C2)vertcat(C2{:}));
                          %For robustness against some outliers (points with perpendicular signatures overestimating noise), use weighted median instead of mean:
                            X = flatten({
                              eDState.noiseEstimation.pixelSDViaRadialNorms4G4K(:,1:eDState.current.k);
                              eDState.noiseEstimation.pixelSDViaRadialNorms4P4K(1:eDState.current.k,:);
                            });
                            W = flatten({ %Assumption: the top of the signature points are not affected by strong perpendicular signatures and thus provide a useful noise estimation base after the current dissection and correction for dimension loss:
                              eDState.noiseEstimation.absW4G4K(:,1:eDState.current.k)/nP; 
                              eDState.noiseEstimation.absW4P4K(1:eDState.current.k,:)/nG; 
                            });
                            [SX,CDF] = ecdfW(X,W);
                            [~,i4Median] = min(abs(CDF-0.5));
                            pixelSDInSignaturePerpendicularDirections = SX(i4Median);
                        end
                    else %we already have a noise distribution built from previous detection iterations; use ir instead of just the radial distribution around the current signature axis:
                      pixelSDInSignaturePerpendicularDirections = eDState.noiseEstimation.preciseNoiseSD(eDState.current.k);
                    end
                  %for sample noise distribution below, define BCertainNonNoiseInPerpendicularDirections2D via cutoff:
                    if(true)
                      %use estiamted pixelSDInSignaturePerpendicularDirections to center the pixelSDViaRadialNorms4? around 1:
                        ZPerpendicular4P = pixelSDViaRadialNorms4P/pixelSDInSignaturePerpendicularDirections; 
                        ZPerpendicular4G = pixelSDViaRadialNorms4G/pixelSDInSignaturePerpendicularDirections;
                      %Compute 2D scores of perpendicularity:
                        Z2D = bsxfun(@plus, ZPerpendicular4G/nP, ZPerpendicular4P/nG)/(1/nP+1/nG);
                      %consider all points free of perpendicular signatures that are around Z2D=1+-3*SD(Z2D):
                        %Find the first/maximal Gaussian mode:
                          %estimate the mode by the maximum of the empirical PDF:
                            [PDF,XI] = ksdensity(min(10,Z2D(:)));
                              %figure; plot(XI,PDF);
                            [~,maxi] = max(PDF);
                            mu = XI(maxi);
                          sdZ2DaroundMu = sqrt(BM.meanW((Z2D(Z2D<=mu)-mu).^2,1,1))/sqrt(1-2/pi);  %exclude points >mode for a more robust SD estimation.
                          maxAllowedRadialSDDistribSDs = 3; 
                          BCertainNonNoiseInPerpendicularDirections2D = (Z2D-mu)/sdZ2DaroundMu > maxAllowedRadialSDDistribSDs;
                    end
                end
            %Sample the noise distribution from perspective of the current signature:
              %select points as noise that are not signifianct in the current signature direction and have less than 3 perpendicular SDs from the signature top signature points around the signature axis (in 3D: cylinder cut around the signature axis):
                if(true)
                  noiseW2D4outOfSignaturePointsb4Dissection = 1-W2D;
                  %only include certain non-signature points that were not neuralized to not overestimate the SD:
                    %minNonSignatureWeight = 0.95; %should be high enough to only select non-dissectd points; mind 0.95, weil sonst Punkte in Effektrichtung mitselektiert werden und deren ParallelDist die SD hochtreibt.
                    %noiseW2D(noiseW2D<minNonSignatureWeight^2) = 0;
                    maxDissectionStrength = 1e-3;
                    noiseW2D4outOfSignaturePointsb4Dissection = noiseW2D4outOfSignaturePointsb4Dissection .* (dissectionStrengths2DOS <= maxDissectionStrength);
                  %exclude points in perpendicular signatures:
                    %noiseW2D(bsxfun(@or, BCertainNonNoiseInPerpendicularDirections4G, BCertainNonNoiseInPerpendicularDirections4P)) = 0;
                    noiseW2D4outOfSignaturePointsb4Dissection(BCertainNonNoiseInPerpendicularDirections2D) = 0;
                  %Get L2Rs (projection norm loss corrected for all previous dissections) before dissection for the selected out-of-signature non-perpendicular points and push them into the noise distrib:
                    %Correct L2Rs for norm loss by previous dissections/projections:
                      correctedL2Rs4outOfSignaturePointsb4Dissection = eDState.current.L2Rs ...
                        ./ bsxfun(@times, euclidNoiseNormUnderestimationFactorsDueToPreviousDissections4G, euclidNoiseNormUnderestimationFactorsDueToPreviousDissections4P) ...
                      ;
                end
              %points significantly in the signature and significantly out of perpendicular signatures are perfect noise estimates after dissection, too:
                if(true)
                  noiseW2DafterDissection = W2D;
                  %exclude strong perpendicular points as before (cylinder cut around the signature axis)
                    noiseW2DafterDissection(BCertainNonNoiseInPerpendicularDirections2D) = 0;
                  %only include certain signature-points that were fully dissectd:
                    %minSignatureWeight = 0.5; %should be high enough to only select fully dissectd points; must not be too high to prevent dim subselecting below the dim that is corrected.
                    %noiseW2DafterDissection(noiseW2DafterDissection<minSignatureWeight^2) = 0;
                    minDissectionStrength = 1-1e-3;
                    noiseW2DafterDissection = noiseW2DafterDissection .* (dissectionStrengths2DOS >= minDissectionStrength);
%                   %decrease weight for pixels with dissections:
%                     noiseW2DafterDissection = noiseW2DafterDissection ./ (1+eDState.noiseEstimation.sumDissectionStrength2DoverAllSignatures);
                  %Get (projection norm loss corrected) L2Rs after dissection for the selected in-signature non-perpendicular points and push them into the noise distrib:
                    correctedL2Rs4inOfSignaturePointsAfterDissection = (eDState.current.L2Rs + eDState.current.shifts) ...
                      ./ bsxfun(@times, underestimationFactorsDueToPreviousPlusCurrentDissections4G, underestimationFactorsDueToPreviousPlusCurrentDissections4P) ...
                    ;
                end
              %update the noise distribution by adding new points, weighted relative to already recorded values for the same pixels:
                if(true)
                  %BNewNoise4outOfSignaturePointsb4Dissection = noiseW2D4outOfSignaturePointsb4Dissection>0 & ~isnan(correctedL2Rs4outOfSignaturePointsb4Dissection);
                  BNewNoise4outOfSignaturePointsb4Dissection = noiseW2D4outOfSignaturePointsb4Dissection>0.5 & ~isnan(correctedL2Rs4outOfSignaturePointsb4Dissection);
                  %<-Failsave: prevent stark overestimation from allowing too much perpendicular signature distance (mode of Z2D not estimated robustly enough):
                    %localSD = sqrt(weightedMean_supportingNaNs((correctedL2Rs4outOfSignaturePointsb4Dissection(BNewNoise4outOfSignaturePointsb4Dissection)-0).^2, noiseW2D4outOfSignaturePointsb4Dissection(BNewNoise4outOfSignaturePointsb4Dissection), 1));
                    localSD = sqrt(weightedMean_supportsNaNs(...  %.nullL2Rs are NaN-initialized.
                      (correctedL2Rs4outOfSignaturePointsb4Dissection(BNewNoise4outOfSignaturePointsb4Dissection)-0).^2 ...
                     ,noiseW2D4outOfSignaturePointsb4Dissection(BNewNoise4outOfSignaturePointsb4Dissection) ...
                    ,1));
                    if(localSD/pixelSDInSignaturePerpendicularDirections > 2); %1.5)
                      SDCM_printStatus(3,'   - Noise estimation: Not using points that are out of the current signature focus for noise estimation, as they seem to have a non-noise signal from other yet undetected signatures.');
                      BNewNoise4outOfSignaturePointsb4Dissection(:) = false;
                    end
                  %BNewNoise4inSignaturePointsAfterDissection = noiseW2DafterDissection>0 & ~isnan(correctedL2Rs4inOfSignaturePointsAfterDissection);
                  BNewNoise4inSignaturePointsAfterDissection = noiseW2DafterDissection>0.5 & ~isnan(correctedL2Rs4inOfSignaturePointsAfterDissection); %nehme nur Tophaelfte des Effekts.
                    if(true)
                      BOldNoiseEstiamtesAvailable = ~isnan(eDState.noiseEstimation.nullL2Rs);
                      eDState.noiseEstimation.nullL2Rs = nansum(cat(3 ...
                        ,BNewNoise4outOfSignaturePointsb4Dissection .* noiseW2D4outOfSignaturePointsb4Dissection .* abs(correctedL2Rs4outOfSignaturePointsb4Dissection) ...
                        ,BNewNoise4inSignaturePointsAfterDissection .* noiseW2DafterDissection .* abs(correctedL2Rs4inOfSignaturePointsAfterDissection) ...
                        ,BOldNoiseEstiamtesAvailable .* eDState.noiseEstimation.noiseW2D .* abs(eDState.noiseEstimation.nullL2Rs) ...
                      ),3) ./ nansum(cat(3 ...
                        ,BNewNoise4outOfSignaturePointsb4Dissection .* noiseW2D4outOfSignaturePointsb4Dissection ...
                        ,BNewNoise4inSignaturePointsAfterDissection .* noiseW2DafterDissection ...
                        ,BOldNoiseEstiamtesAvailable .* eDState.noiseEstimation.noiseW2D ...
                      ),3);
                      eDState.noiseEstimation.noiseW2D = nanmax(cat(3 ... %no not sum weights as this overestimates the true mass for the ttests.
                        ,BNewNoise4outOfSignaturePointsb4Dissection .* noiseW2D4outOfSignaturePointsb4Dissection ... %use the maximum noise probability from perspective of all detected signatures, but never give one pixel a mass >1 (can lead to underestimated p values!)
                        ,BNewNoise4inSignaturePointsAfterDissection .* noiseW2DafterDissection ...
                        ,BOldNoiseEstiamtesAvailable .* eDState.noiseEstimation.noiseW2D  ...
                      ),[],3);
                    end
                end
              %update .preciseNoiseSD (e.g. for computing SNRs and for KS reference nullCDF and threshold):
                eDState.noiseEstimation.preciseNoiseSD(eDState.current.k+1) = sqrt(weightedMean_supportsNaNs(... %.nullL2Rs are NaN-initialized.
                   (eDState.noiseEstimation.nullL2Rs(:)-0).^2 ...
                  ,eDState.noiseEstimation.noiseW2D(:) ...
                ,1));
                eDState.noiseEstimation.halfNormalAbsNoiseMedian(eDState.current.k+1) = eDState.noiseEstimation.preciseNoiseSD(eDState.current.k+1) * sqrt(2)*erfinv(1/2);
              %update noise distribution based reference statistics required for p values for signal strengths during qualification in step 1:
                eDState.noiseEstimation.absMean2D(eDState.current.k+1) = weightedMean_supportsNaNs(... %.nullL2Rs are NaN-initialized.
                   abs(eDState.noiseEstimation.nullL2Rs(:)) ...
                  ,eDState.noiseEstimation.noiseW2D(:) ...
                );
                eDState.noiseEstimation.nsBehindAbsMean2D(eDState.current.k+1) = BM.meanW(eDState.noiseEstimation.noiseW2D(:),~isnan(eDState.noiseEstimation.nullL2Rs(:)),1,true);
                eDState.noiseEstimation.SD4absMean2D(eDState.current.k+1) = sqrt(weightedMean_supportsNaNs(... %this is the SD relative to the abs(signal) niveau of the signature; needed later for signal strength p values.
                   (abs(eDState.noiseEstimation.nullL2Rs(:))-eDState.noiseEstimation.absMean2D(eDState.current.k+1)).^2 ... %sqrt(E[(X-mu)^2])
                  ,eDState.noiseEstimation.noiseW2D(:)...
                ));
            %Compute signal noise probabilities after dissection of the current detected signature for all genes and samples (used by step1 of the next detection iteration, e.g. for qualification):
              if(true)
                %fast approximation via unweighted chi^2 test:
                  %math notes:
                    % - die Radialnormen sqrt(sum(x_i^2)) folgen einer chi-Verteilung mit nG bzw. nP Freiheitsgraden
                    % - die Radialvarianzen sum(x_i^2) folgen einer chi^2-Verteilung mit nG bzw. nP Freiheitsgraden
                  Z = (eDState.current.L2Rs + eDState.current.shifts - 0)/eDState.noiseEstimation.preciseNoiseSD(eDState.current.k+1);
                  chiSquared4G = BM.meanW(Z.^2,1,2,true);
                  nu4G = BM.meanW(~isnan(Z),1,2,true); %make ready for missing values.
                    eDState.noiseEstimation.log10PNoise4G = exactUpperChi2CDF(chiSquared4G,nu4G)/log(10);
                  chiSquared4P = BM.meanW(Z.^2,1,1,true);
                  nu4P = BM.meanW(~isnan(Z),1,1,true); %make ready for missing values.
                    eDState.noiseEstimation.log10PNoise4P = exactUpperChi2CDF(chiSquared4P,nu4P)/log(10);
                    %figure;scatter(log10(eDState.noiseEstimation.log10PNoise4P),exactUpperChi2CDF(chiSquared4P,nG)/log(10));
              end
          end
        %Compute signature statistics: (not yet done for step 2 axes, as not needed for regression; compute here after the noise estimation update and not in step 2 to output a p value estimate for the signal strength even for k==1)
          [ signatureDefinition.step3_regression.signatureAbsMean2D, signatureDefinition.step3_regression.sampleSize4signatureAbsMean2D, signatureDefinition.step3_regression.signatureAbsSD2D, signatureDefinition.step3_regression.log10_p4SignalStrength ...
           ,signatureDefinition.step3_regression.signatureCorrInExtendedFocus, signatureDefinition.step3_regression.log10_p ...
          ] = computeSignatureStatistics(...
             eDState.current.L2Rs ...
            ,signatureDefinition.step3_regression.signedExtendedW4G, signatureDefinition.step3_regression.signedExtendedW4P ...
            ,signatureDefinition.step3_regression.R4G, signatureDefinition.step3_regression.R4P ...
            ,signatureDefinition.step3_regression.signatureSizeByCorrSum4G, signatureDefinition.step3_regression.signatureSizeByCorrSum4P ...
            ,signatureDefinition.step3_regression.log10_p4Correlations ...
            ,struct(... %use the current estimation for the noise distribution (from after the eDState.current.k's dissection; this also allows a retrospect estimate for the signal p value for signature k==1).
                'halfNormalAbsNoiseMedian', eDState.noiseEstimation.halfNormalAbsNoiseMedian(eDState.current.k+1) ...
               ,'preciseNoiseSD', eDState.noiseEstimation.preciseNoiseSD(eDState.current.k+1) ...
               ,'absMean2D', eDState.noiseEstimation.absMean2D(eDState.current.k+1) ...
               ,'nsBehindAbsMean2D', eDState.noiseEstimation.nsBehindAbsMean2D(eDState.current.k+1) ...
               ,'SD4absMean2D', eDState.noiseEstimation.SD4absMean2D(eDState.current.k+1) ...
             ) ... 
            ,BM ...
          );
        %Performance/invalidate candidate initial representatives in the cache for STEP1 that were affected by the current effect:
          if(~inInfo.applicationMode.bEnabled && inInfo.searchStrategy.performance.bCachePotentialInitialRepresentatives(min(end,sL)))
            eDState.current.precompute = eDState.current.precompute_previousIterations; %restore saved after STEP1.
            %Estimate items to remove from dissection strengths:
              if(true)
                signatureDissectionStrengthsOS = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision);
                  signatureDissectionStrengthsOS(...
                    signatureDefinition.step3_regression.BNonzeroRows4signatureEigensignal ...
                   ,signatureDefinition.step3_regression.BNonzeroCols4signatureEigensignal ...
                  ) = signatureDefinition.step3_regression.signatureDissectionStrengthsOS; %we shift the negative eigensignal of the effect to neutralize the effect from the signal.
                affectedI = +find(max(signatureDissectionStrengthsOS,[],2) > inInfo.searchStrategy.performance.maxAllowedDissectionStrengthToKeepInCache(min(end,sL)));
                affectedJ = -find(max(signatureDissectionStrengthsOS,[],1) > inInfo.searchStrategy.performance.maxAllowedDissectionStrengthToKeepInCache(min(end,sL)));
                ImJToRemove = intersect(eDState.current.precompute.ImJ, [affectedI',affectedJ]);
              end
              SDCM_printStatus(2,'   - invalidate %d/%d cached initial representative candidates, since their signal was changed by the current detected effect.\n', length(ImJToRemove), length(eDState.current.precompute.ImJ));
              precomputeCandidatesInParallel(ImJToRemove, 'clearFromCache');
            eDState.current.precompute_previousIterations = eDState.current.precompute; %write in saved slot for next Step1 to restore from.
          end
          %<-dev.note: das Problem ist, dass Fremdneutralisierungen den Effektfokus des Kandidaten beenflussen können, selbst wenn dessen Signal unberührt bleibt => keine 100%ige resume-Fähigkeit, wenn dieser cache aktiviert ist!
        %Quality control: compute signature strengths and correlations in the signature focus with the signal AFTER dissection for dissection control:
          L2RsAfterDissection = eDState.current.L2Rs + eDState.current.shifts;
          %Compute remaining signature strengths (as in step 3, i.e. as projections on the gene respectively sample curves):
            if(true)
              signatureEigensignalOSOrigUnits = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision);
                signatureEigensignalOSOrigUnits(...
                  signatureDefinition.step3_regression.BNonzeroRows4signatureEigensignal...
                 ,signatureDefinition.step3_regression.BNonzeroCols4signatureEigensignal...
                ) = signatureDefinition.step3_regression.signatureEigensignal;

              signatureDefinition.explanations.afterDissection.signatureStrengths4G = 'The projection of genes AFTER dissection onto .sampleAxis, i.e. the Euclidean distance in R^nP of the input genes parallel to the sample axis. (Refined in step 3 by projecting onto the gene-specific point on the regressed sample curve instead of the sample axis.) Can be used to assess the dissection effectiveness.';
                            signatureDefinition.afterDissection.signatureStrengths4G = projectW(...
                                L2RsAfterDissection ... vectorsInSourceSpace...
                               ,bsxfun(@times, sign(signatureDefinition.step3_regression.R4G), signatureEigensignalOSOrigUnits) ... axesInSourceSpace %one sample axis per gene.
                               ,abs(signatureDefinition.step3_regression.signedFocusedW4P) ... weights4sourceDims
                               ,2 ... sourceSpaceMatrixDim
                            ,BM.meanW, BM.euclidW) / BM.euclidW(ones(1,nP),abs(signatureDefinition.step3_regression.signedFocusedW4P),2); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.
              signatureDefinition.explanations.afterDissection.signatureStrengths4P = 'The projections of samples AFTER dissection onto .geneAxis, i.e. the Euclidean distance in R^nG of the input samples parallel to the gene axis. (Refined in step 3 by projecting onto the sample-specific point on the regressed gene curve instead of the gene axis.) Can be used to assess the dissection effectiveness.';
                            signatureDefinition.afterDissection.signatureStrengths4P = projectW(...
                                L2RsAfterDissection ... vectorsInSourceSpace...
                               ,bsxfun(@times, sign(signatureDefinition.step3_regression.R4P), signatureEigensignalOSOrigUnits) ... axesInSourceSpace %one gene axis per sample.
                               ,abs(signatureDefinition.step3_regression.signedFocusedW4G) ... weights4sourceDims
                               ,1 ... sourceSpaceMatrixDim
                            ,BM.meanW, BM.euclidW) / BM.euclidW(ones(nG,1),abs(signatureDefinition.step3_regression.signedFocusedW4G),1); %normalize like for signatureStrengths to become indepent of the projection target space and get a unit "L2Rs per pixel" to be comparable.

              if(inInfo.internal.bDevPlots)
                figure('Position',[679, 403, 1243, 577]); 

                subplot(1,2,1);
                  plot(signatureDefinition.step2_finalSignatureAxes.signatureStrengths4P, signatureDefinition.afterDissection.signatureStrengths4P,'.');
                  xlabel('signatureStrengths4P before dissection');
                  ylabel('signatureStrengths4P after dissection');
                  title('dissection effectiveness / gene prespective');
                  daspect([1 1 1]);
                  box on; grid on;
                subplot(1,2,2);
                  plot(signatureDefinition.step2_finalSignatureAxes.signatureStrengths4G, signatureDefinition.afterDissection.signatureStrengths4G,'.');
                  xlabel('signatureStrengths4G before dissection');
                  ylabel('signatureStrengths4G after dissection');
                  title('dissection effectiveness / sample prespective');
                  daspect([1 1 1]);
                  box on; grid on;
                drawnow;
              end
            end
          %Compute remaining correlations with signature axes based on the L2Rs after dissection (will be part of the exported definition table):
            if(true)
              signatureDefinition.explanations.afterDissection.R4G = 'The correlations of all genes'' log2(ratio)s AFTER dissection with the sample axis (computed analogously to .step2_finalSignatureAxes.R4G). Can be used to assess the dissection effectiveness.';
                          [signatureDefinition.afterDissection.R4G, sampleSizes4R4G] = uncenteredWeightedCorrelation(...
                             signatureDefinition.step2_finalSignatureAxes.sampleAxis_origUnits, L2RsAfterDissection, 2 ...
                            ,abs(signatureDefinition.step3_regression.signedFocusedW4P_withPerpendicularSpace) ...
                            ,BM.meanW ...
                          );
                            [signatureDefinition.afterDissection.P4R4G, log10_p4Correlations4G] = pValues4Correlations(...
                              signatureDefinition.afterDissection.R4G, sampleSizes4R4G ...
                            ,1);
              signatureDefinition.explanations.afterDissection.R4P = 'The correlations of all samples'' log2(ratio)s AFTER dissection with the gene axis (computed analogously to .step2_finalSignatureAxes.R4P). Can be used to assess the dissection effectiveness.';
                          [signatureDefinition.afterDissection.R4P, sampleSizes4R4P] = uncenteredWeightedCorrelation(...
                             signatureDefinition.step2_finalSignatureAxes.geneAxis_origUnits, L2RsAfterDissection, 1 ...
                            ,abs(signatureDefinition.step3_regression.signedFocusedW4G_withPerpendicularSpace) ...
                            ,BM.meanW ...
                          );
                            [signatureDefinition.afterDissection.P4R4P, log10_p4Correlations4P] = pValues4Correlations(...
                              signatureDefinition.afterDissection.R4P, sampleSizes4R4P ...
                            ,2); 
                            signatureDefinition.afterDissection.log10_p4Correlations = min(log10_p4Correlations4G, log10_p4Correlations4P); %these are two indirectly dependent p values (same L2Rs base); the min should be an upper bound of the combined p value; too conservative?!
            end
          %Signal strength statistics:
            if(true)
              %Exclude the user-specified number of initial/lab signatures when calculating the signal-to-explainedSignal-ratio or the signal-to-original ratio, if any:
                originalSignalWithoutLabSignatures = eDState.initialSignal.L2Rs;
                explainedSignalWithoutLabSignatures = eDState.current.explainedSignal;
                if(inInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot > 0)
                  for k=1:min(inInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot,eDState.current.k-1)
                    shifts = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision);
                    shifts(...
                      eDState.signatures(k).step3_regression.BNonzeroRows4signatureEigensignal...
                     ,eDState.signatures(k).step3_regression.BNonzeroCols4signatureEigensignal...
                    ) = -eDState.signatures(k).step3_regression.signatureEigensignal; %we shift the negative eigensignal of the signature to dissect the signature from the signal.

                    originalSignalWithoutLabSignatures = eDState.initialSignal.L2Rs + shifts;
                    explainedSignalWithoutLabSignatures = explainedSignalWithoutLabSignatures - shifts;
                  end
                end
              %Average strength of the signal that is consistent with the signature's correlations:
                if(true)
                  W2D = sqrt(abs(bsxfun(@times, signatureDefinition.step3_regression.signedExtendedW4G, signatureDefinition.step3_regression.signedExtendedW4P)));
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsConsistentOriginalL2Rs = 'Representative original log2(ratio) strength for the current signature (mean of the original log2(ratio)s, using signed weights to focus on the pixels affected by the current signature and add them constructively as long as they are consistent with the signature''s correlation signature)';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsConsistentOriginalL2Rs = BM.meanW(...
                                 abs(originalSignalWithoutLabSignatures(:)) ... %representative L2R that was originally measured for the pixels carrying the signature
                               ,W2D(:));
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsConsistentAlreadyExplainedL2Rs = 'Representative already explained log2(ratio) strength in the current signature focus.';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsConsistentAlreadyExplainedL2Rs = BM.meanW(...
                                 abs(explainedSignalWithoutLabSignatures(:)) ... %representative already explained L2R for the pixels carrying the signature
                               ,W2D(:));
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsConsistentCurrentUnexplainedL2Rs = 'Representative yet unexplained log2(ratio) strength in the current signature focus.';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsConsistentCurrentUnexplainedL2Rs = BM.meanW(...
                                 abs(eDState.current.L2Rs(:)) ... %representative unexplained L2R
                               ,W2D(:));
                               if(signatureDefinition.statistics.signalStrength.focusedMeanAbsConsistentCurrentUnexplainedL2Rs ~= signatureDefinition.step3_regression.signatureAbsMean2D)
                                 warning('code validation: signatureDefinition.statistics.signalStrength.focusedMeanAbsConsistentCurrentUnexplainedL2Rs ~= signatureDefinition.step3_regression.signatureAbsMean2D');
                               end
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsConsistentExplainedByRegressedL2Rs = 'Representative dissected log2(ratio) strength in the current signature focus.';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsConsistentExplainedByRegressedL2Rs = BM.meanW(...
                                 abs(-eDState.current.shifts(:)) ... %representative L2R that gets dissectd
                               ,W2D(:));
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsConsistentRemainingUnexplainedL2Rs = 'Representative remaining unexplained log2(ratio) strength in the current signature focus, i.e. the signal parts that are still superposed to the current signature.';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsConsistentRemainingUnexplainedL2Rs = BM.meanW(...
                                 abs(eDState.current.L2Rs(:) + eDState.current.shifts(:)) ... %representative L2R that gets dissectd
                               ,W2D(:));
                end
              %Average strength of the full signal (including possibly inconsistent parts relative to the signature):
                %<-Note: for the ratioExplainedByRegressedVsOriginal statistic we need to look at the full (absolute) signal (including parts inconsistent to the current signature (but in the same focus), because the consistent signal nearly always matches nearly the focusedMeanAbsExplainedByRegressedL2Rs (since the bimonotoinic signature signal is consistent to the signature's eigenOrder per construction):
                if(true)
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsOriginalL2Rs = 'Representative original log2(ratio) strength for the current signature (mean of the original absolute log2(ratio)s, using absolute weights to focus on the pixels affected by the current signature)';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsOriginalL2Rs = BM.meanW(BM.meanW(...
                                 abs(originalSignalWithoutLabSignatures) ... %representative L2R that was originally measured for the pixels carrying the signature
                               ,abs(signatureDefinition.step3_regression.signedFocusedW4G), 1), abs(signatureDefinition.step3_regression.signedFocusedW4P), 2);
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsAlreadyExplainedL2Rs = 'Representative already explained log2(ratio) strength for the current signature (mean of the sum of all previous signature eigensignal log2(ratio)s, using absolute weights to focus on the pixels affected by the current signature)';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsAlreadyExplainedL2Rs = BM.meanW(BM.meanW(...
                                 abs(explainedSignalWithoutLabSignatures) ... %representative already explained L2R for the pixels carrying the signature
                               ,abs(signatureDefinition.step3_regression.signedFocusedW4G), 1), abs(signatureDefinition.step3_regression.signedFocusedW4P), 2);
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsCurrentUnexplainedL2Rs = 'Representative yet unexplained log2(ratio) strength for the current signature (mean of the current unexplained log2(ratio)s, using absolute weights to focus on the pixels affected by the current signature)';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsCurrentUnexplainedL2Rs = BM.meanW(BM.meanW(...
                                 abs(eDState.current.L2Rs) ... %representative unexplained L2R
                               ,abs(signatureDefinition.step3_regression.signedFocusedW4G), 1), abs(signatureDefinition.step3_regression.signedFocusedW4P), 2);
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsExplainedByRegressedL2Rs = 'Representative dissected log2(ratio) strength for the current signature (mean of the current regressed signature eigensignal log2(ratio)s, using absolute weights to focus on the pixels affected by the current signature)';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsExplainedByRegressedL2Rs = BM.meanW(BM.meanW(...
                                 abs(-eDState.current.shifts) ... %representative L2R that gets dissected
                               ,abs(signatureDefinition.step3_regression.signedFocusedW4G), 1), abs(signatureDefinition.step3_regression.signedFocusedW4P), 2);
                  signatureDefinition.explanations.statistics.signalStrength.focusedMeanAbsRemainingUnexplainedL2Rs = 'Representative remaining unexplained log2(ratio) strength still superposed to the current signature (mean of the current remaining log2(ratio)s, using absolute weights to focus on the pixels affected by the current signature)';
                               signatureDefinition.statistics.signalStrength.focusedMeanAbsRemainingUnexplainedL2Rs = BM.meanW(BM.meanW(...
                                 abs(eDState.current.L2Rs + eDState.current.shifts) ... %representative L2R that gets dissected
                               ,abs(signatureDefinition.step3_regression.signedFocusedW4G), 1), abs(signatureDefinition.step3_regression.signedFocusedW4P), 2);
                end
              %Comparisons:
                signatureDefinition.explanations.statistics.signalStrength.ratioExplainedByRegressedVsOriginal = 'Roughly estimates the ratio of the signal part explained by the current regressed signature eigensignal log2(ratio)s versus the original signal strength. Based on the representative signal strengths (see .focusedMeanAbsExplainedByRegressedL2Rs and .focusedMeanAbsOriginalL2Rs).';
                             signatureDefinition.statistics.signalStrength.ratioExplainedByRegressedVsOriginal = min(1, signatureDefinition.statistics.signalStrength.focusedMeanAbsExplainedByRegressedL2Rs / signatureDefinition.statistics.signalStrength.focusedMeanAbsOriginalL2Rs); %min(1,.) because the focusedMeanAbsExplainedByRegressedL2Rs does not have to integrate noise and can thus be larger than focusedMeanAbsOriginalL2Rs.
                signatureDefinition.explanations.statistics.signalStrength.ratioAlreadyExplainedVsOriginal = 'Roughly estimates the ratio of the signal parts already explained by previously detected signature eigensignals versus the original signal strength. Based on representative signal strengths (see .focusedMeanAbsAlreadyExplainedL2Rs and .focusedMeanAbsOriginalL2Rs).';
                             signatureDefinition.statistics.signalStrength.ratioAlreadyExplainedVsOriginal = signatureDefinition.statistics.signalStrength.focusedMeanAbsAlreadyExplainedL2Rs / signatureDefinition.statistics.signalStrength.focusedMeanAbsOriginalL2Rs;
                signatureDefinition.explanations.statistics.signalStrength.ratioExplainedByRegressedVsCurrentUnexplained = 'Roughly estimates the ratio of the signal part explained by the current signature eigensignal log2(ratio)s versus the remaining still unexplained signal. For already perfectly bi-monotonic empirical data this would be near one. For highly overlapped signatures, very noisy empirical data or non-monotonic empirical data (may happen on application stage if the applied signature does not exist in the current dataset), this ratio approaches zero.';
                             signatureDefinition.statistics.signalStrength.ratioExplainedByRegressedVsCurrentUnexplained = min(1, signatureDefinition.statistics.signalStrength.focusedMeanAbsExplainedByRegressedL2Rs / signatureDefinition.statistics.signalStrength.focusedMeanAbsCurrentUnexplainedL2Rs); %min(1,.) because the focusedMeanAbsExplainedByRegressedL2Rs does not have to integrate noise and can thus be larger than focusedMeanAbsCurrentUnexplainedL2Rs.
            end
        %Some cutoff-based signature statistics (e.g. for defining/deriving classical flat-set gene signatures of for selected genes/samples in signature definition plots):
          if(true)
            %signature membership mask:
              if(true)
                %Select members based on max dissection strengths:
                  signatureDefinition.explanations.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot = 'Minimum dissection strength anywhere in a gene row or sample column to be considered as influenced by the detected signature.';
                               signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot = inInfo.postprocessing.cutoffs4statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot;
                  BSignatureGenes = nanmax(dissectionStrengths2DOS,[],2) >= signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot;
                  BSignatureSamples = nanmax(dissectionStrengths2DOS,[],1) >= signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot;
                %Define those genes as representative that have a dissection strength >= inInfo.postprocessing.cutoffs4statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot:
                  signatureDefinition.explanations.statistics.signatureMembership.nRepCoregGenes = sprintf('Number of correlated genes in the signature that have a dissection strengths >=%0.2f for any sample.', signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembership.nRepCoregGenes = sum(BSignatureGenes & signatureDefinition.step3_regression.R4G>=0);
                  signatureDefinition.explanations.statistics.signatureMembership.nRepAntiregGenes = sprintf('Number of anti-correlated genes in the signature that have a dissection strengths >=%0.2f for any sample.', signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembership.nRepAntiregGenes = sum(BSignatureGenes & signatureDefinition.step3_regression.R4G<=0);
                  signatureDefinition.explanations.statistics.signatureMembership.nNonRepGenes = sprintf('Number of genes having dissection strengths <%0.2f for all samples.', signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembership.nNonRepGenes = nG-signatureDefinition.statistics.signatureMembership.nRepCoregGenes-signatureDefinition.statistics.signatureMembership.nRepAntiregGenes;
                %Define those samples as representative that have a dissection strength >= inInfo.postprocessing.cutoffs4statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot:
                  signatureDefinition.explanations.statistics.signatureMembership.nRepAffectedSamples = sprintf('Number of correlated samples in the signature that have a dissection strengths >=%0.2f for any gene.', signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembership.nRepAffectedSamples = sum(BSignatureSamples & signatureDefinition.step3_regression.R4P>=0);
                  signatureDefinition.explanations.statistics.signatureMembership.nRepAntiAffectedSamples = sprintf('Number of anti-correlated samples in the signature that have a dissection strengths >=%0.2f for any gene.', signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembership.nRepAntiAffectedSamples = sum(BSignatureSamples & signatureDefinition.step3_regression.R4P<=0);
                  signatureDefinition.explanations.statistics.signatureMembership.nNonRepSamples = sprintf('Number of samples having dissection strengths <%0.2f for all genes.', signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembership.nNonRepSamples = nP-signatureDefinition.statistics.signatureMembership.nRepAffectedSamples-signatureDefinition.statistics.signatureMembership.nRepAntiAffectedSamples;
              end
              
            %relCorr based top counts:
              if(true)
                %Define those genes as representative that have a relative correlation > inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G:
                  signatureDefinition.explanations.statistics.relCorrGtTopThreshold.ratio4G = 'Configured minimum ratio of the maximum absolute gene correlation of the signature for genes to be counted as top/bottom gene (green vertical cut lines in the default plot).';
                               signatureDefinition.statistics.relCorrGtTopThreshold.ratio4G = inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4G;
                  signatureDefinition.explanations.statistics.relCorrGtTopThreshold.nTopCoregGenes = sprintf('Number of genes that are correlated >=%0.2f*maximum=%0.2f with the current sample axis.', signatureDefinition.statistics.relCorrGtTopThreshold.ratio4G, signatureDefinition.statistics.relCorrGtTopThreshold.ratio4G*max(abs(signatureDefinition.step3_regression.R4G)));
                               signatureDefinition.statistics.relCorrGtTopThreshold.nTopCoregGenes = sum(signatureDefinition.step3_regression.R4G/max(abs(signatureDefinition.step3_regression.R4G)) >= +signatureDefinition.statistics.relCorrGtTopThreshold.ratio4G);
                  signatureDefinition.explanations.statistics.relCorrGtTopThreshold.nTopAntiregGenes = sprintf('Number of genes that are anti-correlated <=-%0.2f*maximum=%0.2f with the current sample axis.', signatureDefinition.statistics.relCorrGtTopThreshold.ratio4G, -signatureDefinition.statistics.relCorrGtTopThreshold.ratio4G*max(abs(signatureDefinition.step3_regression.R4G)));
                               signatureDefinition.statistics.relCorrGtTopThreshold.nTopAntiregGenes = sum(signatureDefinition.step3_regression.R4G/max(abs(signatureDefinition.step3_regression.R4G)) <= -signatureDefinition.statistics.relCorrGtTopThreshold.ratio4G);
                  signatureDefinition.explanations.statistics.relCorrGtTopThreshold.nTopNonregGenes = sprintf('Number of genes with no high (ratio = %0.2f) relative correlation.', signatureDefinition.statistics.relCorrGtTopThreshold.ratio4G);
                               signatureDefinition.statistics.relCorrGtTopThreshold.nNonTopRegGenes = nG-signatureDefinition.statistics.relCorrGtTopThreshold.nTopCoregGenes-signatureDefinition.statistics.relCorrGtTopThreshold.nTopAntiregGenes;
                %Define those samples as representative that have a relative correlation > inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P:
                  signatureDefinition.explanations.statistics.relCorrGtTopThreshold.ratio4P = 'Configured minimum ratio of the maximum absolute sample correlation of the signature for samples to be counted as top/bottom samples (green vertical cut lines in the default plot).';
                               signatureDefinition.statistics.relCorrGtTopThreshold.ratio4P = inInfo.postprocessing.cutoffs4statistics.relCorrTopThresholds.ratio4P;
                  signatureDefinition.explanations.statistics.relCorrGtTopThreshold.nTopAffectedSamples = sprintf('Number of samples that are correlated >=%0.2f*maximum=%0.2f with the current sample axis.', signatureDefinition.statistics.relCorrGtTopThreshold.ratio4P, signatureDefinition.statistics.relCorrGtTopThreshold.ratio4P*max(abs(signatureDefinition.step3_regression.R4P)));
                               signatureDefinition.statistics.relCorrGtTopThreshold.nTopAffectedSamples = sum(signatureDefinition.step3_regression.R4P/max(abs(signatureDefinition.step3_regression.R4P)) >= +signatureDefinition.statistics.relCorrGtTopThreshold.ratio4P);
                  signatureDefinition.explanations.statistics.relCorrGtTopThreshold.nTopAntiAffectedSamples = sprintf('Number of samples that are anti-correlated <=-%0.2f*maximum=%0.2f with the current sample axis.', signatureDefinition.statistics.relCorrGtTopThreshold.ratio4P, -signatureDefinition.statistics.relCorrGtTopThreshold.ratio4P*max(abs(signatureDefinition.step3_regression.R4P)));
                               signatureDefinition.statistics.relCorrGtTopThreshold.nTopAntiAffectedSamples = sum(signatureDefinition.step3_regression.R4P/max(abs(signatureDefinition.step3_regression.R4P)) <= -signatureDefinition.statistics.relCorrGtTopThreshold.ratio4P);
                  signatureDefinition.explanations.statistics.relCorrGtTopThreshold.nNonTopAffectedSamples = sprintf('Number of samples with no high (ratio = %0.2f) relative correlation.', signatureDefinition.statistics.relCorrGtTopThreshold.ratio4P);
                               signatureDefinition.statistics.relCorrGtTopThreshold.nNonTopAffectedSamples = nP-signatureDefinition.statistics.relCorrGtTopThreshold.nTopAffectedSamples-signatureDefinition.statistics.relCorrGtTopThreshold.nTopAntiAffectedSamples;
              end
            %relCorr based outreach counts:
              if(true)
                %Define those genes as representative that have a relative correlation > inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4G:
                  signatureDefinition.explanations.statistics.relCorrGtNoiseThreshold.ratio4G = 'Configured minimum ratio of the maximum absolute gene correlation of the signature for genes to be counted as participating non-noise gene.';
                               signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4G = inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4G;
                  signatureDefinition.explanations.statistics.relCorrGtNoiseThreshold.nRepCoregGenes = sprintf('Number of genes that are correlated >=%0.2f*maximum=%0.2f with the current sample axis.', signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4G, signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4G*max(abs(signatureDefinition.step3_regression.R4G)));
                               signatureDefinition.statistics.relCorrGtNoiseThreshold.nRepCoregGenes = sum(signatureDefinition.step3_regression.R4G/max(abs(signatureDefinition.step3_regression.R4G)) >= +signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4G);
                  signatureDefinition.explanations.statistics.relCorrGtNoiseThreshold.nRepAntiregGenes = sprintf('Number of genes that are anti-correlated <=-%0.2f*maximum=%0.2f with the current sample axis.', signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4G, -signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4G*max(abs(signatureDefinition.step3_regression.R4G)));
                               signatureDefinition.statistics.relCorrGtNoiseThreshold.nRepAntiregGenes = sum(signatureDefinition.step3_regression.R4G/max(abs(signatureDefinition.step3_regression.R4G)) <= -signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4G);
                  signatureDefinition.explanations.statistics.relCorrGtNoiseThreshold.nNonRepGenes = sprintf('Number of genes with lower relative correlation than ratio=%0.2f, considered as noise wrt. the current signature.', signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4G);
                               signatureDefinition.statistics.relCorrGtNoiseThreshold.nNonRepGenes = nG-signatureDefinition.statistics.relCorrGtNoiseThreshold.nRepCoregGenes-signatureDefinition.statistics.relCorrGtNoiseThreshold.nRepAntiregGenes;
                %Define those samples as representative that have a relative correlation > inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4G:
                  signatureDefinition.explanations.statistics.relCorrGtNoiseThreshold.ratio4P = 'Configured minimum ratio of the maximum absolute sample correlation of the signature for samples to be counted as affected non-noise sample.';
                               signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4P = inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4P;
                  signatureDefinition.explanations.statistics.relCorrGtNoiseThreshold.nRepAffectedSamples = sprintf('Number of samples that are correlated >=%0.2f*maximum=%0.2f with the current sample axis.', signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4P, signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4P*max(abs(signatureDefinition.step3_regression.R4P)));
                               signatureDefinition.statistics.relCorrGtNoiseThreshold.nRepAffectedSamples = sum(signatureDefinition.step3_regression.R4P/max(abs(signatureDefinition.step3_regression.R4P)) >= +signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4P);
                  signatureDefinition.explanations.statistics.relCorrGtNoiseThreshold.nRepAntiAffectedSamples = sprintf('Number of samples that are anti-correlated <=-%0.2f*maximum=%0.2f with the current sample axis.', signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4P, -signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4P*max(abs(signatureDefinition.step3_regression.R4P)));
                               signatureDefinition.statistics.relCorrGtNoiseThreshold.nRepAntiAffectedSamples = sum(signatureDefinition.step3_regression.R4P/max(abs(signatureDefinition.step3_regression.R4P)) <= -signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4P);
                  signatureDefinition.explanations.statistics.relCorrGtNoiseThreshold.nNonRepSamples = sprintf('Number of samples with lower relative correlation than ratio=%0.2f, considered as noise wrt. the current signature.', signatureDefinition.statistics.relCorrGtNoiseThreshold.ratio4P);
                               signatureDefinition.statistics.relCorrGtNoiseThreshold.nNonRepSamples = nP-signatureDefinition.statistics.relCorrGtNoiseThreshold.nRepAffectedSamples-signatureDefinition.statistics.relCorrGtNoiseThreshold.nRepAntiAffectedSamples;
              end

            %signature membership mask and relCorrRatio based:
              if(true)
                %Select members based on max dissection strengths:
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.minRequireddDissectionStrengthToShowInPlot = 'Minimum dissection strength anywhere in a gene row or sample column to be considered as influenced by the detected signature.';
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.minRequireddDissectionStrengthToShowInPlot = inInfo.postprocessing.cutoffs4statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot;
                  BSignatureGenes = nanmax(dissectionStrengths2DOS,[],2) >= signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.minRequireddDissectionStrengthToShowInPlot;
                  BSignatureSamples = nanmax(dissectionStrengths2DOS,[],1) >= signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.minRequireddDissectionStrengthToShowInPlot;
                %Define those genes as representative that have a dissection strength >= signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.minRequireddDissectionStrengthToShowInPlot:
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G = 'Configured minimum ratio of the maximum absolute gene correlation of the signature for genes to be counted as participating non-noise gene.';
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G = inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4G;
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepCoregGenes = sprintf('Number of genes correlated >=%0.2f*maximum=%0.2f in the signature that have a dissection strengths >=%0.2f for any sample.', signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G, signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G*max(abs(signatureDefinition.step3_regression.R4G)), signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepCoregGenes = sum(BSignatureGenes & signatureDefinition.step3_regression.R4G/max(abs(signatureDefinition.step3_regression.R4G)) >= +signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G);
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAntiregGenes = sprintf('Number of genes correlated <=-%0.2f*maximum=%0.2f in the signature that have a dissection strengths >=%0.2f for any sample.', signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G, signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G*max(abs(signatureDefinition.step3_regression.R4G)), signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAntiregGenes = sum(BSignatureGenes & signatureDefinition.step3_regression.R4G/max(abs(signatureDefinition.step3_regression.R4G)) <= -signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G);
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nNonRepGenes = sprintf('Number of genes in the signature with lower relative correlation than ratio=%0.2f*maximum or dissection strengths <=%0.2f for any sample.', signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4G, signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nNonRepGenes = nG-signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepCoregGenes-signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAntiregGenes;
                %Define those samples as representative that have a dissection strength >= signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.minRequireddDissectionStrengthToShowInPlot:
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P = 'Configured minimum ratio of the maximum absolute sample correlation of the signature for samples to be counted as affected non-noise sample.';
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P = inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4P;
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAffectedSamples = sprintf('Number of samples correlated >=%0.2f*maximum=%0.2f in the signature that have a dissection strengths >=%0.2f for any gene.', signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P, signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P*max(abs(signatureDefinition.step3_regression.R4P)), signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAffectedSamples = sum(BSignatureSamples & signatureDefinition.step3_regression.R4P/max(abs(signatureDefinition.step3_regression.R4P)) >= +signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P);
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAntiAffectedSamples = sprintf('Number of samples correlated <=-%0.2f*maximum=%0.2f in the signature that have a dissection strengths >=%0.2f for any gene.', signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P, signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P*max(abs(signatureDefinition.step3_regression.R4P)), signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAntiAffectedSamples = sum(BSignatureSamples & signatureDefinition.step3_regression.R4P/max(abs(signatureDefinition.step3_regression.R4P)) <= -signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P);
                  signatureDefinition.explanations.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nNonRepSamples = sprintf('Number of samples in the signature with lower relative correlation than ratio=%0.2f*maximum or dissection strengths <=%0.2f for any gene.', signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.ratio4P, signatureDefinition.statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot);
                               signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nNonRepSamples = nP-signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAffectedSamples-signatureDefinition.statistics.signatureMembershipAndRelCorrGtNoiseThreshold.nRepAntiAffectedSamples;
              end
              
            %p4R value based significant counts (using inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha):
              if(true)
                %Define those genes as representative that have a p value < inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G:
                  signatureDefinition.explanations.statistics.pLtSignificanceThreshold.alpha4G = 'Configured maximum p value for the correlation of a gene with the sample axis to count it as top correlated / bottom anti-correlated gene of the signature.';
                               signatureDefinition.statistics.pLtSignificanceThreshold.alpha4G = inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G;
                  signatureDefinition.explanations.statistics.pLtSignificanceThreshold.nTopCoregGenes = sprintf('Number of regulated genes that are significantly correlated with the current sample axis (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G);
                               signatureDefinition.statistics.pLtSignificanceThreshold.nTopCoregGenes = sum(signatureDefinition.step3_regression.R4G>0 & signatureDefinition.step3_regression.P4R4G < inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G);
                  signatureDefinition.explanations.statistics.pLtSignificanceThreshold.nTopAntiregGenes = sprintf('Number of regulated genes that are significantly anti-correlated with the current sample axis (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G);
                               signatureDefinition.statistics.pLtSignificanceThreshold.nTopAntiregGenes = sum(signatureDefinition.step3_regression.R4G<0 & signatureDefinition.step3_regression.P4R4G < inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G);
                  signatureDefinition.explanations.statistics.pLtSignificanceThreshold.nTopNonregGenes = sprintf('Number of not significantly (alpha = %0.5f) correlated genes.', inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G);
                               signatureDefinition.statistics.pLtSignificanceThreshold.nNonTopRegGenes = nG-signatureDefinition.statistics.pLtSignificanceThreshold.nTopCoregGenes-signatureDefinition.statistics.pLtSignificanceThreshold.nTopAntiregGenes;
                %Define those samples as representative that have a p value < inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P:
                  signatureDefinition.explanations.statistics.pLtSignificanceThreshold.alpha4P = 'Configured maximum p value for the correlation of a sample with the samples signature to count it as top correlated / bottom anti-correlated sample of the signature.';
                               signatureDefinition.statistics.pLtSignificanceThreshold.alpha4P = inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P;
                  signatureDefinition.explanations.statistics.pLtSignificanceThreshold.nTopAffectedSamples = sprintf('Number of affected samples that are significantly correlated with the current gene axis (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P);
                               signatureDefinition.statistics.pLtSignificanceThreshold.nTopAffectedSamples = sum(signatureDefinition.step3_regression.R4P>0 & signatureDefinition.step3_regression.P4R4P < inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P);
                  signatureDefinition.explanations.statistics.pLtSignificanceThreshold.nTopAntiAffectedSamples = sprintf('Number of affected samples that are significantly anti-correlated with the current gene axis (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P);
                               signatureDefinition.statistics.pLtSignificanceThreshold.nTopAntiAffectedSamples = sum(signatureDefinition.step3_regression.R4P<0 & signatureDefinition.step3_regression.P4R4P < inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P);
                  signatureDefinition.explanations.statistics.pLtSignificanceThreshold.nNonTopAffectedSamples = sprintf('Number of not significantly (alpha = %0.5f) correlated samples.', inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4P);
                               signatureDefinition.statistics.pLtSignificanceThreshold.nNonTopAffectedSamples = nP-signatureDefinition.statistics.pLtSignificanceThreshold.nTopAffectedSamples-signatureDefinition.statistics.pLtSignificanceThreshold.nTopAntiAffectedSamples;
              end
            %p4R value based outreach counts (using inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha):
              if(true)
                %Define those genes as representative that have a p value < inInfo.postprocessing.cutoffs4statistics.significanceThresholds.alpha4G:
                  signatureDefinition.explanations.statistics.pLtNoiseThreshold.alpha4G = 'Configured maximum p value for the correlation of a gene with the sample axis to not count it as noise gene with respect to the signature. Note that the correlation strength should be considered too, especially for large datasets.';
                               signatureDefinition.statistics.pLtNoiseThreshold.alpha4G = inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4G;
                  signatureDefinition.explanations.statistics.pLtNoiseThreshold.nRepCoregGenes = sprintf('Number of regulated genes that are correlated and not considered noise wrt. to the current signature (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4G);
                               signatureDefinition.statistics.pLtNoiseThreshold.nRepCoregGenes = sum(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder)>0 & signatureDefinition.step3_regression.P4R4G < inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4G);
                  signatureDefinition.explanations.statistics.pLtNoiseThreshold.nRepAntiregGenes = sprintf('Number of regulated genes that are anti-correlated and not considered noise wrt. the current signature (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4G);
                               signatureDefinition.statistics.pLtNoiseThreshold.nRepAntiregGenes = sum(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4geneOrder)<0 & signatureDefinition.step3_regression.P4R4G < inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4G);
                  signatureDefinition.explanations.statistics.pLtNoiseThreshold.nNonRepGenes = sprintf('Number of not correlated genes considered as noise wrt. the current signature (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4G);
                               signatureDefinition.statistics.pLtNoiseThreshold.nNonRepGenes = nG-signatureDefinition.statistics.pLtNoiseThreshold.nRepCoregGenes-signatureDefinition.statistics.pLtNoiseThreshold.nRepAntiregGenes;
                %Define those samples as representative that have a p value < inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P:
                  signatureDefinition.explanations.statistics.pLtNoiseThreshold.alpha4P = 'Configured maximum p value for the correlation of a sample with the samples signature to not count it as noise sample with respect to the signature. Note that the correlation strength should be considered too, especially for large datasets.';
                               signatureDefinition.statistics.pLtNoiseThreshold.alpha4P = inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P;
                  signatureDefinition.explanations.statistics.pLtNoiseThreshold.nRepAffectedSamples = sprintf('Number of affected samples that are correlated and not considered noise wrt. the current signature (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P);
                               signatureDefinition.statistics.pLtNoiseThreshold.nRepAffectedSamples = sum(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder)>0 & signatureDefinition.step3_regression.P4R4P < inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P);
                  signatureDefinition.explanations.statistics.pLtNoiseThreshold.nRepAntiAffectedSamples = sprintf('Number of affected samples that are anti-correlated and not considered noise wrt. the current signature (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P);
                               signatureDefinition.statistics.pLtNoiseThreshold.nRepAntiAffectedSamples = sum(signatureDefinition.step3_regression.(signatureDefinition.reference.metric4sampleOrder)<0 & signatureDefinition.step3_regression.P4R4P < inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P);
                  signatureDefinition.explanations.statistics.pLtNoiseThreshold.nNonRepSamples = sprintf('Number of not correlated samplees considered noise wrt. the current signature (alpha = %0.5f).', inInfo.postprocessing.cutoffs4statistics.noiseThresholds.alpha4P);
                               signatureDefinition.statistics.pLtNoiseThreshold.nNonRepSamples = nP-signatureDefinition.statistics.pLtNoiseThreshold.nRepAffectedSamples-signatureDefinition.statistics.pLtNoiseThreshold.nRepAntiAffectedSamples;
              end
          end
      %% EXPORT) save signature definition and invoke plot function:
        SDCM_printStatus(0,'### Save and export signature definition %d...\n', eDState.current.k); drawnow;
          outInfo.reference.timing.b4_startOfExportsAndPlotting(end+1) = now;
        %Export dir for current signature:
          inInfo.export.subdir4signature = fullfile(inInfo.export.rootDir);
          inInfo.export.fileNamePrefix4signature = ''; %defined in plotSignature.m
            if(~isempty(inInfo.export.fileNamePrefix4signature))
              inInfo.export.fileNamePrefix4signature = [inInfo.export.fileNamePrefix4signature, ', '];
            end
        %Signature statistics description:
          if(true)
            signatureDefinition.statistics.sSignatureStatistics = sprintf(...
               '%03d, focused signature = corr=(%2.1fg,%2.1fs,%0.2fr), p=(%1.0ecorr,%1.0esignal), affects=(%d+%dg,%d+%ds), avgSignal=%0.2f=%0.0f%%orig'...
              ,signatureDefinition.reference.k...
              ...
              ,signatureDefinition.step3_regression.signatureSizeByCorrSum4G ...
              ,signatureDefinition.step3_regression.signatureSizeByCorrSum4P ...
              ,signatureDefinition.step3_regression.signatureCorrInExtendedFocus ...
              ,10^signatureDefinition.step3_regression.log10_p4Correlations, 10^signatureDefinition.step3_regression.log10_p4SignalStrength ...
              ...
              ,signatureDefinition.statistics.relCorrGtTopThreshold.nTopCoregGenes ...
              ,signatureDefinition.statistics.relCorrGtTopThreshold.nTopAntiregGenes ...
              ,signatureDefinition.statistics.relCorrGtTopThreshold.nTopAntiAffectedSamples ... %as in the plots on the left show the anti-affected number of samples first.
              ,signatureDefinition.statistics.relCorrGtTopThreshold.nTopAffectedSamples ...
              ...
              ,signatureDefinition.statistics.signalStrength.focusedMeanAbsExplainedByRegressedL2Rs ...
              ,100*signatureDefinition.statistics.signalStrength.ratioExplainedByRegressedVsOriginal ...
            );
          end
          %Create a small info file with the core statistics of the signature:
            if(inInfo.export.infoFile4focusedSignatureStats.bEnabled)
              fileName = sprintf('%s.info', signatureDefinition.statistics.sSignatureStatistics);
              targetDir = inInfo.export.subdir4signature;
                if(ispc() && ~(length(targetDir)>=3 && strcmp(targetDir(1:3),'\\?'))) %add \\?\ prefix to support long paths.
                  %If we already have a network UNC path line \\Server\Share\..., we need to replace \\ by UNC\ to get long path support in Windows:
                    bIsAlreadyANetworkPath = length(targetDir)>=2 && strcmp(targetDir(1:2),'\\');
                    if(bIsAlreadyANetworkPath)
                      targetDir = ['UNC',targetDir(2:end)];
                    end
                  targetDir = ['\\?\',targetDir];
                end
              fid = fopen(fullfile(targetDir,fileName),'w');
              fprintf(fid,'%s',signatureDefinition.statistics.sSignatureStatistics);
              fclose(fid);
            end
        %Complete signature reference info:
          if(true)
            %Append reference information to every single signature definition:
              signatureDefinition.explanations.reference.exportDir = 'Path to the exported data for this signature.';
                           signatureDefinition.reference.exportDir = inInfo.export.subdir4signature;
              signatureDefinition.explanations.reference.fileNamePrefix4signature = 'File name prefix for exported data for this signature.';
                           signatureDefinition.reference.fileNamePrefix4signature = inInfo.export.fileNamePrefix4signature;
            %Save fields needed for resume feature:
              signatureDefinition.explanations.forResume = 'Backup of detection state data required for the resume feature to work.';
                signatureDefinition.reference.forResume.BGeneDidNotQualifyAfterDetectionOfThisSignature = eDState.current.performance.BGeneDidNotQualify; %needed for resume feature.
                signatureDefinition.reference.forResume.BSampleDidNotQualifyAfterDetectionOfThisSignature = eDState.current.performance.BSampleDidNotQualify; %needed for resume feature.
                signatureDefinition.reference.forResume.bNoSignatureDetectedSinceLastBDidNotQualifyReset = eDState.current.performance.bNoSignatureDetectedSinceLastBDidNotQualifyReset; %needed for resume feature.
                signatureDefinition.reference.forResume.reconstrutionChecksums.L2Rs_before = nansum(abs(eDState.current.L2Rs(:)));
                signatureDefinition.reference.forResume.reconstrutionChecksums.shifts = nansum(abs(eDState.current.shifts(:)));
                signatureDefinition.reference.forResume.noiseEstimation = eDState.noiseEstimation; 
            %Save additional data later needed by some plots:
              signatureDefinition.explanations.forPlots = 'Copies of various intermediade results from the detection phase that are needed to create the default plot for this signature.';
                signatureDefinition.forPlots.overview.original.L2Rs = eDState.initialSignal.L2Rs;
                signatureDefinition.forPlots.overview.current.L2Rs = eDState.current.L2Rs;
                signatureDefinition.forPlots.overview.current.stds4L2Rs2D = eDState.current.stds4L2Rs;
                signatureDefinition.forPlots.overview.current.explainedSignal = eDState.current.explainedSignal;
                signatureDefinition.forPlots.overview.current.shifts = eDState.current.shifts;
          end
        %MAT) create slim version of the current signatureDefinition and export it as MAT file:
          if(inInfo.export.matFile4eachSignatureDefinition.bEnabled)
            signatureDef = signatureDefinition;
              %Remove plotting info for the default plot, if not enabled for export:
                if(~inInfo.export.matFile4eachSignatureDefinition.bIncludeDefPlotPanels && isfield(signatureDef.forPlots, 'overview'))
                  signatureDef.forPlots = rmfield(signatureDef.forPlots, 'overview');
                end
                if(~inInfo.export.matFile4eachSignatureDefinition.bIncludeRegressedSignalEstimationSteps && isfield(signatureDef.forPlots, 'dissectionPipelines'))
                  signatureDef.forPlots = rmfield(signatureDef.forPlots, 'dissectionPipelines');
                end
                if(~inInfo.export.matFile4eachSignatureDefinition.bIncludeBimonotonicConvergenceSteps && isfield(signatureDef.forPlots, 'bimonotonicConvergence'))
                  signatureDef.forPlots = rmfield(signatureDef.forPlots, 'bimonotonicConvergence');
                end
              %Hotfix: remove some plotting info for signatureDefinitions with very large memory consumption:
                if(ilv(whos('signatureDef'),@(S)S.bytes)/1024^3>1 && isfield(signatureDef,'forPlots') && isfield(signatureDef.forPlots,'dissectionPipelines'))
                  warning('removing signatureDef.forPlots.dissectionPipelines before saving to disk, since the variable is >1GB');
                  signatureDef.forPlots = rmfield(signatureDef.forPlots,'dissectionPipelines');
                end
              %Save .forPlots separately as it consumes much disk space and slows down resumeFromDisk dramatically:
                if(isstruct(signatureDef.forPlots))
                  forPlots = signatureDef.forPlots;
                  exportFilename = fullfile(inInfo.export.subdir4signature, sprintf('%03d, pipeline state for plots.mat',signatureDefinition.reference.k));
                  SDCM_printStatus(2,' - saving extra info for plots of the detected signature to: %s ...\n', exportFilename);
                  try
                    if(~inInfo.export.matFile4eachSignatureDefinition.bDisablePipelineStateForPlotsFileExport(min(end,sL)))
                      save(exportFilename, 'forPlots');
                    end
                    clear forPlots;
                    signatureDef.forPlots = {'loadFromFile',exportFilename,'.forPlots'};
                  catch ex
                    warning('ERROR writing export to %s\nError details: %s\nThis may be caused by a full hard disk or file system access restrictions. SKIPPING this export and continuing detection... confirm with >>dbcont or inspect', exportFilename, ex.message);
                    keyboard;
                  end
                end
              %Save .noiseEstimation separately as it consumes much disk space and slows down resumeFromDisk dramatically:
                if(isstruct(signatureDef.reference.forResume.noiseEstimation))
                  noiseEstimation = signatureDef.reference.forResume.noiseEstimation;
                  exportFilename = fullfile(inInfo.export.subdir4signature, sprintf('%03d, noise distribution.mat',signatureDefinition.reference.k));
                  SDCM_printStatus(2,' - saving current noise estimation to: %s ...\n', exportFilename);
                  try
                    if(~inInfo.export.matFile4eachSignatureDefinition.bDisableNoiseDistributionFileExport(min(end,sL)))
                      save(exportFilename, 'noiseEstimation');
                    end
                    clear noiseEstimation;
                    signatureDef.reference.forResume.noiseEstimation = {'loadFromFile',exportFilename,'.noiseEstimation'};
                    %Mem saver / remove signatureDefinition.reference.forResume.noiseEstimation if export is enabled and replace it by a file load link:
                      if(inInfo.export.matFile4eachSignatureDefinition.bEnabled && inInfo.export.matFile4eachSignatureDefinition.bRemoveNoiseEstimationResumeDataAfterExport && isstruct(signatureDefinition.reference.forResume.noiseEstimation))
                        signatureDefinition.reference.forResume.noiseEstimation = {'loadFromFile',exportFilename,'.noiseEstimation'};
                      end
                  catch ex
                    warning('ERROR writing export to %s\nError details: %s\nThis may be caused by a full hard disk or file system access restrictions. SKIPPING this export and continuing detection... confirm with >>dbcont or inspect', exportFilename, ex.message);
                    keyboard;
                  end
                end
              exportFilename = fullfile(inInfo.export.subdir4signature, sprintf('%03d, definition.mat',signatureDefinition.reference.k));
              SDCM_printStatus(2,' - saving definition of the detected signature to: %s ...\n', exportFilename);
              try
                save(exportFilename, 'signatureDef');
              catch ex
                warning('ERROR writing export to %s\nError details: %s\nThis may be caused by a full hard disk or file system access restrictions. SKIPPING this export and continuing detection... confirm with >>dbcont or inspect', exportFilename, ex.message);
                keyboard;
              end
            clear signatureDef;
          end
        %MEM) save signature in memory (in the global eDState variable):
          if(true)
            signatureDef = signatureDefinition;
            %Mem saver / remove redundant reference info from signature definition:
              if(isfield(signatureDef.reference,'inInfo'))
                signatureDef.reference = rmfield(signatureDef.reference,'inInfo'); %this is only included to make the signature saved to disc become self-contained.
              end
            %Mem saver / remove plotting info from RAM, if enabled:
              if(inInfo.export.matFile4eachSignatureDefinition.bEnabled && inInfo.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport)
                if(~inInfo.export.matFile4eachSignatureDefinition.bEnabled)
                  warning('The plotting info removed by inInfo.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport cannot be restored, since ~inInfo.export.matFile4eachSignatureDefinition.bEnabled, i.e. no standalone default plots of this signature can be created at a later time.');
                end
                if(isfield(signatureDef,'forPlots')) %make structs array compatible wrt. fields, if we resumed from different save plotsInfo config.
                  signatureDef = rmfield(signatureDef, 'forPlots');
                end
                %remove all plot data (from all signatures in the struct array):
                  if(isfield(eDState.signatures,'forPlots')) %make structs array compatible wrt. fields, if we resumed from different save plotsInfo config.
                    eDState.signatures = rmfield(eDState.signatures, 'forPlots');
                  end
              else
                %interactively late-enable inInfo.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport if we have run into a low memory situation:
                  if(ispc())
                    remainingMem = memory;
                    if(~inInfo.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport && remainingMem.MaxPossibleArrayBytes/1024^2 < 500)
                      warning('<PAUSED>: less than 500MB RAM left; consider enabling "inInfo.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport = true;" before dbcont to avoid out of memory errors');
                      keyboard;
                      if(inInfo.export.matFile4eachSignatureDefinition.bRemoveNoncoreInfoFromRAMandOutputAfterExport && isfield(eDState.signatures,'forPlots'))
                        eDState.signatures = rmfield(eDState.signatures, 'forPlots');
                      end
                    end
                  end
              end
            %Save in global store eDState:
              if(isempty(eDState.signatures))
                eDState.signatures = signatureDef;
              else
                eDState.signatures(eDState.current.k) = signatureDef;
              end
            clear signatureDef;
          end
        %PLOT) if plot config is passed, invoke the plot function:
          if(inInfo.plots.signatureDefinitions.bEnabled)
            %Add information needed for plotting:
              if(inInfo.applicationMode.bEnabled)
                signatureDefinition.reference.maxk = inInfo.applicationMode.maxk_externallyValidatedSignatures;
                if(ischar(signatureDefinition.reference.maxk))
                  assert(strcmp(signatureDefinition.reference.maxk,'DEFAULT'), 'inInfo.applicationMode.maxk_externallyValidatedSignatures; must be numeric or "DEFAULT"');
                  signatureDefinition.reference.maxk = length(previouslyDetectedSignatures);
                end
              end
              %feature .nInitialSignaturesNotToPlot: Hide a user-specified number of initial/brightness/lab signatures from the plots:
                if(inInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot > 0)
                  for k=1:min(inInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot,eDState.current.k-1)
                    shifts = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision);
                    shifts(...
                      eDState.signatures(k).step3_regression.BNonzeroRows4signatureEigensignal...
                     ,eDState.signatures(k).step3_regression.BNonzeroCols4signatureEigensignal...
                    ) = -eDState.signatures(k).step3_regression.signatureEigensignal; %we shift the negative eigensignal of the signature to dissect the signature from the signal.

                    signatureDefinition.forPlots.overview.original.L2Rs = signatureDefinition.forPlots.overview.original.L2Rs + shifts;
                    signatureDefinition.forPlots.overview.current.explainedSignal = signatureDefinition.forPlots.overview.current.explainedSignal - shifts;
                  end
                end
            try
              plotDetectedSignature(signatureDefinition, inInfo); %also passes on the relevant export settings in inInfo.
            catch ex
              warning('ERROR during plotting: %s\nThis may be caused by an incompatibility of your local installation or platform with various plotting requirements. You may want to use "dbstop if caught error" (or execute >> plotDetectedSignature(signatureDefinition, inInfo) now) to debug the error yourself; or ask someone with MATLAB experience to help you. Note that the plotting code is provided AS IS and is not supported.\nInteractive break for manual retry/inspection; use >>dbcont to ignore this plotting error and continue.', ex.message);
              %ex.stack.name
              %ex.stack.line
              keyboard;
            end
          end
      %% DISSECT the signature signal and continue to next iteration:
        %Dissect the signature signal from the current signal and accumulate the explained signal:
          if(~inInfo.applicationMode.bEnabled || ~inInfo.applicationMode.bApplyEveryExternalSignatureToInitialSignal)
            eDState.current.L2Rs = eDState.current.L2Rs + eDState.current.shifts; %dissect signature (shifts equal the negative signature eigensignal).
            eDState.current.explainedSignal = eDState.current.explainedSignal + eDState.current.shifts; %negative of already explained signal parts by all so far detected signatures.
          end
          eDState.signatures(eDState.current.k).reference.forResume.reconstrutionChecksums.L2Rs_after = nansum(abs(eDState.current.L2Rs(:)));
        %Consistency check: The current L2Is must be identical to the original L2Is plus all dissected signature signals:
          if(inInfo.internal.bConsistencyCheckAfterEachDetectionIteration)
            SDCM_printStatus(0,'### State consistency check (disable inInfo.internal.bConsistencyCheckAfterEachDetectionIteration to skip)...\n');
            %Assert accumulative consistency first:
              epsilon = 10^-3;
              if(nanmax(ilsub(abs(...
                  eDState.current.L2Rs ... %current unexplained signal 
               - (eDState.initialSignal.L2Rs + eDState.current.explainedSignal) ...
              ),':')) > epsilon)
                warning('>>> consistency check: detected accumulative error in applied .explainedSignal; entering interactive mode to investigate; use dbcont to ignore and continue.');
                beep;
                keyboard;
              end
            %Reconstruct the current signal from the saved signature baselines:
              if(true)
                reconstructed = struct();
                  reconstructed.current.k = eDState.current.k;
                  reconstructed.current.L2Rs = eDState.initialSignal.L2Rs;
                  reconstructed.current.explainedSignal = zeros(nG,nP);
                    if(inInfo.applicationMode.bEnabled)
                      assert(all(~isnan(inInfo.applicationMode.externalInitialShifts(:))), 'some inInfo.applicationMode.externalInitialShifts were NaN.');
                      reconstructed.current.explainedSignal = reconstructed.current.explainedSignal + inInfo.applicationMode.externalInitialShifts;
                      reconstructed.current.L2Rs = reconstructed.current.L2Rs + inInfo.applicationMode.externalInitialShifts;
                    end
                currentShifts = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision); %preallocate. single is important (must be the same datatype as when the checksums like .signatures(k).reference.forResume.reconstrutionChecksums were created!
                  for k=1:reconstructed.current.k
                    %Check initial L2Rs:
                      if(k==reconstructed.current.k) %performance/only check for last:
                        assert(abs( nansum(abs(reconstructed.current.L2Rs(:))) - eDState.signatures(k).reference.forResume.reconstrutionChecksums.L2Rs_before ) < 10^-2, 'L2Is_before reconstruction checksum mismatch');
                      end
                    %Initialize:
                      currentShifts(:) = 0;
                    %Add dissection shifts (with puned genes ==0):
                      currentShifts(...
                        eDState.signatures(k).step3_regression.BNonzeroRows4signatureEigensignal...
                       ,eDState.signatures(k).step3_regression.BNonzeroCols4signatureEigensignal...
                      ) = -eDState.signatures(k).step3_regression.signatureEigensignal;
                    %Check reconstructed shifts:
                      assert(abs( nansum(abs(currentShifts(:))) - eDState.signatures(k).reference.forResume.reconstrutionChecksums.shifts ) < 10^-2, 'shifts reconstruction checksum mismatch');
                    %Apply shifts and check result:
                    	if(~inInfo.applicationMode.bEnabled || ~inInfo.applicationMode.bApplyEveryExternalSignatureToInitialSignal)
                        reconstructed.current.L2Rs = reconstructed.current.L2Rs + currentShifts;
                        if(...
                          isfield(eDState.signatures(k).reference.forResume.reconstrutionChecksums,'L2Is_after') ... %field not available, if loaded from disk.
                        )
                          assert(abs( nansum(abs(reconstructed.current.L2Rs(:))) - eDState.signatures(k).reference.forResume.reconstrutionChecksums.L2Rs_after ) < 10^-2, 'L2Is_after reconstruction checksum mismatch');
                        end
                        reconstructed.current.explainedSignal = reconstructed.current.explainedSignal + currentShifts;
                      end
                    SDCM_printStatus(2, '    - reconstructed shifts for already detected signature %d/%d\n',k,eDState.current.k);
                  end
              end
            %Assert reconstruction consistency:
              epsilon = 10^-3;              
              if(nanmax(ilsub(abs(...
                 eDState.current.L2Rs ... %current unexplained signal 
               - reconstructed.current.L2Rs ... %must always be reconstructible from the full original input signal by dissecting detected&regressed signature eigensignals
              ),':')) > epsilon)
                warning('>>> consistency check: detected error in reconstructed .explainedSignal; entering interactive mode to investigate; use dbcont to ignore and continue.');
                beep;
                keyboard;
              end
              if(nanmax(ilsub(abs(...
                  eDState.current.L2Rs ... %current unexplained signal 
               - (eDState.initialSignal.L2Rs + reconstructed.current.explainedSignal) ... %must always be reconstructible from the full original input signal by dissecting detected&regressed signature eigensignals
              ),':')) > epsilon)
                warning('>>> consistency check: detected error in reconstructed .explainedSignal; entering interactive mode to investigate; use dbcont to ignore and continue.');
                beep;
                keyboard;
              end
            SDCM_printStatus(1,'   <- state consistent.\n');
          end
        %allow optional interactive breakpoint:
          if(inInfo.internal.bDevEnableInteractiveBreaks)
            interactiveCheckPoint();
          end
        %write/flush the detection log, if enabled:
          if(inInfo.export.bCreateDetectionLogFile)
            diary off; %flush console to log.
            diary(fullfile(inInfo.export.rootDir,sLogFileName)); %continue logging.
          end            
    end
    %Collect outputs:
      if(true)
        outInfo.explanations.reference.initialSignal.L2Is = 'initial log2(intensities) after any enabled preprocessing steps';
                     outInfo.reference.initialSignal.L2Is = eDState.initialSignal.L2Is;
        outInfo.explanations.reference.initialSignal.base4G = 'initial gene baseline relative to which log2(ratio)s are determined (=the gene-wise cohort mean(L2Is), if .preprocessing.bConvertIntoCohortRelativeL2RsAtStart)';
                     outInfo.reference.initialSignal.base4G = eDState.initialSignal.base4G;
        outInfo.explanations.reference.initialSignal.base4P = 'initial sample baseline (zeros(1,nP) in the default; included for implementation symmetry)';
                     outInfo.reference.initialSignal.base4P = eDState.initialSignal.base4P;
        outInfo.explanations.reference.initialSignal.L2Rs = 'initial log2(ratio)s relative to .base4G; this is the start signal for detection and dissection into signatures (M_0 in the paper)';
                     outInfo.reference.initialSignal.L2Rs = eDState.initialSignal.L2Rs;
        outInfo.explanations.reference.inInfo = 'reference of the used detection configuration (see .reference.inInfo.explanations for parameter details)';
                     outInfo.reference.inInfo = inInfo;
        outInfo.explanations.signatures = 'structure array of all detected signatures and their signatures (see .explanations of every single signature for a detailed documentation of their output fields)';
                     outInfo.signatures = eDState.signatures(1:eDState.current.k-1);
        outInfo.explanations.remainingSignal.L2Rs = 'remaining unexplained log2(ratio)s after all signatures that qualified relative to the detection configuration (respectively were externally provided on the application stage) have been dissectd';
                     outInfo.remainingSignal.L2Rs = eDState.current.L2Rs;
      end
  %% Close the log and return results:
    outInfo.reference.timing.d_returningOutput = now;
    SDCM_printStatus(0,'################################################################################################################################################################\n');
    SDCM_printStatus(0,'### FINISHED) Returning %d signatures at %s after %s computation time (days:hours:minutes:seconds).\n', length(outInfo.signatures), datestr(now), datestr(-outInfo.reference.timing.a_startOfPreprocessing+outInfo.reference.timing.d_returningOutput,'dd:HH:MM:SS'));
    if(inInfo.export.nStatusOutputLevel>=2 && ~inInfo.applicationMode.bEnabled)
      for k=1:length(outInfo.signatures)
        SDCM_printStatus(2, sprintf(' %03d) %s\n', k, outInfo.signatures(k).step1_initialRepresentative.sInitialRepresenative));
      end
    end
    SDCM_printStatus(0,'#############\n');
    SDCM_printStatus(0,'#############\n');
    %stop and save detection.log, if enabled:
      if(inInfo.export.bCreateDetectionLogFile)
        diary off;
      end
end

