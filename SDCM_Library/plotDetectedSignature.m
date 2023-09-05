% ABSTRACT
%  Plot detected signature and the associated iteration of signal dissection.
%  panels: original signal + already explained signal parts 
%        = current standardized signal in the detected signature 
% LICENSE
%  Note that the plot code is provided AS IS and may not work with your local 
%  platform and installation as intended. You are free to adapt the plot code 
%  to your needs as long as properly cited, but this code is not supported.
% AUTHOR
%  (C) Michael Grau, 2012-2014

function outInfo = plotDetectedSignature(signatureDefinition, inInfo)
  %%.Parameter logic:
    if(nargin==0) error('required input parameter signatureDef is missing'); end;
    if(nargin<2) inInfo = struct(); end
    assert(isstruct(inInfo) && isscalar(inInfo), 'inInfo must be a scalar structure');
    defConfig = defaultConfig();
    inInfo = mergeStructs(inInfo, defConfig);
      function defaultInfo = defaultConfig()
        %Common descriptions:
          defaultInfo.plots.descriptions.annotations4cols = 'REQUIRED';
          defaultInfo.plots.descriptions.annotations4rows = 'REQUIRED';
          defaultInfo.plots.descriptions.sTitlePrefix = 'DEFAULT';
          defaultInfo.plots.descriptions.sSamplesMnemonic = 'samples';
          defaultInfo.plots.descriptions.sGenesMnemonic = 'genes';
          defaultInfo.plots.descriptions.sSamplesLabel = 'samples';
          defaultInfo.plots.descriptions.sGenesLabel = 'genes';
          defaultInfo.plots.descriptions.sSignalMnemonic4ColorbarLabel = 'log$_2$(ratio)s';
          defaultInfo.plots.descriptions.sSignalMnemonic4SecondaryColorbarLabel = 'log$_2$(ratio)s';
          defaultInfo.plots.descriptions.rowLabels = 'DEFAULT';
          defaultInfo.plots.descriptions.colLabels = 'DEFAULT';
          %Take annotations from rowIDs if provided:
            if(isfield(inInfo,'reference') && isfield(inInfo.reference,'rowIDs'))
              defaultInfo.plots.descriptions.rowLabels = inInfo.reference.rowIDs;
            end
          %Get col labels from annotations4rows for plots, if possible:
            if(isfield(inInfo.plots,'descriptions') && ~isfield(inInfo.plots.descriptions,'colLabels') && isfield(inInfo.plots.descriptions,'annotations4cols'))
              if(isfield(inInfo.plots.descriptions.annotations4cols,'GSM_ID'))
                defaultInfo.plots.descriptions.colLabels = cellfun(...
                   @(gsmid)[strrep(strrep(gsmid,'_','.'),'-','.')]...
                  ,{inInfo.plots.descriptions.annotations4cols.GSM_Id}'...
                ,'UniformOutput',false);
              elseif(isfield(inInfo.plots.descriptions.annotations4cols,'figureAnnotation'))
                defaultInfo.plots.descriptions.colLabels = {inInfo.plots.descriptions.annotations4cols.figureAnnotation}';                
              else
                error('while generating inInfo.plots.descriptions.colLabels from annotations4cols: inInfo.plots.descriptions.annotations4cols.figureAnnotation missing');
              end
          %Take annotations from colIDs if provided:
            elseif(isfield(inInfo,'reference') && isfield(inInfo.reference,'colIDs'))
              defaultInfo.plots.descriptions.colLabels = inInfo.reference.colIDs;
            end

        %Main Signature definition plot:
          defaultInfo.plots.signatureDefinitions.bEnabled = true;
          defaultInfo.plots.signatureDefinitions.bMoveCurrentGlobalSignalOffsetsToExplainedSignal = true;
          %Gene and sample selection to plot:
            defaultInfo.plots.signatureDefinitions.selection.bExcludeIfBelowMinRequireddDissectionStrength = true;
            defaultInfo.plots.signatureDefinitions.selection.bExcludeIfBelowRelCorrNoiseThresholds = true; %false; %
            defaultInfo.plots.signatureDefinitions.selection.bExcludeIfBelowW4GThreshold = false; 
            defaultInfo.plots.signatureDefinitions.selection.bExcludeIfBelowSignatureStrengthQuantileThresholds = false; %added for visualization in 202301
            defaultInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.bEnabled = false;
              defaultInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.BMask4Plot4G = true(signatureDefinition.reference.nG,1);
              defaultInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.BMask4Plot4P = true(1,signatureDefinition.reference.nP);
            defaultInfo.plots.selection.onlyUniqueGenes.bEnabled = false;
              defaultInfo.plots.selection.onlyUniqueGenes.of = 'abs(signatureDefinition.step3_regression.signedExtendedW4G)';
              defaultInfo.plots.selection.onlyUniqueGenes.topN = 40;
          %Plot optics:
            defaultInfo.plots.signatureDefinitions.sTemplate = 'inPipeline'; %default, shows detected or applied signature in pipeline original-explained=current=eigensignal+remaining.
            %Colormap:
              if(true)
                colorResolution=252;
                ramp = ((1:colorResolution/2)/(colorResolution/2))';  
                CM = flipud([[flipud(ramp) zeros(size(ramp)) zeros(size(ramp))]; [0 0 0]; [zeros(size(ramp)) zeros(size(ramp)), ramp]]); %Nature Methods  does not want red/green.
                CM = bsxfun(@times,sum(CM,2).^(1/3),(1-bsxfun(@times,1-sum(CM,2)/3,exp(-3*CM))).^1); %mehr dynamic range.
                defaultInfo.plots.colormap = CM;
                defaultInfo.plots.colormap_withNaNColor = [[.5 .5 .5];CM];
                clear CM;
              end
            defaultInfo.plots.maxRowLabelsPerPlot = 50;
            defaultInfo.plots.maxColLabelsPerPlot = 20;
            defaultInfo.plots.fontSize4GeneLabels = 6.7;
            defaultInfo.plots.color4T4P = [1 0.73 0.55];
            defaultInfo.plots.color4T4G = [0.49 0.85 0];
            defaultInfo.plots.signatureDefinitions.cLim4L2Rs = 'auto';
            defaultInfo.plots.signatureDefinitions.cLim4SDs = [-1.5,1.5];
            defaultInfo.plots.signatureDefinitions.xLim4geneStrengthsOverlay = 'auto';
            defaultInfo.plots.signatureDefinitions.yLim4sampleStrengthsOverlay = 'auto';
        %Signature table export:
          defaultInfo.export.table4eachSignatureDefinition.bEnabled = true;            
        %Simple disssection QC plot:
          defaultInfo.plots.dissectionEffectiveness.bEnabled = false; %true;

        %Export and visibility options:
          defaultInfo.plots.bVisible = 'DEFAULT'; %for all plots.
          defaultInfo.export.rootDir = 'REQUIRED';
            defaultInfo.export.fileNamePrefix4signature = '';
            defaultInfo.export.filePrefix = ''; %for all files.
          defaultInfo.export.plots.bEnabled = false; %true;
            defaultInfo.export.plots.bCloseAfterExport = true; %for all plots.
          defaultInfo.export.plots.exportFormats = {'eps'};
        %Other defaults shared with SDCM (define here again for standalone calls):
          defaultInfo.preprocessing.numericTargetPrecision = 'single';
          defaultInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot = 0; %does not show these initial brightness/cohort/lab signatures in the summary plots of every detected signature.
          defaultInfo.applicationMode.bEnabled = false;
          defaultInfo.applicationMode.externalInitialShifts = 0;
          defaultInfo.applicationMode.bApplyEveryExternalSignatureToInitialSignal = true;
          defaultInfo.applicationMode.kOffset4filenames = 0;
          defaultInfo.applicationMode.maxk_externallyValidatedSignatures = 'DEFAULT';
        %deprec/do not enable:
          defaultInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature = false; 
          defaultInfo.plots.signatureDefinitions.bShowSDL2Rs4SignatureEigensignal = false;
      end
    checkParams();
      function checkParams()
        %Export settings:
          if(~islogical(inInfo.export.plots.bEnabled) || ~isscalar(inInfo.export.plots.bEnabled))
            error('inInfo.export.plots.bEnabled must be a scalar logical.')
          end
          if(inInfo.export.plots.bEnabled || inInfo.export.table4eachSignatureDefinition.bEnabled)
            if(~ischar(inInfo.export.rootDir) || strcmp(inInfo.export.rootDir,'REQUIRED'))
              warning('inInfo.export.rootDir must be configured, if .export.plots.bEnabled; now DISABLING plot and table export');
              inInfo.export.plots.bEnabled = false;
              inInfo.export.table4eachSignatureDefinition.bEnabled = false;
            elseif(~exist(inInfo.export.rootDir,'dir'))
              mkdir(inInfo.export.rootDir);
            end
            if(~ischar(inInfo.export.filePrefix))
              error('inInfo.export.filePrefix must be a string');
            end
            if(~islogical(inInfo.export.plots.bCloseAfterExport) || length(inInfo.export.plots.bCloseAfterExport)~=1)
              error('inInfo.export.plots.bCloseAfterExport must be a logical');
            end
            if(ischar(inInfo.plots.bVisible) && strcmp(inInfo.plots.bVisible,'DEFAULT'))
              inInfo.plots.bVisible = ~inInfo.export.plots.bEnabled;
            end
            if(~islogical(inInfo.plots.bVisible) || ~isscalar(inInfo.plots.bVisible))
              error('inInfo.plots.bVisible must be a scalar logical.')
            end
            if(~iscellstr(inInfo.export.plots.exportFormats) || (inInfo.export.plots.bEnabled&&isempty(inInfo.export.plots.exportFormats)))
              error('inInfo.export.plots.exportFormats must be a cell string of file extensions for figure exports, if plot exports are enabled.')
            end
          end
        %Default export dir:
          if(~isfield(inInfo.export,'subdir4signature'))
            inInfo.export.subdir4signature = inInfo.export.rootDir;
          end 
      end
  %Initialize:
    fprintf('### Create plots for signature %d...\n', signatureDefinition.reference.k); drawnow;
    %Get context from saved reference:
      nG = signatureDefinition.reference.nG;
      nP = signatureDefinition.reference.nP;
      fn4OrderMetric4G = signatureDefinition.reference.metric4geneOrder;
      fn4OrderMetric4P = signatureDefinition.reference.metric4sampleOrder;
    %Get dissection mask:
      dissectionStrengths2DOS = zeros(nG,nP,inInfo.preprocessing.numericTargetPrecision);
        dissectionStrengths2DOS(...
          signatureDefinition.step3_regression.BNonzeroRows4signatureEigensignal...
         ,signatureDefinition.step3_regression.BNonzeroCols4signatureEigensignal...
        ) = signatureDefinition.step3_regression.signatureDissectionStrengthsOS; %we shift the negative eigensignal of the signature to dissect the signature from the signal.      
    %Support NaN color:
      bSupportNanColor = ...
           any(isnan(signatureDefinition.forPlots.overview.original.L2Rs(:)))... %missing values.
        || any(dissectionStrengths2DOS(:)==0) ... %out-of-focus values.
      ;
      fn4Colormap = iif(bSupportNanColor, 'colormap_withNaNColor', 'colormap');
    outInfo = struct();

  %% EPS) Signature definition plot:
    if(inInfo.plots.signatureDefinitions.bEnabled)
      %% Initialize:
        fprintf(' - plotting and exporting signature definition in its detection/application iteration...\n');
        %Plot template:
          bInPipelinePlot = strcmp(inInfo.plots.signatureDefinitions.sTemplate,'inPipeline'); %default, shows detected or applied signature in pipeline original-explained=current=eigensignal+remaining.
          bCenterOnlyPlot = strcmp(inInfo.plots.signatureDefinitions.sTemplate,'centerOnly'); %for paper: only show the center plot for the provided "current signal" and the provided "signature strengths" overlays.
        %Select and order genes and samples to plot:
          %Get signature orders:
            SI = signatureDefinition.step3_regression.eigenSI;
            SJ = signatureDefinition.step3_regression.eigenSJ;
          %Compute selected rows/cols:
            B4G = true(nG,1);
            B4P = true(1,nP);
            if(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowMinRequireddDissectionStrength)
              BSignatureGenes = nanmax(dissectionStrengths2DOS,[],2) >= inInfo.postprocessing.cutoffs4statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot;
              BSignatureSamples = nanmax(dissectionStrengths2DOS,[],1) >= inInfo.postprocessing.cutoffs4statistics.signatureMembership.minRequireddDissectionStrengthToShowInPlot;
              B4G = B4G & BSignatureGenes;
              B4P = B4P & BSignatureSamples;
            end
            if(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowRelCorrNoiseThresholds)
%               BSufficientRelativeCorr4G = abs(signatureDefinition.step2_finalSignatureAxes.R4G)/max(abs(signatureDefinition.step2_finalSignatureAxes.R4G)) ...
%                 >= inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4G ...
%               ;
%               BSufficientRelativeCorr4P = abs(signatureDefinition.step2_finalSignatureAxes.R4P)/max(abs(signatureDefinition.step2_finalSignatureAxes.R4P)) ...
%                 >= inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4P ...
%               ;
              BSufficientRelativeCorr4G = abs(signatureDefinition.step3_regression.R4G)/max(abs(signatureDefinition.step3_regression.R4G)) ...
                >= inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4G ...
              ;
              BSufficientRelativeCorr4P = abs(signatureDefinition.step3_regression.R4P)/max(abs(signatureDefinition.step3_regression.R4P)) ...
                >= inInfo.postprocessing.cutoffs4statistics.relCorrNoiseThresholds.ratio4P ...
              ;
              %<-dev.note: use step3_regression.R4? for plot backwards compatibility if updating correlations in step 3 was enabled.
              B4G = B4G & BSufficientRelativeCorr4G;
              B4P = B4P & BSufficientRelativeCorr4P;
            end
            if(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowW4GThreshold)
              if(inInfo.plots.signatureDefinitions.selection.minRequiredW4G>0)
                BSufficientW4G = abs(signatureDefinition.step3_regression.signedFocusedW4G) ...
                  >= inInfo.plots.signatureDefinitions.selection.minRequiredW4G ... %backwards compat visualization option used in paper.
                ;
                B4G = B4G & BSufficientW4G;
              end
              if(inInfo.plots.signatureDefinitions.selection.minRequiredW4P>0)
                BSufficientW4P = abs(signatureDefinition.step3_regression.signedFocusedW4P) ...
                  >= inInfo.plots.signatureDefinitions.selection.minRequiredW4P ... %backwards compat visualization option used in paper.
                ;
                B4P = B4P & BSufficientW4P;
              end
            end
            if(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowSignatureStrengthQuantileThresholds) %added for visualization in 202301
              X = abs(signatureDefinition.step3_regression.signatureStrengths4G);
              [~,SI4X] = sort(X);
              Q = SI2SR(SI4X)/length(X);
              BSufficientRelativeSignatureStrength4G = Q  >= inInfo.postprocessing.cutoffs4statistics.relSignatureStrengthThresholds.quantile4G;
              
              X = abs(signatureDefinition.step3_regression.signatureStrengths4P);
              [~,SI4X] = sort(X);
              Q = SI2SR(SI4X)/length(X);
              BSufficientRelativeSignatureStrength4P = Q >= inInfo.postprocessing.cutoffs4statistics.relSignatureStrengthThresholds.quantile4P;

              B4G = B4G & BSufficientRelativeSignatureStrength4G;
              B4P = B4P & BSufficientRelativeSignatureStrength4P;
            end
            if(inInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.bEnabled) %added for visualization in 202301
              assert(all(size(inInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.BMask4Plot4G)==size(B4G)), 'input validation: size(inInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.BMask4Plot4G) must be nG*1');
                B4G = B4G & inInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.BMask4Plot4G;
              assert(all(size(inInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.BMask4Plot4P)==size(B4P)), 'input validation: size(inInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.BMask4Plot4P) must be 1*nP');
                B4P = B4P & inInfo.plots.signatureDefinitions.selection.byCallerProvidedSamplesMask.BMask4Plot4P;
            end
            if(true && inInfo.plots.selection.onlyUniqueGenes.bEnabled) %internal/temp for paper: show unique genes by removing duplicates in signature strength order:
              warning('custom paper visualization criteria are active');
              %Paper/print figure: Top n unique genes by gene strengths:
                switch(inInfo.plots.selection.onlyUniqueGenes.of)
                  case 'abs(signatureDefinition.step3_regression.signatureStrengths4G)'
                    selectTopUniqueGenesBy = abs(signatureDefinition.step3_regression.signatureStrengths4G);
                      selectTopUniqueGenesBy(isnan(selectTopUniqueGenesBy)) = 0;
                  case 'abs(signatureDefinition.step3_regression.signedFocusedW4G)'
                    selectTopUniqueGenesBy = abs(signatureDefinition.step3_regression.signedFocusedW4G);
                      selectTopUniqueGenesBy(isnan(selectTopUniqueGenesBy)) = 0;
                  case 'abs(signatureDefinition.step3_regression.R4G)'
                    selectTopUniqueGenesBy = abs(signatureDefinition.step3_regression.R4G);
                      selectTopUniqueGenesBy(isnan(selectTopUniqueGenesBy)) = 0;
                  case 'abs(signatureDefinition.step3_regression.signedExtendedW4G)'
                    selectTopUniqueGenesBy = abs(signatureDefinition.step3_regression.signedExtendedW4G);
                      selectTopUniqueGenesBy(isnan(selectTopUniqueGenesBy)) = 0;
                  otherwise
                    error('inInfo.plots.selection.onlyUniqueGenes.of=%s not implemented', inInfo.plots.selection.onlyUniqueGenes.of);
                end
                [~,SI2] = sort(selectTopUniqueGenesBy,'descend');
              %Exclude duplicates in effect order:
                BIsDuplicate = false(size(B4G));
                alreadyVisitedGeneIDs = [];
                for ii2=1:length(SI2)
                  geneID = inInfo.plots.highlightGenes.measuredGeneIDs(SI2(ii2));
                  BIsDuplicate(SI2(ii2)) = isnan(geneID) || ismember(geneID, alreadyVisitedGeneIDs);
                  alreadyVisitedGeneIDs(end+1) = inInfo.plots.highlightGenes.measuredGeneIDs(SI2(ii2));
                end
              %Select top n unique genes by their gene strengths for plotting the signature "head":
                B4G = B4G & ~BIsDuplicate;
                topNtarget = inInfo.plots.selection.onlyUniqueGenes.topN; %201905; unified plot with RNA-seq vali cohort.
                iiCut = find(cumsum(B4G(SI2))>=topNtarget,1,'first');
                  if(isempty(iiCut))
                    iiCut = length(SI2);
                    warning('only %d<=%d=topNtarget genes match selection criteria', sum(B4G), topNtarget);
                  end
                  assert(isscalar(iiCut), 'not enough probesets fulfill the condition');
                B4G = B4G & ismember((1:nG)',SI2(1:iiCut));
            end
          %Remove deselected:
            SI(~B4G(SI)) = [];
            SJ(~B4P(SJ)) = [];
        %Create figure window:
          %Adaptive figure height logic:
            minFigureHeight = 637;
            stdFigureHeight = 685;
            maxFigureHeight = 1000;
            figureHeight = stdFigureHeight;
            local_maxRowLabelsPerPlot = inInfo.plots.maxRowLabelsPerPlot;
            genesToPlotableRatio = length(SI)/local_maxRowLabelsPerPlot;
              genesToPlotableRatio = min(genesToPlotableRatio, maxFigureHeight/stdFigureHeight);
              genesToPlotableRatio = max(genesToPlotableRatio, minFigureHeight/stdFigureHeight);
            if(genesToPlotableRatio>1) %allow increase of figure height for many genes, if this allows to still plot gene names:
              local_maxRowLabelsPerPlot = ceil(local_maxRowLabelsPerPlot*genesToPlotableRatio + 1);
              figureHeight = min(maxFigureHeight, figureHeight*genesToPlotableRatio);
            elseif(genesToPlotableRatio<1) %verkleinere Hoehe, bei wenig Genzeilen bis maximal auf minFigureHeight.
              figureHeight = max(minFigureHeight, figureHeight*genesToPlotableRatio);
            end
            if(figureHeight==0) warning('figureHeight==0 computed; check that display config includes at least one gene/feature.'); figureHeight=100; end
            local_maxRowLabelsPerPlot = inInfo.plots.maxRowLabelsPerPlot;
          %Workaround: TeX fonts are all wrong, if the first EPS export in the current session of Matlab occurrs/occurred while the figure was invisible:
            global plotSignature_bFirstEPSExportDoneInThisSession;
              if(isempty(plotSignature_bFirstEPSExportDoneInThisSession)) plotSignature_bFirstEPSExportDoneInThisSession = false; end;
            if(~inInfo.export.plots.bEnabled && ~inInfo.plots.bVisible)
              warning('~inInfo.export.plots.bEnabled, but also ~inInfo.plots.bVisible; this will generate hidden figures that are only accesible via handles in the overall output; note that many hidden figures can consume much RAM.');
              inInfo.plots.bVisible = true;
            end
          if(~plotSignature_bFirstEPSExportDoneInThisSession && ~inInfo.plots.bVisible)
            warning('Workaround: Setting figure to Visible before its EPS export, since this is the first EPS export of the current Matlab session and if the exported figure is invisible, Matlab initializeses TeX wrongly, causing weird/wrong fonts...');
            inInfo.plots.bVisible = true;
          end
          signatureDefinition.plots.overview.f = figure(...
             'Name', sprintf('detected signature %d', signatureDefinition.reference.k) ...
            ,'Position',[73 49 1849 figureHeight] ...
            ,'Visible', iif(inInfo.plots.bVisible,'on','off')...
            ,'Color','w' ...
          );
            plotSignature_bFirstEPSExportDoneInThisSession = true;
        %Hard-coded plot config:
          bSDs = false;
          %Always show the top+bottom gene labels:
            bPrintGeneLabels = true;
            inInfo.plotConfig.nMaxTopRowLabels = floor(local_maxRowLabelsPerPlot/2);
            inInfo.plotConfig.nMaxBottomRowLabels = floor(local_maxRowLabelsPerPlot/2);
            bPrintSampleLabels = length(SJ)<=inInfo.plots.maxColLabelsPerPlot && length(inInfo.plots.descriptions.colLabels)==nP;
          %Color limits and apply color map:
            cLim4L2Rs = inInfo.plots.signatureDefinitions.cLim4L2Rs;
              if(ischar(cLim4L2Rs)&&strcmp(cLim4L2Rs,'auto'))
                cLim4L2Rs = 0.8* [-3,3]*sqrt(weightedMean_supportsNaNs(...
                  (signatureDefinition.forPlots.overview.current.L2Rs(:)-0).^2 ... %sqrt(E[(X-mu)^2])
                ,ilsub(abs(bsxfun(@times,signatureDefinition.step3_regression.signedFocusedW4G, signatureDefinition.step3_regression.signedFocusedW4P)),':')));
                %,ilsub(abs(bsxfun(@times,signatureDefinition.step2_finalSignatureAxes.signedFocusedW4G,signatureDefinition.step2_finalSignatureAxes.signedFocusedW4P)),':')));
                assert(~any(isnan(cLim4L2Rs)), 'could not automatically derive color bar limits; define inInfo.plots.signatureDefinitions.cLim4L2Rs on caller level.');
              end
            cLim4SDs  = inInfo.plots.signatureDefinitions.cLim4SDs;
            bAdditionalColorbarOnTheLeft = bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4SignatureEigensignal;
            if(~bAdditionalColorbarOnTheLeft)
              if(bSDs||(inInfo.plots.signatureDefinitions.bShowSDL2Rs4SignatureEigensignal&&inInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature))
                sColorbarLabel = sprintf('%s (standardized for $\\mathbf{M_{%d}^S},\\mathbf{E_{%d}^%s}$%s)',inInfo.plots.descriptions.sSignalMnemonic4ColorbarLabel,signatureDefinition.reference.k,signatureDefinition.reference.k, 'S', iif(bSupportNanColor,', missing values are gray',''));
              elseif(inInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature)
                sColorbarLabel = sprintf('%s (standardized for $\\mathbf{M_{%d}^S}}$%s)',inInfo.plots.descriptions.sSignalMnemonic4ColorbarLabel,signatureDefinition.reference.k', iif(bSupportNanColor&&bInPipelinePlot,', out-of-focus and missing values are gray',''));
              elseif(inInfo.plots.signatureDefinitions.bShowSDL2Rs4SignatureEigensignal)
                sColorbarLabel = sprintf('%s (standardized for $\\mathbf{E_{%d}^%s}$%s)',inInfo.plots.descriptions.sSignalMnemonic4ColorbarLabel,signatureDefinition.reference.k, 'S', iif(bSupportNanColor&&bInPipelinePlot,', out-of-focus and missing values are gray',''));
              else
                sColorbarLabel = sprintf('%s%s',inInfo.plots.descriptions.sSignalMnemonic4ColorbarLabel, iif(bSupportNanColor&&bInPipelinePlot,' (out-of-focus and missing values are gray)',''));
              end
            else %simply print the units of the right remaining signal panel, whereto the colorbar gets assigned (we have a separate colorbar for alternate units and the panel titles show the user which unit applies):
              if(bSDs)
                sColorbarLabel = sprintf('standardized %s',inInfo.plots.descriptions.sSignalMnemonic4ColorbarLabel, iif(bSupportNanColor&&bInPipelinePlot,' (out-of-focus and missing values are gray)',''));
                sLeftColorbarLabel = sprintf('%s',inInfo.plots.descriptions.sSignalMnemonic4SecondaryColorbarLabel);
              else
                sColorbarLabel = sprintf('%s%s',inInfo.plots.descriptions.sSignalMnemonic4ColorbarLabel, iif(bSupportNanColor,' (missing values are gray)',''));
                sLeftColorbarLabel = sprintf('standardized %s',inInfo.plots.descriptions.sSignalMnemonic4SecondaryColorbarLabel);
              end
            end
            colormap(inInfo.plots.(fn4Colormap)); %sdet for all subplots in the figure
          %Configure Tex prologue and epilogue:
            latexFontSize = 13;
            %sTexPre  = '\makebox[4in][c]{\textsf{';
            %sTexPre  = '\makebox[4in][c]{\sffamily{';
            sTexPre  = '\makebox[4in][c]{\textrm{';
            sTexPost = '}}';  
      %% Plot data matrices:
        %!Subplots:
          %Original log2(ratio)s:
            if(bInPipelinePlot)
              signatureDefinition.plots.overview.aOriginalSignal = subplot(1,5,1,'Layer','bottom','FontSize',latexFontSize-3,'LineWidth',2); hold on;

              data2Plot = ilsub(...
                ilv(signatureDefinition.forPlots.overview.original.L2Rs...
                  ,@(M)iif(~bSDs,M,@()M./signatureDefinition.forPlots.overview.current.stds4L2Rs2D)...
                )...
              ,{SI,SJ});
                if(bSupportNanColor)
                  if(bSDs)
                    cLim=cLim4SDs;
                  else
                    cLim=cLim4L2Rs;
                  end
                  lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(inInfo.plots.(fn4Colormap),1);
                  data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                end
              imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                if(bSDs)
                  caxis(cLim4SDs);
                else
                  caxis(cLim4L2Rs);
                end

              h=title(signatureDefinition.plots.overview.aOriginalSignal,{
                sprintf('%s\\textbf{a)} original %ssignal $\\mathbf{M_%d%s}(I_{%d},J_{%d})$%s',sTexPre...
                  ,iif(bSDs,'standardized ','') ...
                  ,0+double(min(signatureDefinition.reference.k-1,inInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot)) ...
                  ,iif(bSDs,'^S','') ...
                  ,signatureDefinition.reference.k, signatureDefinition.reference.k ...
                ,sTexPost)
                sprintf('%sin %s signature order%s',sTexPre,iif(~inInfo.applicationMode.bEnabled,'the detected','applied'),sTexPost)
              },'Interpreter','LaTex','FontSize',latexFontSize); 
            end
          %Already explained signal / previous shifts:
            if(bInPipelinePlot)
              signatureDefinition.plots.overview.aExplainedSignal = subplot(1,5,2,'Layer','bottom','FontSize',latexFontSize-3,'LineWidth',2); hold on;

              %signalOffsets = bsxfun(@plus, signatureDefinition.step3_regression.base4G, signatureDefinition.step3_regression.base4P);
              signalOffsets = 0; 
              if(bSDs)
                data2Plot = -signatureDefinition.forPlots.overview.current.explainedSignal./signatureDefinition.forPlots.overview.current.stds4L2Rs2D;
              else
                if(inInfo.plots.signatureDefinitions.bMoveCurrentGlobalSignalOffsetsToExplainedSignal) %otherwise the bi-monotonicity can be hidden by ugly overlapping global offset stripes, if they are of the magnitude of the signal.
                  data2Plot = bsxfun(@plus,-signatureDefinition.forPlots.overview.current.explainedSignal,+signalOffsets);
                else
                  data2Plot = -signatureDefinition.forPlots.overview.current.explainedSignal;
                end
              end
                data2Plot = data2Plot(SI,SJ); %show in eigenOrder.
                if(bSupportNanColor)
                  if(bSDs)
                    cLim=cLim4SDs;
                  else
                    cLim=cLim4L2Rs;
                  end
                  lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(inInfo.plots.(fn4Colormap),1);
                  data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                end
              imagesc(data2Plot);
                if(bSDs)
                  caxis(cLim4SDs);
                else
                  caxis(cLim4L2Rs);
                end

              if(inInfo.applicationMode.bEnabled)
                if(inInfo.applicationMode.bApplyEveryExternalSignatureToInitialSignal)
                  h=title([{
                    sprintf('%s\\textbf{b)} %ssignal already dissected for%s',sTexPre,iif(bSDs,'standardized ',''),sTexPost);
                    sprintf('%sunvalidated lab signatures%s',sTexPre,sTexPost);
                  };iif(inInfo.plots.signatureDefinitions.bMoveCurrentGlobalSignalOffsetsToExplainedSignal,{
                    sprintf('%s(plus any global offsets)%s',sTexPre,sTexPost);
                  },{})],'Interpreter','LaTex','FontSize',latexFontSize); 
                else
                  h=title([{
                    sprintf('%s\\textbf{b)} %ssignal already explained by%s',sTexPre,iif(bSDs,'standardized ',''),sTexPost);
                    sprintf('%sprevious signatures $\\sum_{k=%d}^{%d}\\mathbf{E_k%s}(I_{%d},J_{%d})$%s'...
                      ,sTexPre...
                      ,1+double(inInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot), signatureDefinition.reference.k-1 ...
                        ,iif(bSDs,'^S','')...
                        ,signatureDefinition.reference.k,signatureDefinition.reference.k ...
                      ,sTexPost...
                    );
                  };iif(inInfo.plots.signatureDefinitions.bMoveCurrentGlobalSignalOffsetsToExplainedSignal,{
                    sprintf('%s(plus any global offsets)%s',sTexPre,sTexPost);
                  },{})],'Interpreter','LaTex','FontSize',latexFontSize); 
                end
              else
                h=title([{
                  sprintf('%s\\textbf{b)} %ssignal already explained by%s',sTexPre,iif(bSDs,'standardized ',''),sTexPost);
                  sprintf('%sprevious signatures $\\sum_{k=%d}^{%d}\\mathbf{E_k%s}(I_{%d},J_{%d})$%s'...
                    ,sTexPre...
                    ,1+double(inInfo.plots.signatureDefinitions.nInitialSignaturesNotToPlot), signatureDefinition.reference.k-1 ...
                      ,iif(bSDs,'^S','')...
                      ,signatureDefinition.reference.k,signatureDefinition.reference.k ...
                    ,sTexPost...
                  );
                };iif(inInfo.plots.signatureDefinitions.bMoveCurrentGlobalSignalOffsetsToExplainedSignal,{
                  sprintf('%s(plus any global offsets)%s',sTexPre,sTexPost);
                },{})],'Interpreter','LaTex','FontSize',latexFontSize); 
              end
            end
          %current log2(ratio)s:
            if(true)
              if(bInPipelinePlot)
                signatureDefinition.plots.overview.aCurrentL2Rs = subplot(1,5,3,'Layer','bottom','FontSize',latexFontSize-3,'LineWidth',2); hold on;
              end
              if(bCenterOnlyPlot)
                signatureDefinition.plots.overview.aCurrentL2Rs = axes('Layer','bottom','FontSize',latexFontSize-3,'LineWidth',2); hold on;
                set(signatureDefinition.plots.overview.f, 'Position',  [73 1 739 1003]);
              end
              if(bInPipelinePlot||bCenterOnlyPlot)
                %signalOffsets = bsxfun(@plus, signatureDefinition.step3_regression.base4G, signatureDefinition.step3_regression.base4P);
                signalOffsets = 0; 
                if(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature)
                  data2Plot = bsxfun(@plus,signatureDefinition.forPlots.overview.current.L2Rs,-signalOffsets)./signatureDefinition.forPlots.overview.current.stds4L2Rs2D;
                else
                  if(inInfo.plots.signatureDefinitions.bMoveCurrentGlobalSignalOffsetsToExplainedSignal) %otherwise the bi-monotonicity can be hidden by ugly overlapping global offset stripes, if they are of the magnitude of the signal.
                    data2Plot = bsxfun(@plus,signatureDefinition.forPlots.overview.current.L2Rs,-signalOffsets);
                  else
                    data2Plot = signatureDefinition.forPlots.overview.current.L2Rs;
                  end
                end
                  data2Plot = data2Plot(SI,SJ); %show in eigenOrder.
                  if(bSupportNanColor)
                    if(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature)
                      cLim=cLim4SDs;
                    else
                      cLim=cLim4L2Rs;
                    end
                    lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(inInfo.plots.(fn4Colormap),1);
                    data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                  end
              end
              imagesc(data2Plot); 
                if(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature)
                  caxis(cLim4SDs);
                else
                  caxis(cLim4L2Rs);
                end

              if(bInPipelinePlot)
                h=title({
                  sprintf('%s\\textbf{c)} current %ssignal $\\mathbf{M_{%d}%s}(I_{%d},J_{%d})$ in %s signature order%s',sTexPre...
                    ,iif(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature,'standardized ','')...
                    ,iif(inInfo.applicationMode.bEnabled && inInfo.applicationMode.bApplyEveryExternalSignatureToInitialSignal, 0, signatureDefinition.reference.k-1) ... %M_1 is the current signal for the iteration that yields M_2 after dissecting signature 2.
                    ,iif(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature,'^S','')...
                    ,signatureDefinition.reference.k...
                    ,signatureDefinition.reference.k...
                    ,iif(~inInfo.applicationMode.bEnabled, 'detected', 'applied') ...
                  ,sTexPost);
                  sprintf('%s($p_{\\mathrm{corr}}=%s, p_{\\mathrm{signal}}=%s$)%s'...
                    ,sTexPre...
                    ,ils(signatureDefinition.step3_regression.log10_p4Correlations,0,'1',-Inf,'0_{\mathrm{underflow}}','otherwise',sprintf(...
                       '%0.1f\\cdot10^{%d}' ...
                      ,ilv(signatureDefinition.step3_regression.log10_p4Correlations,@(log10_p)10^(log10_p-(ceil(log10_p)-1))) ...
                      ,ilv(signatureDefinition.step3_regression.log10_p4Correlations,@(log10_p)ceil(log10_p)-1) ...
                     ))...
                    ,ils(signatureDefinition.step3_regression.log10_p4SignalStrength,0,'1',-Inf,'0_{\mathrm{underflow}}','otherwise',sprintf(...
                       '%0.1f\\cdot10^{%d}' ...
                      ,ilv(signatureDefinition.step3_regression.log10_p4SignalStrength,@(log10_p)10^(log10_p-(ceil(log10_p)-1))) ...
                      ,ilv(signatureDefinition.step3_regression.log10_p4SignalStrength,@(log10_p)ceil(log10_p)-1) ...
                     ))...
                    ,sTexPost...
                  );
                  ' ';' ';' ';' ';' ';' ';' '
                },'Interpreter','LaTex','FontSize',latexFontSize,'VerticalAlignment','middle'); 
%                     ,ils(signatureDefinition.step2_finalSignatureAxes.log10_p4Correlations,0,'1',-Inf,'0_{\mathrm{underflow}}','otherwise',sprintf(...
%                        '%0.1f\\cdot10^{%d}' ...
%                       ,ilv(signatureDefinition.step2_finalSignatureAxes.log10_p4Correlations,@(log10_p)10^(log10_p-(ceil(log10_p)-1))) ...
%                       ,ilv(signatureDefinition.step2_finalSignatureAxes.log10_p4Correlations,@(log10_p)ceil(log10_p)-1) ...
%                      ))...
%                     ,ils(signatureDefinition.step2_finalSignatureAxes.log10_p4SignalStrength,0,'1',-Inf,'0_{\mathrm{underflow}}','otherwise',sprintf(...
%                        '%0.1f\\cdot10^{%d}' ...
%                       ,ilv(signatureDefinition.step2_finalSignatureAxes.log10_p4SignalStrength,@(log10_p)10^(log10_p-(ceil(log10_p)-1))) ...
%                       ,ilv(signatureDefinition.step2_finalSignatureAxes.log10_p4SignalStrength,@(log10_p)ceil(log10_p)-1) ...
%                      ))...
%                     ,sTexPost...
              end
              if(bCenterOnlyPlot)
                h=title(signatureDefinition.plots.overview.aCurrentL2Rs,{
                  sprintf('%s%ssignal $\\mathbf{M_{%d}%s}(I_{%d},J_{%d})$ in the %s signature order%s',sTexPre ...
                    ,[iif(inInfo.applicationMode.bEnabled,'validation ',''), iif(bSDs,'standardized ','')] ...
                    ,iif(inInfo.applicationMode.bEnabled && inInfo.applicationMode.bApplyEveryExternalSignatureToInitialSignal, 0, signatureDefinition.reference.k-1) ... %M_1 is the current signal for the iteration that yields M_2 after dissecting signature 2.
                    ,iif(bSDs,'^S','') ...
                    ,signatureDefinition.reference.k,signatureDefinition.reference.k ...
                    ,iif(~inInfo.applicationMode.bEnabled, 'detected', 'applied') ...
                  ,sTexPost)
                  %sprintf('%sin %s signature%s',sTexPre,iif(~inInfo.applicationMode.bEnabled,'the detected','applied'),sTexPost)
                  ' '
                  ' '
                },'Interpreter','LaTex','FontSize',latexFontSize);                   
              end
              %Write labels: sample and genes representative t cutoffs and numbers of the partitions:
                if(true)
                  h=xlabel(signatureDefinition.plots.overview.aCurrentL2Rs, [iif(bCenterOnlyPlot,cell(0,1),{});{
                    sprintf('%s$%d_{\\mathrm{left}}$+$%d_{\\mathrm{hidden}}$+$%d_{\\mathrm{right}}$ %s%s'...
                      ,sTexPre...
                      ,sum(signatureDefinition.step3_regression.(fn4OrderMetric4P)(SJ)<0) ... %202301: corrected the left/right flip of the unequal sign.
                      ,length(signatureDefinition.step3_regression.eigenSJ)-length(SJ) ... 
                      ,sum(signatureDefinition.step3_regression.(fn4OrderMetric4P)(SJ)>=0) ... 
                      ,inInfo.plots.descriptions.sSamplesMnemonic ...
                    ,sTexPost);
                  }],'Interpreter','LaTex','FontSize',latexFontSize,'Color',inInfo.plots.color4T4P*0.505,'VerticalAlignment','top');
                  positionBackup.H = findobj(signatureDefinition.plots.overview.f,'Type','axes');
                  positionBackup.Pos = get(positionBackup.H,'Position');
                    if(~iscell(positionBackup.Pos)) positionBackup.Pos = {positionBackup.Pos}; end
                  h=ylabel(signatureDefinition.plots.overview.aCurrentL2Rs, [{
                    sprintf('%s$%d_{\\mathrm{bottom}}$+$%d_{\\mathrm{hidden}}$+$%d_{\\mathrm{top}}$ %s%s' ...
                      ,sTexPre...
                      ,sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)>=0) ... 
                      ,length(signatureDefinition.step3_regression.eigenSI)-length(SI) ... 
                      ,sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)<0) ...
                      ,inInfo.plots.descriptions.sGenesMnemonic ...
                    ,sTexPost);
                  };iif(bCenterOnlyPlot,cell(6,1),{})],'Interpreter','LaTex','FontSize',latexFontSize,'Color',inInfo.plots.color4T4G*0.505); 
                  set(positionBackup.H,{'Position'},positionBackup.Pos); %hotfix against changing axes position on ylabel.
                end
            end
            %Setup center axes:
              if(bCenterOnlyPlot)
                A = signatureDefinition.plots.overview.aCurrentL2Rs; 
                l=1;
                set(A(l),'Position',[.263 .175 .55 .73]);
                set(A(l),'YDir','rev'); %'normal');
                box(A(l),'on');
                  set(A(l),'LineWidth',1);
                  set(A(l),'Layer','top');
                  set(A(l),'XColor',[0 0 1/2]/2);
                  set(A(l),'YColor',[0 0 1/2]/2);
                  set(A(l),'TickDir','out');
                  set(A(l),'TickLength',[0,0]);
                %already set individually above/caxis(A(l),cLim4L2Rs);
                axis(A(l),'tight');
                %Position Y axis
                  set(A(l),'YAxisLocation','left');
              end
          %Eigensignal * membership mask:
            if(bInPipelinePlot)
              signatureDefinition.plots.overview.aShifts = subplot(1,5,4,'Layer','bottom','FontSize',latexFontSize-3,'LineWidth',2); hold on;
              %signalOffsets = bsxfun(@plus, signatureDefinition.step3_regression.base4G, signatureDefinition.step3_regression.base4P);
              signalOffsets = 0; 
              data2Plot = ilsub(...
                ilv(-signatureDefinition.forPlots.overview.current.shifts... %minus since we show the signature's eigensignal and not the shifts for signal dissection
                  ,@(M)iif(~(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4SignatureEigensignal), M, bsxfun(@plus,M,-signalOffsets)./signatureDefinition.forPlots.overview.current.stds4L2Rs2D) ...
                )...
              ,{SI,SJ});
                if(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4CurrentSignalInSignature)
                  cLim=cLim4SDs;
                else
                  cLim=cLim4L2Rs;
                end
                if(bSupportNanColor)
                  lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(inInfo.plots.(fn4Colormap),1);
                  data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                end
                %Gray shading based on dissectionStrengths2DOS:
                  if(true)
                    I4CM = round(1 + (size(inInfo.plots.(fn4Colormap),1)-1) * max(0,min(1, (data2Plot-cLim(1))/(cLim(2)-cLim(1))) ));
                    RGB = reshape(inInfo.plots.(fn4Colormap)(I4CM,:), [size(data2Plot),3]);
                    RGB = bsxfun(@plus...
                     ,bsxfun(@times, dissectionStrengths2DOS(SI,SJ), RGB) ...
                     ,bsxfun(@times, 1-dissectionStrengths2DOS(SI,SJ), cat(3, .5, .5, .5)) ...
                    );
                  end
              %imagesc(data2Plot); 
              imagesc(RGB); 
                if(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4SignatureEigensignal)
                  caxis(signatureDefinition.plots.overview.aShifts,cLim4SDs);
                else
                  caxis(signatureDefinition.plots.overview.aShifts,cLim4L2Rs);
                end

              h=title({
                sprintf('%s\\textbf{d)} new explaining bi-monotonic %s%s',sTexPre,iif(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4SignatureEigensignal,'standardized ',''),sTexPost)
                sprintf('%ssignature eigensignal $\\mathbf{E_{%d}%s}(I_{%d},J_{%d})$%s'...
                  ,sTexPre...
                  ,signatureDefinition.reference.k...
                    ,iif(bSDs||inInfo.plots.signatureDefinitions.bShowSDL2Rs4SignatureEigensignal,'^S','')...
                    ,signatureDefinition.reference.k,signatureDefinition.reference.k...
                  ,sTexPost...
                );
                ...' ';' ';
              },'Interpreter','LaTex','FontSize',latexFontSize); 
            end
          %Remaining signal:
            if(bInPipelinePlot)
              signatureDefinition.plots.overview.aRemainingSignal = subplot(1,5,5,'Layer','bottom','FontSize',latexFontSize-3,'LineWidth',2); hold on;

              %signalOffsets = bsxfun(@plus, signatureDefinition.step3_regression.base4G, signatureDefinition.step3_regression.base4P);
              signalOffsets = 0; 
              data2Plot = ilsub(...
                ilv(signatureDefinition.forPlots.overview.current.L2Rs + signatureDefinition.forPlots.overview.current.shifts ...
                  ,@(M)iif(~bSDs, M, bsxfun(@plus,M,-signalOffsets) ./ signatureDefinition.forPlots.overview.current.stds4L2Rs2D) ...
                )...
              ,{SI,SJ});
                if(bSupportNanColor)
                  if(bSDs)
                    cLim=cLim4SDs;
                  else
                    cLim=cLim4L2Rs;
                  end
                  lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(inInfo.plots.(fn4Colormap),1);
                  data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                end
              imagesc(data2Plot);
                if(bSDs)
                  caxis(signatureDefinition.plots.overview.aRemainingSignal, cLim4SDs);
                else
                  caxis(signatureDefinition.plots.overview.aRemainingSignal, cLim4L2Rs);
                end

              h=title({
                sprintf('%s\\textbf{e)} remaining %ssignal $\\mathbf{M_{%d}%s}(I_{%d},J_{%d})$%s'...
                  ,sTexPre...
                  ,iif(bSDs,'standardized ','')...
                  ,signatureDefinition.reference.k ... %M_1 is the current signal for the iteration that yields M_2 after dissecting signature 2.
                    ,iif(bSDs,'^S',''),signatureDefinition.reference.k,signatureDefinition.reference.k ...
                  ,sTexPost...
                );
                sprintf('%s$=\\mathbf{M_{%d}%s}(I_{%d},J_{%d})-\\mathbf{E_{%d}%s}(I_{%d},J_{%d})$%s'...
                  ,sTexPre...
                  ,signatureDefinition.reference.k-1 ... %M_1 is the current signal for the iteration that yields M_2 after dissecting signature 2.
                    ,iif(bSDs,'^S',''),signatureDefinition.reference.k,signatureDefinition.reference.k ...
                  ,signatureDefinition.reference.k...
                    ,iif(bSDs,'^S','')...
                    ,signatureDefinition.reference.k,signatureDefinition.reference.k...
                  ,sTexPost...
                );
              },'Interpreter','LaTex','FontSize',latexFontSize); 
            end
          %Add colorbar for all subplots on the right:
            if(true)
              if(bInPipelinePlot)
                peerAxes = signatureDefinition.plots.overview.aRemainingSignal;
              elseif(bCenterOnlyPlot)
                peerAxes = signatureDefinition.plots.overview.aCurrentL2Rs;
              end
              backupP = get(peerAxes,'Position');
              cbh = colorbar('peer',peerAxes,'FontSize',latexFontSize-3,'LineWidth',2);
              set(peerAxes, 'Position', backupP);
              set(cbh,'Position',ilv(get(cbh,'Position'),@(P)illa(P,1,(1*P(1)+1*1)/2)));
              h=ylabel(cbh...
                ,sprintf('%s%s%s',sTexPre,sColorbarLabel,sTexPost) ... %print standardized, if the right remaining signal was displayed in standardized units and if we have an extra colorbar.
                ,'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90 ...
              ); 
                set(h,'Position',ilv(get(h,'Position'),@(P)illa(P,1,1.51*P(1))));
              %remove the gray color from the color bar:
                if(bSupportNanColor) 
                  hImg = findobj(cbh,'Type','Image');
                  if(~isempty(hImg)) %no longer possible in R2015b
                    CData = get(hImg,'CData');
                    if(CData(end)==size(inInfo.plots.(fn4Colormap),1))
                      CData(1)=2;
                    end
                    set(findobj(cbh,'Type','Image'),'CData',CData);
                  end
                end
            end
            if(bCenterOnlyPlot)
              %set(cbh,'Position',ilv(get(cbh,'Position'),@(P)[P(1),P(2)+P(4)*0/3,P(3)*0.8,P(4)/4]),'LineWidth',2);
              set(cbh,'Position',ilv(get(cbh,'Position'),@(P)[P(1),P(2)+P(4)*3/4,P(3)*0.8,P(4)/4]),'LineWidth',1);
            end
        %Labels:
          %Setup axes positions and write '-' '=' '=' '+' between subplots:
            if(bInPipelinePlot)
              %Get axes:
                A = findobj(gcf,'Type','axes');
                A = A(ilnth(2,@sort,cellfun(@(P)P(1),get(A,'Position'))));
                pDiff = 0.025;
                height = ilsub(get(A(1),'Position'),4);
                  height = height*0.91;
              for l=1:5
                set(A(l),'Position',ilv(get(A(l),'Position'),@(P)[...
                  P(1)-pDiff/2-iif(l<3,pDiff,0)-iif(l==3,pDiff*1/4,0)+iif(l>3,pDiff/2,0)...
                 ,P(2)-pDiff*1/4 ...
                 ,P(3)+pDiff...
                 ,height ...
                ])); %widen axes / no space in-between
                set(A(l),'YDir','rev'); %'normal');
                box(A(l),'on');
                  set(A(l),'LineWidth',1);
                  set(A(l),'Layer','top');
                  set(A(l),'XColor',[0 0 1/2]/2);
                  set(A(l),'YColor',[0 0 1/2]/2);
                  set(A(l),'TickDir','out');
                  set(A(l),'TickLength',[0,0]);
                %already set individually above/caxis(A(l),cLim4L2Rs);
                axis(A(l),'tight');
                %Position Y axis
                  if(l==1||l==3)
                    set(A(l),'YAxisLocation','left');
                  else
                    set(A(l),'YAxisLocation','right');
                  end
                %Hide redundant x/y tick labels:
                  if(l~=1)
                    set(A(l),'YTickLabel',{});
                  end
              end
              %decrease axes but center plot and write '-' '=' '=' '+':
                shrinkTo = 0.85;
                for l=1:5
                  if(l~=3)
                    set(A(l),'Position',ilv(get(A(l),'Position'),@(P)[P(1)+P(3)*(1-shrinkTo)/2,P(2)+P(4)*(1-shrinkTo)/2,P(3)*shrinkTo,P(4)*shrinkTo]));
                  end
                  %Write -==+ texts between plots:
                    offsetRatio = 2.5*pDiff;
                    switch(l)
                      case 2
                        halfWayBetweenSubplotsInSubplotUnits = (-sum(ilsub(get(A(l-1),'Position'),[1,3])) + ilsub(get(A(l),'Position'),1))/2 * abs(diff(xlim(A(l))))/ilsub(get(A(l),'Position'),3);
                        %text( 0-length(SJJ)*1.15*offsetRatio, length(SII)/2+1/2, '$$\mathbf{-}$$'...
                        %h = text( ilsub(xlim(A(l)),1)-halfWayBetweenSubplotsInSubplotUnits, length(SII)/2+1/2, '$$\mathbf{-}$$'...
                        %  ,'Interpreter','LaTex','Color',[0 0 0.8],'FontSize',31,'HorizontalAlignment','center','VerticalAlignment','middle'...
                        %  ,'Parent',A(l)...
                        %);
                        h = text( ilsub(xlim(A(l)),1)-halfWayBetweenSubplotsInSubplotUnits*0.91, length(SI)/2+1/2, '-'...
                          ,'Color',[0.73 0.73 0.73],'FontSize',55,'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','middle'...
                          ,'Parent',A(l)...
                        );
                        %text(length(SJJ)+length(SJJ)*1.09*offsetRatio, length(SII)/2+1/2, '$$\mathbf{=}$$'...
                        %text(  ilsub(xlim(A(l)),2)+halfWayBetweenSubplotsInSubplotUnits, length(SII)/2+1/2, '$$\mathbf{=}$$'...
                        %  ,'Interpreter','LaTex','Color',[0 0 0.8],'FontSize',31,'HorizontalAlignment','center','VerticalAlignment','middle'...
                        %  ,'Parent',A(l)...
                        %);
                        h = text(  ilsub(xlim(A(l)),2)+halfWayBetweenSubplotsInSubplotUnits*0.61, length(SI)/2+1/2, '='...
                          ,'Color',[0.73 0.73 0.73],'FontSize',37,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle'...
                          ,'Parent',A(l)...
                        );
                      case 4
                        halfWayBetweenSubplotsInSubplotUnits = (-sum(ilsub(get(A(l-1),'Position'),[1,3])) + ilsub(get(A(l),'Position'),1))/2 * abs(diff(xlim(A(l))))/ilsub(get(A(l),'Position'),3);
                        %text( 0-length(SJJ)*0.61*offsetRatio, length(SII)/2+1/2, '$$\mathbf{=}$$'...
                        %text( ilsub(xlim(A(l)),1)-halfWayBetweenSubplotsInSubplotUnits*0.391, length(SII)/2+1/2, '$$\mathbf{=}$$'...
                        %  ,'Interpreter','LaTex','Color',[0 0 0.5],'FontSize',31,'HorizontalAlignment','center','VerticalAlignment','middle'...
                        %  ,'Parent',A(l)...
                        %);
                        h = text( ilsub(xlim(A(l)),1)-halfWayBetweenSubplotsInSubplotUnits*0.391, length(SI)/2+1/2, '='...
                          ,'Color',[0.73 0.73 0.73],'FontSize',37,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle'...
                          ,'Parent',A(l)...
                        );
                        %text(length(SJJ)+length(SJJ)*1.39*offsetRatio, length(SII)/2+1/2, '$$\mathbf{+}$$'...
                        %text( ilsub(xlim(A(l)),2)+1*halfWayBetweenSubplotsInSubplotUnits*0.85, length(SII)/2+1/2, '$$\mathbf{+}$$'...
                        %  ,'Interpreter','LaTex','Color',[0 0 0.5],'FontSize',31,'HorizontalAlignment','center','VerticalAlignment','middle'...
                        %  ,'Parent',A(l)...
                        %);
                        h = text( ilsub(xlim(A(l)),2)+1*halfWayBetweenSubplotsInSubplotUnits*0.85, length(SI)/2+1/2, '+'...
                          ,'Color',[0.73 0.73 0.73],'FontSize',55,'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','middle'...
                          ,'Parent',A(l)...
                        );
                    end
                end
                set(cbh,'Position',ilv(get(A(1),'Position'),get(cbh,'Position'),@(P1,P2)illa(P2,[2,4],P1([2,4]))));
            end
          %Add another colorbar on the left, if different color scales were used for L2Rs and SDs:
            if(bAdditionalColorbarOnTheLeft && bInPipelinePlot)
              if(bSDs) %we already plotted the colorbar for standardized units on the right => show a colorbar for original units below the first panel:
                %drawnow;
                backupP = get(signatureDefinition.plots.overview.aOriginalSignal,'Position');
                cbh = colorbar('peer',signatureDefinition.plots.overview.aOriginalSignal,'FontSize',latexFontSize-3,'location','SouthOutside');
                set(signatureDefinition.plots.overview.aOriginalSignal, 'Position', backupP);
                %set(cbh,'Position',ilv(get(cbh,'Position'),@(P)illa(P,1,(1*P(1)+1*1)/2)));
                h=xlabel(cbh...
                  ,sprintf('%s%s%s',sTexPre,sLeftColorbarLabel,sTexPost) ... %print standardized, if the right remaining signal was displayed in standardized units and if we have an extra colorbar.
                  ,'Interpreter','LaTex','FontSize',latexFontSize ...
                ); 
                %remove the gray color from the color bar:
                  if(bSupportNanColor) 
                    hImg = findobj(cbh,'Type','Image');
                    if(~isempty(hImg)) %no longer possible in R2015b
                      CData = get(hImg,'CData');
                      if(CData(end)==size(inInfo.plots.(fn4Colormap),1))
                        CData(1)=2;
                      end
                      set(findobj(cbh,'Type','Image'),'CData',CData);
                    else
                      warning('TODO: >=R2015b: The new graphics engine no longer allows access to the colorbar image data. Find a workaround to not show the NaN color.');
                    end
                  end
              else %we already plotted the colorbar for original units on the right => plot a colorbar for standardized units under the eigensignal panel:
                %drawnow;
                backupP = get(signatureDefinition.plots.overview.aShifts,'Position');
                cbh = colorbar('peer',signatureDefinition.plots.overview.aShifts,'FontSize',latexFontSize-3,'location','SouthOutside');
                set(signatureDefinition.plots.overview.aShifts, 'Position', backupP);
                %set(cbh,'Position',ilv(get(cbh,'Position'),@(P)illa(P,1,(1*P(1)+1*1)/2)));
                h=xlabel(cbh...
                  ,sprintf('%s%s%s',sTexPre,sLeftColorbarLabel,sTexPost) ... %print standardized, if the right remaining signal was displayed in standardized units and if we have an extra colorbar.
                  ,'Interpreter','LaTex','FontSize',latexFontSize ...
                ); 
                %remove the gray color from the color bar:
                  if(bSupportNanColor) 
                    hImg = findobj(cbh,'Type','Image');
                    if(~isempty(hImg)) %no longer possible in R2015b
                      CData = get(hImg,'CData');
                      if(CData(end)==size(inInfo.plots.(fn4Colormap),1))
                        CData(1)=2;
                      end
                      set(findobj(cbh,'Type','Image'),'CData',CData);
                    else
                      warning('TODO: >=R2015b: The new graphics engine no longer allows access to the colorbar image data. Find a workaround to not show the NaN color.');
                    end
                  end
              end
            end
          %Signature title on the left:
            if(strcmp(inInfo.plots.descriptions.sTitlePrefix,'DEFAULT'))
              inInfo.plots.descriptions.sTitlePrefix = iif(~inInfo.applicationMode.bEnabled, 'detected signature', 'applied signature');
            end
            if(bInPipelinePlot)
              peerAxes = signatureDefinition.plots.overview.aOriginalSignal;
            elseif(bCenterOnlyPlot)
              peerAxes = signatureDefinition.plots.overview.aCurrentL2Rs;
            end
            if(~inInfo.applicationMode.bEnabled)
              h = text(ilsub(xlim(peerAxes),1)-abs(diff(xlim(peerAxes)))/13, mean(ylim(peerAxes))...
                ,[{sprintf('\\textrm{%s %d}',inInfo.plots.descriptions.sTitlePrefix, signatureDefinition.reference.k)};
                  {' ';' '}; iif(bPrintGeneLabels,{' '},{});
                 ]...
                ,'Parent',peerAxes,'Interpreter','LaTex'...
                ,'Color',[0 0 0.5],'FontWeight','bold','FontSize',latexFontSize+13,'Rotation',90 ...
                ,'HorizontalAlignment','center','VerticalAlignment','bottom' ...
              );
            else
              h=text(ilsub(xlim(peerAxes),1)-abs(diff(xlim(peerAxes)))/13, mean(ylim(peerAxes))...
                ,[{sprintf('\\textrm{%s %d/%d}',inInfo.plots.descriptions.sTitlePrefix, signatureDefinition.reference.k...
                           ,signatureDefinition.reference.maxk ..., iif(~isnan(signatureDefinition.reference.kRef), sprintf(' (\\#%d)',signatureDefinition.reference.kRef),'')...
                  )};
                  {' ';' '}; iif(bPrintGeneLabels,{' '},{});
                 ]...
                ,'Parent',peerAxes,'Interpreter','LaTex'...
                ,'Color',[0 0 0.5],'FontWeight','bold','FontSize',latexFontSize+13,'Rotation',90 ...
                ,'HorizontalAlignment','center','VerticalAlignment','bottom' ...
              );
            end
          %print gene labels, if <local_maxRowLabelsPerPlot:
            if(bInPipelinePlot)
              peerAxes = signatureDefinition.plots.overview.aOriginalSignal;
            elseif(bCenterOnlyPlot)
              peerAxes = signatureDefinition.plots.overview.aCurrentL2Rs;
            end
            if(bPrintGeneLabels)
              hold(peerAxes,'on');
              set(peerAxes,'YTick',[]);
              nSII = length(SI); nSJJ = length(SJ);
              if(~isfield(inInfo.plots,'highlightGenes'))
                shrinkFactor = 1; %0.97;
                inInfo.plotConfig.topBottomLabelsTrapezoidColor = [.5 .5 .5];
                %Compute the number of prinable rows (from available place and available top/bottom genes):
                  nMaxTopRowLabels = min(inInfo.plotConfig.nMaxTopRowLabels, sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)<0));
                  nMaxBottomRowLabels = min(inInfo.plotConfig.nMaxBottomRowLabels, sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)>=0)); 
                  %Verteile gewonnen Platz um:
                    if(nMaxTopRowLabels<inInfo.plotConfig.nMaxTopRowLabels)
                      nMaxBottomRowLabels = min(inInfo.plotConfig.nMaxBottomRowLabels + (inInfo.plotConfig.nMaxTopRowLabels-nMaxTopRowLabels), sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)>=0)); 
                    elseif(nMaxBottomRowLabels<inInfo.plotConfig.nMaxBottomRowLabels)
                      nMaxTopRowLabels = min(inInfo.plotConfig.nMaxTopRowLabels + (inInfo.plotConfig.nMaxBottomRowLabels-nMaxBottomRowLabels), sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)<0)); 
                    end
                %Get row labels:
                  rowLabels = inInfo.plots.descriptions.rowLabels;
                  %for test scenarios, highlight the simulated signature number that occurs most often:
                    bIsTestScenarioLabelling = any(cellfun(@(s)~isempty(strfind(s,'#1+')),rowLabels));
                    if(bIsTestScenarioLabelling)
                      topLabels = rowLabels(SI(1:min(nSII,nMaxTopRowLabels)));
                      bottomLabels = rowLabels(SI(nSII:-1:max([1,nSII-nMaxBottomRowLabels+1])));
                      topLabels = cellfun(@(s)strsplit(s,','), topLabels, 'UniformOutput',false);
                      bottomLabels = cellfun(@(s)strsplit(s,','), bottomLabels, 'UniformOutput',false);
                      topLabels = horzcat(topLabels{:});
                      bottomLabels = horzcat(bottomLabels{:});
                      K = 50;
                      countsPerEffect = nan(1,K);
                      for k=1:K;
                        countsPerEffect(k) = max(...
                          sum(strcmp(sprintf('#%d+',k),topLabels))-sum(strcmp(sprintf('#%d-',k),topLabels)) + sum(strcmp(sprintf('#%d-',k),bottomLabels))-sum(strcmp(sprintf('#%d+',k),bottomLabels)) ...
                         ,sum(strcmp(sprintf('#%d-',k),topLabels))-sum(strcmp(sprintf('#%d+',k),topLabels)) + sum(strcmp(sprintf('#%d+',k),bottomLabels))-sum(strcmp(sprintf('#%d-',k),bottomLabels)) ...
                        );
                      end
                      [maxCount,maxk] = max(countsPerEffect + (1:length(countsPerEffect))/1000); %+(1:length)/1000 to highlight the last effect in case of equal counts.
                      %if(maxCount>=5)
                      if(maxCount / (min(nSII,nMaxTopRowLabels) + (nSII-max([1,nSII-nMaxBottomRowLabels+1]))) >= 0.4) %at least 50% of displayed top genes must be from the same simulated effect to highlight it.
                        rowLabels = strrep(rowLabels,sprintf('#%d+',maxk),sprintf('\\color[rgb]{0 0 0.8}#%d+\\color[rgb]{0 0 0}',maxk));
                        rowLabels = strrep(rowLabels,sprintf('#%d-',maxk),sprintf('\\color[rgb]{0 0 0.8}#%d-\\color[rgb]{0 0 0}',maxk));
                      end
                    end
                %xOffset = iif(sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)<0)>nMaxTopRowLabels, 3.55, 0.75);
                xOffset = iif(2/(nMaxTopRowLabels+nMaxBottomRowLabels)*floor(shrinkFactor*nSII/2)>1, 3.61, 0.85);
                fcnRowLabelYOffsetFromRowIndex = @(i)iif(isinf(nMaxTopRowLabels),i,1+(i-1)*max(1,2/(nMaxTopRowLabels+nMaxBottomRowLabels)*shrinkFactor*floor(nSII/2)));
                  for i=1:min(nSII,nMaxTopRowLabels)
                    h=text(0.5 - nSJJ/67*xOffset,fcnRowLabelYOffsetFromRowIndex(i),[rowLabels{SI(i)}] ...
                      ,'HorizontalAlignment','right', 'VerticalAlignment','middle'...
                      ,'FontSize',inInfo.plots.fontSize4GeneLabels, 'FontName','Helvetica' ...
                      ,'Interpreter','Tex'...
                      ...,'BackgroundColor', 'none'..., 'Margin',1, 'Rotation', 0 ...
                      ,'Parent',peerAxes...
                    );
                  end
                %xOffset = iif(sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)>=0)>nMaxBottomRowLabels, 3.55, 0.75);
                xOffset = iif(2/(nMaxTopRowLabels+nMaxBottomRowLabels)*floor(shrinkFactor*nSII/2)>1, 3.61, 0.85);
                fcnRowLabelYOffsetFromRowIndex = @(i)iif(isinf(nMaxBottomRowLabels),i,nSII-(nSII-i)*max(1,2/(nMaxTopRowLabels+nMaxBottomRowLabels)*floor(shrinkFactor*nSII/2)));
                  for i=nSII:-1:max([1,nSII-nMaxBottomRowLabels+1])
                    h=text(0.5 - nSJJ/67*xOffset,fcnRowLabelYOffsetFromRowIndex(i),[rowLabels{SI(i)}] ...
                      ,'HorizontalAlignment','right', 'VerticalAlignment','middle'...
                      ,'FontSize',inInfo.plots.fontSize4GeneLabels, 'FontName','Helvetica' ...
                      ,'Interpreter','Tex'...
                      ...,'BackgroundColor', 'none'..., 'Margin',1, 'Rotation', 0 ...
                      ,'Parent',peerAxes...
                    );
                  end
                %If we only show the top/bottom gene labels, draw trapezoid patches to indicate that:
                  xLimBackup = xlim(peerAxes);
                  yLimBackup = ylim(peerAxes);
                  if(nMaxTopRowLabels>0 && 2/(nMaxTopRowLabels+nMaxBottomRowLabels)*floor(shrinkFactor*nSII/2)>1)
                    %xOffset = iif(floor(nSII/2)>nMaxTopRowLabels, nSJJ/100*1, 0);
                    hTopLabelsPatch = patch(...
                       ...[0.5, 0.5, 0-1.45*nSJJ/50-xOffset, 0-1.51*nSJJ/50-xOffset, 0-1.51*nSJJ/50-xOffset]...
                       0.5 + nSJJ/67*[0, 0, -0.2, -2.5, -3, -3]...
                      ,0.5+[0, nMaxTopRowLabels, nMaxTopRowLabels, floor(nSII*nMaxTopRowLabels/(nMaxTopRowLabels+nMaxBottomRowLabels)-1)-nSII/200, floor(nSII*nMaxTopRowLabels/(nMaxTopRowLabels+nMaxBottomRowLabels)-1)-nSII/200, -nSII/133]...
                      ,inInfo.plotConfig.topBottomLabelsTrapezoidColor...
                      ,'Parent',peerAxes,'EdgeColor','none','Clipping','off'...
                    );
                  end
                  if(nMaxBottomRowLabels>0 && 2/(nMaxTopRowLabels+nMaxBottomRowLabels)*floor(shrinkFactor*nSII/2)>1) % && sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)>=0)>nMaxBottomRowLabels)
                    %xOffset = iif(ceil(nSII/2)<nSII-nMaxBottomRowLabels+1, nSJJ/100*1, 0);
                    hBottomLabelsPatch = patch(...
                       ...[0.5, 0.5, 0-1.45*nSJJ/50-xOffset, 0-1.51*nSJJ/50-xOffset, 0-1.51*nSJJ/50-xOffset]...
                       0.5 + nSJJ/67*[0, 0, -0.2, -2.5, -3, -3]...
                      ,0.5+[nSII, nSII-nMaxBottomRowLabels, nSII-nMaxBottomRowLabels, nSII+1-floor(nSII*nMaxBottomRowLabels/(nMaxTopRowLabels+nMaxBottomRowLabels)+1)+nSII/200, nSII+1-floor(nSII*nMaxBottomRowLabels/(nMaxTopRowLabels+nMaxBottomRowLabels)+1)+nSII/200, nSII+nSII/133]...
                      ,inInfo.plotConfig.topBottomLabelsTrapezoidColor...
                      ,'Parent',peerAxes,'EdgeColor','none','Clipping','off'...
                    );
                  end
                  xlim(peerAxes,xLimBackup);
                  ylim(peerAxes,yLimBackup);
              end
              if(isfield(inInfo.plots,'highlightGenes')) %internal/temp for paper.
                xLimBackup = xlim(peerAxes);
                yLimBackup = ylim(peerAxes);

                %Schreibe alle markierten Gene auf (Ueberlappungen werden im postprocessing fuer das Paper entfernt)
                  alreadyPrintedGeneIDs = [];
                  nMaxTopRowLabels = sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)<0);
                  for i=1:nMaxTopRowLabels
                    if(ismember(inInfo.plots.highlightGenes.measuredGeneIDs(SI(i)), alreadyPrintedGeneIDs)) continue; end
                    if(~ismember(inInfo.plots.highlightGenes.measuredGeneIDs(SI(i)), inInfo.plots.highlightGenes.eligibleForHighlight)) continue; end
                    sGene = [inInfo.plots.highlightGenes.measuredGeneSymbols{SI(i)}];
                      if(isempty(sGene) || ~ischar(sGene)) continue; end
                    h = text(0.5 - nSJJ/67*0.85, i, sGene ...
                      ,'HorizontalAlignment','right', 'VerticalAlignment','middle'...
                      ,'FontSize',inInfo.plots.fontSize4GeneLabels, 'FontName','Helvetica' ...
                      ,'Interpreter','Tex'...
                      ,'Parent',peerAxes...
                    );
                    alreadyPrintedGeneIDs(end+1) = inInfo.plots.highlightGenes.measuredGeneIDs(SI(i));
                  end

                  alreadyPrintedGeneIDs = [];
                  nMaxBottomRowLabels = sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)>=0); 
                  for i=nSII:-1:max([1,nSII-nMaxBottomRowLabels+1])
                    if(ismember(inInfo.plots.highlightGenes.measuredGeneIDs(SI(i)), alreadyPrintedGeneIDs)) continue; end
                    if(~ismember(inInfo.plots.highlightGenes.measuredGeneIDs(SI(i)), inInfo.plots.highlightGenes.eligibleForHighlight)) continue; end
                    sGene = [inInfo.plots.highlightGenes.measuredGeneSymbols{SI(i)}];
                      if(isempty(sGene) || ~ischar(sGene)) continue; end
                    h = text(0.5 - nSJJ/67*0.85, i, sGene  ...
                      ,'HorizontalAlignment','right', 'VerticalAlignment','middle'...
                      ,'FontSize',inInfo.plots.fontSize4GeneLabels, 'FontName','Helvetica' ...
                      ,'Interpreter','Tex'...
                      ,'Parent',peerAxes...
                    );
                    alreadyPrintedGeneIDs(end+1) = inInfo.plots.highlightGenes.measuredGeneIDs(SI(i));
                  end

                %set better figure size:
                  %set(signatureDefinition.plots.overview.f,'Position',[73 1 100uu0*min(1,(nSJJ)/(498)) 1000*min(1,log(nSII)/log(500))]);
                  if(nSJJ>0 && nSII>0)
                    set(signatureDefinition.plots.overview.f,'Position',[73 1 1000*min(1,(nSJJ)/(498)) 1000*log(20)/log(500)*min(1,nSII/20)]);
                  else
                    warning('no genes/samples fulfil the configured plot filter');
                  end
                  xlim(peerAxes,xLimBackup);
                  ylim(peerAxes,yLimBackup);
              end
            else
              %Prevent Matlabs "x 10^4" type axes display by writing out all numbers:
                %YTicks = get(peerAxes,'YTick');
                YTicks = floor(ylim(peerAxes));
                  YTicks = round(linspace(YTicks(1),YTicks(2),iif(length(SI)>inInfo.plots.maxRowLabelsPerPlot/2, min(length(SI),5), length(SI)+1)));
                  YTicks(YTicks==0) = 1;
                  YTicks = unique(YTicks);
                set(peerAxes,'YTickMode','manual');
                set(peerAxes,'YTick',YTicks);
                if(strcmp(get(peerAxes,'YDir'),'reverse'))
                  set(peerAxes,'YTickLabel',cellfun(@(n)num2str(n), num2cell(fliplr(YTicks)), 'UniformOutput',false));
                else
                  set(peerAxes,'YTickLabel',cellfun(@(n)num2str(n), num2cell(YTicks), 'UniformOutput',false));
                end
              h=ylabel(peerAxes,sprintf('%s%s%s',sTexPre,inInfo.plots.descriptions.sGenesLabel,sTexPost),'Interpreter','LaTex','FontSize',latexFontSize);
            end
            %Show number of genes on the right of the remaining signal:
              YTicks = floor(ylim(peerAxes));
                YTicks = round(linspace(YTicks(1),YTicks(2),iif(length(SI)>inInfo.plots.maxRowLabelsPerPlot/2, min(length(SI),5), length(SI)+1)));
                YTicks(YTicks==0) = 1;
                YTicks = unique(YTicks);
              if(bInPipelinePlot)
                set(signatureDefinition.plots.overview.aRemainingSignal,'YTickMode','manual');
                set(signatureDefinition.plots.overview.aRemainingSignal,'YTick',YTicks);
                if(strcmp(get(peerAxes,'YDir'),'reverse'))
                  set(signatureDefinition.plots.overview.aRemainingSignal,'YTickLabel',cellfun(@(n)num2str(n), num2cell(fliplr(YTicks)), 'UniformOutput',false));
                else
                  set(signatureDefinition.plots.overview.aRemainingSignal,'YTickLabel',cellfun(@(n)num2str(n), num2cell(YTicks), 'UniformOutput',false));
                end
                set(signatureDefinition.plots.overview.aRemainingSignal,'YAxisLocation','right');
              end
          %print sample labels, if <local_maxColLabelsPerPlot
            if(bInPipelinePlot)
              peerAxes = signatureDefinition.plots.overview.aOriginalSignal;
            elseif(bCenterOnlyPlot)
              peerAxes = signatureDefinition.plots.overview.aCurrentL2Rs;
            end
            if(bPrintSampleLabels)
              hold(peerAxes,'on');
              set(peerAxes,'XTick',[]);
              for j=1:length(SJ)
                h=text(j,nSII+nSII/50,sprintf('%d) %s',j,inInfo.plots.descriptions.colLabels{SJ(j)}) ...
                  ,'HorizontalAlignment','right', 'VerticalAlignment','middle'...
                  ,'FontSize',inInfo.plots.fontSize4GeneLabels, 'FontName','Helvetica' ...
                  ,'Interpreter','Tex'...
                  ...,'BackgroundColor', 'none'..., 'Margin',1...
                  ,'Rotation', 90 ...
                  ,'Parent',peerAxes...
                );
              end
              %Show number of samples on all axes expect the first:
                %XTicks = get(peerAxes,'XTick');%
                XTicks = floor(xlim(peerAxes));
                  XTicks = round(linspace(XTicks(1),XTicks(2),iif(length(SJ)>inInfo.plots.maxColLabelsPerPlot/2, min(length(SJ),5), length(SJ)+1)));
                  XTicks(XTicks==0) = 1;
                  XTicks = unique(XTicks);
                if(bInPipelinePlot)
                  set(signatureDefinition.plots.overview.aExplainedSignal,'XTickMode','manual');
                  set(signatureDefinition.plots.overview.aExplainedSignal,'XTick',XTicks);
                  set(signatureDefinition.plots.overview.aExplainedSignal,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                  set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTickMode','manual');
                  set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTick',XTicks);
                  set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                  set(signatureDefinition.plots.overview.aShifts,'XTickMode','manual');
                  set(signatureDefinition.plots.overview.aShifts,'XTick',XTicks);
                  set(signatureDefinition.plots.overview.aShifts,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                  set(signatureDefinition.plots.overview.aRemainingSignal,'XTickMode','manual');
                  set(signatureDefinition.plots.overview.aRemainingSignal,'XTick',XTicks);
                  set(signatureDefinition.plots.overview.aRemainingSignal,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                end
            else
              %Prevent Matlabs "x 10^4" type axes display by writing out all numbers:
                %XTicks = get(peerAxes,'XTick');
                XTicks = floor(xlim(peerAxes));
                  XTicks = round(linspace(XTicks(1),XTicks(2),iif(length(SJ)>inInfo.plots.maxColLabelsPerPlot/2, min(length(SJ),5), length(SJ)+1)));
                  XTicks(XTicks==0) = 1;
                  XTicks = unique(XTicks);
                set(peerAxes,'XTickMode','manual');
                set(peerAxes,'XTick',XTicks);
                set(peerAxes,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                if(bInPipelinePlot)
                  set(signatureDefinition.plots.overview.aExplainedSignal,'XTickMode','manual');
                  set(signatureDefinition.plots.overview.aExplainedSignal,'XTick',XTicks);
                  set(signatureDefinition.plots.overview.aExplainedSignal,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                  if(~(isfield(inInfo.plots,'displayPredefinedClasses4aCurrentL2Rs')&&isfield(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs,'samples')))
                    set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTickMode','manual');
                    set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTick',XTicks);
                    set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                  else %leve space for the color underbars:
                    set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTickMode','manual');
                    set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTick',[]);
                    set(signatureDefinition.plots.overview.aCurrentL2Rs,'XTickLabel',{});
                  end
                  if(~(isfield(inInfo.plots,'displayPredefinedClasses4aShifts')&&isfield(inInfo.plots.displayPredefinedClasses4aShifts,'samples')))
                    set(signatureDefinition.plots.overview.aShifts,'XTickMode','manual');
                    set(signatureDefinition.plots.overview.aShifts,'XTick',XTicks);
                    set(signatureDefinition.plots.overview.aShifts,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                  else %leve space for the color underbars:
                    set(signatureDefinition.plots.overview.aShifts,'XTickMode','manual');
                    set(signatureDefinition.plots.overview.aShifts,'XTick',[]);
                    set(signatureDefinition.plots.overview.aShifts,'XTickLabel',{});
                  end
                  set(signatureDefinition.plots.overview.aRemainingSignal,'XTickMode','manual');
                  set(signatureDefinition.plots.overview.aRemainingSignal,'XTick',XTicks);
                  set(signatureDefinition.plots.overview.aRemainingSignal,'XTickLabel',cellfun(@(n)num2str(n), num2cell(XTicks), 'UniformOutput',false));
                end
              %define X axis descriptions for all subplots in the center:
                if(bInPipelinePlot)
                  h=xlabel(peerAxes ...
                    ,sprintf('%s%s%s',sTexPre,inInfo.plots.descriptions.sSamplesLabel,sTexPost)...
                  ,'Interpreter','LaTex','FontSize',latexFontSize); 
                end
            end
        %Overlays:
          %Ueberlagertes Einplotten der Anker in den zentralen Plot:
            if(true) %gene strengths (or configured gene ordering metric):
              %delete(aOverlay4RowTStats)
              aOverlay4RowTStats = axes(...
                'Parent',signatureDefinition.plots.overview.f...
               ,'Position',get(signatureDefinition.plots.overview.aCurrentL2Rs,'Position')...
               ,'FontSize',get(signatureDefinition.plots.overview.aCurrentL2Rs,'FontSize')...
               ,'Color','none'...
               ,'TickDir','in' ...
               ,'TickLength',[0.01,0.01] ...
               ,'LineWidth',0.01 ...
               ,'FontSize',latexFontSize-3 ...
               ,'XColor',inInfo.plots.color4T4G*0.505 ...
               ,'YTick',[] ...
               ,'YTickLabel',{} ...
               ,'YDir','rev' ... %to be compatible with gene order of heatmap.
               ,'XAxisLocation','top' ...
               ,'LineWidth',2 ...
              );
              hold(aOverlay4RowTStats,'on');
              if(~inInfo.applicationMode.bEnabled)
                B = signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)>=0; 
                  hTGenesPlot = plot(aOverlay4RowTStats...
                    ,signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI(B)), ilsub(1:length(SI),B) ...
                    ,'>','MarkerSize',4.3,'MarkerFaceColor',inInfo.plots.color4T4G*0.97,'LineWidth',0.5,'Color',inInfo.plots.color4T4G*0.97 ...
                  );
                B = signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)<0; 
                  hTGenesPlot = plot(aOverlay4RowTStats...
                    ,signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI(B)), ilsub(1:length(SI),B) ...
                    ,'<','MarkerSize',4.3,'MarkerFaceColor',inInfo.plots.color4T4G*0.97,'LineWidth',0.5,'Color',inInfo.plots.color4T4G*0.97 ...
                  );
              else
                set(aOverlay4RowTStats,'XTickLabel',{});
              end
                        inInfo.plots.signatureDefinitions.xLim4geneStrengthsOverlay = 'auto';
              if(ischar(inInfo.plots.signatureDefinitions.xLim4geneStrengthsOverlay)&&strcmp(inInfo.plots.signatureDefinitions.xLim4geneStrengthsOverlay,'auto'))
                axis(aOverlay4RowTStats,'tight')
                xlim(aOverlay4RowTStats,[-1.03,1.03]*max(abs(xlim(aOverlay4RowTStats))));
              else
                xlim(aOverlay4RowTStats,inInfo.plots.signatureDefinitions.xLim4geneStrengthsOverlay);
              end
              set(aOverlay4RowTStats,'XTick',ilv(get(aOverlay4RowTStats,'XTick'),@(X)iif(mod(length(X),2)==0 || length(X)<=5, X, X(2:2:end))));
              ylim(aOverlay4RowTStats,ylim(signatureDefinition.plots.overview.aCurrentL2Rs));
              %Write xlabel: info about gene t axis and sample t cutoffs:
                switch(fn4OrderMetric4G)
                  %Usual sample axis base:
                    case 'R4G'
                      h=xlabel({
                        iif(~inInfo.applicationMode.bEnabled, sprintf('%scorrelations $r^g$ with $a^g$%s',sTexPre,sTexPost), sprintf('%spreviously learnt $r^g$%s',sTexPre,sTexPost))
                      },'Interpreter','LaTex','FontSize',latexFontSize);                     
                      xlim(aOverlay4RowTStats,[-1,1]);
                    case 'R4GtimesAbsSampleAxis'
                      h=xlabel({
                        sprintf('%s$r^g$ times the gene axis $|a^g|$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize); 
                    case 'R4GtimesAbsSdSampleAxis'
                      h=xlabel({
                        sprintf('%s$r^g$ times the standardized gene axis $|a^{g,std}|$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize); 
                  %Usual gene axis base:
                    case 'geneAxis_origUnits'
                      h=xlabel({
                        sprintf('%sgene axis $a^g$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize); 
                    case 'geneAxis_sdUnits'
                      h=xlabel({
                        sprintf('%sstandardized gene axis $a^g$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize); 
                  %Compound for gene axis base:
                    case 'geneAxisTimesAbsR4G'
                      h=xlabel({
                        sprintf('%sgene axis $a^g$ times $|r^g|$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize); 
                    case 'sdSampleAxisTimesAbsR4G'
                      h=xlabel({
                        sprintf('%sstandardized gene axis $a^{g,std}$ times |$r^g|$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize); 
                    case 'signatureStrengths4G'
                      h=xlabel({
                        sprintf('%ssignature strengths for %s%s',sTexPre,inInfo.plots.descriptions.sGenesMnemonic,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize);                     
                    case 'geneStrengths4originalSignal'
                      h=xlabel({
                        sprintf('%ssignature gene strengths for %s%s',sTexPre,inInfo.plots.descriptions.sGenesMnemonic,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize);                     
                    case 'signatureStrengths4GtimesAbsR4G'
                      h=xlabel({
                        sprintf('%ssignature strengths for %s times $|r^g|$%s',sTexPre,inInfo.plots.descriptions.sGenesMnemonic,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize);                     
                  otherwise 
                    warning('overlay plot not implemented for fn4OrderMetric4G=%s', fn4OrderMetric4G);
                    h=xlabel({
                      sprintf('%sunknown gene order metric %s%s',sTexPre,strrep(fn4OrderMetric4G,'_','\_'),sTexPost)
                    },'Interpreter','LaTex','FontSize',latexFontSize); 
                end
              %Plotte t-Nulllinie:
                if(bCenterOnlyPlot)
                  yWhereXIsRepT4G = sum(signatureDefinition.step3_regression.(fn4OrderMetric4G)(SI)<0) + 1/2; %+1/2 to show line in-between genes (important if there are only few lines in the plot)
                    line(xlim(aOverlay4RowTStats),[yWhereXIsRepT4G,yWhereXIsRepT4G],'LineStyle','-','Color',inInfo.plots.color4T4G*0.97,'LineWidth',1.5,'Parent',aOverlay4RowTStats);
                end
              %Plotte Trennlinien zwischen Signatur-Genhaelften entsprechend aktiver Genselektion:
                if(bInPipelinePlot)
                  fn4Cutlines = '';
                    if(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowMinRequireddDissectionStrength && inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowRelCorrNoiseThresholds)
                      fn4Cutlines = 'signatureMembershipAndRelCorrGtNoiseThreshold';
                    elseif(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowMinRequireddDissectionStrength)
                      fn4Cutlines = 'signatureMembership';
                    elseif(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowRelCorrNoiseThresholds)
                      fn4Cutlines = 'relCorrGtNoiseThreshold';
                    end
                  if(~isempty(fn4Cutlines))
                    yWhereXIsRepT4G = signatureDefinition.statistics.(fn4Cutlines).nRepAntiregGenes + 1/2; %+1/2 to show line in-between genes (important if there are only few lines in the plot)
                      if(~isempty(yWhereXIsRepT4G) && yWhereXIsRepT4G>1 && yWhereXIsRepT4G<length(SI))
                        line(xlim(aOverlay4RowTStats),[yWhereXIsRepT4G,yWhereXIsRepT4G],'LineStyle','-','Color',inInfo.plots.color4T4G*0.97,'LineWidth',1.5,'Parent',aOverlay4RowTStats);
                      end
                    yWhereXIsRepT4G = length(SI)+1 - signatureDefinition.statistics.(fn4Cutlines).nRepCoregGenes - 1/2; %-1/2 to show line in-between genes (important if there are only few lines in the plot)
                      if(~isempty(yWhereXIsRepT4G) && yWhereXIsRepT4G>1 && yWhereXIsRepT4G<length(SI))
                        line(xlim(aOverlay4RowTStats),[yWhereXIsRepT4G,yWhereXIsRepT4G],'LineStyle','-','Color',inInfo.plots.color4T4G*0.97,'LineWidth',1.5,'Parent',aOverlay4RowTStats);
                      end
                  else
                    warning('cannot draw cut lines between signature gene halves as this has not been implemented yet for the the current gene selection configuration');
                  end
                end
            end
            if(true) %sample strenghts (or configured sample ordering metric):
              %delete(aOverlay4ColTStats)
              aOverlay4ColTStats = axes(...
                'Parent',signatureDefinition.plots.overview.f...
               ,'Position',get(signatureDefinition.plots.overview.aCurrentL2Rs,'Position')...
               ,'FontSize',get(signatureDefinition.plots.overview.aCurrentL2Rs,'FontSize')...
               ,'Color','none'...
               ,'TickDir','in' ...
               ,'TickLength',[0.01,0.01] ...
               ,'LineWidth',0.01 ...
               ,'FontSize',latexFontSize-3 ...
               ,'YColor',inInfo.plots.color4T4P*0.505 ...
               ,'XTickLabel',{} ...
               ,'YAxisLocation','right' ...
               ,'LineWidth',2 ...
              );
              hold(aOverlay4ColTStats,'on');
              B = signatureDefinition.step3_regression.(fn4OrderMetric4P)(SJ)>=0; 
                hTSamplesPlot = plot(aOverlay4ColTStats...
                  ,ilsub(1:length(SJ),B), ilsub(signatureDefinition.step3_regression.(fn4OrderMetric4P)(SJ),B) ...
                  ,'v','MarkerSize',4.3,'MarkerFaceColor',inInfo.plots.color4T4P*0.883,'LineWidth',0.5,'Color',inInfo.plots.color4T4P*0.883 ...
                );
              B = signatureDefinition.step3_regression.(fn4OrderMetric4P)(SJ)<0; 
                hTSamplesPlot = plot(aOverlay4ColTStats...
                  ,ilsub(1:length(SJ),B), ilsub(signatureDefinition.step3_regression.(fn4OrderMetric4P)(SJ),B) ...
                  ,'^','MarkerSize',4.3,'MarkerFaceColor',inInfo.plots.color4T4P*0.883,'LineWidth',0.5,'Color',inInfo.plots.color4T4P*0.883 ...
                );
              if(ischar(inInfo.plots.signatureDefinitions.yLim4sampleStrengthsOverlay)&&strcmp(inInfo.plots.signatureDefinitions.yLim4sampleStrengthsOverlay,'auto'))
                axis(aOverlay4ColTStats,'tight');
                ylim(aOverlay4ColTStats,[-1.03,1.03]*max(abs(ylim(aOverlay4ColTStats))));
              else
                ylim(aOverlay4ColTStats,inInfo.plots.signatureDefinitions.yLim4sampleStrengthsOverlay);
              end
                set(aOverlay4ColTStats,'YTick',ilv(get(aOverlay4ColTStats,'YTick'),@(X)iif(mod(length(X),2)==0 || length(X)<=5, X, X(2:2:end))));
              xlim(aOverlay4ColTStats,xlim(signatureDefinition.plots.overview.aCurrentL2Rs));
              %Write ylabel: info about sample t axis and gene t cutoffs:
                switch(fn4OrderMetric4P)
                  %Usual sample axis base:
                    case 'sampleAxis_origUnits'
                      h=ylabel({
                        ''
                        sprintf('%ssample axis $a^s$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                    case 'sampleAxis_sdUnits'
                      h=ylabel({
                        ''
                        sprintf('%sstandardized sample axis $a^s$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                  %Usual gene axis base:
                    case 'R4P'
                      h=ylabel({
                        sprintf('%scorrelations $r^s$ with $a^s$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                      ylim(aOverlay4ColTStats,[-1,1]);
                    case 'R4PtimesAbsGeneAxis'
                      h=ylabel({
                        sprintf('%s$r^s$ times the sample axis $|a^s|$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                    case 'R4PtimesAbsSdGeneAxis'
                      h=ylabel({
                        sprintf('%s$r^s$ times the standardized sample axis $|a^s|$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                  %Compound for gene axis base:
                    case 'sampleAxisTimesAbsR4P'
                      h=ylabel({
                        sprintf('%ssample axis $a^s$ times $|r^s|$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                    case 'sdGeneAxisTimesAbsR4P'
                      h=ylabel({
                        sprintf('%sstandardized sample axis $a^{s,std}$ times $|r^s|$%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                    case 'signatureStrengths4P'
                      h=ylabel({
                        sprintf('%ssignature strengths in %s%s',sTexPre,inInfo.plots.descriptions.sSamplesMnemonic,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                    case 'sampleStrengths4originalSignal'
                      h=ylabel({
                        sprintf('%ssignature sample strengths in %s%s',sTexPre,inInfo.plots.descriptions.sSamplesMnemonic,sTexPost)
                        sprintf('%sin the validation cohort%s',sTexPre,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                    case 'signatureStrengths4PtimesAbsR4P'
                      h=ylabel({
                        sprintf('%ssignature strengths for %s times $|r^s|$%s',sTexPre,inInfo.plots.descriptions.sSamplesMnemonic,sTexPost)
                      },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                  otherwise 
                    warning('overlay plot not implemented for fn4OrderMetric4P=%s', fn4OrderMetric4P);
                    h=ylabel({
                      sprintf('%sunknown sample order metric: %s%s',sTexPre,strrep(fn4OrderMetric4P,'_','\_'),sTexPost)
                    },'Interpreter','LaTex','FontSize',latexFontSize,'Rotation',-90);
                end
                  set(h,'Position',ilv(get(h,'Position'),@(P)illa(P,1,P(1)+length(SJ)*8.5/100)));
              %Plotte t-Nulllinie:
                if(bCenterOnlyPlot)
%                     line(xlim(aOverlay4ColTStats),[0,0],'LineStyle',':','Color',inInfo.plots.color4T4P*0.883,'LineWidth',1.5,'Parent',aOverlay4ColTStats);
                  xWhereYIsRepT4P = sum(signatureDefinition.step3_regression.(fn4OrderMetric4P)(SJ)<0) + 1/2; %+1/2 to show line in-between genes (important if there are only few lines in the plot)
                    line([xWhereYIsRepT4P,xWhereYIsRepT4P],ylim(aOverlay4ColTStats),'LineStyle','-','Color',inInfo.plots.color4T4P*0.883,'LineWidth',1.5,'Parent',aOverlay4ColTStats);
                end
              %Plotte Trennlinien zwischen Signatur-Samplehaelften entsprechend aktiver Genselektion:
                if(bInPipelinePlot)
                  fn4Cutlines = '';
                    if(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowMinRequireddDissectionStrength && inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowRelCorrNoiseThresholds)
                      fn4Cutlines = 'signatureMembershipAndRelCorrGtNoiseThreshold';
                    elseif(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowMinRequireddDissectionStrength)
                      fn4Cutlines = 'signatureMembership';
                    elseif(inInfo.plots.signatureDefinitions.selection.bExcludeIfBelowRelCorrNoiseThresholds)
                      fn4Cutlines = 'relCorrGtNoiseThreshold';
                    end
                  if(~isempty(fn4Cutlines))
                    xWhereYIsRepT4P = signatureDefinition.statistics.(fn4Cutlines).nRepAntiAffectedSamples + 1/2; %+1/2 to show line in-between genes (important if there are only few lines in the plot)
                      if(~isempty(xWhereYIsRepT4P) && xWhereYIsRepT4P>1 && xWhereYIsRepT4P<length(SJ))
                        line([xWhereYIsRepT4P,xWhereYIsRepT4P],ylim(aOverlay4ColTStats),'LineStyle','-','Color',inInfo.plots.color4T4P*0.883,'LineWidth',1.5,'Parent',aOverlay4ColTStats);
                      end
                    xWhereYIsRepT4P = length(SJ)+1 - signatureDefinition.statistics.(fn4Cutlines).nRepAffectedSamples - 1/2; %-1/2 to show line in-between genes (important if there are only few lines in the plot)
                      if(~isempty(xWhereYIsRepT4P) && xWhereYIsRepT4P>1 && xWhereYIsRepT4P<length(SJ))
                        line([xWhereYIsRepT4P,xWhereYIsRepT4P],ylim(aOverlay4ColTStats),'LineStyle','-','Color',inInfo.plots.color4T4P*0.883,'LineWidth',1.5,'Parent',aOverlay4ColTStats);
                      end
                  else
                    warning('cannot draw cut lines between signature gene halves as this has not been implemented yet for the the current gene selection configuration');
                  end
                end
            end
            %<-Add zoom handler to synchronize zoom with the underlying signatureDefinition.plots.overview.aCurrentL2Rs and aOverlay4RowTStats:
              addPanAndZoomHandler(gcf, aOverlay4ColTStats, aOverlay4RowTStats, signatureDefinition.plots.overview.aCurrentL2Rs);
              zoom(gcf,'on');
          %Add color underbar and legend for sample classes if provided:
            if(isfield(inInfo.plots,'displayPredefinedClasses4aCurrentL2Rs')&&isfield(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs,'samples') ...
             || bCenterOnlyPlot && isfield(inInfo.plots,'displayPredefinedClasses4aShifts')&&isfield(inInfo.plots.displayPredefinedClasses4aShifts,'samples') ...
            ) %<-Note: usually, not displayPredefinedClasses4aCurrentL2Rs, but displayPredefinedClasses4aShifts below is configured for the fourth heatmap, not the third.
              if(~(isfield(inInfo.plots,'displayPredefinedClasses4aCurrentL2Rs')&&isfield(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs,'samples')) ...
                && bCenterOnlyPlot && isfield(inInfo.plots,'displayPredefinedClasses4aShifts')&&isfield(inInfo.plots.displayPredefinedClasses4aShifts,'samples') ...
              )
                inInfo.plots.displayPredefinedClasses4aCurrentL2Rs = inInfo.plots.displayPredefinedClasses4aShifts;
              end
              for cil=1:length(inInfo.plots.displayPredefinedClasses4aShifts)
                assert(length(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples)==nP, 'length(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples)~=nP');
                aUnderbarTarget = signatureDefinition.plots.overview.aCurrentL2Rs;
                %202301: support for value->color function handles:
                  if(isfield(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil),'fcn4SampleValue2Color') && ~isempty(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).fcn4SampleValue2Color))
                    assert(isa(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).fcn4SampleValue2Color, 'function_handle'), 'if specified, inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(%d).fcn4SampleValue2Color must be a function_handle',cil);
                    orderedSampleClassColors = arrayfun(...
                      inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).fcn4SampleValue2Color...
                     ,inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples(SJ)' ...
                    ,'Uni',false);
                      orderedSampleClassColors = vertcat(orderedSampleClassColors{:});
                      assert(size(orderedSampleClassColors,2)==3, 'expected to get RGB triples from user function inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(%d).fcn4SampleValue2Color',cil);
                    %For the legend, generate some (min,medium,max) pseudo-categories:
                      uniqueSampleClasses = {'min','median','max'};
                      colors4SampleClasses = [
                        inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).fcn4SampleValue2Color(nanmin(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples(SJ)'));
                        inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).fcn4SampleValue2Color(nanmedian(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples(SJ)'));
                        inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).fcn4SampleValue2Color(nanmax(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples(SJ)'));
                      ];
                %Default: Categorical coloring, using optional .CM4SampleClasses:
                  else
                    uniqueSampleClasses = unique(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples);
                    %Define class colors:
                      CM4SampleClasses = [0 .75 1; 1 0 1; 1 .75 0; 0.25 1 0; 0.25 0.75 1; 0.75 0.25 1; 1 0 0.5; 0.5 1 0; 1 0.5 0];
                      %colors4SampleClasses = lines(length(uniqueSampleClasses));
                      if(size(CM4SampleClasses,1)<length(uniqueSampleClasses))
                        CM4SampleClasses = [CM4SampleClasses; lines(length(uniqueSampleClasses))];
                      end
                      CM4SampleClasses = min(1,sqrt(0.2+CM4SampleClasses));
                      if(isfield(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil),'CM4SampleClasses'))
                        CM4SampleClasses = inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).CM4SampleClasses;
                      end
                      colors4SampleClasses = CM4SampleClasses(1:length(uniqueSampleClasses),:);
                    orderedSampleClassColors = joinAndFill(uniqueSampleClasses(:), colors4SampleClasses, inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples(SJ)');
                      orderedSampleClassColors = reshape(orderedSampleClassColors,1,length(SJ),3);
                  end
                xRange = xlim(aUnderbarTarget);
                  xRange(1) = xRange(1)+1/2;
                  xRange(end) = xRange(end)-1/2;
                YLim = ylim(aUnderbarTarget);
                yRange = YLim(2)+[1/500,21.7/500*iif(bInPipelinePlot,maxFigureHeight/figureHeight,1)]*abs(-YLim(1)+YLim(2));
                  yRange = yRange + (cil-1+2.5)*abs(diff(yRange));
                hSampleClasses = image(...
                   'CData',repmat(orderedSampleClassColors,16,1,1) ... %prevent anti-aliasing in PDF viewers.
                  ,'Parent',aUnderbarTarget ...
                  ,'XData',xRange ...
                  ,'YData',yRange ...
                  ,'Clipping','off'...
                );
                  ylim(aUnderbarTarget, YLim+1); %workaround, as ylim thinks that it was not changed by the image command, but it was
                  ylim(aUnderbarTarget, YLim); %resote ylim, if Matlab moved it due to the added image.
                if(~isfield(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil),'legendHeader4samples'))
                  %Backwards compatibility (renamed a field name after a typo...):
                    if(isfield(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil),'legentHeader4samples')) 
                      inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).legendHeader4samples = inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).legentHeader4samples;
                  %Default:
                    else
                      inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).legendHeader4samples = 'original cohorts:';
                    end
                end
                bMoveFullLegendsBelowTheHeatmaps = true; %new flavor 202301. TODO: connect to public options.
                if(~bMoveFullLegendsBelowTheHeatmaps) %old standard:
                  hLegend4SampleClasses = text(...
                     iif(mod(cil,2)==0,min(xRange)-abs(diff(xlim(aUnderbarTarget)))/25, max(xRange)+abs(diff(xlim(aUnderbarTarget)))/25) ...
                    ,mean(yRange)...
                    ,[{['\color[rgb]{1 1 1}',sprintf('%d) ',cil),inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).legendHeader4samples]}...
                      ,cellfun(...
                          @(s,c)sprintf('\\color[rgb]{%0.2f %0.2f %0.2f}%s',c(1),c(2),c(3),s)...
                         ,uniqueSampleClasses...
                         ,num2cell(colors4SampleClasses,2)'...
                       ,'UniformOutput',false)...
                     ] ...
                    ,'Parent',aUnderbarTarget ...
                    ,'HorizontalAlignment',iif(mod(cil,2)==0,'right','left'), 'VerticalAlignment',iif(mod(cil,3)==1,'middle',iif(mod(cil,3)==2,'middle','top')) ...
                    ,'FontSize',latexFontSize-5, 'FontName','Helvetica', 'FontWeight','bold' ...
                    ,'Interpreter','Tex'...
                    ,'Clipping','off'...
                    ,'BackgroundColor',[.67 .67 .67] ... [.5 .5 .5] ...
                    ,'Margin',2.5 ...
                  );
                end
                if(bMoveFullLegendsBelowTheHeatmaps)
                  hLegend4SampleClasses_abbrevOnly = text(...
                     ...iif(mod(cil,2)==0,min(xRange)-abs(diff(xlim(aUnderbarTarget)))/25, max(xRange)+abs(diff(xlim(aUnderbarTarget)))/25) ...
                     max(xRange)+abs(diff(xlim(aUnderbarTarget)))/250 ...
                    ,mean(yRange)...
                    ,[{['\color[rgb]{0 0 0}',sprintf('%d) ',cil),inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).legendHeader4samples]}...
                     ] ...
                    ,'Parent',aUnderbarTarget ...
                    ,'HorizontalAlignment',iif(mod(cil,2)==0,'right','left'), 'VerticalAlignment',iif(mod(cil,3)==1,'middle',iif(mod(cil,3)==2,'middle','top')) ...
                    ,'FontSize',latexFontSize-5, 'FontName','Helvetica', 'FontWeight','bold' ...
                    ,'Interpreter','Tex'...
                    ,'Clipping','off'...
                    ,'BackgroundColor','none' ...[.67 .67 .67] ... [.5 .5 .5] ...
                    ,'Margin',2.5 ...
                  );
                  hLegend4SampleClasses = text(...
                     min(xRange) + (cil-length(inInfo.plots.displayPredefinedClasses4aShifts)/2) * abs(diff(xlim(aUnderbarTarget)))/3 ...
                    ,YLim(2)+[1/500,21.7/500*iif(bInPipelinePlot,maxFigureHeight/figureHeight,1)]*abs(-YLim(1)+YLim(2)) + (length(inInfo.plots.displayPredefinedClasses4aShifts)-1+1+2.5)*abs(diff(yRange)) ...
                    ,[{['\color[rgb]{1 1 1}',sprintf('%d) ',cil),inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).legendHeader4samples]}...
                      ,cellfun(...
                          @(s,c)sprintf('\\color[rgb]{%0.2f %0.2f %0.2f}%s',c(1),c(2),c(3),s)...
                         ,uniqueSampleClasses...
                         ,num2cell(colors4SampleClasses,2)'...
                       ,'UniformOutput',false)...
                     ] ...
                    ,'Parent',aUnderbarTarget ...
                    ,'HorizontalAlignment','left', 'VerticalAlignment','top' ...
                    ,'FontSize',latexFontSize-5, 'FontName','Helvetica', 'FontWeight','bold' ...
                    ,'Interpreter','Tex'...
                    ,'Clipping','off'...
                    ,'BackgroundColor',[.67 .67 .67] ... [.5 .5 .5] ...
                    ,'Margin',2.5 ...
                  );
                end
              end
            end
            if(bInPipelinePlot && isfield(inInfo.plots,'displayPredefinedClasses4aShifts')&&isfield(inInfo.plots.displayPredefinedClasses4aShifts,'samples'))
              aUnderbarTarget = signatureDefinition.plots.overview.aShifts;
              XLim = xlim(aUnderbarTarget);
              for cil=1:length(inInfo.plots.displayPredefinedClasses4aShifts)
                assert(length(inInfo.plots.displayPredefinedClasses4aShifts(cil).samples)==nP, 'length(inInfo.plots.displayPredefinedClasses4aShifts(%d).samples)~=nP', cil);
                %202301: support for value->color function handles:
                  if(isfield(inInfo.plots.displayPredefinedClasses4aShifts(cil),'fcn4SampleValue2Color') && ~isempty(inInfo.plots.displayPredefinedClasses4aShifts(cil).fcn4SampleValue2Color))
                    assert(isa(inInfo.plots.displayPredefinedClasses4aShifts(cil).fcn4SampleValue2Color, 'function_handle'), 'if specified, inInfo.plots.displayPredefinedClasses4aShifts(%d).fcn4SampleValue2Color must be a function_handle',cil);
                    orderedSampleClassColors = arrayfun(...
                      inInfo.plots.displayPredefinedClasses4aShifts(cil).fcn4SampleValue2Color ...
                     ,inInfo.plots.displayPredefinedClasses4aShifts(cil).samples(SJ)' ...
                    ,'Uni',false);
                      orderedSampleClassColors = vertcat(orderedSampleClassColors{:});
                      assert(size(orderedSampleClassColors,2)==3, 'expected to get RGB triples from user function inInfo.plots.displayPredefinedClasses4aShifts(%d).fcn4SampleValue2Color',cil);
                      orderedSampleClassColors = reshape(orderedSampleClassColors,1,length(SJ),3);
                    %For the legend, generate some (min,medium,max) pseudo-categories:
                      uniqueSampleClasses = {'5%tile','25%tile','median','75%tile','95%tile'};
                      X = inInfo.plots.displayPredefinedClasses4aShifts(cil).samples;
                        X = X(~isnan(X));
                      colors4SampleClasses = [
                        inInfo.plots.displayPredefinedClasses4aShifts(cil).fcn4SampleValue2Color(quantile(X,0.05));
                        inInfo.plots.displayPredefinedClasses4aShifts(cil).fcn4SampleValue2Color(quantile(X,0.25));
                        inInfo.plots.displayPredefinedClasses4aShifts(cil).fcn4SampleValue2Color(quantile(X,0.50));
                        inInfo.plots.displayPredefinedClasses4aShifts(cil).fcn4SampleValue2Color(quantile(X,0.75));
                        inInfo.plots.displayPredefinedClasses4aShifts(cil).fcn4SampleValue2Color(quantile(X,0.95));
                      ];
                %Default: Categorical coloring, using optional .CM4SampleClasses:
                  else
                    uniqueSampleClasses = unique(inInfo.plots.displayPredefinedClasses4aShifts(cil).samples);
                    %Define class colors:
                      CM4SampleClasses = [0 .75 1; 1 0 1; 1 .75 0; 0.25 1 0; 0.25 0.75 1; 0.75 0.25 1; 1 0 0.5; 0.5 1 0; 1 0.5 0];
                      %colors4SampleClasses = lines(length(uniqueSampleClasses));
                      if(size(CM4SampleClasses,1)<length(uniqueSampleClasses))
                        CM4SampleClasses = [CM4SampleClasses; lines(length(uniqueSampleClasses))];
                      end
                      CM4SampleClasses = min(1,sqrt(0.2+CM4SampleClasses));
                      if(isfield(inInfo.plots.displayPredefinedClasses4aShifts(cil),'CM4SampleClasses'))
                        CM4SampleClasses = inInfo.plots.displayPredefinedClasses4aShifts(cil).CM4SampleClasses;
                      end
                      colors4SampleClasses = CM4SampleClasses(1:length(uniqueSampleClasses),:);
                    orderedSampleClassColors = joinAndFill(uniqueSampleClasses(:), colors4SampleClasses, inInfo.plots.displayPredefinedClasses4aShifts(cil).samples(SJ)');
                      orderedSampleClassColors = reshape(orderedSampleClassColors,1,length(SJ),3);
                  end
                xRange = XLim;
                  xRange(1) = xRange(1)+1/2;
                  xRange(end) = xRange(end)-1/2;
                YLim = ylim(aUnderbarTarget);
                %yRange = YLim(2)+[2/500,15.7/500*maxFigureHeight/figureHeight]*abs(-YLim(1)+YLim(2));
                yRange = YLim(2)+[2/500,12.7/500*maxFigureHeight/figureHeight]*abs(-YLim(1)+YLim(2));
                  yRange = yRange + (cil-1)*abs(diff(yRange));
                hSampleClasses = image(...
                   'CData',repmat(orderedSampleClassColors,16,1,1) ... %prevent anti-aliasing in PDF viewers.
                  ,'Parent',aUnderbarTarget ...
                  ,'XData',xRange ...
                  ,'YData',yRange ...
                  ,'Clipping','off'...
                );
                  ylim(aUnderbarTarget, YLim+1); %workaround, as ylim thinks that it was not changed by the image command, but it was
                  ylim(aUnderbarTarget, YLim); %restore ylim, if Matlab moved it due to the added image.
                if(~isfield(inInfo.plots.displayPredefinedClasses4aShifts(cil),'legendHeader4samples'))
                  inInfo.plots.displayPredefinedClasses4aShifts(cil).legendHeader4samples = 'original cohorts:';
                end
                bMoveFullLegendsBelowTheHeatmaps = true; %new flavor 202301. TODO: connect to public options.
                if(~bMoveFullLegendsBelowTheHeatmaps) %old standard:
                  hLegend4SampleClasses = text(...
                     iif(mod(cil,2)==1,min(xRange)-abs(diff(xlim(aUnderbarTarget)))/25, max(xRange)+abs(diff(xlim(aUnderbarTarget)))/25) ...
                    ,mean(yRange)...
                    ,[{['\color[rgb]{1 1 1}',inInfo.plots.displayPredefinedClasses4aShifts(cil).legendHeader4samples]}...
                     ,cellfun(...
                        @(s,c)sprintf('\\color[rgb]{%0.2f %0.2f %0.2f}%s',c(1),c(2),c(3),s)...
                       ,uniqueSampleClasses...
                       ,num2cell(colors4SampleClasses,2)'...
                      ,'UniformOutput',false)...
                     ] ...
                    ,'Parent',aUnderbarTarget ...
                    ...,'HorizontalAlignment',iif(mod(cil,2)==1,'right','left'), 'VerticalAlignment',iif(mod(cil,3)==1,'bottom',iif(mod(cil,3)==2,'middle','top')) ...
                    ,'HorizontalAlignment',iif(mod(cil,2)==1,'right','left'), 'VerticalAlignment',iif(mod(cil,3)==1,'middle',iif(mod(cil,3)==2,'middle','top')) ...
                    ,'FontSize',latexFontSize-5, 'FontName','Helvetica', 'FontWeight','bold' ...
                    ,'Interpreter','Tex'...
                    ,'Clipping','off'...
                    ,'BackgroundColor',[.67 .67 .67] ... [.5 .5 .5] ...
                    ,'Margin',2.5 ...
                  );
                end
                if(bMoveFullLegendsBelowTheHeatmaps)
                  %delete(hLegend4SampleClasses_abbrevOnly)
                  hLegend4SampleClasses_abbrevOnly = text(...
                     ...iif(mod(cil,2)==0,min(xRange)-abs(diff(xlim(aUnderbarTarget)))/25, max(xRange)+abs(diff(xlim(aUnderbarTarget)))/25) ...
                     max(xRange)+abs(diff(xlim(aUnderbarTarget)))/250 ...
                    ,mean(yRange)...
                    ,[{['\color[rgb]{0 0 0}',sprintf('%d) ',cil),inInfo.plots.displayPredefinedClasses4aShifts(cil).legendHeader4samples]}...
                     ] ...
                    ,'Parent',aUnderbarTarget ...
                    ...,'HorizontalAlignment',iif(mod(cil,2)==0,'right','left'), 'VerticalAlignment',iif(mod(cil,3)==1,'middle',iif(mod(cil,3)==2,'middle','top')) ...
                    ,'HorizontalAlignment','left', 'VerticalAlignment','middle' ...
                    ,'FontSize',latexFontSize-5, 'FontName','Helvetica', 'FontWeight','bold' ...
                    ,'Interpreter','Tex'...
                    ,'Clipping','off'...
                    ,'BackgroundColor','none'...[.67 .67 .67] ... [.5 .5 .5] ...
                    ,'Margin',2.5 ...
                  );
                  %delete(hLegend4SampleClasses);
                  hLegend4SampleClasses = text(...
                     min(xRange) + (cil-1) * abs(diff(xlim(aUnderbarTarget)))/2.5 ...
                    ...,mean(YLim(2)+[1/500,21.7/500*iif(bInPipelinePlot,maxFigureHeight/figureHeight,1)]*abs(-YLim(1)+YLim(2))) + (0+2.5)*abs(diff(yRange)) ...
                    ,YLim(1) - abs(diff(YLim))*0.33 ...
                    ,[{['\color[rgb]{1 1 1}',sprintf('%d) ',cil),inInfo.plots.displayPredefinedClasses4aShifts(cil).legendHeader4samples]}...
                      ,cellfun(...
                          @(s,c)sprintf('\\color[rgb]{%0.2f %0.2f %0.2f}%s',c(1),c(2),c(3),s)...
                         ,uniqueSampleClasses...
                         ,num2cell(colors4SampleClasses,2)'...
                       ,'UniformOutput',false)...
                     ] ...
                    ,'Parent',aUnderbarTarget ...
                    ,'HorizontalAlignment','left', 'VerticalAlignment','top' ...
                    ,'FontSize',latexFontSize-5, 'FontName','Helvetica', 'FontWeight','bold' ...
                    ,'Interpreter','Tex'...
                    ,'Clipping','off'...
                    ,'BackgroundColor',[.67 .67 .67] ... [.5 .5 .5] ...
                    ,'Margin',2.5 ...
                  );
                end
                xlim(XLim); %prevent auto-changing XLim due to text objects beyond the clipping bounds and auto-zoom-out.
              end
              if(length(inInfo.plots.displayPredefinedClasses4aShifts)>1)
                h=xlabel(signatureDefinition.plots.overview.aShifts...
                  ,[repmat({' '},ceil(1.33*(length(inInfo.plots.displayPredefinedClasses4aShifts)-1)),1);get(get(signatureDefinition.plots.overview.aShifts,'XLabel'),'String')] ...
                  ,'Interpreter','LaTex','FontSize',latexFontSize,'Color',inInfo.plots.color4T4P*0.505,'VerticalAlignment','top'...
                );              
              end
            end
      %% Export overview plot:
        if(inInfo.export.plots.bEnabled)
          saveFormat = inInfo.export.plots.exportFormats; %= {'eps'};
          userDefinedPaperType = '1:1';
          exportDPI = 300;
          if(~ispc()) %workaround against truncated heatmaps in eps exports for Linux.
            drawnow;
          end
          saveFig(...
            signatureDefinition.plots.overview.f...
           ,inInfo.export.subdir4signature ...
           ,sprintf('%ssignature overview',iif(~isempty(inInfo.export.fileNamePrefix4signature),inInfo.export.fileNamePrefix4signature,sprintf('%03d, ',signatureDefinition.reference.k)))...
           ,saveFormat...
           ,userDefinedPaperType, exportDPI ...
          );
            fprintf('   <-saved plot: %s\n', fullfile(inInfo.export.subdir4signature, [sprintf('%soverview',iif(~isempty(inInfo.export.fileNamePrefix4signature),inInfo.export.fileNamePrefix4signature,sprintf('%03d, ',signatureDefinition.reference.k))), '.', saveFormat{1}]));
          if(inInfo.export.plots.bCloseAfterExport)
            close(signatureDefinition.plots.overview.f);
          else
            set(signatureDefinition.plots.overview.f, 'Visible','on');
          end
        end
    end
  %% EPS) Disssection QC plot:
    if(inInfo.plots.dissectionEffectiveness.bEnabled && ~inInfo.applicationMode.bEnabled)
      if(~inInfo.export.plots.bEnabled && ~inInfo.plots.bVisible)
        warning('~inInfo.export.plots.bEnabled, but also ~inInfo.plots.bVisible => setting plot to visible.');
        inInfo.plots.bVisible = true;
      end
      signatureDefinition.plots.dissectionSignatureOnTValues.f = figure('Position',[679, 403, 1243, 577], 'Visible', iif(inInfo.plots.bVisible,'on','off'), 'Color','w');
        subplot(1,2,1);
          plot(signatureDefinition.step3_regression.signatureStrengths4G, signatureDefinition.afterDissection.signatureStrengths4G,'.');
          xlabel('signatureStrengths4G before dissection');
          ylabel('signatureStrengths4G after dissection');
          title(sprintf('detected signature %d\ndissection effectiveness\ngene prespective', signatureDefinition.reference.k));
          daspect([1 1 1]);
          box on; grid on;
        subplot(1,2,2);
          plot(signatureDefinition.step3_regression.signatureStrengths4P, signatureDefinition.afterDissection.signatureStrengths4P,'.');
          xlabel('signatureStrengths4P before dissection');
          ylabel('signatureStrengths4P after dissection');
          title(sprintf('detected signature %d\ndissection effectiveness\nsample prespective', signatureDefinition.reference.k));
          daspect([1 1 1]);
          box on; grid on;
      %Export dissection plot:
        if(inInfo.export.plots.bEnabled)
          saveFormat = {'eps'};
          saveFig(...
            signatureDefinition.plots.dissectionSignatureOnTValues.f...
           ,inInfo.export.subdir4signature...
           ,sprintf('%sdissection effectiveness',iif(~isempty(inInfo.export.fileNamePrefix4signature),inInfo.export.fileNamePrefix4signature,sprintf('%03d, ',signatureDefinition.reference.k)))...
           ,saveFormat...
           ,userDefinedPaperType, exportDPI...
          );
            fprintf('   <-saved plot: %s\n', fullfile(inInfo.export.subdir4signature, [sprintf('%sdissection effectiveness',iif(~isempty(inInfo.export.fileNamePrefix4signature),inInfo.export.fileNamePrefix4signature,sprintf('%03d, ',signatureDefinition.reference.k))), '.', saveFormat{1}]));
          if(inInfo.export.plots.bCloseAfterExport && ishandle(signatureDefinition.plots.overview.f))
            close(signatureDefinition.plots.overview.f);
          end
        end
    end
  %% XLSX) Export signature including axes, correlations, weights, all before and after dissection:
    if(inInfo.export.table4eachSignatureDefinition.bEnabled)
      %Get descriptions:
        rowLabels = cellfun(@(s)strrep(s,'\color[rgb]{0 0 0.67}',''),inInfo.plots.descriptions.rowLabels,'UniformOutput',false);
        columnLabels = inInfo.plots.descriptions.colLabels';
      %Gene perspective:
        if(true)
          saveAsXLSXInInfo4GenesWorksheet = struct();
            saveAsXLSXInInfo4GenesWorksheet.sWorksheetName = sprintf('signature %s order %03d', inInfo.plots.descriptions.sGenesMnemonic, signatureDefinition.reference.k);
            saveAsXLSXInInfo4GenesWorksheet.formatting.colors.background = {1,':',[.8 .8 .8]};
            saveAsXLSXInInfo4GenesWorksheet.formatting.colors.bMACExcelColorBug = true;
            saveAsXLSXInInfo4GenesWorksheet.formatting.fontAlignment = {':','2:end','center'};
            saveAsXLSXInInfo4GenesWorksheet.formatting.fontBold = {1,':'};
            saveAsXLSXInInfo4GenesWorksheet.formatting.wrapText = {1,':'};
            saveAsXLSXInInfo4GenesWorksheet.display.freezePanes.rowsUntil = 1;
            saveAsXLSXInInfo4GenesWorksheet.display.columnWidths = nan(1,9+iif(inInfo.plots.signatureDefinitions.bEnabled,1,0)-iif(~inInfo.applicationMode.bEnabled,0,2));
              saveAsXLSXInInfo4GenesWorksheet.display.columnWidths(3:end) = 13;
            saveAsXLSXInInfo4GenesWorksheet.postprocess.createColumnGroups.columnTreePaths = [...
              repmat({{}},1,5+iif(inInfo.plots.signatureDefinitions.bEnabled,1,0)) ...
             ,repmat({{iif(~inInfo.applicationMode.bEnabled,'before signature dissection','signature application')}},1,2) ...
             ,iif(~inInfo.applicationMode.bEnabled,repmat({{'after signature dissection'}},1,2),{}) ...
            ];
          flatTable4GenesWorksheet = [
            [{ 'row #', inInfo.plots.descriptions.sGenesLabel...'row ID'...
             ,'|b^g_l\hat> (signature gene axis)' ...
             ,'|w^g_l\hat> (signature gene focus)' ...
             ,'|v^g_l\hat> (extended gene focus)' ...
             },iif(inInfo.plots.signatureDefinitions.bEnabled, {...
              'selected by plot filter'...
             },{}),{
              iif(~inInfo.applicationMode.bEnabled...
                ,'|u^g> (gene strengths for genes in M_(k-1))'...
                ,'|u^g> (gene strengths for genes in the provided input signal)'...
              )...
             ,iif(~inInfo.applicationMode.bEnabled...
                ,'|r^g> (correlations with the signature''s sample axis for M_(k-1))'...
                ,'|r^g> (correlations with the signature''s sample axis for the provided input signal)'...
              )...
              ...,'|p^g> (p values for correlations)' ...
             },iif(~inInfo.applicationMode.bEnabled, {...
              '|u^g> (signature strengths for genes in the remaining signal M_k)'...
             ,'|r^g> (correlations with the signature''s sample axis for the remaining signal M_k)' ...
             ...,'|p^g> (p values for correlations)' ...
             },{})];
            [[ num2cell(1:nG)', rowLabels...
             ,num2cell(signatureDefinition.step2_finalSignatureAxes.geneAxis_origUnits)...
             ,num2cell(abs(signatureDefinition.step3_regression.signedFocusedW4G))...
             ,num2cell(abs(signatureDefinition.step3_regression.signedExtendedW4G))...
            ],iif(inInfo.plots.signatureDefinitions.bEnabled, [...
              num2cell(ismember((1:nG)',SI)) ... 
            ],cell(nG,0)), [ ...
              num2cell(signatureDefinition.step3_regression.signatureStrengths4G)...
             ,num2cell(signatureDefinition.step3_regression.R4G)...
             ...,num2cell(signatureDefinition.step3_regression.P4R4G)...
            ],iif(~inInfo.applicationMode.bEnabled, [...
              num2cell(signatureDefinition.afterDissection.signatureStrengths4G) ... 
             ,num2cell(signatureDefinition.afterDissection.R4G)...
             ...,num2cell(signatureDefinition.afterDissection.P4R4G)...
            ],cell(nG,0))] ...
          ];
%             [[ num2cell(1:nG)', rowLabels...
%              ,num2cell(signatureDefinition.step2_finalSignatureAxes.geneAxis_origUnits)...
%              ,num2cell(abs(signatureDefinition.step2_finalSignatureAxes.signedFocusedW4G))...
%              ,num2cell(abs(signatureDefinition.step2_finalSignatureAxes.signedExtendedW4G))...
%             ],iif(inInfo.plots.signatureDefinitions.bEnabled, [...
%               num2cell(ismember((1:nG)',SI)) ... 
%             ],cell(nG,0)), [ ...
%               num2cell(signatureDefinition.step3_regression.signatureStrengths4G)...
%              ,num2cell(signatureDefinition.step2_finalSignatureAxes.R4G)...
%              ...,num2cell(signatureDefinition.step2_finalSignatureAxes.P4R4G)...
%             ],iif(~inInfo.applicationMode.bEnabled, [...
%               num2cell(signatureDefinition.afterDissection.signatureStrengths4G) ... 
%              ,num2cell(signatureDefinition.afterDissection.R4G)...
%              ...,num2cell(signatureDefinition.afterDissection.P4R4G)...
%             ],cell(nG,0))] ...

            saveAsXLSXInInfo4GenesWorksheet.dataFormats = {
               '2:end', '2:end', '0.00'
               '2:end',find(cellfun(@(s)~isempty(strfind(s,'|p')),flatTable4GenesWorksheet(1,:))),'[>=0.0001]0.0000;0.0E+00'
            };
            %Sort in eigenOrder:
              flatTable4GenesWorksheet(2:end,:) = flatTable4GenesWorksheet(signatureDefinition.step3_regression.eigenSI+1,:);
            %Filter config:
              saveAsXLSXInInfo4GenesWorksheet.autoFilter.cols = 1:size(flatTable4GenesWorksheet,2);
            saveAsXLSXInInfo4GenesWorksheet.display.rowHeights = nan(size(flatTable4GenesWorksheet,1),1);
              saveAsXLSXInInfo4GenesWorksheet.display.rowHeights(1) = 85;
        end
      %Sample perspective:
        if(true)
          saveAsXLSXInInfo4SamplesWorksheet = struct();
            saveAsXLSXInInfo4SamplesWorksheet.sWorksheetName = sprintf('signature %s order %03d', inInfo.plots.descriptions.sSamplesMnemonic, signatureDefinition.reference.k);
            saveAsXLSXInInfo4SamplesWorksheet.formatting.colors.background = {1,':',[.8 .8 .8]};
            saveAsXLSXInInfo4SamplesWorksheet.formatting.colors.bMACExcelColorBug = true;
            saveAsXLSXInInfo4SamplesWorksheet.formatting.fontAlignment = {':','2:end','center'};
            saveAsXLSXInInfo4SamplesWorksheet.formatting.fontBold = {1,':'};
            saveAsXLSXInInfo4SamplesWorksheet.formatting.wrapText = {1,':'};
            saveAsXLSXInInfo4SamplesWorksheet.display.freezePanes.rowsUntil = 1;
            saveAsXLSXInInfo4SamplesWorksheet.display.columnWidths = nan(1,9+iif(inInfo.plots.signatureDefinitions.bEnabled,1,0)-iif(~inInfo.applicationMode.bEnabled,0,2));
              saveAsXLSXInInfo4SamplesWorksheet.display.columnWidths(3:end) = 13;
            saveAsXLSXInInfo4SamplesWorksheet.saveToXLSXFilename = fullfile(...
               inInfo.export.subdir4signature...
              ,[iif(~isempty(inInfo.export.fileNamePrefix4signature),inInfo.export.fileNamePrefix4signature,sprintf('%03d, ',signatureDefinition.reference.k)),'definition','.xlsx']... %nochmals mit signature name since Excel cannot open multiple files with the same name.
            );
            saveAsXLSXInInfo4SamplesWorksheet.bOverwrite = true;
            saveAsXLSXInInfo4SamplesWorksheet.postprocess.createColumnGroups.columnTreePaths = [...
              repmat({{}},1,5+iif(inInfo.plots.signatureDefinitions.bEnabled,1,0)) ...
             ,repmat({{iif(~inInfo.applicationMode.bEnabled,'before signature dissection','signature application')}},1,2) ...
             ,iif(~inInfo.applicationMode.bEnabled,repmat({{'after signature dissection'}},1,2),{}) ...
            ];
            saveAsXLSXInInfo4SamplesWorksheet.postprocess.lastSelection = {saveAsXLSXInInfo4GenesWorksheet.sWorksheetName, 2, 3};
          additionalSampleColumns.headers = {};
          additionalSampleColumns.data = cell(length(columnLabels),0);
          %Add color underbars as additional columns, if available:
            if(isfield(inInfo.plots,'displayPredefinedClasses4aCurrentL2Rs'))
              for cil=1:length(inInfo.plots.displayPredefinedClasses4aCurrentL2Rs)
                additionalSampleColumns.headers{end+1} = inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).legendHeader4samples;
                additionalSampleColumns.data = [additionalSampleColumns.data, inInfo.plots.displayPredefinedClasses4aCurrentL2Rs(cil).samples(:)];
                saveAsXLSXInInfo4SamplesWorksheet.postprocess.createColumnGroups.columnTreePaths{end+1} = {'clinical data'};
                saveAsXLSXInInfo4SamplesWorksheet.display.columnWidths(end+1) = NaN;
              end
            end
            if(bInPipelinePlot && isfield(inInfo.plots,'displayPredefinedClasses4aShifts'))
              for cil=1:length(inInfo.plots.displayPredefinedClasses4aShifts)
                additionalSampleColumns.headers{end+1} = inInfo.plots.displayPredefinedClasses4aShifts(cil).legendHeader4samples;
                C = inInfo.plots.displayPredefinedClasses4aShifts(cil).samples(:);
                  if(isnumeric(C)) C = num2cell(C); end
                additionalSampleColumns.data = [additionalSampleColumns.data, C];
                saveAsXLSXInInfo4SamplesWorksheet.postprocess.createColumnGroups.columnTreePaths{end+1} = {'clinical data'};
                saveAsXLSXInInfo4SamplesWorksheet.display.columnWidths(end+1) = NaN;
              end
            end
          flatTable4SamplesWorksheet = [
            [{'col #', inInfo.plots.descriptions.sSamplesLabel...
             ,'|b^s_l\hat> (signature sample axis)' ...
             ,'|w^s_l\hat> (signature sample focus)' ...
             ,'|v^s_l\hat> (extended sample focus)' ...
             },iif(inInfo.plots.signatureDefinitions.bEnabled, {...
              'selected by plot filter'...
             },{}),{
              iif(~inInfo.applicationMode.bEnabled...
                ,'|u^s> (signature strengths for samples in M_(k-1))'...
                ,'|u^s> (signature strengths for samples in the provided input signal)'...
              )...
             ,iif(~inInfo.applicationMode.bEnabled...
                ,'|r^s> (correlations with the signature''s gene axis for M_(k-1))'...
                ,'|r^s> (correlations with the signature''s gene axis for the provided input signal)'...
              )...
             ...,'|p^s> (p values for correlations)' ...
            },iif(~inInfo.applicationMode.bEnabled, {...
              '|u^s> (signature strengths for samples in M_k)'...
             ,'|r^s> (correlations with the signature''s gene axis for M_k)'...
             ...,'|p^s> (p values for correlations)' ...
            },{}),additionalSampleColumns.headers];
            [[ num2cell(1:nP)', columnLabels...
             ,num2cell(signatureDefinition.step2_finalSignatureAxes.sampleAxis_origUnits)' ...
             ,num2cell(abs(signatureDefinition.step3_regression.signedFocusedW4P))' ...
             ,num2cell(abs(signatureDefinition.step3_regression.signedExtendedW4P))' ...
            ],iif(inInfo.plots.signatureDefinitions.bEnabled, [...
              num2cell(ismember((1:nP)',SJ)) ... 
            ],cell(nG,0)),[ ...
              num2cell(signatureDefinition.step3_regression.signatureStrengths4P)' ...
             ,num2cell(signatureDefinition.step3_regression.R4P)' ...
             ...,num2cell(signatureDefinition.step3_regression.P4R4P)' ...
            ],iif(~inInfo.applicationMode.bEnabled, [...
              num2cell(signatureDefinition.afterDissection.signatureStrengths4P)' ... 
             ,num2cell(signatureDefinition.afterDissection.R4P)' ...
             ...,num2cell(signatureDefinition.afterDissection.P4R4P)' ...
            ],cell(nP,0)) ...
             ,additionalSampleColumns.data ...
            ] ...
          ];
%             [[ num2cell(1:nP)', columnLabels...
%              ,num2cell(signatureDefinition.step2_finalSignatureAxes.sampleAxis_origUnits)' ...
%              ,num2cell(abs(signatureDefinition.step2_finalSignatureAxes.signedFocusedW4P))' ...
%              ,num2cell(abs(signatureDefinition.step2_finalSignatureAxes.signedExtendedW4P))' ...
%             ],iif(inInfo.plots.signatureDefinitions.bEnabled, [...
%               num2cell(ismember((1:nP)',SJ)) ... 
%             ],cell(nG,0)),[ ...
%               num2cell(signatureDefinition.step3_regression.signatureStrengths4P)' ...
%              ,num2cell(signatureDefinition.step2_finalSignatureAxes.R4P)' ...
%              ...,num2cell(signatureDefinition.step2_finalSignatureAxes.P4R4P)' ...
%             ],iif(~inInfo.applicationMode.bEnabled, [...
%               num2cell(signatureDefinition.afterDissection.signatureStrengths4P)' ... 
%              ,num2cell(signatureDefinition.afterDissection.R4P)' ...
%              ...,num2cell(signatureDefinition.afterDissection.P4R4P)' ...
%             ],cell(nP,0)) ...
%              ,additionalSampleColumns.data ...
%             ] ...

            saveAsXLSXInInfo4SamplesWorksheet.dataFormats = {
               '2:end', '2:end', '0.00'; 
               '2:end',find(cellfun(@(s)~isempty(strfind(s,'|p')),flatTable4GenesWorksheet(1,:))),'[>=0.0001]0.0000;0.0E+00'
            };
            %Sort in eigenOrder:
              flatTable4SamplesWorksheet(2:end,:) = flatTable4SamplesWorksheet(signatureDefinition.step3_regression.eigenSJ+1,:);
            %Filter config:
              saveAsXLSXInInfo4SamplesWorksheet.autoFilter.cols = 1:size(flatTable4SamplesWorksheet,2);
            saveAsXLSXInInfo4SamplesWorksheet.display.rowHeights = nan(size(flatTable4SamplesWorksheet,1),1);
              saveAsXLSXInInfo4SamplesWorksheet.display.rowHeights(1) = 85;
        end
      %Save tables in RAM:
        outInfo.excelExport.flatTable4GenesWorksheet = flatTable4GenesWorksheet;
        outInfo.excelExport.saveAsXLSXInInfo4GenesWorksheet = saveAsXLSXInInfo4GenesWorksheet;
        outInfo.excelExport.flatTable4SamplesWorksheet = flatTable4SamplesWorksheet;
        outInfo.excelExport.saveAsXLSXInInfo4SamplesWorksheet = saveAsXLSXInInfo4SamplesWorksheet;
      %Export to disk:
        if(ispc() && ~isempty(which('saveAsXLSX')))
          saveAsXLSXInInfo4GenesWorksheet.bAsBackgroundJob = true;
          saveAsXLSX(flatTable4GenesWorksheet,saveAsXLSXInInfo4GenesWorksheet, flatTable4SamplesWorksheet,saveAsXLSXInInfo4SamplesWorksheet);
            fprintf('   <-scheduled export of table: %s\n', saveAsXLSXInInfo4SamplesWorksheet.saveToXLSXFilename);        
        else
          T4genes = num2cell(flatTable4GenesWorksheet,1);
            T4genes = table(T4genes{:});
            writetable(T4genes, strrep(saveAsXLSXInInfo4SamplesWorksheet.saveToXLSXFilename,'.xlsx','-genes.txt'),'WriteVariableNames',false,'Delimiter','tab');
              fprintf('   <-exported table: %s\n', strrep(saveAsXLSXInInfo4SamplesWorksheet.saveToXLSXFilename,'.xlsx','-genes.txt'));
          T4samples = num2cell(flatTable4SamplesWorksheet,1);
            T4samples = table(T4samples{:});
            writetable(T4samples, strrep(saveAsXLSXInInfo4SamplesWorksheet.saveToXLSXFilename,'.xlsx','-samples.txt'),'WriteVariableNames',false,'Delimiter','tab');
              fprintf('   <-exported table: %s\n', strrep(saveAsXLSXInInfo4SamplesWorksheet.saveToXLSXFilename,'.xlsx','-samples.txt'));
        end
    end  
end

%Helper function: Keep axes limits in sync after the user panned/zoomed the top overlay axes:
  function addPanAndZoomHandler(f, aOverlay4ColTStats, aOverlay4RowTStats, lc_aCurrentL2Rs)
    zoomState = struct();
      zoomState.orig.aOverlay4ColT.xlim = xlim(aOverlay4ColTStats);
      zoomState.orig.aOverlay4ColT.ylim = ylim(aOverlay4ColTStats);
      zoomState.orig.aOverlay4RowT.xlim = xlim(aOverlay4RowTStats);
      zoomState.orig.aOverlay4RowT.ylim = ylim(aOverlay4RowTStats);
      zoomState.orig.aCurrentL2Rs.xlim = xlim(lc_aCurrentL2Rs);
      zoomState.orig.aCurrentL2Rs.ylim = ylim(lc_aCurrentL2Rs);
      zoomState.trafos.colT2L2RX = @(X)(X-zoomState.orig.aOverlay4ColT.xlim(1))/(-zoomState.orig.aOverlay4ColT.xlim(1)+zoomState.orig.aOverlay4ColT.xlim(2))*(-zoomState.orig.aCurrentL2Rs.xlim(1)+zoomState.orig.aCurrentL2Rs.xlim(2))+zoomState.orig.aCurrentL2Rs.xlim(1);
      zoomState.trafos.colT2L2RY = @(Y)fliplr(-(Y-zoomState.orig.aOverlay4ColT.ylim(1))/(-zoomState.orig.aOverlay4ColT.ylim(1)+zoomState.orig.aOverlay4ColT.ylim(2))*(-zoomState.orig.aCurrentL2Rs.ylim(1)+zoomState.orig.aCurrentL2Rs.ylim(2))+zoomState.orig.aCurrentL2Rs.ylim(2));
      zoomState.trafos.colT2rowTX = @(X)(X-zoomState.orig.aOverlay4ColT.xlim(1))/(-zoomState.orig.aOverlay4ColT.xlim(1)+zoomState.orig.aOverlay4ColT.xlim(2))*(-zoomState.orig.aOverlay4RowT.xlim(1)+zoomState.orig.aOverlay4RowT.xlim(2))+zoomState.orig.aOverlay4RowT.xlim(1);
      zoomState.trafos.colT2rowTY = @(Y)fliplr(-(Y-zoomState.orig.aOverlay4ColT.ylim(1))/(-zoomState.orig.aOverlay4ColT.ylim(1)+zoomState.orig.aOverlay4ColT.ylim(2))*(-zoomState.orig.aOverlay4RowT.ylim(1)+zoomState.orig.aOverlay4RowT.ylim(2))+zoomState.orig.aOverlay4RowT.ylim(2));
    set(zoom(f),'ActionPostCallback',@(f,e)updateAxesLimits(aOverlay4ColTStats, aOverlay4RowTStats, lc_aCurrentL2Rs, zoomState));
    set(pan(f),'ActionPostCallback',@(f,e)updateAxesLimits(aOverlay4ColTStats, aOverlay4RowTStats, lc_aCurrentL2Rs, zoomState));
  end
    function updateAxesLimits(aOverlay4ColTStats, aOverlay4RowTStats, lc_aCurrentL2Rs, zoomState)
      set(aOverlay4RowTStats,'XLim',zoomState.trafos.colT2rowTX(get(aOverlay4ColTStats,'XLim')), 'YLim',zoomState.trafos.colT2rowTY(get(aOverlay4ColTStats,'YLim')));
      set(lc_aCurrentL2Rs   ,'XLim',zoomState.trafos.colT2L2RX(get(aOverlay4ColTStats,'XLim')),  'YLim',zoomState.trafos.colT2L2RY(get(aOverlay4ColTStats,'YLim')));
    end



    
