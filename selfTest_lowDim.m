%ABSTRACT
%  Demo and self test of SDCM for a low-dimensional 3D test signal
%  simulating four distinct interactions.
%SYNTAX
%  Just run this script (press F5).

%% Initialize:
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
        error('Cannot find SDCM library functions subfolder at %s. Please add the library functions manually to your Matlab path. Pressing F5 and then changing the directory to this self test file should also be sufficient.', sPathToSelf);
      end
    end

%% Simulate a 3D signal containing four distinct signatures:
  simulation3D = struct();
  simulation3D.noiseSD = 0.1; %0.2;

  %a simple linear interaction with normally distributed sample strengths:
    simulation3D.signatures(1).geneAxis = [-0.2 0.2 1]';
    simulation3D.signatures(1).geneAxis = simulation3D.signatures(1).geneAxis/norm(simulation3D.signatures(1).geneAxis);
    simulation3D.signatures(1).data = bsxfun(@times, randn(1,1000), simulation3D.signatures(1).geneAxis);
    simulation3D.signatures(1).data = simulation3D.signatures(1).data + simulation3D.noiseSD*randn(size(simulation3D.signatures(1).data));
  %an interaction that is quadratic with respect to one gene:
    simulation3D.signatures(2).geneAxis = [-0.8 -0.2 0.5]';
    simulation3D.signatures(2).geneAxis = simulation3D.signatures(2).geneAxis /norm(simulation3D.signatures(2).geneAxis);
    simulation3D.signatures(2).data = bsxfun(@times, max(-3,min(3,randn(1,500))), simulation3D.signatures(2).geneAxis);
    simulation3D.signatures(2).data(3,:) = sign(simulation3D.signatures(2).data(3,:)).*abs(simulation3D.signatures(2).data(3,:)).^2;
    simulation3D.signatures(2).data = simulation3D.signatures(2).data + simulation3D.noiseSD*randn(size(simulation3D.signatures(2).data));
  %an interaction using the logistic function to mimick a saturation signature:
    simulation3D.signatures(3).geneAxis = [-0.5 0.75 0.2]';
    simulation3D.signatures(3).geneAxis = simulation3D.signatures(3).geneAxis /norm(simulation3D.signatures(3).geneAxis);
    S = randn(1,500);
    simulation3D.signatures(3).data = bsxfun(@times, S, simulation3D.signatures(3).geneAxis);
    simulation3D.signatures(3).data(1,:) = 3*(-0.5 + 1./(1+exp(-4*(simulation3D.signatures(3).data(1,:)-0)))); 
    simulation3D.signatures(3).data = simulation3D.signatures(3).data + simulation3D.noiseSD*randn(size(simulation3D.signatures(3).data));
  %a on/off interaction with an offset mimicking an activation threshold:
    simulation3D.signatures(4).geneAxis = [-.4 -.7 -.7]'; 
    simulation3D.signatures(4).geneAxis = simulation3D.signatures(4).geneAxis/norm(simulation3D.signatures(4).geneAxis);
    simulation3D.signatures(4).data = bsxfun(@times, randn(1,500), simulation3D.signatures(4).geneAxis);
    B = simulation3D.signatures(4).data(1,:) > quantile(simulation3D.signatures(4).data(1,:),0.2);
    simulation3D.signatures(4).data(:,B) = 0; %one-sided offset signature.
    simulation3D.signatures(4).data = simulation3D.signatures(4).data + simulation3D.noiseSD*randn(size(simulation3D.signatures(4).data));

  %Combine and create indices
    simulation3D.combinedData = [];
    simulation3D.sampleIndices4simulatedSignatures = {};
    sampleCursor = 0;
    for k=1:4
      simulation3D.combinedData = [simulation3D.combinedData, simulation3D.signatures(k).data];
      simulation3D.signatures(k).ISamplesInCombinedData = sampleCursor+(1:size(simulation3D.signatures(k).data,2)); %needed for plotting.
      sampleCursor=sampleCursor+size(simulation3D.signatures(k).data,2);
    end
    for k=1:4
      simulation3D.signatures(k).BSamplesInCombinedData = ismember(1:size(simulation3D.combinedData,2),simulation3D.signatures(k).ISamplesInCombinedData);  %needed for plotting.
    end
    simulation3D.nG = size(simulation3D.combinedData,1);
    simulation3D.nP = size(simulation3D.combinedData,2);
  %% Plot simulated data:
    %Plot config:
      CM = [1 .5 .5; .5 1 .5; .5 .5 1; 1 .5 1; .5 1 1; 1 1 .5; .5 .5 .5];
      fcnAdjustBrightness = @(color,zeroIsBlackOneIsWhite)iif(zeroIsBlackOneIsWhite<0.5, color*zeroIsBlackOneIsWhite, color*(1-2*(zeroIsBlackOneIsWhite-0.5))+[1 1 1]*2*(zeroIsBlackOneIsWhite-0.5));            
      bShowTicksAndLabels = false;
      lim = [-2,2]; %[-2.23,2.23]; %[-2,2];
      nMarkerSize = 11.5; %10.3; %13; %10.3;
      commonView = [11.5, 37];
      nCameraViewAngle = 9.1;
      nCameraTarget = [-0.0607003 -0.176287 -0.13055];
    %Initial data:
      f = figure('Position',[1 469 565 535], 'Name', 'all points, uni color','Color','w');
      a = axes('FontSize',13,'Projection','perspective');
      if(true)
        hold on;
        %export point-per-point (plot the points in z order for EPS export usign 'SortMethod' 'childorder' to prevent EPS artefacts)
          [~,SI4Z] = sort(simulation3D.combinedData(2,:),'descend');
          for ii=1:length(SI4Z)
            kSignature = ils(SI4Z(ii), @(i)simulation3D.signatures(1).BSamplesInCombinedData(i),1, @(i)simulation3D.signatures(2).BSamplesInCombinedData(i),2, @(i)simulation3D.signatures(3).BSamplesInCombinedData(i),3, @(i)simulation3D.signatures(4).BSamplesInCombinedData(i),4);
            plot3(simulation3D.combinedData(1,SI4Z(ii)), simulation3D.combinedData(2,SI4Z(ii)), simulation3D.combinedData(3,SI4Z(ii)), 'o' ...
              ,'MarkerEdgeColor', fcnAdjustBrightness(CM(kSignature,:),0.5-0.07) ...
              ,'MarkerFaceColor', fcnAdjustBrightness(CM(kSignature,:),0.5+0.07) ...
            ,'MarkerSize', nMarkerSize/2.5,'LineWidth',1);
          end
        daspect([1 1 1]); pbaspect([1 1 1]); 
        view(commonView);
        axis tight; box on; grid off;
        if(bShowTicksAndLabels)
          set(gca,'FontSize',13,'FontName','Cambria','FontWeight','Bold');
          xlabel('gene x'); ylabel('gene y'); zlabel('gene z');
        else
          %set(gca,'XTick',[],'YTick',[],'ZTick',[]);
          set(gca,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{});
        end
        set(gca,'LineWidth',1.51,'Color','w');
        set(gca,'XColor',[0 0 0], 'YColor',[0 0 0], 'ZColor',[0 0 0], 'Clipping','off', 'Box','on');%'off');
        set(gca,'CameraViewAngle',nCameraViewAngle);
        set(gca,'CameraTarget',nCameraTarget);
      end
      xlim(lim);ylim(lim);zlim(lim);
    
%% Dissect the signal:
  configuration = SDCM_defaultConfig(simulation3D.nG, simulation3D.nP);
    configuration.plots.bVisible = true; %show plots for user.          
    configuration.export.nStatusOutputLevel = 1; %disable most status messages.#
    configuration.internal.bAllowDissectionOfTinyDatasets = true; %disable warnings for small datasets.
  detection = SDCM(simulation3D.combinedData, configuration);

%% Plot detected signatures as gene curves: 
  fprintf('\nNote: this self test for randomly generated input data can be considered passed, if generated plots resemble expected results at \\SDCM_selfTest_expectedResults\\selfTest_lowDim_expectedResult.png\n');
  %Plot config:
    CM = [1 .5 .5; .5 1 .5; .5 .5 1; 1 .5 1; .5 1 1; 1 1 .5; .5 .5 .5];
    fcnAdjustBrightness = @(color,zeroIsBlackOneIsWhite)iif(zeroIsBlackOneIsWhite<0.5, color*zeroIsBlackOneIsWhite, color*(1-2*(zeroIsBlackOneIsWhite-0.5))+[1 1 1]*2*(zeroIsBlackOneIsWhite-0.5));            
    bShowTicksAndLabels = false;
    lim = [-2,2]; %[-2.23,2.23]; %[-2,2];
    nMarkerSize = 11.5; %10.3; %13; %10.3;
    commonView = [11.5, 37];
    nCameraViewAngle = 9.1;
    nCameraTarget = [-0.0607003 -0.176287 -0.13055];
  %full signal, all signature curves:
    if(true)
      f = figure('Position',[1 1 565 667], 'Name', 'all points, uni color','Color','w');
      nSubplots = 1;
      a=axes('FontSize',13,'Projection','perspective');
        if(true)
          hold on;
          %Exportiere Punkte einzeln in z-Order f�r EPS-Export mit 'SortMethod'='childorder' (um Effektkurve nicht zu zerst�ckeln)
            [~,SI4Z] = sort(simulation3D.combinedData(2,:),'descend');
            for ii=1:length(SI4Z)
              kSignature = ils(SI4Z(ii), @(i)simulation3D.signatures(1).BSamplesInCombinedData(i),1, @(i)simulation3D.signatures(2).BSamplesInCombinedData(i),2, @(i)simulation3D.signatures(3).BSamplesInCombinedData(i),3, @(i)simulation3D.signatures(4).BSamplesInCombinedData(i),4);
              plot3(simulation3D.combinedData(1,SI4Z(ii)), simulation3D.combinedData(2,SI4Z(ii)), simulation3D.combinedData(3,SI4Z(ii)), 'o' ...
                ,'MarkerEdgeColor', fcnAdjustBrightness(CM(kSignature,:),0.5-0.07) ...
                ,'MarkerFaceColor', fcnAdjustBrightness(CM(kSignature,:),0.5+0.07) ...
              ,'MarkerSize', nMarkerSize/2.5,'LineWidth',1);
            end

          daspect([1 1 1]); pbaspect([1 1 1]); 
          view([10.3,29.5]); %view([11.5,27.1]); %%view([25,31]);
          view(commonView);
          axis tight; box on; grid off;
          if(bShowTicksAndLabels)
            set(gca,'FontSize',13,'FontName','Cambria','FontWeight','Bold');
            xlabel('gene x'); ylabel('gene y'); zlabel('gene z');
          else
            %set(gca,'XTick',[],'YTick',[],'ZTick',[]);
            set(gca,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{});
          end
          set(gca,'LineWidth',1.51,'Color','w');
          set(gca,'XColor',[0 0 0], 'YColor',[0 0 0], 'ZColor',[0 0 0], 'Clipping','off', 'Box','on');%'off');
          %set(gca,'Box','off','XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
          set(gca,'CameraViewAngle',nCameraViewAngle);
          set(gca,'CameraTarget',nCameraTarget);
        end
        xlim(lim);ylim(lim);zlim(lim);
      %Show the signature curve in gene space:
        if(length(detection.signatures)>=1)
          firstEigensignal = zeros(size(simulation3D.combinedData));
          firstEigensignal(detection.signatures(1).step3_regression.BNonzeroRows4signatureEigensignal, detection.signatures(1).step3_regression.BNonzeroCols4signatureEigensignal) = ...
            detection.signatures(1).step3_regression.signatureEigensignal ./ detection.signatures(1).step3_regression.signatureDissectionStrengthsOS  ... %even in case of dissectionStrength<1, show the full regressed eigensignal.
          ;
          signatureCurveInGeneSpace = firstEigensignal(:, detection.signatures(1).step3_regression.eigenSJ);
            signatureCurveInGeneSpace(:,all(signatureCurveInGeneSpace==0,1)) = []; %exclude those samples no weight in this effect.
          hLine = plot3(signatureCurveInGeneSpace(1,:), signatureCurveInGeneSpace(2,:), signatureCurveInGeneSpace(3,:)...
            ,'-', 'Color', [1,0.91,0], 'LineWidth', 7.9 ...                    
          );
            hC = get(gca,'Children');
              hC(ismember(hC,hLine)) = [];
            set(gca,'Children',[hC;hLine]);
            set(gca,'SortMethod','childorder'); %!make sure that the yellow line is not split by all points in the plot.
        end
        if(length(detection.signatures)>=2)
          secondEigensignal = zeros(size(simulation3D.combinedData));
          secondEigensignal(detection.signatures(2).step3_regression.BNonzeroRows4signatureEigensignal, detection.signatures(2).step3_regression.BNonzeroCols4signatureEigensignal) = ...
            detection.signatures(2).step3_regression.signatureEigensignal ./ detection.signatures(2).step3_regression.signatureDissectionStrengthsOS ... %even in case of dissectionStrength<1, show the full regressed eigensignal.
          ;
          signatureCurveInGeneSpace = secondEigensignal(:, detection.signatures(2).step3_regression.eigenSJ);
            signatureCurveInGeneSpace(:,all(signatureCurveInGeneSpace==0,1)) = [];
          hLine = plot3(signatureCurveInGeneSpace(1,:), signatureCurveInGeneSpace(2,:), signatureCurveInGeneSpace(3,:) ...
            ,'-', 'Color', [1,0.91,0], 'LineWidth', 7.9 ...                    
          );                
            hC = get(gca,'Children');
              hC(ismember(hC,hLine)) = [];
            set(gca,'Children',[hC;hLine]);
            set(gca,'SortMethod','childorder'); %!make sure that the yellow line is not split by all points in the plot.
        end
        if(length(detection.signatures)>=3)
          thirdEigensignal = zeros(size(simulation3D.combinedData));
          thirdEigensignal(detection.signatures(3).step3_regression.BNonzeroRows4signatureEigensignal, detection.signatures(3).step3_regression.BNonzeroCols4signatureEigensignal) = ...
            detection.signatures(3).step3_regression.signatureEigensignal ./ detection.signatures(3).step3_regression.signatureDissectionStrengthsOS  ... %even in case of dissectionStrength<1, show the full regressed eigensignal.
          ;
          signatureCurveInGeneSpace = thirdEigensignal(:, detection.signatures(3).step3_regression.eigenSJ);
            signatureCurveInGeneSpace(:,all(signatureCurveInGeneSpace==0,1)) = [];
          hLine = plot3(signatureCurveInGeneSpace(1,:), signatureCurveInGeneSpace(2,:), signatureCurveInGeneSpace(3,:) ...
            ,'-','Color', [1,0.91,0], 'LineWidth', 7.9 ...                    
          );
            hC = get(gca,'Children');
              hC(ismember(hC,hLine)) = [];
            set(gca,'Children',[hC;hLine]);
            set(gca,'SortMethod','childorder'); %!make sure that the yellow line is not split by all points in the plot.
        end
        if(length(detection.signatures)>=4)
          fourthEigensignal = zeros(size(simulation3D.combinedData));
          fourthEigensignal(detection.signatures(4).step3_regression.BNonzeroRows4signatureEigensignal, detection.signatures(4).step3_regression.BNonzeroCols4signatureEigensignal) = ...
            detection.signatures(4).step3_regression.signatureEigensignal ./ detection.signatures(4).step3_regression.signatureDissectionStrengthsOS ... %even in case of dissectionStrength<1, show the full regressed eigensignal.
          ;
          signatureCurveInGeneSpace = fourthEigensignal(:, detection.signatures(4).step3_regression.eigenSJ);
            signatureCurveInGeneSpace(:,all(signatureCurveInGeneSpace==0,1)) = [];
          hLine = plot3(signatureCurveInGeneSpace(1,:), signatureCurveInGeneSpace(2,:), signatureCurveInGeneSpace(3,:) ...
            ,'-', 'Color', [1,0.91,0], 'LineWidth', 7.9 ...
          );
            hC = get(gca,'Children');
              hC(ismember(hC,hLine)) = [];
            set(gca,'Children',[hC;hLine]);
            set(gca,'SortMethod','childorder'); %!make sure that the yellow line is not split by all points in the plot.
        end
        title(a,'SDCM');
    end
  %full signal, PCA:
    if(true)
      %All points, uni color:
        f = figure('Position',[577 1 565 667], 'Name', 'all points, uni color','Color','w');
        nSubplots = 1;
        a=axes('FontSize',13,'Projection','perspective');
          if(true)
            hold on;
            %Exportiere Punkte einzeln in z-Order f�r EPS-Export mit 'SortMethod'='childorder' (um Effektkurve nicht zu zerst�ckeln)
              [~,SI4Z] = sort(simulation3D.combinedData(2,:),'descend');
              for ii=1:length(SI4Z)
                kSignature = ils(SI4Z(ii), @(i)simulation3D.signatures(1).BSamplesInCombinedData(i),1, @(i)simulation3D.signatures(2).BSamplesInCombinedData(i),2, @(i)simulation3D.signatures(3).BSamplesInCombinedData(i),3, @(i)simulation3D.signatures(4).BSamplesInCombinedData(i),4);
                plot3(simulation3D.combinedData(1,SI4Z(ii)), simulation3D.combinedData(2,SI4Z(ii)), simulation3D.combinedData(3,SI4Z(ii)), 'o' ...
                  ,'MarkerEdgeColor', fcnAdjustBrightness(CM(kSignature,:),0.5-0.07) ...
                  ,'MarkerFaceColor', fcnAdjustBrightness(CM(kSignature,:),0.5+0.07) ...
                ,'MarkerSize', nMarkerSize/2.5,'LineWidth',1);
              end
            daspect([1 1 1]); pbaspect([1 1 1]); 
            view([10.3,29.5]); %view([11.5,27.1]); %%view([25,31]);
            view(commonView);
            axis tight; box on; grid off;
            if(bShowTicksAndLabels)
              set(gca,'FontSize',13,'FontName','Cambria','FontWeight','Bold');
              xlabel('gene x'); ylabel('gene y'); zlabel('gene z');
            else
              %set(gca,'XTick',[],'YTick',[],'ZTick',[]);
              set(gca,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{});
            end
            set(gca,'LineWidth',1.51,'Color','w');
            set(gca,'XColor',[0 0 0], 'YColor',[0 0 0], 'ZColor',[0 0 0], 'Clipping','off', 'Box','on');%'off');
            %set(gca,'Box','off','XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
            set(gca,'CameraViewAngle',nCameraViewAngle);
            set(gca,'CameraTarget',nCameraTarget);
          end
          xlim(lim);ylim(lim);zlim(lim);
      %PCA results:
        [principalComponents, ~, explainedVariances] = pca(simulation3D.combinedData','Centered',false,'Economy',false);
        for d=1:3
          direction = principalComponents(:,d);
          line(...
             [-direction(1),direction(1)]*sqrt(explainedVariances(d))*7 ...
            ,[-direction(2),direction(2)]*sqrt(explainedVariances(d))*7 ...
            ,[-direction(3),direction(3)]*sqrt(explainedVariances(d))*7 ...
            ,'LineWidth', 7.9 ...
            ,'Color', [1,0.91,0] ...
            ,'LineStyle','-' ...
          );
        end
        title('PCA');
    end

