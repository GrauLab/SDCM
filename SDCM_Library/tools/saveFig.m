%ABSTRACT
% Save figure to file helper function.
%SYNTAX
% exportedFiles = saveFig(f, rootDir, filePrefix, formats...
%   ,pngPaperType, exportDPI, userDefinedPaperOrientation...
%   ,enforceFigPosForExportUnderParallelization...
% )
%AUTHOR
% (C) Michael Grau, 2011
% Library function for SDCM standalone deployment. 

function exportedFiles = saveFig(f, targetDir, filePrefix, exportFormats...
  ,userDefinedPaperType, exportDPI, userDefinedPaperOrientation...
  ,enforceFigPosForExportUnderParallelizationXorScaleFactor...
)
  %% Initialize
    exportedFiles = {};
    %Default params:
      if(nargin<4 || isempty(exportFormats)) exportFormats={'eps'}; end;
        if(ischar(exportFormats)) exportFormats={exportFormats}; end %support single format like 'eps'
      if(nargin<5 || isempty(userDefinedPaperType)) userDefinedPaperType='1:1'; end; %'A3';
      if(nargin<6 || isempty(exportDPI)) exportDPI = 400; end;
      if(nargin<7 || isempty(userDefinedPaperOrientation)) userDefinedPaperOrientation='portrait'; end;
    %Create target folder, if needed:
      if(~exist(targetDir,'dir'))
        if(ispc() && ~(length(targetDir)>=3 && strcmp(targetDir(1:3),'\\?'))) %add \\?\ prefix to support long paths.
          %If we already have a network UNC path line \\Server\Share\..., we need to replace \\ by UNC\ to get long path support in Windows:
            bIsAlreadyANetworkPath = length(targetDir)>=2 && strcmp(targetDir(1:2),'\\');
            if(bIsAlreadyANetworkPath)
              targetDir = ['UNC',targetDir(2:end)];
            end
          targetDir = ['\\?\',targetDir];
        end
        if(~exist(targetDir,'dir'))
          try
            mkdir(targetDir);
          catch ex
            warning('ERROR: Cannot create non-existing directory [%s]. If it is a UNC path, maybe the network is not ready? Please create it manually before using >>dbcont.');
            keyboard;
          end
        end
      end
    bExportUnderParallelizationDebug = false;
      if(bExportUnderParallelizationDebug)
        disp(sprintf('\nDEBUGGING figure export on worker side for "%s":',fullfile(targetDir, [filePrefix, '.', exportFormats{1}])));
        disp(sprintf('<- screen size: %s', mat2str(get(0,'ScreenSize')))); 
      end
    %Get information about screen size:
    screenSize = get(0,'ScreenSize');
      screenAspectRatio = screenSize(3)/screenSize(4);
    scaleFactor = 1;
    if(nargin>=8)
      if(isscalar(enforceFigPosForExportUnderParallelizationXorScaleFactor))
        scaleFactor = enforceFigPosForExportUnderParallelizationXorScaleFactor;
        enforceFigPosForExportUnderParallelizationXorScaleFactor = [];
      else
        enforceFigPosForExportUnderParallelization = enforceFigPosForExportUnderParallelizationXorScaleFactor;
      end
    end
  %% Backup figure's original properties:
    backupUnits = get(f,'Units');
    backupPos = get(f,'Position');      
    if(bExportUnderParallelizationDebug)
      disp(sprintf('<- figure position before drawnow: %s (units: %s)',mat2str(backupPos), backupUnits));
    end
    backupRenderer = get(f, 'Renderer');
    backupVisibility = get(f, 'Visible');
  %% figure positioning logic to get correct aspect ratios while exporting under parallelization (use normalized units):
    %get enforceFigPosForExportUnderParallelization:
      set(f,'Units','centimeters');
      if(nargin<8 || isempty(enforceFigPosForExportUnderParallelizationXorScaleFactor))
        if(~isWorker())
          enforceFigPosForExportUnderParallelization = get(f,'Position');
          if(bExportUnderParallelizationDebug)
            disp(sprintf('-> using current figure position as enforceFigPosForExportUnderParallelization = %s', mat2str(enforceFigPosForExportUnderParallelization)));
          end
        else %if we are a worker instance, the figure position is unreliable 
             % => calculate the outerBox(all axes) aspect ratio and infer enforceFigPosForExportUnderParallelization from these:
          HAxes = findobj(f,'Type','axes');
            if(bExportUnderParallelizationDebug)
              disp(sprintf('<- axes OuterPosition(s) inside of figure: ')); 
                disp(cellfun(...
                   @mat2str...
                  ,iif(length(HAxes)>1,get(HAxes,'OuterPosition'),{get(HAxes,'OuterPosition')})...
                  ,'UniformOutput',false...
                ));
            end
          %Note: don't go to normalized units here (may disrupt aspect ratio); instead /max() below: 
            backupUnits4Axes = get(HAxes,'Units');
              if(~iscell(backupUnits4Axes)) backupUnits4Axes={backupUnits4Axes}; end
              set(HAxes,'Units','centimeters');
            PAxes = get(HAxes,'OuterPosition');
              set(f,'ActivePositionProperty', 'OuterPosition'); %further resize events should keep the OuterPosition.
              if(~iscell(PAxes)) PAxes={PAxes}; end
            set(HAxes,{'Units'},backupUnits4Axes(:));
          %calculate axes aspect ratio:
            PAxes = vertcat(PAxes{:});
            outerAxesPosition = nan(1,4);
              outerAxesPosition(1) = min(PAxes(:,1));
              outerAxesPosition(2) = min(PAxes(:,2));
              outerAxesPosition(3) = max(PAxes(:,1)+PAxes(:,3))-0*outerAxesPosition(1); %(1) nicht abziehen, denn 0* funktioniert bei GO plots aus workers heraus!/ (evtl. weil innerFigure-Verschiebung irrelevant?)
              outerAxesPosition(4) = max(PAxes(:,2)+PAxes(:,4))-0*outerAxesPosition(2); %(2) nicht abziehen, denn 0* funktioniert bei GO plots aus workers heraus!/ (evtl. weil innerFigure-Verschiebung irrelevant?)
              if(bExportUnderParallelizationDebug)
                disp(sprintf('   -> outerBox of child axes (in unchanged units): %s', mat2str(outerAxesPosition)));
              end
          %Use outer axes dimensions as target size for the figure (and the paperPosition and the export width&height later):
            enforceFigPosForExportUnderParallelization = outerAxesPosition;
          if(bExportUnderParallelizationDebug)
            disp(sprintf('-> determined enforceFigPosForExportUnderParallelization in worker based on axes = %s', mat2str(enforceFigPosForExportUnderParallelization)));
          end
        end
      else %user-specified enforceFigPosForExportUnderParallelization:
        if(any(enforceFigPosForExportUnderParallelization>25)) %user has input it in pixel units => converet to centimeters
          enforceFigPosForExportUnderParallelization = enforceFigPosForExportUnderParallelization/get(0,'ScreenPixelsPerInch')*2.54;
        end
      end
    %define targetAspectRatio and set enforceFigPosForExportUnderParallelization to figure as centimeters:
      targetAspectRatio = enforceFigPosForExportUnderParallelization(3)/enforceFigPosForExportUnderParallelization(4);
        if(bExportUnderParallelizationDebug)
          disp(sprintf('-> targetAspectRatio = %f', targetAspectRatio ));
        end
      %allow for maximal figure sized wrt. current screen resolution before window size cropping occurs.
        enforceFigPosForExportUnderParallelization(1) = 0; 
        enforceFigPosForExportUnderParallelization(2) = 0;
      %feature scaleFactor:
        enforceFigPosForExportUnderParallelization = enforceFigPosForExportUnderParallelization*scaleFactor;
        %<-Note: this is useful for for easy relative font size manipuilation in png exports.
      %APPLY figure position
        sicherheitsRandFaktor = 1; %Sicherheitsrand.
        set(f,'Position', sicherheitsRandFaktor.*enforceFigPosForExportUnderParallelization);
      %redraw (and thereby invoke any resizeFcn) (important to do this explicitly in workers for correct export under parallelization):
        if(isWorker())
          set(f,'Visible','off');
          drawnow;
          set(f,'Visible','on'); %important for correct export under parallelization.
          drawnow; %important for correct export under parallelization.
        end
        if(bExportUnderParallelizationDebug)
          disp(sprintf('<- actual figure position (after drawnow, if worker): %s (units: %s)', mat2str(get(f,'Position')),get(f,'Units')));
          disp(sprintf('-> axes OuterPosition(s) inside of figure are now: ')); 
            disp(cellfun(...
               @mat2str...
              ,iif(@iscell...
                ,get(findobj(f,'Type','axes'),'OuterPosition')...
                ,{get(findobj(f,'Type','axes'),'OuterPosition')}...
               )...
              ,'UniformOutput',false...
            ));
        end
    %fix for changing aspect ratio due to too large figure windows:
      actualAspectRatio = get(f,'Position');
      actualAspectRatio = actualAspectRatio(3)/actualAspectRatio(4);
      if(~isWorker() && max(actualAspectRatio/targetAspectRatio,targetAspectRatio/actualAspectRatio)>=1.03)
        %Just warn:
          warning([...
                'Concerning export of figure "%s":'...
               ,'\nThe target aspect ratio = %f was not retained after resizing for printout; instead the resulting aspect ratio was = %f.'...
               ,'\n<-This can happen, if the figure window became larger than the screen size due to scaling (the operating system crops the window size in this case).' ...
               ,'\n=>The aspect ratio and expected relative font sizes may be off. To prevent this, design the figure for a lower screen resolution or use a lower scaleFactor.' ...
             ]...
            ,fullfile(targetDir, [filePrefix, '.', exportFormats{1}]) ...
            ,targetAspectRatio ...
            ,actualAspectRatio ...
          );
      end
  %% Paper setup:
    if(~strcmp(userDefinedPaperType,'1:1'))
      set(f,'PaperType',userDefinedPaperType);
      set(f,'PaperOrientation', userDefinedPaperOrientation);
    else
      set(f,'PaperType','<custom>');
      if(enforceFigPosForExportUnderParallelization(3)>enforceFigPosForExportUnderParallelization(4))
        %set(f,'PaperOrientation', 'landscape');
        %<-Workaround: funktioniert nicht, da Bild von Matlab nicht mitgedreht wird (nur die Seite), also lieber PDF manuell drehen im Viewer:
          set(f,'PaperOrientation', 'portrait');
      else
        set(f,'PaperOrientation', 'portrait');
      end
      set(f,'PaperSize',enforceFigPosForExportUnderParallelization(3:4));
    end
    set(f,'PaperPositionMode','auto');
    set(f,'InvertHardcopy','on');
    %Set the paper position to be compatible with the screen aspect ratio so that we can use the same normalized units on both output devices:
      set(f,'PaperUnits','centimeters');
      set(f,'PaperPositionMode','manual');        
      paperSize = get(f,'PaperSize');
      %Set the paper position to the box of the figure's actual size (1:1) and center the box on the current paper type:
        manualCenteredPaperPosition = [0,0,enforceFigPosForExportUnderParallelization(3:4)];
          if(manualCenteredPaperPosition(3)<paperSize(1)+1e-4)
            manualCenteredPaperPosition(1) = max(0,(paperSize(1)-manualCenteredPaperPosition(3))/2);
          else
            if(any(ismember({'eps','pdf'},exportFormats)))
              warning('[saveFig(%s)]: the figure width is %f cm, but the paper width is only %f cm; use paper type "1:1" or select a paper type with increased size instead of the current paper type "%s" to avoid the exported figure being cropped by its bounding box.', fullfile(targetDir,[filePrefix,'.*']), manualCenteredPaperPosition(3), paperSize(1), userDefinedPaperType);
            end
          end
          if(manualCenteredPaperPosition(4)<paperSize(2)+1e-4)
            manualCenteredPaperPosition(2) = max(0,(paperSize(2)-manualCenteredPaperPosition(4))/2);
          else
            if(any(ismember({'eps','pdf'},exportFormats)))
              warning('[saveFig(%s)]: the figure height is %f cm, but the paper height is only %f cm; use paper type "1:1" or select a paper type with increased size instead of the current paper type "%s" to avoid the exported figure being cropped by its bounding box.', fullfile(targetDir,[filePrefix,'.*']), manualCenteredPaperPosition(4), paperSize(2), userDefinedPaperType);
            end
          end
        set(f,'PaperPosition',manualCenteredPaperPosition);
    if(bExportUnderParallelizationDebug)
      disp(sprintf('<- paper size = %s, paper position = %s', mat2str(get(f,'PaperSize')), mat2str(get(f,'PaperPosition')))); 
    end
  %% define figure export setup:
    %Note: display all available export settings via: hgexport('factorystyle') 
    Exportsetup = hgexport('factorystyle');
      %define the figure's actual size as exort width&height to exactly match/fill the paperPosition box:
        Exportsetup.Units = 'centimeters';
        Exportsetup.Width = manualCenteredPaperPosition(3);
        Exportsetup.Height = manualCenteredPaperPosition(4);
      Exportsetup.Color = 'rgb';
      %Exportsetup.Color = 'cmyk';
      Exportsetup.Background = 'w';
      Exportsetup.FixedFontSize = 3; %4;
      Exportsetup.ScaledFontSize = 100;
      Exportsetup.FontMode = 'none'; %'scaled';
      Exportsetup.FontSizeMin = 3; %4;
      Exportsetup.FixedLineWidth = 1;
      Exportsetup.ScaledLineWidth = 'auto';
      Exportsetup.LineMode = 'none';
      Exportsetup.LineWidthMin = 0.25;
      Exportsetup.PSLevel = 2;
      Exportsetup.Renderer = 'painters';
      Exportsetup.Resolution = exportDPI; %inInfo.plots.export.dpi;
      Exportsetup.Bounds = 'loose'; %'tight' %loose bounds make sure that the graphics/axes have borders like on screen.
  %% export all selected formats:
    %#1/workaround for error 'Printing of uicontrols is not supported on this platform.' when exporting a figure from within a worker:
      H_uiControls = findall(f,'Type','uicontrol');
      bRestoreUIControls = false;
      if(~isempty(H_uiControls) && isWorker())
        %Include the parents, if they are not axes/figures:
          HParents = get(H_uiControls, 'Parent');
            HParents = vertcat(HParents{:});
            HParents(strcmp('axes',get(HParents,'Type'))) = [];
            HParents(strcmp('figure',get(HParents,'Type'))) = [];
          H_uiControls = [H_uiControls; HParents];
        backupVisibleStates = get(H_uiControls, 'Visible');
        set(H_uiControls, 'Visible','off');
        bRestoreUIControls = true;
      end
    %long path support on Windows:
      if(ispc() && ~(length(targetDir)>=3 && strcmp(targetDir(1:3),'\\?')))
        %If we already have a network UNC path line \\Server\Share\..., we need to replace \\ by UNC\ to get long path support in Windows:
          bIsAlreadyANetworkPath = length(targetDir)>=2 && strcmp(targetDir(1:2),'\\');
          if(bIsAlreadyANetworkPath)
            targetDir = ['UNC',targetDir(2:end)];
          end
        targetDir = ['\\?\',targetDir];
      end
    Exs = {};
    for format=sort(exportFormats(:)'); format=format{1};
      try
        %Restore renderer property for next export format, if it was changed in a previous iteration:
          set(f,'Renderer',backupRenderer);
        switch(format)
          case 'eps'
            exportedFiles{end+1} = fullfile(targetDir,[filePrefix,'.eps']);
            %!Only the painters renderer can export in vectorized format; enforce it here:
              set(f...
                ,'RendererMode','manual'...
                ,'Renderer','painters'...
              );
            %Apply export setup:
              Exportsetup.Format = 'eps';
              Exportsetup.Renderer = 'painters';
              setappdata(f,'Exportsetup',Exportsetup);
            %Launch export:                                                                     
              if(~isempty(get(f,'ResizeFcn')))
                try 
                  ilx(get(f,'ResizeFcn'), f, []);
                catch ex
                  warning('Error while executing resize function of figure: %s', ex.message);
                end
              end
              bSaved = false; nRetries = 0;
                while(nRetries <= 7)
                  try 
                    hgexport(...
                       f ...
                      ,exportedFiles{end}...
                      ,Exportsetup ...
                    );
                    bSaved = true;
                    break;
                  catch subex %Workaround: for unknown reson sometimes we get "File not found or permission denied" saving to the network, but not so after a short pause:
                    warning('Could not export plot to [%s]; retrying in 25 seconds...\nError details: %s\n', exportedFiles{end}, subex.message);
                    beep;
                    drawnow;                    
                    pause(25);
                    nRetries = nRetries + 1;
                    continue;
                  end
                end
              if(~bSaved) error('Could not save %s after %d retries.', exportedFiles{end}, nRetries); end
          case 'pdf'
            exportedFiles{end+1} = fullfile(targetDir,[filePrefix,'.pdf']);
            %!Only the painters renderer can export in vectorized format; enforce it here:
              set(f...
                ,'RendererMode','manual'...
                ,'Renderer','painters'...
              );
            %Apply export setup:
              Exportsetup.Format = 'pdf';
              Exportsetup.Renderer = 'painters';
              setappdata(f,'Exportsetup',Exportsetup);
            %Workaround: Matlab exports images with one patch per pixel, but PDF viewers usually anti-alias nearby patches causing wrong/unsharp colors.
            %=> temporarily replace all images by higher-res ones:
              bSubsampleImagesBeforeExport = true;
              if(bSubsampleImagesBeforeExport)
                imageH = findall(f,'Type','image');
                subsample = 8;
                backupImageData = {};
                for iih=1:length(imageH)
                  ih=imageH(iih);
                  backupImageData{iih} = {get(ih,'XData'), get(ih,'YData'), get(ih,'CData')};
                  if(size(backupImageData{iih}{3},1)==1 || size(backupImageData{iih}{3},2)==1) continue; end; %dont subsample images with only one pixel row/col like color bars.
                  if(size(backupImageData{iih}{3},1)>=100 && size(backupImageData{iih}{3},2)>=100) continue; end; %dont subsample large images
                  xAxis = linspace(backupImageData{iih}{1}(1), backupImageData{iih}{1}(2), max(2,size(get(ih,'CData'),2)));
                    dX = diff(xAxis); dX(end+1) = dX(end);
                  xiAxis = cellfun(@(i)linspace(xAxis(i)-dX(i)/2,xAxis(i)+dX(i)/2,subsample+1),num2cell(1:length(xAxis)),'UniformOutput',false);
                    xiAxis = unique(horzcat(xiAxis{:}));
                  yAxis = linspace(backupImageData{iih}{2}(1), backupImageData{iih}{2}(2), max(2,size(get(ih,'CData'),1)))';
                    dY = diff(yAxis); dY(end+1) = dY(end);
                  yiAxis = cellfun(@(i)linspace(yAxis(i)-dY(i)/2,yAxis(i)+dY(i)/2,subsample+1),num2cell(1:length(yAxis)),'UniformOutput',false);
                    yiAxis = unique(horzcat(yiAxis{:}))';
                  subsampledCData = interp2(...
                    xAxis, yAxis ...
                   ,ilv(get(ih,'CData'),@(M)iif(size(get(ih,'CData'),2)==1,[M,M],M),@(M)iif(size(get(ih,'CData'),1)==1,[M;M],M)) ...
                   ,xiAxis, yiAxis ...
                   ,'nearest' ... %we want to create a sharp image.
                   ,NaN ... %as interp2 cannot extrapolate, we pad border half border pixels below:
                  );
                  %Pad border half pixels:
                    subsampledCData(1:subsample/2,:) = subsampledCData(subsample/2+1:subsample,:);
                    subsampledCData(end-subsample/2+1:end,:) = subsampledCData(end-subsample+1:end-subsample/2,:);
                    subsampledCData(:,1:subsample/2) = subsampledCData(:,subsample/2+1:subsample);
                    subsampledCData(:,end-subsample/2+1:end) = subsampledCData(:,end-subsample+1:end-subsample/2);
                  set(ih,'CData',subsampledCData,'XData',xiAxis,'YData',yiAxis);
                end
              end
            %Launch export:                                                                     
              if(~isempty(get(f,'ResizeFcn')))
                try 
                  ilx(get(f,'ResizeFcn'), [], []);
                catch ex
                  warning('Error while executing resize function of figure: %s', ex.message);
                end
              end
              bSaved = false; nRetries = 0;
                while(nRetries <= 7)
                  try 
                    hgexport(...
                       f ...
                      ,exportedFiles{end}...
                      ,Exportsetup ...
                    );
                    bSaved = true;
                    break;
                  catch subex %Workaround: for unknown reson sometimes we get "File not found or permission denied" saving to the network, but not so after a short pause:
                    pause(10);
                    nRetries = nRetries + 1;
                    continue;
                  end
                end
              if(~bSaved) error('Could not save %s after %d retries.', exportedFiles{end}, nRetries); end
            %Undo image subsampling:
              if(bSubsampleImagesBeforeExport)
                for iih=1:length(imageH)
                  ih=imageH(iih);
                  set(ih,'CData',backupImageData{iih}{3},'XData',backupImageData{iih}{1},'YData',backupImageData{iih}{2});
                end
              end
          case 'png'
            exportedFiles{end+1} = fullfile(targetDir,[filePrefix,'.png']);
            %Apply export setup:
              Exportsetup.Format = 'png';
              Exportsetup.Renderer = 'auto'; %'zbuffer';
              Exportsetup.Resolution = exportDPI;
              if(isWorker() && Exportsetup.Resolution>get(0,'ScreenPixelsPerInch'))
                warning('As Matlab exports from a worker to png format via a temporary ps file that will always be converted at 72dpi, the desired resolution cannot be reached; use eps format instead and convert via GhostScript later.');
              end
              setappdata(f,'Exportsetup',Exportsetup);
            %Launch export:
              bSaved = false; nRetries = 0;
                while(nRetries <= 7)
                  try 
                    hgexport(...
                       f ...
                      ,exportedFiles{end}...
                      ,Exportsetup ...
                    );
%                     %Only the print command rescales such that the fonts get smaller:
%                       print(...
%                         f ...
%                        ,exportedFiles{end} ...
%                        ,['-d','png'] ...
%                        ,['-r',num2str(exportDPI)] ...
%                       );
                    bSaved = true;
                    break;
                  catch subex %Workaround: for unknown reson sometimes we get "File not found or permission denied" saving to the network, but not so after a short pause:
                    warning('Could not save %s; retrying in 10 seconds... Details: %s', exportedFiles{end}, subex.message);
                    pause(10);
                    nRetries = nRetries + 1;
                    continue;
                  end
                end
              if(~bSaved) error('Could not save %s after %d retries.', exportedFiles{end}, nRetries); end
          case 'fig'
            exportedFiles{end+1} = fullfile(targetDir,[filePrefix,'.fig']);
            %Launch export:
              hgsave(f, exportedFiles{end});
          otherwise
            error(['image exporting format "',format,'" not implemented']);
        end
      catch ex; Exs{end+1}=ex;
      end
    end
  %% Restore figure's original properties:
    %setappdata(f,'Exportsetup',[]);
    set(f, 'Renderer', backupRenderer);
    set(f, 'Units', backupUnits); 
    set(f, 'Position', backupPos);       
    set(f, 'Visible', backupVisibility);
    %#1/restore:
      if(bRestoreUIControls)
        set(H_uiControls,{'Visible'},backupVisibleStates);
      end
    %Throw error, if any:
      if(~isempty(Exs))
        Exs = Exs{1};
        throw(Exs);
      end
end

