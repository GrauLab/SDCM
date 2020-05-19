%ABSTRACT
% Library function for SDCM. Singal smoothening in the space of signature 
% strengths based on 2D Fourier transforms.

function [smoothSignal, fcn_smoothSignal4orderMetricCoordinates] = smoothenInSignatureStrengthsSpace(...
   orderMetric4G_orig, orderMetric4P_orig... %order metric (in the default, the signature gene and sample strengths, i.e. their projetions on the respective signature axes)
  ,signal ...
  ,inInfo ...
  ,bSupressStatusmessages ... %no need to show status messages for every bimonotonic regression iteration; show only once per high level signature curve regression iteration.
  ,fullSignal_orig, signatureDefinition... 
) 
  %% Initialize:
    nG = size(signal,1);
    nP = size(signal,2);      
    fcnNumericTargetPrecision = iif(strcmp(inInfo.preprocessing.numericTargetPrecision,'single'),@single,@double);
    global eDState;
    sL = min(eDState.current.multiPassesLevel, inInfo.searchStrategy.nDefinedPasses);
  %% Rescale (and downscale) the signal into the space using equidistance signal gene and sample strengths as axes/grid intervals:
    %Span the target space of signature strengths for rescaling (Note: smoothening with a constant window over equidistant bins of gene and sample strengths corresponds to a smoothening with an adaptive smoothening window size in the original matrix index space)
      if(true)
        %get the reduced dimensions of the target space (downscaling together with rescaling):
          %Note: for fast FFT we need powers of 2.
          nGp2 = inInfo.dissection.signatureSpace.resolution4G;
            if(isa(nGp2,'function_handle')) nGp2 = nGp2(nG); end;
          nPp2 = inInfo.dissection.signatureSpace.resolution4P;
            if(isa(nPp2,'function_handle')) nPp2 = nPp2(nP); end;

        %axes and meshgrid:
          %Ensure that there is enough extrapolation (at least 1 sigma):
            min4G_orig = min(orderMetric4G_orig);
            max4G_orig = max(orderMetric4G_orig);
            min4P_orig = min(orderMetric4P_orig);
            max4P_orig = max(orderMetric4P_orig);
            bSigma4GSpecifiedInAbsoluteUnits = ~isnan(inInfo.dissection.gaussianKernel.sigma4G);
              if(bSigma4GSpecifiedInAbsoluteUnits)
                sigma4G = inInfo.dissection.gaussianKernel.sigma4G;
              else
                sigma4G = (-min4G_orig+max4G_orig)/nGp2*inInfo.dissection.gaussianKernel.eqdResolutionInSigma4G;
              end
            bSigma4PSpecifiedInAbsoluteUnits = ~isnan(inInfo.dissection.gaussianKernel.sigma4P);
              if(bSigma4PSpecifiedInAbsoluteUnits)
                sigma4P = inInfo.dissection.gaussianKernel.sigma4P;
              else
                sigma4P = (-min4P_orig+max4P_orig)/nPp2*inInfo.dissection.gaussianKernel.eqdResolutionInSigma4P;
              end
            extrapFactor4G = max(0.05, 3*sigma4G/(-min4G_orig+max4G_orig));
            extrapFactor4P = max(0.05, 3*sigma4P/(-min4P_orig+max4P_orig));

          min4G = min4G_orig-(-min4G_orig+max4G_orig)*extrapFactor4G;
          max4G = max4G_orig+(-min4G_orig+max4G_orig)*extrapFactor4G;
          orderMetric4G_eqd = linspace(min4G,max4G,nGp2)';

          min4P = min4P_orig-(-min4P_orig+max4P_orig)*extrapFactor4P;
          max4P = max4P_orig+(-min4P_orig+max4P_orig)*extrapFactor4P;
          orderMetric4P_eqd = linspace(min4P,max4P,nPp2);

          %[orderMetric4P_eqd_mesh,orderMetric4G_eqd_mesh] = meshgrid(orderMetric4P_eqd,orderMetric4G_eqd);
          [orderMetric4G_eqd_mesh,orderMetric4P_eqd_mesh] = ndgrid(orderMetric4G_eqd,orderMetric4P_eqd);
      end
        if(true) %relatively fast implementation using accumarray:
          %Compute eqd target bins for every orig pixel:
            IInEqd4OrigG = interp1(orderMetric4G_eqd, (1:nGp2)', orderMetric4G_orig,'nearest','extrap');
            IInEqd4OrigP = interp1(orderMetric4P_eqd, 1:nPp2,    orderMetric4P_orig,'nearest','extrap');
            IInEqd4OrigG_mesh = repmat(IInEqd4OrigG,1,nP);
            IInEqd4OrigP_mesh = repmat(IInEqd4OrigP,nG,1);
            targetBins4SourcePixel = [IInEqd4OrigG_mesh(:),IInEqd4OrigP_mesh(:)];
          %Accumulate source pixel in non-empty eqd target bins:
            bSignalHasNaNs = any(isnan(signal(:)));
            fcnAccumulation = iif(bSignalHasNaNs, @nanmean, @mean);
            signal_eqd = accumarray(targetBins4SourcePixel, signal(:), [nGp2,nPp2], fcnAccumulation, fcnNumericTargetPrecision(NaN));
              %figure; imagesc(orderMetric4G_eqd,orderMetric4P_eqd,signal_eqd); caxis([-5,5]);
          %fill empty eqd target bins via scattered nearest interpolation:
            %NOrigInEqdBins = accumarray(targetBins4SourcePixel, ones(size(targetBins4SourcePixel,1),1), [nGp2,nPp2]);
            NOrigInEqdBins = accumarray(targetBins4SourcePixel, ~isnan(signal(:)), [nGp2,nPp2]);
              %figure; imagesc(orderMetric4G_eqd,orderMetric4P_eqd,NOrigInEqdBins); caxis([-5,5]);
            BNonEmpty = NOrigInEqdBins>0;
            %Focus the subgrid with nonzero corresponding original values (and extrapolate cells in this subgrid with no orig values, if needed):
              signal_eqd_nonEmptyRowsAndCols = signal_eqd(any(BNonEmpty,2),any(BNonEmpty,1));
              if(sum(sum(isnan(signal_eqd_nonEmptyRowsAndCols)))/numel(signal_eqd_nonEmptyRowsAndCols)>0.5)
                warning('Filling %d/%d grid cells in the signautre strength-spanned space wherein no original value lies; there are original values in other cell of the same grid rows/cols, however. Maybe the signal has a high NaNs ratio?',sum(sum(isnan(signal_eqd_nonEmptyRowsAndCols))),numel(signal_eqd_nonEmptyRowsAndCols));
              end
              if(any(any(isnan(signal_eqd_nonEmptyRowsAndCols))))
                signal_eqd_nonEmptyRowsAndCols = imputeNaNsViaNearest2DNeighbours(signal_eqd_nonEmptyRowsAndCols);
              end
              %<-Note: can happen with the zero signatureing: all zero rows/cols may extent to areas that would otherwise be NaN and excluded on caller level.
            %Interpolate values for eqd grid cells without orig values in their interval:
              if(~inInfo.internal.bEnablePostProductionCodeOptimizations) %dev.note: replaced by symmetrized version below for better reproducability; see below.
                fcn_signal_eqd = griddedInterpolant(...
                  {orderMetric4G_eqd(any(BNonEmpty,2)), orderMetric4P_eqd(any(BNonEmpty,1))} ...
                 ,signal_eqd_nonEmptyRowsAndCols ...
                 ,'nearest'... %interpolation method.
                 ,'nearest'... %extrapolation method.
                );
                signal_eqd = fcn_signal_eqd({orderMetric4G_eqd,orderMetric4P_eqd}); 
              end
              %dev.note:fft2 and ifft2 can yield differences in the order of 1e-15 for differnet machines/CPUs; These differences
              %         propagate to orderMetric4G_eqd/orderMetric4P_eqd axes/grids. We want to interpolate from the
              %         BNonEmpty subgrid to the full grid, but for any (nonEmpty,(n+1)*empty,nonEmpty) submask, 'nearest'
              %         is not well-defined for the center empty cell, as it is exactly between the the nonEmpty grid points.
              %         The 1e-15 deviations from fft2/ifft2 between machins can escalate to changes on the order of signal values
              %         in this way. (nearest can effectively takte the value from the left nonEmpty grid cell in one machine and
              %         from the right nonEmpty grid cell in the other machine...).
              %         To prevent such unwanted difficulties for reproducability, we need to symmetrize, i.e. for any
              %         empty cell exactly between two nonEmpty grid cells we must use the mean of the left and the right value
              %         instead. To facilitate this, we interpolate with 'nearest' for +dx/4 and -dx/4 and use the mean
              %         of resulting matrices (interpolated grid cells that are clearly nearer to one target cell will not change
              %         their value by this). 
                if(inInfo.internal.bEnablePostProductionCodeOptimizations)
                  dx = mean(diff(orderMetric4G_eqd));
                  dy = mean(diff(orderMetric4P_eqd));
                  fcn_signal_eqd_left_top = griddedInterpolant(...
                    {orderMetric4G_eqd(any(BNonEmpty,2))-dx/4, orderMetric4P_eqd(any(BNonEmpty,1))-dy/4} ...
                   ,signal_eqd_nonEmptyRowsAndCols ...
                   ,'nearest'... %interpolation method.
                   ,'nearest'... %extrapolation method.
                  );
                  fcn_signal_eqd_left_bottom = griddedInterpolant(...
                    {orderMetric4G_eqd(any(BNonEmpty,2))-dx/4, orderMetric4P_eqd(any(BNonEmpty,1))+dy/4} ...
                   ,signal_eqd_nonEmptyRowsAndCols ...
                   ,'nearest'... %interpolation method.
                   ,'nearest'... %extrapolation method.
                  );
                  fcn_signal_eqd_right_top = griddedInterpolant(...
                    {orderMetric4G_eqd(any(BNonEmpty,2))+dx/4, orderMetric4P_eqd(any(BNonEmpty,1))-dy/4} ...
                   ,signal_eqd_nonEmptyRowsAndCols ...
                   ,'nearest'... %interpolation method.
                   ,'nearest'... %extrapolation method.
                  );
                  fcn_signal_eqd_right_bottom = griddedInterpolant(...
                    {orderMetric4G_eqd(any(BNonEmpty,2))+dx/4, orderMetric4P_eqd(any(BNonEmpty,1))+dy/4} ...
                   ,signal_eqd_nonEmptyRowsAndCols ...
                   ,'nearest'... %interpolation method.
                   ,'nearest'... %extrapolation method.
                  );
                  signal_eqd = (...
                      fcn_signal_eqd_left_top({orderMetric4G_eqd,orderMetric4P_eqd}) ...
                    + fcn_signal_eqd_left_bottom({orderMetric4G_eqd,orderMetric4P_eqd}) ...
                    + fcn_signal_eqd_right_top({orderMetric4G_eqd,orderMetric4P_eqd}) ...
                    + fcn_signal_eqd_right_bottom({orderMetric4G_eqd,orderMetric4P_eqd}) ...
                  )/4;
                  %<-Note: Confimed that this symmetrization restores reprocudability between R2012b/i7.ffts and R2016a/Xeon.ffts (although fft2/ifft2 produce differences in the order of 1e-15 leading to likewise differences in grid cell positions.)
                end
            %figure; imagesc(orderMetric4G_eqd,orderMetric4P_eqd,signal_eqd); caxis([-5,5]);
        end
      %figure; imagesc(orderMetric4G_eqd,orderMetric4P_eqd,signal_eqd); caxis([-5,5]);
  %% Smoothening/convolution in the target space via 2D Fourier transform:
    if(~bSupressStatusmessages)
      SDCM_printStatus(2,'   - Fourier-transform the 2D signal, then fold it with the smoothing kernel by multiplication and finally use inverse Fourier transforms to get back to the original signal space.\n');
    end
    %FFT: Fourier transform the signal:
      %pad circularly (with sign flip) onto doubled size to avoid border-effects of the Fourier transform due to non-periodicity:
        bUsePadding = true;
        if(bUsePadding)
          %first2lastColAbsL2RsRatio = fcnMean(abs(Z_orig_dsct(:,1)))/fcnMean(abs(Z_orig_dsct(:,end))); %match the absolute L2Rs strength to avoid steps due to circular padding.
          first2lastColAbsL2RsRatio = 1; %!do not try to repair the jump here as in some situations this ratio estimation is so off that is confuses the fourier transform and then leads to overcompensating shifts.
          padded_signal_eqd = [-signal_eqd(:,end/2+1:end)*first2lastColAbsL2RsRatio, signal_eqd, -signal_eqd(:,1:end/2)/first2lastColAbsL2RsRatio];
          %first2lastRowAbsL2RsRatio = fcnMean(abs(Z_orig_dsct(1,:)))/fcnMean(abs(Z_orig_dsct(end,:))); %match the absolute L2Rs strength to avoid steps due to circular padding.
          first2lastRowAbsL2RsRatio = 1; %do not try to repair the jump here as in some situations this ratio estimation is so off that is confuses the fourier transform and then leads to overcompensating shifts.
          padded_signal_eqd = [-padded_signal_eqd(end/2+1:end,:)*first2lastRowAbsL2RsRatio; padded_signal_eqd; -padded_signal_eqd(1:end/2,:)/first2lastRowAbsL2RsRatio];
        end
      fourierTransformed_padded_signal_eqd = fft2(padded_signal_eqd);
    %Gaussian folding kernel:
      if(true)
        %Check for too low resolution / too small sigma:
          if(sigma4G/mean(diff(orderMetric4G_eqd))<4)
            warning('sigmaG/mean(diff(orderMetric4G_eqd))=%0.2f < 4 pixels => increase inInfo.dissection.signatureSpace.resolution4G or inInfo.dissection.sigma*; now increasing to 4 pixels per sigma.', sigma4G/mean(diff(orderMetric4G_eqd)));
            sigma4G = mean(diff(orderMetric4G_eqd))*4;
          end
          if(sigma4P/mean(diff(orderMetric4P_eqd))<4)
            warning('sigmaP/mean(diff(orderMetric4P_eqd))=%0.2f < 4 pixels => increase inInfo.dissection.signatureSpace.resolution4P or inInfo.dissection.sigma*; now increasing to 4 pixels per sigma.', sigma4P/mean(diff(orderMetric4P_eqd)));
            sigma4P = mean(diff(orderMetric4P_eqd))*4;
          end
        %Create kernel in equidistant t space:
          gaussianKernel_eqd = exp(-(((orderMetric4P_eqd_mesh-mean(orderMetric4P_eqd))/sigma4P).^2+((orderMetric4G_eqd_mesh-mean(orderMetric4G_eqd))/sigma4G).^2)/2); 
            gaussianKernel_eqd = gaussianKernel_eqd/sum(sum(gaussianKernel_eqd)); %normalize.
        %->FFT: Fourier transform the Gaussian kernel:
          %pad with zeros onto double size (needed to be size-compatible with padded signal)
            if(bUsePadding)
              %paddedZ_kern = padarray(Z_kern,size(Z_kern)/2,0);
              padded_gaussianKernel_eqd = [
                zeros(size(gaussianKernel_eqd,1)/2,2*size(gaussianKernel_eqd,2));
                [zeros(size(gaussianKernel_eqd,1),size(gaussianKernel_eqd,2)/2), gaussianKernel_eqd, zeros(size(gaussianKernel_eqd,1),size(gaussianKernel_eqd,2)/2)];
                zeros(size(gaussianKernel_eqd,1)/2,2*size(gaussianKernel_eqd,2));
              ];
            else
              padded_gaussianKernel_eqd = gaussianKernel_eqd;
            end
          fourierTransformed_padded_gaussianKernel_eqd = fft2(padded_gaussianKernel_eqd); %mit Null forsetzen
      end
    %!Multiply in Fourier space to quickly calculate the convolution in original (eqd strengths-)space:
      fourierTransformed_convolution = fourierTransformed_padded_signal_eqd.*fourierTransformed_padded_gaussianKernel_eqd;
    %inverse Fourier transform back to original eqd t space:
      padded_convolution_eqd = ifft2(fourierTransformed_convolution);
      %<-Remove roundoff complexity:
        %paddedZ_eqd_conv(imag(paddedZ_eqd_conv)<10^-4) = ilsub(real(paddedZ_eqd_conv),imag(paddedZ_eqd_conv)<10^-4);
        %assert(all(all(imag(paddedZ_eqd_conv)==0)),'paddedZ_eqd_conv is complex');
        padded_convolution_eqd = real(padded_convolution_eqd); %Performance. (previous assertion was always true).
      %->unpad:
        if(bUsePadding)
          convolution_eqd = padded_convolution_eqd(size(signal_eqd,1)/2+1:end-size(signal_eqd,1)/2, size(signal_eqd,2)/2+1:end-size(signal_eqd,2)/2);
        else
          convolution_eqd = padded_convolution_eqd;
        end
  %% interpolate the smoothed 2D baseline signal from in the eqd signature strengths space back to the original signal/matrix index space:
    %baseline for space spanned by the equidistant order metrics:
      signatureEigensignal_eqd = fcnNumericTargetPrecision(convolution_eqd); %single to save RAM.
    %interpolate signal onto non-eqd original space (repeat commands as if we were reconstructing the baseSignal4Samples from the final signature definition output):
      fcn_smoothSignal4orderMetricCoordinates = griddedInterpolant(...
        {orderMetric4G_eqd, orderMetric4P_eqd} ...
       ,signatureEigensignal_eqd ...
       ,'linear'... %interpolation method.
       ,'none'... %extrapolation method.
      );
      smoothSignal = fcn_smoothSignal4orderMetricCoordinates({orderMetric4G_orig, orderMetric4P_orig});
    %Dev plot:
      if(inInfo.internal.bDevPlots)
        baselineEstimation_f = figure('Position',[73 547 1849 421]);
          s=10;
          cLim = 0.8*[-1,1]*max(abs([min(fullSignal_orig(:)),max(fullSignal_orig(:))]));
          %Colormap:
            CM = inInfo.plots.colormap;
            %CM = sqrt(inInfo.plots.colormap);
          %Indices f�r orig:
            r4GThreshold = 0; %0.1;
            if(length(signatureDefinition.step3_regression.R4G)==length(orderMetric4G_orig))
              SII = find(abs(signatureDefinition.step3_regression.R4G(signatureDefinition.step3_regression.eigenSI))>r4GThreshold);
            else
              SII = 1:length(orderMetric4G_orig);
            end
          baselineEstimation_aCurrentL2Rs = subplot(1,s,1, 'LineWidth',2);
            imagesc(fullSignal_orig(SII,:)); %orderMetric4G_orig(SII),orderMetric4P_orig,
            colormap(CM); caxis(cLim);
            title('Z_orig','Interpreter','none');
            set(baselineEstimation_aCurrentL2Rs,'Position',get(baselineEstimation_aCurrentL2Rs,'Position').*[0.5,1,1,1]);
          if(true)
            subplot(1,s,2, 'LineWidth',2);
              imagesc(orderMetric4P_eqd,orderMetric4G_eqd,signal_eqd);
              colormap(CM); caxis(cLim);
              title(sprintf('signal_eqd\n(resolution: %d*%d)',nGp2,nPp2),'Interpreter','none');
            subplot(1,s,3, 'LineWidth',2);
              imagesc(padded_signal_eqd); %.*paddedFourierWeights_eqd
              colormap(CM); caxis(cLim);
              title('paddedZ_eqd','Interpreter','none');
            subplot(1,s,4, 'LineWidth',2);
              imagesc(gaussianKernel_eqd); caxis([0,max(max(gaussianKernel_eqd))]);
              title('(unpadded) Z_kern','Interpreter','none');
            subplot(1,s,5, 'LineWidth',2);
              imagesc(padded_convolution_eqd); 
              title(sprintf('convolution result\npaddedZ_eqd_conv'),'Interpreter','none');
              colormap(CM); caxis(cLim);
            subplot(1,s,6, 'LineWidth',2);
              imagesc(orderMetric4P_eqd,orderMetric4G_eqd,convolution_eqd); 
              colormap(CM); caxis(cLim);
              title(sprintf('unpadded convolution result\nZ_eqd_conv'),'Interpreter','none');
          end
          subplot(1,s,7, 'LineWidth',2);
            imagesc(orderMetric4P_orig,orderMetric4G_orig(SII),signal(SII,:)); 
            colormap(CM); caxis(cLim);
            title(sprintf('smooth signature in original axes\nsignatureEigensignal_orig_smooth'),'Interpreter','none');
          subplot(1,s,8, 'LineWidth',2);
            %imagesc(orderMetric4P_eqd,orderMetric4G_eqd,bimonotonic_convolution_eqd); 
            %imagesc(signatureSignalStrengths4G_eqd,signatureSignalStrengths4P_eqd,bimonotonic_signatureEigensignal_subsampled); 
            %imagesc(orderMetric4P_orig,orderMetric4G_orig(SII),signatureEigensignal_orig(SII,:)); 
            imagesc(smoothSignal(SII,:)); %leave gene numbers on axes.
            colormap(CM); caxis(cLim);
            title(sprintf('after Fourier2D smoothing\nsignatureEigensignal_orig_prepared_smooth'),'Interpreter','none');
          baselineEstimation_aRemainingL2Rs = subplot(1,s,[9,10], 'LineWidth',2);
            imagesc(fullSignal_orig(SII,:)-smoothSignal(SII,:)); 
            colormap(CM); caxis(cLim); cbh=colorbar;
            title(sprintf('remaining signal\nZ_orig-signatureEigensignal_orig_prepared_smooth'),'Interpreter','none');
            set(baselineEstimation_aRemainingL2Rs,'Position',ilv(get(baselineEstimation_aRemainingL2Rs,'Position'),@(P)[(2*P(1)+1*1)/3,P(2:end)]));
          %Equalize axes widths:
            Hs = setdiff(findobj(baselineEstimation_f,'Type','axes'),[cbh]);
            Poss = get(Hs,'Position');
            Poss = vertcat(Poss{:});
            Poss(:,3) = mean(Poss(:,3));
            Poss = num2cell(Poss,2);
            set(Hs,{'Position'},Poss);
          %Manual/Zoom to top genes/samples:
            zoom(baselineEstimation_f,'on');
            if(false)
              set([baselineEstimation_aCurrentL2Rs, baselineEstimation_aRemainingL2Rs], 'YLim',[1,100])

              set([baselineEstimation_aCurrentL2Rs, baselineEstimation_aRemainingL2Rs],'YLim',[nG-100,nG])
            end
      end
end

