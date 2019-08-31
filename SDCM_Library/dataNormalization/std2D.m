%ABSTRACT
% Library function for SDCM. Standardize in 2D (both rows and cols) to get rid of tumor load and gene affinity effects as preparation for smooth neutralization (form of whitening).
  function SD2D4L2Rs = std2D(L2Rs,fcnVar,epsilon)
    %Initialize:
      %if(nargin<2) fcnVar = @(X,dim)nanvar(X,0,dim); end
      %if(nargin<2) @(X,dim)fcnMean((X-0).^2,dim); end %precenteredVariance (the power of the signal)
      if(nargin<3) epsilon = sqrt(fcnVar(L2Rs(:),1))/100; end
      stdL2Rs = L2Rs;
      nSDIterations = 0;
      maxSDIterations = 500;

    %While not standardized:
      while(true); nSDIterations = nSDIterations + 1;
        if(nSDIterations > maxSDIterations) %Note: did not converge in 100 iterations for very high rates of NaNs in L2Rs because the epsilon was also tiny (1e-6); this can be ignored, but print a warning:
          warning('2D SDs did not converge after %d iterations; remaining variance differences <= %1.1e; now accepting AS IS.', maxSDIterations, max(d1,d2));
          break;
        end
        %calculate variances:
          columnVariances = fcnVar(stdL2Rs,1);
          rowVariances = fcnVar(stdL2Rs,2);
          %define the 2D variances as the sqrt of the corresponding row and column variances:
            twoDimVariances = sqrt(bsxfun(@times, rowVariances, columnVariances));
          %flat vectors support: if some genes are flat already, their variances are zero. For these cases, divide by 1 instead by 0 to avoid NaNs.
            twoDimVariances(twoDimVariances==0) = 1; 
          %missing values support: for all_NaN vectors, define the variance as one:
            twoDimVariances(isnan(twoDimVariances)) = 1; 
        %smooth tranition to convergence target via sqrt (do not go directly the full step to prevent numeric oscillations in some cases):
          stdL2Rs = stdL2Rs./sqrt(twoDimVariances); 
        %check convergence:
          d1 = nanmax(abs(nanmax(twoDimVariances)-nanmin(twoDimVariances)));
          d2 = nanmax(abs(nanmax(twoDimVariances,[],2)-nanmin(twoDimVariances,[],2)));
          if(max(d1,d2)<epsilon)
            SDCM_printStatus(3,'      <- combined 2D SDs converged after %d iterations.\n', nSDIterations); drawnow;
            break;
          else
            SDCM_printStatus(4,'      <- remaining max distance from a row respectively column SD of one: %f\n', max(d1,d2)); drawnow;
          end
      end

    %Output:
      stdL2Rs(stdL2Rs==0) = 1; %for zero L2Rs, divide by 1 instead by 0 to avoid NaNs.
      SD2D4L2Rs = L2Rs./stdL2Rs;
      SD2D4L2Rs(SD2D4L2Rs==0) = 1; %for zero L2Rs, divide by 1 instead by 0 to avoid NaNs.
  end
  %<-Note: geometrically viewed, std(X)==1 via fcnVar(X)=precenteredVariance(X)==1 is equivalent to euclid(X)==1. This means all column vectors 
  %        are projected on the R^nG unit sphere and simultanously all row vectors are projected on the R^nP unit sphere. 
  %        This results in "as much angular space dOmega as possible per effect" and removes any influence from length scales.
