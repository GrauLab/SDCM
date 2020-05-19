%ABSTRACT
%  Library function for SDCM. 
%  Iteratively enforce bi-monotonicity based on GPAV regression; parallelized.

function [bimonotonicX, signatureDefinition] = bimonotonicRegression(...
   X ... %input signal to regress
  ,slopeSigns4G, slopeSigns4P ... %regress increasing or decreasing?
  ,orderMetric4G, orderMetric4P ... %signature gene/sample strengths used for interpolation of regression results for those genes/samples that are not in the performance subspace.
  ,weights2D ... %pixel weights used for regression
  ,performanceSubI4G, performanceSubJ4P ... %only regress for selected rows/cols and interpolate others
  ,signatureDefinition ... 
  ,inInfo ...
) 
  %% Initialize:
    global eDState;
    %initialize output and get dimensions:
      bimonotonicX = X;
      nG = size(X,1);
      nP = size(X,2);
    %some consistency checks:
      assert(isempty(slopeSigns4G) || all(size(slopeSigns4G)==[nG,1]),'slopeSigns4G is not of size [nG,1]');
      assert(isempty(slopeSigns4P) || all(size(slopeSigns4P)==[1,nP]),'slopeSigns4P is not of size [1,nP]');
      assert(all(size(weights2D)==[nG,nP]),'weights2D is not of size [nG,nP]');
      assert(all(weights2D(:)>=0 & weights2D(:)<=1), 'weights2D must be in [0,1]');
    %Slope signs (detect automatically by regressing in both directions):
      if(isempty(slopeSigns4G) && isempty(slopeSigns4P)) 
        slopeSigns4G = nan(nG,1); %determine all monotonicity directions robustly, adaptive to the data (will only cause a factor 2 in the first iteration, as the code stores robustly identified directions).
        slopeSigns4P = nan(1,nP); %determine all monotonicity directions robustly, adaptive to the data (will only cause a factor 2 in the first iteration, as the code stores robustly identified directions).
      end
      bFindOptimalMonotonicDirectionInEachIteration = false; %This is a factor two slower per iteration, but may converge in much fewer iterations.
      BRecheckSlopeSigns4G = false(nG,1);
      BRecheckSlopeSigns4P = false(1,nP);
    %GPAV version selector:
      numericTargetPrecision = class(X);
      bMexVersionAvailable = false && ~isempty(which('GPAV_static_mex'));
      %if(~bMexVersionAvailable)
      %  warning('Compiled GPAV_static_mex not found; falling back to the platform-independent interpreted GPAV version. For performance reasons it is recommended to mex/compile the static version of the GPAV algorithm with Matlab for your local platform and bitness.');
      %end
      if(bMexVersionAvailable && ~strcmp(numericTargetPrecision,'double'))
        warning('converting single->double; TODO//performance: create mex version of GPAV for single precision arithmetic');
        bimonotonicX = double(bimonotonicX);
        slopeSigns4G = double(slopeSigns4G);
        slopeSigns4P = double(slopeSigns4P);
        weights2D = double(weights2D);
      end
    %Performance subspace:
      if(isempty(performanceSubI4G)) performanceSubI4G=(1:nG)'; end
      if(isempty(performanceSubJ4P)) performanceSubJ4P=1:nP; end
      performanceSubI4G = unique(performanceSubI4G);
      performanceSubJ4P = unique(performanceSubJ4P);
      %Remove full NaN rows/cols from the performance subspace (will stay NaNs anyway)
        performanceSubI4G = setdiff(performanceSubI4G, find(all(isnan(X),2)));
        performanceSubJ4P = setdiff(performanceSubJ4P, find(all(isnan(X),1)));
      %Remove exact zeros rows/cols from the performance subspace:
        performanceSubI4G = setdiff(performanceSubI4G, find(all(X==0,2)));
        performanceSubJ4P = setdiff(performanceSubJ4P, find(all(X==0,1)));
      %Remove zero-weighted rows/cols from the performance subspace:
        performanceSubI4G = setdiff(performanceSubI4G, find(all(weights2D==0,2)));
        performanceSubJ4P = setdiff(performanceSubJ4P, find(all(weights2D==0,1)));

      assert(length(performanceSubI4G)>=2, 'code validation: at least two genes must be in the regression performance subspace');
      assert(length(performanceSubJ4P)>=2, 'code validation: at least two samples must be in the regression performance subspace');
      %the corresponding adjacency matrix for GPAV:
        persistent cache4adjecency; %do not recompute the adjacency matrix if dimensions stay the same.
        if(true)
          nG_perfSubI4G = length(performanceSubI4G);
            assert(all(diff(performanceSubI4G)>0), 'performanceSubI4G must be a unique and monotonically increasing subselection of 1:nG');
            %adj4G = sparse( triu(ones(size(Z_eqd_conv,1)),1)-triu(ones(size(Z_eqd_conv,1)),2) );
            adjI4G = 1:nG_perfSubI4G-1;
            adjJ4G = 2:nG_perfSubI4G;
            adj4G = sparse(adjI4G, adjJ4G, ones(length(adjI4G),1), nG_perfSubI4G, nG_perfSubI4G);

            adj4GRows4N = nan(nG_perfSubI4G,1);
            if(bMexVersionAvailable)
              if(isempty(cache4adjecency) || ~isstruct(cache4adjecency) || ~isfield(cache4adjecency,'nP_perfSubI4G') || cache4adjecency.nP_perfSubI4G~=nG_perfSubI4G)
                for j=1:nG_perfSubI4G
                  [rowI, colI] = find(adj4G(:,j));
                  if(isempty(rowI)) continue; end
                  if(isscalar(rowI))
                    adj4GRows4N(j) = rowI;
                    continue;
                  end
                  error('cannot prebuild adjI4N for adjacency matrices with more than one predecessor (as is the case for column j=%d)', j);
                end
                adj4GRows4N = int32(adj4GRows4N);
                cache4adjecency.adj4GRows4N = adj4GRows4N;
                cache4adjecency.nP_perfSubI4G=nG_perfSubI4G;
              else
                adj4GRows4N = cache4adjecency.adj4GRows4N;
              end
            end

          nP_perfSubJ4P = length(performanceSubJ4P);
            assert(all(diff(performanceSubJ4P)>0), 'performanceSubJ4P must be a unique and monotonically increasing subselection of 1:nP');
            %adj4P = sparse( triu(ones(size(Z_eqd_conv,2)),1)-triu(ones(size(Z_eqd_conv,2)),2) );
            adjI4P = 1:nP_perfSubJ4P-1;
            adjJ4P = 2:nP_perfSubJ4P;
            adj4P = sparse(adjI4P, adjJ4P, ones(length(adjI4P),1), nP_perfSubJ4P, nP_perfSubJ4P);

            adj4PRows4N = nan(nP_perfSubJ4P,1);
            if(bMexVersionAvailable)
              if(isempty(cache4adjecency) || ~isstruct(cache4adjecency) || ~isfield(cache4adjecency,'nP_perfSubJ4P') || cache4adjecency.nP_perfSubJ4P~=nP_perfSubJ4P)
                for j=1:nP_perfSubJ4P
                  [rowI, colI] = find(adj4P(:,j));
                  if(isempty(rowI)) continue; end
                  if(isscalar(rowI))
                    adj4PRows4N(j) = rowI;
                    continue;
                  end
                  error('cannot prebuild adjI4N for adjacency matrices with more than one predecessor (as is the case for column j=%d)', j);
                end
                adj4PRows4N = int32(adj4PRows4N);
                cache4adjecency.adj4PRows4N = adj4PRows4N;
                cache4adjecency.nP_perfSubJ4P=nP_perfSubJ4P;
              else
                adj4PRows4N = cache4adjecency.adj4PRows4N;
              end
            end
        end
      %the corresponding interpolation space and interpolation flags:
        BInterpolationSpace4G = ~ismember((1:nG)',performanceSubI4G);
        BInterpolationSpace4P = ~ismember((1:nP), performanceSubJ4P);
        bUseLinearInterpolationForNaNIntervals = true;
        B4G_distinct = [true; diff(orderMetric4G)>eps(max(orderMetric4G))]; %the support grid for interpolation can only have distinct values.
        B4P_distinct = [true, diff(orderMetric4P)>eps(max(orderMetric4P))]; %the support grid for interpolation can only have distinct values.
    %Weights:
      %Weights for gene/sample 1D regressions:
        W4G = weights2D;
        W4P = weights2D;
        %Note: We need 2D weights, because if a column j with high W4P(j) has only few rows with high W4G(i), these rows will be overwritten by values with high W4G(i) from other columns and this value will then be spread in the next iteration because of high W4P(j). However, since W4P(j) was high in the first place, no rows of column j should be overwritten, especially zeros should stay zeros! This can only be realized by a 2D weights concept.
      %Weights for the convergence check:
        W2D_initial = weights2D;
        %do not count exact zero vectors in the convergence test:
          BExactZeros4G = all(X==0,2) | all(weights2D==0,2);
          BExactZeros4P = all(X==0,1) | all(weights2D==0,1);
          W2D_initial(BExactZeros4G,:) = 0; %0.5; %sum(BInterpolationSpace4P)/nP;
          W2D_initial(:,BExactZeros4P) = 0; %0.5; %sum(BInterpolationSpace4G)/nG;
        %do not use the interpolation space in the convergence test (only truly regressed parts that spawn the other pixels by interpolation):
          W2D_initial(BInterpolationSpace4G,:) = 0;
          W2D_initial(:,BInterpolationSpace4P) = 0;
        %Setup the convergence epsilon:
          sd4epsilon = eDState.noiseEstimation.sd4epsilon;
          dL2RsConvergenceEpsilon = sd4epsilon*inInfo.dissection.bimonotonicRegression.epsilon4convergenceInSDs;
  %% Bimonotonic regression loop
    bBiMonotonic = false;
    nBimonotonicConvergenceLoopIterations = 0;
    while(~bBiMonotonic) 
      nBimonotonicConvergenceLoopIterations=nBimonotonicConvergenceLoopIterations+1;
      %some consistency checks:
        assert(~any(any(isnan(W4G(:,performanceSubJ4P)))), 'code validation: at least one W4G was NaN');
        assert(~any(any(isnan(W4P(performanceSubI4G,:)))), 'code validation: at least one W4P was NaN');
        assert(~any(isinf(W4G(:))), 'code validation: at least one W4G was Inf');
        assert(~any(isinf(W4P(:))), 'code validation: at least one W4P was Inf');
        assert(all(size(slopeSigns4G)==[nG,1]), 'slopeSigns4G must be of size [nG,1]');
        assert(all(size(slopeSigns4P)==[1,nP]), 'slopeSigns4P must be of size [1,nP]');
      %% monotonic regression of all genes along the signature sample order:
        SDCM_printStatus(3,'      -> monotonic regression along sample dimension for every gene (iteration %d)...\n', nBimonotonicConvergenceLoopIterations);        
        %Performance: the workload between rows can vary strongly. Sometimes this causes one worker to lag behind. => randomly permute the rows und unpermute later to evently distribute the workload:
          permI4workloadDistribution = randperm(nG);
            permI4workloadDistribution(~ismember(permI4workloadDistribution,performanceSubI4G)) = [];
          bimonotonic_convolution_eqd_permSubI4G = bimonotonicX(permI4workloadDistribution,performanceSubJ4P);
          slopeSigns4G_permSubI4G = slopeSigns4G(permI4workloadDistribution,:);

        Z_eqd_conv_enforcedMonotonicity4G_permSubI4G = nan(length(permI4workloadDistribution),length(performanceSubJ4P));
        W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G = W4P(permI4workloadDistribution,performanceSubJ4P);
        parfor i=1:length(permI4workloadDistribution)
          increasingRow = []; W4increasingRow = []; decreasingRow = []; W4decreasingRow = []; %initialize to prevent temporary variable compilation warnings.
          vz0 = sign(slopeSigns4G_permSubI4G(i));
            if(vz0~=+1 && vz0~=-1) vz0=NaN; end
          vz = +1;
            if(isnan(vz0) || vz==vz0)
              if(bMexVersionAvailable)
                [increasingRow, W4increasingRow] = GPAV_static_mex(...
                  int32(nP_perfSubJ4P) ...
                 ,adj4PRows4N ...
                 ,vz*bimonotonic_convolution_eqd_permSubI4G(i,:) ...
                 ,W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) ...
                 ,true ...
                );
              else
                [increasingRow, W4increasingRow] = GPAV(...
                  adj4P ...
                 ,vz*bimonotonic_convolution_eqd_permSubI4G(i,:) ...
                 ,W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) ...
                 ,[],true ...
                );
              end
            end
          vz = -1;
            if(isnan(vz0) || vz==vz0)
              if(bMexVersionAvailable)
                [decreasingRow, W4decreasingRow] = GPAV_static_mex(...
                  int32(nP_perfSubJ4P) ...
                 ,adj4PRows4N ...
                 ,vz*bimonotonic_convolution_eqd_permSubI4G(i,:) ...
                 ,W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) ...
                 ,true ...
                );
              else
                [decreasingRow, W4decreasingRow] = GPAV(...
                  adj4P ...
                 ,vz*bimonotonic_convolution_eqd_permSubI4G(i,:) ...
                 ,W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) ...
                 ,[],true ...
                );
              end
            end
          if(isnan(vz0))
            sdRatio = nanstd(increasingRow)/(nanstd(decreasingRow)+eps);
            if(sdRatio>1) vz0=+1; else vz0=-1; end
            if(sdRatio>2 || sdRatio<1/2) slopeSigns4G_permSubI4G(i) = vz0; end %Performance: remember robustly identified directions of monotonicity for later iterations.
          end 
          if(vz0==+1)
            vz=+1;
            Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) = vz*increasingRow;
            W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) = W4increasingRow;
          elseif(vz0==-1)
            vz=-1;
            Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) = vz*decreasingRow;
            W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) = W4decreasingRow;
          else
            error('code validation: vz0=%f unexpected.', vz0);
          end
          if(bUseLinearInterpolationForNaNIntervals) %GPAV can handle NaNs, but extrapolates into NaN areas constantly instead of linearly => interpolate linearly instead:
            BWasNaN = isnan(bimonotonic_convolution_eqd_permSubI4G(i,:));
            if(any(BWasNaN))
              pflc_GPAVresult = Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:);
                pflc_GPAVresult(1,BWasNaN) = NaN;
                BValid = ~isnan(pflc_GPAVresult);
                BInterpolationBase4P = B4P_distinct(performanceSubJ4P)&BValid;
                if(sum(BInterpolationBase4P)>2) %otherwise leave NaNs for the next iteration (all NaN cols/rows cannot be interpolated to non-NaN special case).
                  fcn_linearIP = griddedInterpolant(...
                    {ilsub(orderMetric4P(performanceSubJ4P),BInterpolationBase4P)} ...
                   ,pflc_GPAVresult(1,BInterpolationBase4P) ...
                   ,'linear'... %interpolation method.
                   ,'nearest'... %extrapolation method.
                  );
                    Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) = fcn_linearIP({orderMetric4P(performanceSubJ4P)}); 
                end
              pflc_GPAVresult = W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:);
                pflc_GPAVresult(1,BWasNaN) = NaN;
                BValid = ~isnan(pflc_GPAVresult);
                BInterpolationBase4P = B4P_distinct(performanceSubJ4P)&BValid;
                if(sum(BInterpolationBase4P)>2) %otherwise leave NaNs for the next iteration (all NaN cols/rows cannot be interpolated to non-NaN special case).
                  fcn_linearIP = griddedInterpolant(...
                    {ilsub(orderMetric4P(performanceSubJ4P),BInterpolationBase4P)} ...
                   ,pflc_GPAVresult(1,BInterpolationBase4P) ...
                   ,'linear'... %interpolation method.
                   ,'nearest'... %extrapolation method.
                  );
                    W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G(i,:) = fcn_linearIP({orderMetric4P(performanceSubJ4P)}); 
                end
            end
          end
        end
        Z_eqd_conv_enforcedMonotonicity4G = nan(size(bimonotonicX));
          Z_eqd_conv_enforcedMonotonicity4G(permI4workloadDistribution,performanceSubJ4P) = Z_eqd_conv_enforcedMonotonicity4G_permSubI4G;
        W4Z_eqd_conv_enforcedMonotonicity4G = W4P; %monotonically regressing ->4G is based on sample weights and shifts them to other pixels in the same row => use W4P as base weights for pixels outside of the performance subspace.
          W4Z_eqd_conv_enforcedMonotonicity4G(permI4workloadDistribution,performanceSubJ4P) = W4Z_eqd_conv_enforcedMonotonicity4G_permSubI4G;
        slopeSigns4G(permI4workloadDistribution,:) = slopeSigns4G_permSubI4G;
      %% monotonic regression for all samples along the signature gene order:
        SDCM_printStatus(3,'      -> monotonic regression along gene dimension for every sample (iteration %d)...\n', nBimonotonicConvergenceLoopIterations);
        %Performance: the workload between rows can vary strongly. Sometimes this causes one worker to lag behind. => randomly permute the rows und unpermute later to evently distribute the workload:
          permJ4workloadDistribution = randperm(nP);
            permJ4workloadDistribution(~ismember(permJ4workloadDistribution,performanceSubJ4P)) = [];
          bimonotonic_convolution_eqd_permSubJ4P = bimonotonicX(performanceSubI4G,permJ4workloadDistribution);
          slopeSigns4P_permSubJ4P = slopeSigns4P(:,permJ4workloadDistribution);

        Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P = nan(length(performanceSubI4G),length(permJ4workloadDistribution));
        W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P = W4G(performanceSubI4G,permJ4workloadDistribution);
        parfor j=1:length(permJ4workloadDistribution)
          increasingCol = []; W4increasingCol = []; decreasingCol = []; W4decreasingCol = []; %initialize to prevent temporary variable compilation warnings.
          vz0 = sign(slopeSigns4P_permSubJ4P(j));
            if(vz0~=+1 && vz0~=-1) vz0=NaN; end
          vz = +1;
            if(isnan(vz0) || vz==vz0)
              if(bMexVersionAvailable)
                [increasingCol, W4increasingCol] = GPAV_static_mex(...
                  int32(nG_perfSubI4G) ...
                 ,adj4GRows4N ...
                 ,vz*bimonotonic_convolution_eqd_permSubJ4P(:,j)' ...
                 ,W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j)' ...
                 ,true ...
                );
              else
                [increasingCol, W4increasingCol] = GPAV(...
                  adj4G ...
                 ,vz*bimonotonic_convolution_eqd_permSubJ4P(:,j)' ...
                 ,W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j)' ...
                 ,[],true ...
                );
              end
            end
          vz = -1;
            if(isnan(vz0) || vz==vz0)
              if(bMexVersionAvailable)
                [decreasingCol, W4decreasingCol] = GPAV_static_mex(...
                  int32(nG_perfSubI4G) ...
                 ,adj4GRows4N ...
                 ,vz*bimonotonic_convolution_eqd_permSubJ4P(:,j)' ...
                 ,W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j)' ...
                 ,true ...
                );
              else
                [decreasingCol, W4decreasingCol] = GPAV(...
                  adj4G ...
                 ,vz*bimonotonic_convolution_eqd_permSubJ4P(:,j)' ...
                 ,W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j)' ...
                 ,[],true ...
                );
              end
            end
          if(isnan(vz0))
            sdRatio = nanstd(increasingCol)/(nanstd(decreasingCol)+eps);
            if(sdRatio>1) vz0=+1; else vz0=-1; end
            if(sdRatio>2 || sdRatio<1/2) slopeSigns4P_permSubJ4P(j) = vz0; end %Performance: remember robustly identified directions of monotonicity for later iterations.
          end 
          if(vz0==+1)
            vz=+1;
            Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j) = vz*increasingCol;
            W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j) = W4increasingCol;
          elseif(vz0==-1)
            vz=-1;
            Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j) = vz*decreasingCol;
            W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j) = W4decreasingCol;
          else
            error('code validation: vz0=%f unexpected.', vz0);
          end
          if(bUseLinearInterpolationForNaNIntervals) %GPAV can handle NaNs, but extrapolates into NaN areas constantly instead of linearly => interpolate linearly instead:
            BWasNaN = isnan(bimonotonic_convolution_eqd_permSubJ4P(:,j));
            if(any(BWasNaN))
              pflc_GPAVresult = Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j);
                pflc_GPAVresult(BWasNaN,1) = NaN;
                BValid = ~isnan(pflc_GPAVresult);
                BInterpolationBase4G = B4G_distinct(performanceSubI4G)&BValid;
                if(sum(BInterpolationBase4G)>2) %otherwise leave NaNs for the next iteration (all NaN cols/rows cannot be interpolated to non-NaN special case).
                  fcn_linearIP = griddedInterpolant(...
                    {ilsub(orderMetric4G(performanceSubI4G),BInterpolationBase4G)} ...
                   ,pflc_GPAVresult(BInterpolationBase4G,1) ...
                   ,'linear'... %interpolation method.
                   ,'nearest'... %extrapolation method.
                  );
                  Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j) = fcn_linearIP({orderMetric4G(performanceSubI4G)}); 
                end
              pflc_GPAVresult = W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j);
                pflc_GPAVresult(BWasNaN,1) = NaN;
                BValid = ~isnan(pflc_GPAVresult);
                BInterpolationBase4G = B4G_distinct(performanceSubI4G)&BValid;
                if(sum(BInterpolationBase4G)>2) %otherwise leave NaNs for the next iteration (all NaN cols/rows cannot be interpolated to non-NaN special case).
                  fcn_linearIP = griddedInterpolant(...
                    {ilsub(orderMetric4G(performanceSubI4G),BInterpolationBase4G)} ...
                   ,pflc_GPAVresult(BInterpolationBase4G,1) ...
                   ,'linear'... %interpolation method.
                   ,'nearest'... %extrapolation method.
                  );
                  W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P(:,j) = fcn_linearIP({orderMetric4G(performanceSubI4G)}); 
                end
            end
          end
        end
        Z_eqd_conv_enforcedMonotonicity4P = nan(size(bimonotonicX));
          Z_eqd_conv_enforcedMonotonicity4P(performanceSubI4G,permJ4workloadDistribution) = Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P;
        W4Z_eqd_conv_enforcedMonotonicity4P = W4G; %monotonically regressing ->4P is based on gene weights and shifts them to other pixels in the same col => use W4G as base weights for pixels outside of the performance subspace.
          W4Z_eqd_conv_enforcedMonotonicity4P(performanceSubI4G,permJ4workloadDistribution) = W4Z_eqd_conv_enforcedMonotonicity4P_permSubJ4P;
        slopeSigns4P(:,permJ4workloadDistribution) = slopeSigns4P_permSubJ4P;
      %% smooth symmetric convergence:
        %store some infos for plots, if requested:
          if(isa(inInfo.export.matFile4eachSignatureDefinition.bIncludeBimonotonicConvergenceSteps,'function_handle'))
            current_bSaveIntermediateResultsForDevPlots = inInfo.export.matFile4eachSignatureDefinition.bIncludeBimonotonicConvergenceSteps(eDState.current.k);
          else
            current_bSaveIntermediateResultsForDevPlots = inInfo.export.matFile4eachSignatureDefinition.bIncludeBimonotonicConvergenceSteps;
          end
          if(current_bSaveIntermediateResultsForDevPlots) %<-save additional data maybe needed by some plots later:
            signatureDefinition.forPlots.bimonotonicConvergence.unpaddedEqdConvoluted.slopeSigns4P = slopeSigns4P;
            signatureDefinition.forPlots.bimonotonicConvergence.unpaddedEqdConvoluted.slopeSigns4G = slopeSigns4G;
            signatureDefinition.forPlots.bimonotonicConvergence.unpaddedEqdConvoluted.weights2D = weights2D;
            %signatureDefinition.forPlots.bimonotonicConvergence.unpaddedEqdConvoluted.weights4P = weights4P;
            %signatureDefinition.forPlots.bimonotonicConvergence.unpaddedEqdConvoluted.weights4G = weights4G;
            signatureDefinition.forPlots.bimonotonicConvergence.iterations(nBimonotonicConvergenceLoopIterations).initial = bimonotonicX;
            signatureDefinition.forPlots.bimonotonicConvergence.iterations(nBimonotonicConvergenceLoopIterations).monotonic4G = Z_eqd_conv_enforcedMonotonicity4G;
            signatureDefinition.forPlots.bimonotonicConvergence.iterations(nBimonotonicConvergenceLoopIterations).monotonic4P = Z_eqd_conv_enforcedMonotonicity4P;
            signatureDefinition.forPlots.bimonotonicConvergence.iterations(nBimonotonicConvergenceLoopIterations).final = (Z_eqd_conv_enforcedMonotonicity4G+Z_eqd_conv_enforcedMonotonicity4P)/2;
          end
        %Initialize:
          bimonotonicX_previous = bimonotonicX; 
          W2D = zeros(nG,nP);
          bimonotonicX = zeros(nG,nP);
        %Add results from regressing genes along ordered sample columns:
          accuW4G = (W4P*nP+W4Z_eqd_conv_enforcedMonotonicity4P*nG)/(nP+nG);
            BValid4G = ~isnan(accuW4G) ... %can happen is the weights are all-zeros for the whole column.
              & ~isnan(Z_eqd_conv_enforcedMonotonicity4G)... 
            ;
          bimonotonicX(BValid4G) = bimonotonicX(BValid4G) + 0.5*accuW4G(BValid4G).*Z_eqd_conv_enforcedMonotonicity4G(BValid4G);
          W2D(BValid4G) = W2D(BValid4G) + 0.5*accuW4G(BValid4G);
        %Add results from regressing samples along ordered gene rows:
          accuW4P = (W4G*nG+W4Z_eqd_conv_enforcedMonotonicity4G*nP)/(nG+nP);
            BValid4P = ~isnan(accuW4P) ... %can happen if the weights are all-zeros for the whole column.
              & ~isnan(Z_eqd_conv_enforcedMonotonicity4P) ...
            ;
          bimonotonicX(BValid4P) = bimonotonicX(BValid4P) + 0.5*accuW4P(BValid4P).*Z_eqd_conv_enforcedMonotonicity4P(BValid4P);
          W2D(BValid4P) = W2D(BValid4P) + 0.5*accuW4P(BValid4P);
        %estimate signal from regression only, where weights>0 are available (this serves as basis for linear interpolation below).
          BRegressionInfoAvailable = W2D>0;
          regressedX = bimonotonicX;
          regressedX(BRegressionInfoAvailable) = bimonotonicX(BRegressionInfoAvailable)./W2D(BRegressionInfoAvailable); %scale up to 100%.
        %Infer the bimonotonic signal of the signature for pixels outside of samples and genes with high regression weight:
          %<-Note: Outside of the signature focus (low correlations) we have regression weights W2D(i,j)<1 or maybe even ==0. 
          %  Hence, we need to infer the signature signal there by other means than regression.
          %smoothen the signal (adaptive smoothening window depending on the signature strengths via rescaling & Fourier 2D transforms):
            if(true)
              [~,fnc_ip_smoothening] = smoothenInSignatureStrengthsSpace(...
                orderMetric4G, orderMetric4P ...
               ,bimonotonicX_previous ...
               ,inInfo, true ... %inInfo, bSupressStatusmessages
               ,bimonotonicX_previous, signatureDefinition... 
              );
              smoothenedX = fnc_ip_smoothening({orderMetric4G, orderMetric4P});
                %!<-funktioniert aber nicht im 3D, wenn/weil in der Mitte ein starkes inkonsistentes Signal ist, wo der Nulldurchgang des aktuellen Effekts sein sollte => hier braucht man lineare Interpolation 
                %   => nehme Minimum von beidem.
            end
          %linearly interpolate the signal:
            if(true)
              X4G = orderMetric4G(performanceSubI4G);
                BDistinct4G = [true; diff(X4G)>10*eps(max(X4G))]; %distinct values only.
              X4P = orderMetric4P(performanceSubJ4P);
                BDistinct4P = [true, diff(X4P)>10*eps(max(X4P))]; %distinct values only.
              B4G = BDistinct4G & any(BRegressionInfoAvailable(performanceSubI4G,:),2);
              B4P = BDistinct4P & any(BRegressionInfoAvailable(:,performanceSubJ4P),1);
              bCanInterpolate = sum(B4G)>=2 && sum(B4P)>=2;
              if(bCanInterpolate)
                fnc_ip_linearInterpolation = griddedInterpolant(...
                  {X4G(B4G), X4P(B4P)} ...
                 ,regressedX(performanceSubI4G(B4G),performanceSubJ4P(B4P)) ...
                 ,iif(bUseLinearInterpolationForNaNIntervals,'linear','nearest')... %interpolation method.
                 ,'nearest'... %extrapolation method.
                );
                linearlyInterpolatedX = fnc_ip_linearInterpolation({orderMetric4G, orderMetric4P});
              end
                %<-hat im 3D alleine (ohne smoothenInSignatureStrengthsSpace) gut funktioniert
            end
          %Define the final inferred signal as the minimum of smoothenedX and linearlyInterpolatedX 
            %<-Note: linearlyInterpolatedX usually infers an signature better near the zero-passing, if we lack samples with regression weight near zero. Fourier2D tends to behave more like a nearest-extrapolation instead of a linear extrapolation there.
            if(bCanInterpolate)
              interpolatedX = smoothenedX;
              BLinearIPisSmaller = abs(linearlyInterpolatedX)<abs(interpolatedX);
              interpolatedX(BLinearIPisSmaller) = linearlyInterpolatedX(BLinearIPisSmaller);
            else %special case: for a 1D matrix there is no 2D interpolation.
              interpolatedX = smoothenedX;
            end 
        %whereever regression does not provide enough info to define 100% of the signal, fill it with the interpolated signal:
          %bimonotonicX = bimonotonicX + interpolatedX.*(1-W2D); %every pixel now has 100%, i.e. the weighted average of both 1D regression directions and the interpolation results is summed completely now.
          %bimonotonicX(BRegressionInfoAvailable) = bimonotonicX(BRegressionInfoAvailable) ./ W2D(BRegressionInfoAvailable);
          %bimonotonicX = bimonotonicX.*W2D + interpolatedX.*(1-W2D);
          bimonotonicX = regressedX.*W2D + interpolatedX.*(1-W2D);
        %in the next iteration, use symmetric 2D regression weights for both 1D regression directions again:
          W4G_previous = W4G; 
          W4P_previous = W4P; 
          W4G = W2D;
          W4P = W2D;
        %re-check signs (test both decreasing and increasing 1D GPAV-regressions), whenever the average diff sign for a gene or sample flipped:
          BRecheckSlopeSigns4P = sign(nanmean(diff(Z_eqd_conv_enforcedMonotonicity4G,1,1),1))~=sign(nanmean(diff(Z_eqd_conv_enforcedMonotonicity4P,1,1),1)); %need to recheck in all henceforth iterations to not fix a potentially wrong direction.
            slopeSigns4P(1,BRecheckSlopeSigns4P) = NaN;
          BRecheckSlopeSigns4G = sign(nanmean(diff(Z_eqd_conv_enforcedMonotonicity4G,1,2),2))~=sign(nanmean(diff(Z_eqd_conv_enforcedMonotonicity4P,1,2),2)); %need to recheck in all henceforth iterations to not fix a potentially wrong direction.
            slopeSigns4G(BRecheckSlopeSigns4G) = NaN;
          if(bFindOptimalMonotonicDirectionInEachIteration)
            slopeSigns4P(:) = NaN;
            slopeSigns4G(:) = NaN;
          end
      %% check for convergence:
        SQ = (bimonotonicX_previous-bimonotonicX).^2;
        dL2RsConvergence = sqrt(nanvar(SQ(:),W2D_initial(:)));
        bBiMonotonic = dL2RsConvergence <= dL2RsConvergenceEpsilon;
        %if bi-monotonic within the preformance subspace, affitionally check convergence in the full space (based on interpolation):
          if(bBiMonotonic) 
            SQ = (bimonotonicX_previous-bimonotonicX).^2;
            dL2RsConvergence = max(dL2RsConvergence, sqrt(nanvar(SQ(:),1)));
            bBiMonotonic = dL2RsConvergence <= dL2RsConvergenceEpsilon;
          end
        SDCM_printStatus(2-double(~bMexVersionAvailable),'      <- left-over average inconsistency after bimonotonic regression iteration %d: %1.1e*SD (%s, SD=%1.1e)...\n'...
          ,nBimonotonicConvergenceLoopIterations...
          ,dL2RsConvergence/sd4epsilon...
          ,iif(bBiMonotonic, sprintf('<=%1.1e*SD => converged.',inInfo.dissection.bimonotonicRegression.epsilon4convergenceInSDs), sprintf('> %1.1e*SD => not converged, yet',inInfo.dissection.bimonotonicRegression.epsilon4convergenceInSDs)) ...
          ,sd4epsilon...
        );
        %Dev plot:
          if(inInfo.internal.bDevPlots)
            sTitle = sprintf('bimonotonic regression iteration %d', nBimonotonicConvergenceLoopIterations);
            fDevPlot = figure(...
              'Name', sTitle ...
             ,'Position',[73 1 1841 1021]...
             ,'Color','w' ...
            );
              colormap(sqrt(inInfo.plots.colormap_withNaNColor));
              Hs = [];
              cbHs = [];
            %Original, current, next, delta:
              if(true)
                Hs(end+1) = subplot(2,6,1); h1=Hs(end);
                  data2Plot = X; %original
                    cLim = 1.5*[-1,1]*sqrt(nanvar(X(:),abs(X(:))));
                    caxis(cLim);
                    if(true)
                      lowerBound4Data = cLim(1) + 1.5*abs(diff(cLim))/size(inInfo.plots.colormap_withNaNColor,1);
                      data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                    end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('signal in eigenorder\nbefore regression'));
                %figure title:
                  h=text(...
                     ilsub(xlim(h1),1)-abs(diff(xlim(h1)))/20, min(ylim(h1))+abs(diff(ylim(h1)))*0.61 ...
                    ,[{sTitle};{' ';' ';' '}]...
                    ,'Parent',h1 ...
                    ,'Color',[0 0 0.5],'FontWeight','bold','FontSize',19,'Rotation',90 ...
                    ,'HorizontalAlignment','right','VerticalAlignment','bottom' ...
                  );

                Hs(end+1) = subplot(2,6,2);
                  data2Plot = bsxfun(@plus,ismember(1:nG,performanceSubI4G)',+ismember(1:nP,performanceSubJ4P));
                  imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis([0.1,2.1]);
                    %cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('performance subspace mask\n(incomplete rows/cols are interpolated)'));
%                 subplot(2,6,2);
%                   data2Plot = bsxfun(@plus,ismember(1:nG,performanceSubI4G)',+ismember(1:nP,performanceSubJ4P));
%                   imagesc(bsxfun(@or, BInterpolationSpace4G, BInterpolationSpace4P)); %cannot use robustBaseline here as it is cluastering-dependent.
%                     caxis([0,1]);
%                     %cbHs(end+1) = colorbar('location','SouthOutside');
%                   title(sprintf('interpolation mask'));

                Hs(end+1) = subplot(2,6,3);
                  data2Plot = bimonotonicX_previous; %current/before
                    caxis(cLim);
                    if(true)
                      lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                      data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                    end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('start of bimonotonic\nregression iteration %d',nBimonotonicConvergenceLoopIterations));

                Hs(end+1) = subplot(2,6,4);
                  data2Plot = Z_eqd_conv_enforcedMonotonicity4G-Z_eqd_conv_enforcedMonotonicity4P; %differences
                    cLim2 = 3*[-1,1]*nanstd(data2Plot(:));
                    if(all(~isnan(cLim2)))
                      caxis(cLim2);
                    end
                    if(true)
                      lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                      data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                    end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim2);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('differences between\nmonotonic genes\nand monotonic samples'));
                  
%                 Hs(end+1) = subplot(2,6,4);
%                   data2Plot = fnc_ip_convolution_eqd({orderMetric4G, orderMetric4P}); %naive Fourier smoothing+eqdBimonotonicity based.
%                     caxis(cLim);
%                     if(true)
%                       lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
%                       data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
%                     end
%                     imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
%                     caxis(cLim);
%                     cbHs(end+1) = colorbar('location','SouthOutside');
%                   title(sprintf('interpolation base'));


                Hs(end+1) = subplot(2,6,5);
                  data2Plot = bimonotonicX; %common and next
                    caxis(cLim);
                    if(true)
                      lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                      data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                    end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('result of bimonotonic regression\niteration %d (weighted mean ~\nnewW4P*regressionResults4P\n+newW4G*regressionResults4G)',nBimonotonicConvergenceLoopIterations));
                  
                Hs(end+1) = subplot(2,6,6);
                  data2Plot = -bimonotonicX_previous+bimonotonicX; %Veränderung
                    cLim3 = 3*[-1,1]*nanstd(data2Plot(:));
                    caxis(cLim3);
                    if(true)
                      lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                      data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                    end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim3);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('result minus before (average\nabs(signal) delta = %1.1e)%s',dL2RsConvergence,iif(bBiMonotonic,sprintf('\nCONVERGED'),'')));
              end
            %Weights before and after regression and regression results:
              if(true)
                Hs(end+1) = subplot(2,6,7);
                  data2Plot = log10(W4P_previous); %weights used for regression of genes.
                    %cLim4Weights = [0,max(3*nanstd(W4P_previous(:)),3*nanstd(W4G_previous(:)))];
                    %cLim4Weights = [min([W4P_previous(:);W4G_previous(:)]),max([W4P_previous(:);W4G_previous(:)])];
%                     cLim4Weights = [
%                        ( mean(log10([W4P_previous(W4P_previous>0);W4G_previous(W4G_previous>0);W4Z_eqd_conv_enforcedMonotonicity4G(W4Z_eqd_conv_enforcedMonotonicity4G>0);W4Z_eqd_conv_enforcedMonotonicity4P(W4Z_eqd_conv_enforcedMonotonicity4P>0)])) ...
%                        -3*nanstd(log10([W4P_previous(W4P_previous>0);W4G_previous(W4G_previous>0);W4Z_eqd_conv_enforcedMonotonicity4G(W4Z_eqd_conv_enforcedMonotonicity4G>0);W4Z_eqd_conv_enforcedMonotonicity4P(W4Z_eqd_conv_enforcedMonotonicity4P>0)])) )...
%                       ,max(log10([W4P_previous(W4P_previous>0);W4G_previous(W4G_previous>0);W4Z_eqd_conv_enforcedMonotonicity4G(W4Z_eqd_conv_enforcedMonotonicity4G>0);W4Z_eqd_conv_enforcedMonotonicity4P(W4Z_eqd_conv_enforcedMonotonicity4P>0)])) ...
%                     ];
                    %cLim4Weights = [log10(1+min([W4G_previous(:);W4P_previous(:);W4G(:);W4P(:)])), log10(max([W4G_previous(:);W4P_previous(:);W4G(:);W4P(:)]))]
                    %cLim4Weights = [log10(mean([W4G_previous(:);W4P_previous(:);W4G(:);W4P(:)])/5), log10(max([W4G_previous(:);W4P_previous(:);W4G(:);W4P(:)]))];
                    cLim4Weights = [log10(mean([W4G_previous(:);W4P_previous(:);W4G(:);W4P(:)])/5), log10(max([W4G_previous(:);W4P_previous(:);W4G(:);W4P(:)]))];
                    %if(cLim4Weights(1)==0) cLim4Weights(1) = min([W4P_previous(W4P_previous>0);W4G_previous(W4G_previous>0)])/2; end
                    %cLim4Weights = log10(cLim4Weights);
                    caxis(cLim4Weights);
                      if(true)
                        lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                        data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                      end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim4Weights);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('log_{10}(W4P before regression)\n(used for regressing genes\nand newW4P*nP/(nG+nP)\nfor regressionResults4P)'));

                Hs(end+1) = subplot(2,6,8);
                  data2Plot = Z_eqd_conv_enforcedMonotonicity4G; %monotonic genes
                    caxis(cLim);
                    if(true)
                      lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                      data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                    end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('monotonic genes\n(regressionResults4G)'));

                Hs(end+1) = subplot(2,6,9);
                  data2Plot = log10(W4Z_eqd_conv_enforcedMonotonicity4G); %weights used for regression of genes.
                    caxis(cLim4Weights);
                      if(true)
                        lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                        data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                      end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim4Weights);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('log_{10}(weights from genes regression)\n(newW4G*nP/(nG+nP)\nfor regressionResults4G)')); %W4Z_eqd_conv_enforcedMonotonicity4G

                    
                    
                Hs(end+1) = subplot(2,6,10);
                  data2Plot = log10(W4G_previous); %before
                    caxis(cLim4Weights);
                      if(true)
                        lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                        data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                      end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim4Weights);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('log_{10}(W4G before regression)\n(used for regressing samples\nand newW4G*nG/(nG+nP)\nfor regressionResults4G)'));

                Hs(end+1) = subplot(2,6,11);
                  data2Plot = Z_eqd_conv_enforcedMonotonicity4P; %monotonic samples
                    caxis(cLim);
                    if(true)
                      lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                      data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                    end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('monotonic samples\n(regressionResults4P)'));

                Hs(end+1) = subplot(2,6,12);
                  data2Plot = log10(W4Z_eqd_conv_enforcedMonotonicity4P);
                    caxis(cLim4Weights);
                      if(true)
                        lowerBound4Data = ilsub(caxis,1) + 1.5*abs(diff(caxis))/size(inInfo.plots.colormap_withNaNColor,1);
                        data2Plot(data2Plot<lowerBound4Data) = lowerBound4Data; %show underflow expressions in the bluest color and not in NaN-gray.
                      end
                    imagesc(data2Plot); %cannot use robustBaseline here as it is cluastering-dependent.
                    caxis(cLim4Weights);
                    cbHs(end+1) = colorbar('location','SouthOutside');
                  title(sprintf('log_{10}(weights from samples regression)\n(newW4P*nG/(nG+nP)\nfor regressionResults4P)')); %W4Z_eqd_conv_enforcedMonotonicity4P
              end
            %Optics:
              set(cbHs,{'Position'},cellfun(@(pos)[pos(1)+pos(3)/5,pos(2)-2.5*pos(4),pos(3)*3/5,pos(4)*2/3],get(cbHs,'Position'),'UniformOutput',false),'FontSize',8);
            %Manual/Zoom to top genes/samples:
              zoom(fDevPlot,'on');
              if(false)
                set(Hs,'YLim',[1,signatureDefinition.step3_regression.signatureSizeByCorrSum4G])

                set(Hs,'YLim',[nG-signatureDefinition.step3_regression.signatureSizeByCorrSum4G,nG])
              end
            drawnow;
            keyboard;
            close(fDevPlot); %save RAM (do not keep multiple figures open containing multiple full matrix subplots...)
          end
    end
    %restore the numeric precision type of the input, if needed:
      if(~isa(bimonotonicX, numericTargetPrecision))
        fcnNumericTargetPrecision = iif(strcmp(numericTargetPrecision,'single'),@single,@double);
        bimonotonicX = fcnNumericTargetPrecision(bimonotonicX);
        assert(strcmp(class(bimonotonicX), numericTargetPrecision), 'code validation: numeric precision of the output does not match the input precision.');
      end
end

