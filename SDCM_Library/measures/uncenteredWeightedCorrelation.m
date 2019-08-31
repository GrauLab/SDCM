%ABSTRACT
%  Library function for SDCM. 
%  - Computes the uncentered, possibly weighted correlations between num2cell(V1s,corrDim) and num2cell(V2s,corrDim) elements. 
%  - all vectors listed in batchDim=iif(corrDim==2,1,2) are processed separately
%  - R will be of size size(V1s,batchDim)*size(V2s,batchDim); DOFs and Ps likewise.
%  - parallelized for multiple vectors in V1s and V2s.
%  - Supports NaNs.
%NOTE
%  The usual Pearcon correlation centers the inputs via e.g. V1=-fcnMeanW(V1,abs(V1),1)+V1; We use uncentered correlations, 
%  since the inputs are already log2(ratio)s and the meaning of "zero==no regulation wrt. cohort average" should 
%  not be changed. This allows "high positive only" V1s and V2s to get a high correlation score, as intended, even 
%  if samples within the highly positive cluster are not correlated (i.e. not on a line in gene space). Centered correlation
%  would lose information here.

function [Rs, DOFs, Ps] = uncenteredWeightedCorrelation(V1s,V2s,corrDim,Ws, fcnMeanW,bJustDiagonal,BsPerformanceV2sSubspace)
  %%.Initialize and check parameters:
    %Get correlation dimension and check size consistency:
      if(nargin<3) %corrDim
        if(isrow(V1s) && isrow(V2s)) corrDim=2;
        elseif(iscolumn(V1s) && iscolumn(V2s)) corrDim=1;
        else
          error('corrDim must be specified for non-1D inputs');
        end
      else
        assert(isscalar(corrDim),'corrDim must be scalar');
      end
      if(nargin<6) bJustDiagonal=false; end
      batchDim = iif(corrDim==2,1,2);
      nElemsInV1s = size(V1s,batchDim);
      bSingleThreaded = nElemsInV1s==1; %do not invoke workers for a single V1 element (avoid unnecessary overhead).
      nElemsInV2s = size(V2s,batchDim);
      bSufficientSamplesInCorrDim = size(V1s,corrDim)>1 && size(V2s,corrDim)>1;
      assert(size(V1s,corrDim)==size(V2s,corrDim), 'V1s and V2s have an inconsistent number of elements in the correlation dimension %d', corrDim);
    %Special empty case (nothing to do):
      if(nElemsInV1s==0 || nElemsInV2s==0 || ~bSufficientSamplesInCorrDim)
        if(~bJustDiagonal)
          if(corrDim==2) %each gene gets correlated with the gene anchor
            if(~bSufficientSamplesInCorrDim)
              warning('V1s and V2s only have %d samples in the correlation dimension %d; this is not enough to compute correlations; now returning %d*%d correlations of zero with p values of one', size(V1s,corrDim), corrDim, nElemsInV2s,nElemsInV1s);
            end
            Rs = zeros(nElemsInV2s,nElemsInV1s);
            DOFs = ones(nElemsInV2s,nElemsInV1s);
            Ps = ones(nElemsInV2s,nElemsInV1s);
          elseif(corrDim==1) %each sample gets correlated with the sample anchor
            if(~bSufficientSamplesInCorrDim)
              warning('V1s and V2s only have %d samples in the correlation dimension %d; this is not enough to compute correlations; now returning %d*%d correlations of zero with p values of one', size(V1s,corrDim), corrDim, nElemsInV1s,nElemsInV2s);
            end
            Rs = zeros(nElemsInV1s,nElemsInV2s);
            DOFs = ones(nElemsInV1s,nElemsInV2s);
            Ps = ones(nElemsInV1s,nElemsInV2s);
            %currently Rs(k,:) is the kth candidate sample anchor V1s(:,k) correlated with all other sample cols in V2s, i.e. the sample correlations are already in the second dimension where they should be.
          end
        else %bJustDiagonal: do not compute all correlation pairs (i,j), but only compare (i,i) elements of V1s and V2s:
          if(corrDim==2) %each gene gets correlated with the gene anchor
            if(~bSufficientSamplesInCorrDim)
              warning('V1s and V2s only have %d samples in the correlation dimension %d; this is not enough to compute correlations; now returning %d*%d correlations of zero with p values of one', size(V1s,corrDim), corrDim, nElemsInV2s,nElemsInV1s);
            end
            Rs = zeros(nElemsInV2s,1);
            DOFs = ones(nElemsInV2s,1);
            Ps = ones(nElemsInV2s,1);
          elseif(corrDim==1) %each sample gets correlated with the sample anchor
            if(~bSufficientSamplesInCorrDim)
              warning('V1s and V2s only have %d samples in the correlation dimension %d; this is not enough to compute correlations; now returning %d*%d correlations of zero with p values of one', size(V1s,corrDim), corrDim, nElemsInV1s,nElemsInV2s);
            end
            Rs = zeros(1,nElemsInV2s);
            DOFs = ones(1,nElemsInV2s);
            Ps = ones(1,nElemsInV2s);
            %<-Rs(k,:) is the kth candidate sample anchor V1s(:,k) correlated with all other sample cols in V2s, i.e. the sample correlations are already in the second dimension where they should be.
          end
        end
        return;
      end
    %Weighting schemes:
      bW2D4V2sAsFunctionHandle = nargin>=4 && isa(Ws,'function_handle');
        if(bW2D4V2sAsFunctionHandle) assert(nargin(Ws)==2, 'if W2D is provided as function handle, it must be of format @(K1,K2) and return W2D with size(V1s,corrDim) elements in corrDim and size(V2s,batchDim) elements in batchDim'); end
      bWeightedCorrelation = nargin>=4 && ~isscalar(Ws) || bW2D4V2sAsFunctionHandle;
      if(~bWeightedCorrelation)
        if(Ws~=1) warning('scalar Ws=%f is interpreted as computing unweighted correlations', Ws); end
        Ws=1;
        bWeights4V1s = false;
        bW2D4V2s = false;
      end
      if(bWeightedCorrelation)
        if(bW2D4V2sAsFunctionHandle) %Check weights function handle for test case (1,:):
          WTest = Ws(1,1:nElemsInV2s);
          assert(all(size(WTest)==size(V2s)), 'provided weight function does not return a weights for the (1,:) test case in size(V2s)');
          assert(all(WTest(:)>=0), 'provided weight function does not return a weights for the (1,:) test case that are all positive'); %also checks than weights are not NaN which is assumed below.
          bW2D4V2s = false;
        else %Check 1D or 2D weights array for correct size:
          assert(size(Ws,corrDim)==size(V1s,corrDim), 'provided weight vector does not have the same number of elements in the correlation dimension %d as V1s', corrDim);
          assert(all(Ws(:)>=0), 'weights must be positive'); %also checks than weights are not NaN which is assumed below.
          bWeights4V1s = size(Ws,batchDim)==nElemsInV1s;
          bW2D4V2s = ndims(Ws)==3 && size(Ws,3)==nElemsInV2s;
          assert(size(Ws,batchDim)==1 || bWeights4V1s || bW2D4V2s, 'provided weights must be a single vector for usage with all elements in V1s or specify a weight vector for every element in V1s or additionally for every element in V2s in the third dimension, i.e. be of size(V1s) or [size(V1s),size(V2s,batchDim)]; in the 3D case i.e. as nG*nElems*nP 3D array for corrDim==1 and nElems*nP*nG 3D array for corrDim==2');
        end
      end
    if(nargin<5)
      fcnMeanW = @weightedMean_supportsNaNs;
    end
    if(nargin<6)
      bJustDiagonal = false; %Performance option if only V1s(i,:) with V2s(i,:) need to be correlated.
    end
      if(bJustDiagonal)
        assert(nElemsInV1s==nElemsInV2s, 'if bJustDiagonal, the number of provided elements in V1s and V2s in the batchDim=%d must be equal', batchDim);
      end

  %% Performance recursion feature: only correlate in a predetermined V2s subspace:
    if(nargin>=7 && ~isempty(BsPerformanceV2sSubspace) && ~(all(BsPerformanceV2sSubspace{1})&&all(BsPerformanceV2sSubspace{2})))
      assert(~bW2D4V2sAsFunctionHandle, 'bW2D4V2sAsFunctionHandle not implemented with usage of BsPerformanceV2sSubspace, yet');
      if(corrDim==2) %each gene gets correlated with the sample axis
        %Extract data within the performance subspace:
          %Performance: exclude cols with weight==0 for all elements in V1s from the computation as they will not change anything anyway:
            BsPerformanceV2sSubspace{2} = BsPerformanceV2sSubspace{2} | all(all(Ws==0, batchDim),3);

          V1s = V1s(:,BsPerformanceV2sSubspace{2});
          if(bWeightedCorrelation)
            if(bW2D4V2s)
              Ws = Ws(:,BsPerformanceV2sSubspace{2},BsPerformanceV2sSubspace{1});
            else
              Ws = Ws(:,BsPerformanceV2sSubspace{2});
            end
          end
          V2s = V2s(:,BsPerformanceV2sSubspace{2});

          assert(size(BsPerformanceV2sSubspace{1},1)==nElemsInV2s, 'code validation: size(BsPerformanceV2sSubspace{1},1)~=nElemsInV2s');
          V2s = V2s(BsPerformanceV2sSubspace{1},:);
          if(bJustDiagonal)
            V1s = V1s(BsPerformanceV2sSubspace{1},:);
            if(bWeights4V1s)
              Ws = Ws(BsPerformanceV2sSubspace{1},:);
            end
          end
        %Correlate in the performance subspace:
          if(nargout>=3)
            [subRs,subDOFs,subPs] = uncenteredWeightedCorrelation(V1s,V2s,corrDim,Ws, fcnMeanW,bJustDiagonal);              
          else
            [subRs,subDOFs] = uncenteredWeightedCorrelation(V1s,V2s,corrDim,Ws, fcnMeanW,bJustDiagonal);              
          end
        %Reshape outputs to the original space size and fill with zeros:
          if(~bJustDiagonal)
            Rs = zeros(nElemsInV2s,nElemsInV1s);
            DOFs = ones(nElemsInV2s,nElemsInV1s);
            Ps = ones(nElemsInV2s,nElemsInV1s);
            Rs(BsPerformanceV2sSubspace{1},:) = subRs;
            DOFs(BsPerformanceV2sSubspace{1},:) = subDOFs;
            if(nargout>=3)
              Ps(BsPerformanceV2sSubspace{1},:) = subPs;
            end
          else
            Rs = zeros(nElemsInV2s,1);
            DOFs = ones(nElemsInV2s,1);
            Ps = ones(nElemsInV2s,1);
            Rs(BsPerformanceV2sSubspace{1},1) = subRs; Rs=Rs(:,1);
            DOFs(BsPerformanceV2sSubspace{1},1) = subDOFs; DOFs=DOFs(:,1);
            if(nargout>=3)
              Ps(BsPerformanceV2sSubspace{1},1) = subPs; Ps=Ps(:,1);
            end
          end
      elseif(corrDim==1) %each sample gets correlated with the gene axis
        %Extract data within the performance subspace:
          %Performance: exclude rows with weight==0 for all elements in V1s from the computation as they will not change anything anyway:
            BsPerformanceV2sSubspace{1} = BsPerformanceV2sSubspace{1} | all(all(Ws==0, batchDim),3);
            
          V1s = V1s(BsPerformanceV2sSubspace{1},:);
          if(bWeightedCorrelation)
            if(bW2D4V2s)
              Ws = Ws(BsPerformanceV2sSubspace{1},:,BsPerformanceV2sSubspace{2});
            else
              Ws = Ws(BsPerformanceV2sSubspace{1},:);
            end
          end
          V2s = V2s(BsPerformanceV2sSubspace{1},:);

          assert(size(BsPerformanceV2sSubspace{2},2)==nElemsInV2s, 'code validation: size(BsPerformanceV2sSubspace{2},2)~=nElemsInV2s');
          V2s = V2s(:,BsPerformanceV2sSubspace{2});
          if(bJustDiagonal)
            V1s = V1s(:,BsPerformanceV2sSubspace{2});
            if(bWeights4V1s) %bugfix March2015 for diagonal corr of GEP and CNs for CHOP with weights for every gene.
              Ws = Ws(:,BsPerformanceV2sSubspace{2});
            end
          end
        %Correlate in the performance subspace:
          if(nargout>=3)
            [subRs,subDOFs,subPs] = uncenteredWeightedCorrelation(V1s,V2s,corrDim,Ws, fcnMeanW,bJustDiagonal);
          else
            [subRs,subDOFs] = uncenteredWeightedCorrelation(V1s,V2s,corrDim,Ws, fcnMeanW,bJustDiagonal);
          end
        %Reshape outputs to the original space size and fill with zeros:
          if(~bJustDiagonal)
            Rs = zeros(nElemsInV1s,nElemsInV2s);
            DOFs = ones(nElemsInV1s,nElemsInV2s);
            Ps = ones(nElemsInV1s,nElemsInV2s);
            Rs(:,BsPerformanceV2sSubspace{2}) = subRs;
            DOFs(:,BsPerformanceV2sSubspace{2}) = subDOFs;
            if(nargout>=3)
              Ps(:,BsPerformanceV2sSubspace{2}) = subPs;
            end
          else
            Rs = zeros(1,nElemsInV2s);
            DOFs = ones(1,nElemsInV2s);
            Ps = ones(1,nElemsInV2s);
            Rs(1,BsPerformanceV2sSubspace{2}) = subRs; Rs=Rs(1,:);
            DOFs(1,BsPerformanceV2sSubspace{2}) = subDOFs; DOFs=DOFs(1,:);
            if(nargout>=3)
              Ps(1,BsPerformanceV2sSubspace{2}) = subPs; Ps=Ps(1,:);
            end
          end
      end
      return;
    end

  %% Process all elements in V1s:
    BV2sIsNotNaN = ~isnan(V2s); %performance/precompute; needed below for NaN handling.
      if(all(all(BV2sIsNotNaN))) BV2sIsNotNaN=1; end %performance.
    bSingleThreaded = bSingleThreaded && ~bW2D4V2sAsFunctionHandle; 
    if(bSingleThreaded)
      assert(~bW2D4V2sAsFunctionHandle, 'bW2D4V2sAsFunctionHandle not implemented for bSingleThreaded, yet');
      %Unwrap weights for V2s:
        assert(size(Ws,batchDim)==1, 'size(Ws,%d) must be nElemsInV1s==1', batchDim);
        if(bW2D4V2s)
          if(corrDim==1)
            Ws = shiftdim(Ws(:,1,:),2)'; 
          elseif(corrDim==2)
            Ws = shiftdim(Ws(1,:,:),1)';
          end
          assert(all(size(Ws)==size(V2s)), '2D weights for V2s must be specified with the dual dimension as third dimension of Ws: i.e. as nG*1*nP 3D array for corrDim==1 and 1*nP*nG 3D array for corrDim==2');
        else
          assert(size(Ws,batchDim)==1 && length(size(Ws))<=2, 'weights for each element in V2s must be specified as nG*nElemsInV1s 2D array for corrDim==1 and nElemsInV1s*nP 2D array for corrDim==2');
        end          
      %Covariances:
        SQs = bsxfun(@times, V1s, V2s);
        Cov = fcnMeanW(SQs, Ws.^2, corrDim); %Note: both possible NaN sources (V1s and V2s) are passed to fcnMeanW here; therefore no special NaN handling is required here.
      %Estimate degrees of freedoms (number of elements behind the computed correlations; required for p value estimations):
        if(nargout>=2)
          DOFs = bsxfun(@times,Ws,double(~isnan(SQs))); %copy NaNs from V1 and V2s.
          DOFs = bsxfun(@rdivide,DOFs,max(DOFs,[],corrDim)+eps); %can only lose information by (relative) underweighting, but not gain by overweighting; +eps to 0/0->0 instead of NaN.
          DOFs = sum(DOFs, corrDim);
        end
      %V1 variance:
        %precomputed above/BV2sIsNotNaN = ~isnan(V2s); %we must compute the nanvar of V1 only for those elements of V2s(i,:) that are not NaN (otherwise we do not compute a correlation and r>1 may result)
        varV1 = fcnMeanW((V1s-0).^2, bsxfun(@times,Ws.^2,BV2sIsNotNaN), corrDim); %uncentered
      %V2s variances
        BV1IsNotNaN = ~isnan(V1s); %we must compute the nanvar of V2 only for those elements of V1 that are not NaN (otherwise we do not compute a correlation and r>1 may result)
          if(all(BV1IsNotNaN)) BV1IsNotNaN = 1; end %performance.
        varV2s = fcnMeanW((V2s-0).^2, bsxfun(@times,Ws.^2,BV1IsNotNaN), corrDim); %uncentered
      %Correlations ov V1 to V2s:
        R = Cov ./ sqrt(varV1.*varV2s);
        if(any(abs(R)>1.01))
          error('code validation: Correlations abs(R)>1 computed; if provided by the caller, check that the precomputed variances for V2s respect the actual passed weights!');
        end
        R(abs(R)>1) = sign(R(abs(R)>1)); %numeric protection; some single precision operations might yield slightly >1 values; truncate them here.
        R(isnan(R)|isinf(R)) = 0; %isinf(R) can happen for pruned/flat genes. isnan(R) can happen for full NaN rows/columns.
      %reshape/deliver Rs for the rows respectively cols that were correlated:
        if(corrDim==2) %each gene gets correlated with the gene anchor => nG*1 vector
          Rs = R(:); %can use (:) since there is only one element in V1s, so all r are between this element and the various elements of V2s
        elseif(corrDim==1) %each sample gets correlated with the sample anchor => 1*nP vector
          Rs = R(:)'; %can use (:) since there is only one element in V1s, so all r are between this element and the various elements of V2s
        end
    end
    if(~bSingleThreaded)
      %Initialize:
        pflc_nargout = nargout;
        %unweighted case: initialize variables for slicing (needed to prevend index out of bounds error from parfor side (cannot have Ws(i,:) in parfor even if disabled with bWeightedCorrelation))
          if(~bWeightedCorrelation) 
            if(corrDim==1)
              Ws = ones(1,nElemsInV1s); %use scalar weight 1 for each component of the vectors being correlated.
            end
            if(corrDim==2)
              Ws = ones(nElemsInV1s,1); %use scalar weight 1 for each component of the vectors being correlated.
            end
            bW2D4V2s = false;
          end 
        %performance: permute such that the V1s elements to process are in the first dimension for slicing: (else the V1s variable is not sliced when transported to the parallel workers; will be back-transposed inside the parfor loop)
          if(~bW2D4V2sAsFunctionHandle)
            dimsOfWs = length(size(Ws));
            %Consistency check:
              if(corrDim==1) %Ws is a nG*nElemsInV1s*nP or a nG*nElemsInV1s array
                if(dimsOfWs==2)
                  assert(all(size(Ws)==[size(V1s,1),nElemsInV1s])||all(size(Ws)==[1,nElemsInV1s])||all(size(Ws)==[size(V1s,1),1]), '2D weights for corrDim==1 must be of size [size(V1s,1),nElemsInV1s] or [1,nElemsInV1s] (scalar weights) or [size(V1s,1),nElemsInV1s] (common weights for each column that is to be correlated)');
                end
                if(dimsOfWs==3)
                  assert(all(size(Ws)==[size(V1s,1),nElemsInV1s,size(V2s,2)]), '3D weights for corrDim==1 must be of size [size(V1s,1),nElemsInV1s,size(V2s,2)]');
                end
              end 
              if(corrDim==2) %Ws is a nElemsInV1s*nP*nG or a nElemsInV1s*nP array
                if(dimsOfWs==2)
                  assert(all(size(Ws)==[nElemsInV1s,size(V1s,2)])||all(size(Ws)==[nElemsInV1s,1])||all(size(Ws)==[1,size(V1s,2)]), '2D weights for corrDim==2 must be of size [nElemsInV1s,size(V1s,2)] or [nElemsInV1s,1] (scalar weights) or [1,size(V1s,2)] (common weights for each row that is to be correlated)');
                end
                if(dimsOfWs==3)
                  assert(all(size(Ws)==[nElemsInV1s,size(V1s,2),size(V2s,1)]), '3D weights for corrDim==2 must be of size [nElemsInV1s,size(V1s,2),size(V2s,1)]');
                end
              end
            if(corrDim==1) %Ws is a nG*nElemsInV1s*nP or a nG*nElemsInV1s array
              Ws = num2cell(Ws,[1,3:dimsOfWs]); %package nElemsInV1s dim into cell array.
            end 
            if(corrDim==2) %Ws is a nElemsInV1s*nP*nG or a nElemsInV1s*nP array
              Ws = num2cell(Ws,[2,3:dimsOfWs]); %package nElemsInV1s dim into cell array.
            end
          end
          if(corrDim==1) V1s = V1s'; end
      if(~bJustDiagonal)
        Rs = nan(nElemsInV1s,nElemsInV2s);
        DOFs = nan(nElemsInV1s,nElemsInV2s);
        Ps = nan(nElemsInV1s,nElemsInV2s);
        %<-Note: these are the corrDim==1 dimensions; for corrDim==2 we also compute in these dimensions (for parfor slicing), but transpose after parfor.
      else
        Rs = nan(nElemsInV1s,1);
        DOFs = nan(nElemsInV1s,1);
        Ps = nan(nElemsInV1s,1);
        %<-Note: these are the corrDim==1 dimensions; for corrDim==2 we also compute in these dimensions (for parfor slicing), but transpose after parfor.
      end
      parfor i=1:nElemsInV1s
        %Get V1:
          V1 = V1s(i,:);
            if(corrDim==1) V1 = V1'; end %backtransposing; necessary for efficient slicing; see above.
        %Get weights:
          W = []; %prevent temporary variable compilation warings.
          if(~bW2D4V2sAsFunctionHandle)
            if(length(Ws)==1) %common weights for each vector being correlated
              W = Ws{1};
            elseif(length(Ws)==nElemsInV1s)
              W = Ws{i};
            else
              error('code validation: Ws not of length nElemsInV1s or one');
            end
          end
          if(bW2D4V2sAsFunctionHandle)
            W = Ws(i,1:nElemsInV2s);
            assert(all(size(W)==size(V2s)), '2D weights for V2s returned from the provided weights function handle must be of size(V2s)');
          elseif(bW2D4V2s)
            if(corrDim==1)
              W = shiftdim(W(:,1,:),2)';
            elseif(corrDim==2)
              W = shiftdim(W(1,:,:),1)';
            end
            assert(all(size(W)==size(V2s)), '2D weights for V2s must be specified with the dual dimension as third dimension of Ws: i.e. as nG*nElemsInV1s*nP 3D array for corrDim==1 and nElemsInV1s*nP*nG 3D array for corrDim==2');
          else
            assert(size(W,batchDim)==1 && length(size(W))<=2, 'weights for each element in V2s must be specified as nG*nElemsInV1s 2D array for corrDim==1 and nElemsInV1s*nP 2D array for corrDim==2');
          end
        %Covariances:
          SQs = []; %for the compiler as it cannot assert assingment before usage.
          V2s_local = []; %for the compiler as it cannot assert assingment before usage.
          if(bJustDiagonal)
            V2s_local = V2s; %Workaround: need a local copy, because indexing (i,:) or (:,i) depending on the corrDim confuses Matlabs parfor/slicing compiler.
          end
          if(~bJustDiagonal)
            SQs = bsxfun(@times, V1, V2s);
          else %Performance option, if only V1s(i,:) with V2s(i,:) need to be correlated:
            if(corrDim==2)
              %SQs = bsxfun(@times, V1, V2s(i,:));
              SQs = bsxfun(@times, V1, V2s_local(i,:));
            elseif(corrDim==1)
              %SQs = bsxfun(@times, V1, V2s(:,i));
              SQs = bsxfun(@times, V1, V2s_local(:,i));
            end
          end
          Cov = fcnMeanW(SQs, W.^2, corrDim); %Note: both possible NaN sources (V1s and V2s) are passed to fcnMeanW here; therefore no special NaN handling is required here.
        %Estimate degrees of freedoms (number of elements behind the computed correlations; required for p value estimations):
          if(pflc_nargout>=2)
            NsBehindRs = bsxfun(@times,W,double(~isnan(SQs))); %copy NaNs from V1 and V2s.
            NsBehindRs = bsxfun(@rdivide,NsBehindRs,max(NsBehindRs,[],corrDim)+eps); %can only lose information by (relative) underweighting, but not gain by overweighting; +eps to 0/0->0 instead of NaN.
            NsBehindRs = sum(NsBehindRs, corrDim);
            DOFs(i,:) = NsBehindRs;
          end
        %V1 variance:
          varV1 = []; %for the compiler as it cannot assert assingment before usage.
          if(~bJustDiagonal || (corrDim==2 && isrow(BV2sIsNotNaN)) || (corrDim==1 && iscolumn(BV2sIsNotNaN)) )
            varV1 = fcnMeanW((V1-0).^2, bsxfun(@times,W.^2,BV2sIsNotNaN), corrDim); %uncentered
          else
            if(corrDim==2)
              varV1 = fcnMeanW((V1-0).^2, bsxfun(@times,W.^2,BV2sIsNotNaN(i,:)), corrDim); %uncentered
            elseif(corrDim==1)
              varV1 = fcnMeanW((V1-0).^2, bsxfun(@times,W.^2,BV2sIsNotNaN(:,i)), corrDim); %uncentered
            end
          end
        %V2s variances:
          varV2s = []; %for the compiler as it cannot assert assingment before usage.
          BV1IsNotNaN = ~isnan(V1); %we must compute the nanvar of V2 only for those elements of V1 that are not NaN (otherwise we do not compute a correlation and r>1 may result)
            if(all(BV1IsNotNaN)) BV1IsNotNaN = 1; end %performance.
          if(~bJustDiagonal)
            varV2s = fcnMeanW((V2s-0).^2, bsxfun(@times,W.^2,BV1IsNotNaN), corrDim);
          else %Performance option if only V1s(i,:) with V2s(i,:) need to be correlated:
            if(corrDim==2)
              varV2s = fcnMeanW((V2s_local(i,:)-0).^2, bsxfun(@times,W.^2,BV1IsNotNaN), corrDim); %indirect indexing to prevent index out of bounds errors for wrong-dim auto runtime slicing (only for bJustDiagonal / not perf critical)
            elseif(corrDim==1)
              varV2s = fcnMeanW((V2s_local(:,i)-0).^2, bsxfun(@times,W.^2,BV1IsNotNaN), corrDim); %indirect indexing to prevent index out of bounds errors for wrong-dim auto runtime slicing (only for bJustDiagonal / not perf critical)
            end
          end
        %Correlations ov V1 to V2s:
          R = Cov ./ sqrt(varV1.*varV2s); %*fcnVariance(V2s,0,2));
          if(any(abs(R)>1.01))
            error('code validation: correlations abs(R(i=%d))>1 computed; if provided by the caller, check that the precomputed variances for V2s respect the actual passed weights!', i);
          end
          R(abs(R)>1) = sign(R(abs(R)>1)); %numeric protection; some single precision operations might yield slightly >1 values; truncate them here.
          R(isnan(R)|isinf(R)) = 0; %isinf(R) can happen for pruned/flat genes. isnan(R) can happen for full NaN rows/columns.
        %->output: this iteration computed the rs between the ith element of V1s and all elements of V2s.
          Rs(i,:) = R; 
      end
      %reshape/deliver Rs for the rows respectively cols that were correlated:
        if(corrDim==2) %each gene gets correlated with the gene anchor
          Rs = Rs'; %currently Rs(k,:) is the kth candidate gene anchor V1s(k,:) correlated with all other gene rows in V2s, but we want the gene correlations in the first dimension.
          DOFs = DOFs';
          %If only the diagonal was requested/computed:
            if(bJustDiagonal)
              %already a vector/Rs = diag(Rs); %return the correlations as a column vector, i.e. along the dimension of the vectors that were compared/correlated.
              Rs = Rs'; %return the correlations as a column vector, i.e. along the dimension of the vectors that were compared/correlated.
              DOFs = DOFs'; %since NsBehindRs is scalar in the bJustDiagonal case, DOFs is already a row vector; just transpose:
            end
        elseif(corrDim==1) %each sample gets correlated with the sample anchor
          %currently Rs(k,:) is the kth candidate sample anchor V1s(:,k) correlated with all other sample cols in V2s, i.e. the sample correlations are already in the second dimension where they should be.
          %If only the diagonal was requested/computed:
            if(bJustDiagonal)
              %already a vector/Rs = diag(Rs)'; %return the correlations as a row vector, i.e. along the dimension of the vectors that were compared/correlated.
              Rs = Rs'; %return the correlations as a row vector, i.e. along the dimension of the vectors that were compared/correlated.
              DOFs = DOFs'; %since NsBehindRs is scalar in the bJustDiagonal case, DOFs is already a row vector.
            end
        end
    end

  %% Compute p values for the correlations, if requested:
    if(nargout>=3)
      Ps = pValues4Correlations(Rs, DOFs, batchDim); %use a deterministic computation of p values (an approximation in case of weighted correlations that shows good correlation with p values determined by permutation tests)
    end
end
