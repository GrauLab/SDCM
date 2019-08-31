%ABSTRACT
% Library function for SDCM. Weighted projections:

  function projectionsOnNormalizedAxesAsTargetSpaceVectors = projectW(vectorsInSourceSpace, axesInSourceSpace, weights4sourceDims, sourceSpaceMatrixDim, fcnMeanW, euclidW)
    if(nargin<5) fcnMeanW = @weightedMean_supportsNaNs; end
    if(nargin<6) euclidW = @(X,W,dim)fcnMeanW(X.^2,W.^2,dim,true).^(1/2); end %=sqrt(<x|x>_w)=sqrt(<w.x|w.x>)
    
    norm4axesInSourceSpace = euclidW(axesInSourceSpace, abs(weights4sourceDims), sourceSpaceMatrixDim); %axesInSourceSpace normalized within the weights4sourceDims-defined subspace. = sqrt(<axesInSourceSpace|axesInSourceSpace>_weights4sourceDims) = sqrt(<weights4sourceDims.*axesInSourceSpace|weights4sourceDims.*axesInSourceSpace>)
      norm4axesInSourceSpace(norm4axesInSourceSpace==0) = Inf; %make 0/0=NaN => 0/Inf=0, i.e. projection on zero vector is defined as zero. 
      
    projectionsOnNormalizedAxesAsTargetSpaceVectors = fcnMeanW(bsxfun(@times...
      ,vectorsInSourceSpace...
      ,bsxfun(@rdivide, axesInSourceSpace, norm4axesInSourceSpace)...
    ), abs(weights4sourceDims).^2, sourceSpaceMatrixDim, true); %<x|y>_w = <w.*x|w.*y>
  end
