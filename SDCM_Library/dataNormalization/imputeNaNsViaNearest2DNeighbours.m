%ABSTRACT
% Library function for SDCM. Impute NaNs in a 2D matrix via nearest neighbors.

function X = imputeNaNsViaNearest2DNeighbours(X)
  BNans = isnan(X);
    if(all(all(BNans)))
      error('code validation: all(all(isnan(Z_orig))). Either the input data consists only of NaNs or there is an error.');
    end
  nNans=Inf;
  while(nNans>0)
    [I,J] = find(BNans);
    imputeSources = nan(length(I),8);
      leftTopLinI = sub2ind(size(X),    illa(I-1,@(K)K==0,NaN),               illa(J-1,@(K)K==0,NaN));
        B = ~isnan(leftTopLinI);  
        imputeSources(B,1) = X(leftTopLinI(B));
      topLinI = sub2ind(size(X),        I,                                    illa(J-1,@(K)K==0,NaN));
        B = ~isnan(topLinI);      
        imputeSources(B,2) = X(topLinI(B));
      topRightLinI = sub2ind(size(X),   illa(I+1,@(K)K==size(X,1)+1,NaN),illa(J-1,@(K)K==0,NaN));
        B = ~isnan(topRightLinI); 
        imputeSources(B,3) = X(topRightLinI(B));
      rightLinI = sub2ind(size(X),      illa(I+1,@(K)K==size(X,1)+1,NaN),J);
        B = ~isnan(rightLinI);    
        imputeSources(B,4) = X(rightLinI(B));
      bottomRightLinI = sub2ind(size(X),illa(I+1,@(K)K==size(X,1)+1,NaN),illa(J+1,@(K)K==size(X,2)+1,NaN));
        B = ~isnan(bottomRightLinI); 
        imputeSources(B,5) = X(bottomRightLinI(B));
      bottomLinI = sub2ind(size(X),     I,                                    illa(J+1,@(K)K==size(X,2)+1,NaN));
        B = ~isnan(bottomLinI); 
        imputeSources(B,6) = X(bottomLinI(B));
      leftBottomLinI = sub2ind(size(X), illa(I-1,@(K)K==0,NaN),               illa(J+1,@(K)K==size(X,2)+1,NaN));
        B = ~isnan(leftBottomLinI); 
        imputeSources(B,7) = X(leftBottomLinI(B));
      leftLinI = sub2ind(size(X),       illa(I-1,@(K)K==0,NaN),               J);
        B = ~isnan(leftLinI); 
        imputeSources(B,8) = X(leftLinI(B));
    %imputeValues = nanmean(imputeSources,2);
    imputeValues = nanmedian(imputeSources,2); %for pixels near signal plateau borders, let the majority of neighbours decide rather than the average.
      X(sub2ind(size(X),I,J)) = imputeValues;
    BNans = isnan(X);
      last_nNans=nNans;
      nNans = sum(sum(BNans));
      if(last_nNans==nNans)
        error('code validation: could not reduce the number of NaNs.');
      end
  end
end

