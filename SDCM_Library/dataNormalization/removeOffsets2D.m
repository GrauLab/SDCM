%ABSTRACT
% Library function for SDCM. remove global 2D offsets:
  function [L2RsWithoutOffsets, base4G, base4P] = removeOffsets2D(L2Rs,fcnMean,epsilon)
    nG = size(L2Rs,1); nP = size(L2Rs,2);
    L2RsWithoutOffsets = L2Rs;
    base4G = zeros(nG,1);
    base4P = zeros(1,nP);
    
    nIterations=0; 
    while(true) nIterations=nIterations+1;
      remainingBase4G = fcnMean(L2RsWithoutOffsets,2);
      remainingBase4P = fcnMean(L2RsWithoutOffsets,1);
      base4G = base4G + remainingBase4G*0.5; %not a full step to prevent oscillations; use exact half steps so that the result is symmetric around zero.
      base4P = base4P + remainingBase4P*0.5; %not a full step to prevent oscillations; use exact half steps so that the result is symmetric around zero.
      L2RsWithoutOffsets = bsxfun(@plus, -base4G, bsxfun(@plus, -base4P, L2Rs));
      remainingMaxOffset = max(max(abs(remainingBase4G)),max(abs(remainingBase4P)));
      if(remainingMaxOffset < epsilon)
        break;
      end
      if(nIterations>100)
        warning('removing offsets did not converge below epsilon=%0.1e after 100 iterations; now accepting AS IS', epsilon); %never happened; remove?
        break;
      end
    end
  end

