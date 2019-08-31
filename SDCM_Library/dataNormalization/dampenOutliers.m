%ABSTRACT
% Library function for SDCM. Dampens outliers by restricting the maximum neighbour signal ratio at the sort borders:
  function X = dampenOutliers(X, maxAllowedNeighboursRatio, topInXOrderToCheck, dim, precalc_SI4numStab)
    %Initialize:
      bMatrixProcessing = size(X,1)>1 && size(X,2)>1;
      if(isempty(X)) %special case, nothing to do.
        return;
      end
      if(nargin<4)
        if(bMatrixProcessing)
          error('dampenOutliers was called with a matrix input argument and therefore the dampening dimension must be specified.');
        end
        if(size(X,1)<=1)
          dim=2;
        elseif(size(X,2)<=1)
          dim=1;
        end
      end
    %recursive parallel matrix processing:
      if(bMatrixProcessing)
        %Set NaNs to zero temporarily (sort them to the middle):
          BNaNs = isnan(X);
          X(BNaNs) = 0;
        %Let Matlab do the dorting in parallel:
          [~,SI] = sort(X,dim,'ascend');
        %For each row, dampen in parallel:
          switch(dim)
            case 1
              parfor j=1:size(X,2)
                X(:,j) = dampenOutliers(X(:,j), maxAllowedNeighboursRatio, topInXOrderToCheck, dim, SI(:,j));
              end
            case 2
              parfor i=1:size(X,1)
                X(i,:) = dampenOutliers(X(i,:), maxAllowedNeighboursRatio, topInXOrderToCheck, dim, SI(i,:));
              end
          end
        %Recover NaNs:
          X(BNaNs) = NaN;
        return;
      end
    %single vector processing:
      %Temporarily center in dampening dimension (to have nothing near zero at the sort border to avoid extreme ratios):
        base = nanmean(X,dim);
        X = -base + X;
      %Sort if not already precalculated:
        if(nargin<5 || isempty(precalc_SI4numStab))
          X_zeroNaNs = X; X_zeroNaNs(isnan(X_zeroNaNs)) = 0; %temporarily set NaNs to zero (put to the middle of the sort order)
          [~,SI4numStab] = sort(X_zeroNaNs);
        else
          SI4numStab = precalc_SI4numStab;
        end
      %Dampen upper sort end:
        bChanged = true;
        while(bChanged) bChanged = false;
          for topn=1:min(length(X)-1,topInXOrderToCheck)
            bOutlier = abs(X(SI4numStab(topn))) > abs(maxAllowedNeighboursRatio*abs(X(SI4numStab(topn+1))));
            if(bOutlier)
              %warning('numeric stability: reduced X(SI4numStab(%d)) from %f to %f as it was more than inInfo.anchorSearch.outlierRemoval.maxAllowedNeighboursRatio times the next-smaller neighbour in anchor sort order.', topn, X(SI4numStab(topn)), maxAllowedNeighboursRatio*X(SI4numStab(topn+1)));
              X(SI4numStab(topn)) = maxAllowedNeighboursRatio*X(SI4numStab(topn+1));
              bChanged = true;
            end
          end
        end
      %Dampen lower sort end:
        SI4numStab = flipud(SI4numStab(:));
        bChanged = true;
        while(bChanged) bChanged = false;
          for topn=1:min(length(X)-1,topInXOrderToCheck)
            bOutlier = abs(X(SI4numStab(topn))) > abs(maxAllowedNeighboursRatio*abs(X(SI4numStab(topn+1))));
            if(bOutlier)
              %warning('numeric stability: reduced X(SI4numStab(%d)) from %f to %f as it was more than inInfo.anchorSearch.outlierRemoval.maxAllowedNeighboursRatio times the next-smaller neighbour in anchor sort order.', topn, X(SI4numStab(topn)), maxAllowedNeighboursRatio*X(SI4numStab(topn+1)));
              X(SI4numStab(topn)) = maxAllowedNeighboursRatio*X(SI4numStab(topn+1));
              bChanged = true;
            end
          end
        end
      %Restore baseline:
        X = X + base;
  end
