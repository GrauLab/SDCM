%ABSTRACT
%  Monotonic regression using GPAV algorithm.
%SYNTAX
%  Yfit = GPAV(adj,Y,W, Blocks0)
% 
%  Input:
%   * adj adjacency matrix of observations which is sparse, upper triangular 
%     logical matrix, size [N x N]
%   * Y  response vector of observations, size [N x 1]
%   * W  vector of observation weights, size [N x 1]. By default W=ones(N,1)
%   * Blocks  initial blocks which are treated by GPAV. If
%             not specified, initial blocks are observations
%  Here N is the number of observations.
%
%  Output:
%   * Yfit fitted response which is nondecreasing 
%         in all explanatory variables, size [N x 1]
%   * BlocksF cell array containing the partitioning of indices 1..N into
%   blocks. Some cells are [0] which should be ignored.
%AUTHORS: 
%  code based on:
%    Burdakov O., Grimvall A., Sysoev O.
%    Reference: O. Burdakov, O. Sysoev, A. Grimvall and M. Hussian. 
%               "An O(n2) algorithm for isotonic regression."
%                In: G. Di Pillo and M. Roma (Eds) Large-Scale Nonlinear Optimization.
%                Series: Nonconvex Optimization and Its Applications,
%                Springer-Verlag, (2006) 83, pp. 25-33
%    download source: http://users.mai.liu.se/olebu87/SOFTWARE/GPAV/
%  adapted by 
%    Michael Grau, Dez 2013.
%    (added checks&assertions, added output parameters for weights, added 
%    missing values support)

function [YFit, WFit, Y4monotonicBlocks, W4monotonicBlocks, obs2block, Blocks] = GPAV(adj,Y,W,Blocks,bQuiet)
  %Parameter checks:
    %Check that adj is an upper triangular matrix (else the algorithm crashes):
      [I,J]=ind2sub(size(adj),find(adj));
      assert(all(J>=I), 'The observations should be enumerated so that adj(i,j)=1 implies i<j, i.e. adj is upper triangular.');
    assert(all(size(Y)==size(W)), 'size(Y)~=size(W)');
    assert(all(W>=0), 'weights must be positive');
  %Default params:
    if(nargin<4 || isempty(Blocks))
      N=size(adj,1);
      for q = 1:N
        Blocks{q} = q;
      end
      Y4monotonicBlocks=Y;
      W4monotonicBlocks=W;
      obs2block=1:N;
    else
      N=size(Blocks,2);
      obs2block=zeros(N,1);
      W4monotonicBlocks=zeros(N,1);
      Y4monotonicBlocks=zeros(N,1);
      for q = 1:N
        obs2block(Blocks{q}) = [q];
        W4monotonicBlocks(q)=sum(W(Blocks{q}));
        Y4monotonicBlocks(q)=sum(Y(Blocks{q}).*W(Blocks{q}))/W4monotonicBlocks(q);
      end
    end
    if(nargin<5) bQuiet=false; end
  %If the last value is NaN, replace it by the nearest non-NaN value (constant extrapolation):
    if(isnan(Y(end)))
      %warning('the last input value in regression order was NaN; it will now be imputed with the nearest non-NaN value')
      iLastNonNan = find(~isnan(Y),1,'last');
      if(~isempty(iLastNonNan))
        Y(end) = Y(iLastNonNan);
        W(end) = W(iLastNonNan);
      end
      %<-Note: due to the processing order, a NaN at the beginning of Y will automatically be imputed by the next-nearest non-NaN value, but from the end all NaN values will stay, if the last Nan is not set to the nearest non-NaN at the beginning.
    end

  %Main GPAV loop:
    for j=1:N
      if(Blocks{j}==0)
        warning('code validation: Blocks{j=%d} was already 0', j);
        continue;
      end
      while(true)
        %Get all predecessors as per the adjacency matrix:
          [rows4Predecessors, dummy] = find(adj(:,Blocks{j})); %die Vorgaenger (rowI) von den Punkten Blocks{j} bzgl. aller Kanten im Graph.
          rows4Predecessors = obs2block(rows4Predecessors); %uebersetze Punktindices zu den Indices ihrer Bloecke.
            rows4Predecessors(rows4Predecessors>=j) = []; %exclude j itelf and any later Nan Blocks (can be selected by the row below via the | isnan rule.
        %Get most inconsistent predecessor block:
          minRequiredAmonotonicity = 0;
          i4maxPredecessor = int32(0);
          for k=1:length(rows4Predecessors)
            iCandidate = rows4Predecessors(k);
            derivative = Y4monotonicBlocks(iCandidate) - Y4monotonicBlocks(j);
            bIsNanStep = isnan(derivative);
            if(bIsNanStep)
              if(i4maxPredecessor==0||iCandidate<i4maxPredecessor) %&& iCandidate<j) %NaN predecessors (from zero weight intervals) are inconsistent per definition and have to be regressed/inferred/overwritten with the value of neighboring blocks, but only after all nonNan predecessor blocks have been processed (considered to be more inconsistent with the monotony).
                i4maxPredecessor = rows4Predecessors(k);
              end
            else
              if(derivative > minRequiredAmonotonicity) %zeigt Blockindex iCandidate zum aktuellen Punkt inkonsistente Monotonie? und ist er inkonsistenter als ein anderer ggf. schon gefundener inkonsistenter Vorgaengerblock?
                minRequiredAmonotonicity = max(minRequiredAmonotonicity, derivative); %taking the highest of these is equivalent to selecting the absolutely highest predecessor, since Y4monotonicBlocks(j) is a constant.
                i4maxPredecessor = rows4Predecessors(k);
              end
            end
          end
          if(i4maxPredecessor==0) break; end %if there are no more inconsistent predecessor-blocks.
          i = i4maxPredecessor;
        %Merge values:
          %berechne Mittelwert des aktuellen Blocks mit dem inkonsistentesten Block eines Vorgaengers:
            %vA(j) = (xW(i)*vA(i) + xW(j)*vA(j))/(xW(j) + xW(i));
            %->Nan-robust version:
              new_vA_j_s1 = W4monotonicBlocks(i)*Y4monotonicBlocks(i);
              new_vA_j_s2 = W4monotonicBlocks(j)*Y4monotonicBlocks(j);
              new_vA_j_sum = 0;
              new_vA_j_norm = 0;
                if(~isnan(new_vA_j_s1)) 
                  new_vA_j_sum = new_vA_j_sum + new_vA_j_s1; 
                  new_vA_j_norm = new_vA_j_norm + W4monotonicBlocks(i); 
                end
                if(~isnan(new_vA_j_s2)) 
                  new_vA_j_sum = new_vA_j_sum + new_vA_j_s2; 
                  new_vA_j_norm = new_vA_j_norm + W4monotonicBlocks(j); 
                end
              Y4monotonicBlocks(j) = new_vA_j_sum/new_vA_j_norm;
              %<-Note: if both weights are zero, a NaN block will be formed, i.e. it will later be overwritten by the value of the nearest non-Nan block.
          %verschiebe das Gewicht des inkonsistenten Vorgaengerblocks in den aktuellen.
            %weights(j) = weights(j) + weights(i); 
            %->Nan-robust version:
              W4monotonicBlocks(j) = new_vA_j_norm;
        %Add block i into j
          Blocks{j} = [Blocks{i}, Blocks{j}]; %verschiebe die Punkteindices des inkonsistenten Vorgaengerblocks in den aktuellen.
          Blocks{i} = []; %leere den inkonsistenten Vorgaengerblock. Die Reihenfolge der Punkte und das Problem muessen so sein, dass adj upper triag ist; dann ist immer i<j und Blocks{i} kann in spaeteren j--Iterationen keine Probleme verursachen.
          obs2block(Blocks{j})=j; %ergaenze Punkteindices->Blockindices-Uebersetzungstabelle
        Y4monotonicBlocks(i)=0; %leere den Wert des inkonsistenten Vorgaengerblocks.
      end
      %status output:
        if(~bQuiet && mod(j,ceil(N/10))==0)
          SDCM_printStatus('%s%02d%%..',iif(j==ceil(N/10),sprintf('\nGPAV: '),''),round(j/N*100));
        end
    end
    if(~bQuiet)SDCM_printStatus('\n');end

  %Output:
    %Block values:
      YFit = Y4monotonicBlocks(obs2block); %=Y4monotonicBlocks(obs2block(1:N)), i.e. for every observation 1:N, get the Y of its block.
        if(any(isnan(YFit)))
          assert(any(isnan(Y)) || any(isnan(W)) || all(W==0), 'code validation: YFit contains NaNs although the input was NaN-free');
        end
    %Block weights:
      %WFit = W4monotonicBlocks(obs2block); %=W(obs2block(1:N)), i.e. for every observation 1:N, get the W of its block.
      %normalize weights per sample (instead of assigning all samples in a block the block weight):
        [uniqueBlockI,~,uniqueBlockII] = unique(obs2block);
        N4uniqueBlockI = histc(uniqueBlockII,1:length(uniqueBlockI));
        blockSizes = N4uniqueBlockI(uniqueBlockII); blockSizes = reshape(blockSizes,size(YFit));
          assert(all(blockSizes(obs2block)>0), 'code validation: block sizes of used blocks must be > 0');
        WFit = W4monotonicBlocks(obs2block)./blockSizes(obs2block); %=W(obs2block(1:N)), i.e. for every observation 1:N get the Y of its block.
          assert(all(WFit>=0), 'code validation: not all WFit>=0');
end

