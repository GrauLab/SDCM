% ABSTRACT
%  Takes a source ID vector and a source data table, where each row corresponds to 
%  the source ID at the same position.
%  Takes a target ID vector, which determines the target data table layout and content.
%  Fills in data from the source table, where IDs match.
% SYNTAX
%  targetData = joinAndFill(sourceIDs, sourceData, targetIDs, fcnAggregateSourceDataRows, fcnComparer)
% INPUT
%  sourceIDs: column vector of numeric or cellstr sourceIDs
%  targetIDs: column vector of numeric or cellstr targetIDs of the same type
%  sourceData: numeric or cell matrix of source data rows corresponding to source IDs
%  fcnAggregateSourceDataRows: Optional. In case of an n-to-1 match, this function is given
%         all matching sourceData rows and must return one aggregated targetData row.
%         Must also be able to process an empty sourceData submatrix and return a
%         fill value. Is also called with submatrices with only one row (i.e. can 
%         also be used to transform data).
%  fcnComparer: optional custom compare function of signature
%     [bEqual, bSourceRow_gte_TargetRow, bIsEmptyMatch] = cmp(sourceIDsRow, targetIDsRow)
%  or [bEqual, bSourceRow_gte_TargetRow, bIsEmptyMatch] = cmp(sourceIDsRow, targetIDsRow, sDomain)
%     where sDomain is in {'sourceVsTarget', 'targetVsTarget', 'sourceVsSource'} (only needed if
%     the source and target IODs have different number types.)
% OUTPUT
%  targetData: numeric or cell matrix of source data rows corresponding to target IDs
% EXAMPLE:
%  joinAndFill(1:26, cellfun(@char,num2cell('a'+(0:25)),'UniformOutput',false), [1,7])
% AUTHOR
%  (C) Michael Grau, March 2011.
function targetData = joinAndFill(sourceIDs, sourceData, targetIDs, fcnAggregateSourceDataRows, fcnComparer, customSortCols)
  %% Parameter logic:
    bTransposedSourceIDs = false;
    bTransposedSourceData = false;
    bTransposedTargetIDs = false;
    persistent timeOfLastCodeChangeWarning;
      if(isempty(timeOfLastCodeChangeWarning)) timeOfLastCodeChangeWarning=datenummx(clock)-1; end
    %old/sourceIDs and targetIDs must be column vectors; transpose dimensions, if needed:
      %<-backwards compatibility/still allow auto-transpose if we otherwise would only have one ID row and multiple ID cols (possible, but useless joinAndFill call):
        if(size(sourceIDs,1)==1 && size(sourceIDs,2)>size(sourceIDs,1) && size(sourceIDs,2)==size(sourceData,2))
          warning('deprecated: row vector inputs for sourceIDs and sourceData are still accepted and auto-transposed, but joinAndFill will require sourceIDs to be a nDataRecords*nIDColumns array and sourceData to be a nDataRecords*nDataColumns array in the future (auto-transposing is ambiguous for multi ID columns)');
          sourceIDs = sourceIDs';
            bTransposedSourceIDs = ~bTransposedSourceIDs;
          sourceData = sourceData'; %also transpose sourceData, assuming the caller input consistent rows for both; otherwise it will be transposed back below with an additional warning.
            bTransposedSourceData = ~bTransposedSourceData;
          %In this case, also allow depracated syntax of targetIDs as row vector:
            if(size(targetIDs,1)==size(sourceIDs,2))
              targetIDs = targetIDs';
                bTransposedTargetIDs = ~bTransposedTargetIDs;
            end
        end
if(-timeOfLastCodeChangeWarning+datenummx(clock) > 1/24/60/60*10) %no not spam the console with warnings that are only relevant for old code; just 10 times per minute.
      if(size(sourceIDs,2)>size(sourceIDs,1)) 
%2017-01/das ist ambiguous bei z.B. leren [0,1] cell strings als sourceIDs (oder auch bei multicol sourceIDs) => überlasse das dem caller!
%         sourceIDs = sourceIDs'; 
%         sourceData = sourceData';
%         bTransposedSourceInputs = true;
  warning('2017-01 change in joinAndFill: auto-transposing of non-column vectors as sourceIDs and sourceData switched off, as this was ambiguous to cases with empty sourceIDs or a few multicolumn sourceID rows. If this warning is displayed before a subsequent error in old code, check the orientation of joinAndFill inputs on caller level.')
      end
      if(size(targetIDs,2)>size(targetIDs,1)) 
%2017-01/das ist ambiguous bei z.B. leren [0,1] cell strings als sourceIDs (oder auch bei multicol sourceIDs) => überlasse das dem caller!
%           targetIDs=targetIDs'; 
%           bTransposedTargetIDs = true;
  warning('2017-01 change in joinAndFill: auto-transposing of non-column vectors as targetIDs switched off, as this was ambiguous to cases with empty sourceIDs or a few multicolumn sourceID rows. If this warning is displayed before a subsequent error in old code, check the orientation of joinAndFill inputs on caller level.')
      end
  timeOfLastCodeChangeWarning = datenummx(clock);
end
    %IDs must have the same number of columns, if no custom fcnComparer was specified:
      if(size(sourceIDs,2)~=size(targetIDs,2) && nargin<5)
        %Convenience syntax: Try to auto-transpose with warning to correct this:
          if(size(sourceIDs,1)==size(targetIDs,2))
            warning('ID columns count mismatch: size(sourceIDs,2) did not match size(targetIDs,2), but it matched after transposing sourceIDs');
            sourceIDs = sourceIDs';
              bTransposedSourceIDs = ~bTransposedSourceIDs;
            sourceData = sourceData'; %also transpose sourceData, assuming the caller input consistent rows for both; otherwise it will be transposed back below with an additional warning.
              bTransposedSourceData = ~bTransposedSourceData;
          elseif(size(sourceIDs,2)==size(targetIDs,1))
            warning('ID columns count mismatch: size(sourceIDs,2) did not match size(targetIDs,2), but it matched after transposing targetIDs');
            targetIDs = targetIDs';
              bTransposedTargetIDs = ~bTransposedTargetIDs;
          elseif(size(sourceIDs,1)==0 && iscell(sourceIDs)) %Unfortunately, in Matlab {'X'}(false) = {}, but {'X';'Y'}([false;false])=cell(0,1) => in the empty case, replace by empty array of matching column count:
            warning('ID columns count mismatch: size(sourceIDs,2) did not match size(targetIDs,2), but it matched after replacing the empty array for sourceIDs with an empty array of appropriate size. Often this is the case due to subindexing a scalar source cell array {''X''} via {''X''}(false) resulting in {}, instead of {''X''}(false,1) resulting in cell(0,1). Consider this fix on caller level.');
            sourceIDs = cell(0,size(targetIDs,2));
          else
            error('ID columns count mismatch: sourceIDs and targetIDs must be compatible, i.e. have the same number of columns, if no custom fcnComparer is specified');
          end
      end
      assert(size(sourceIDs,2)==size(targetIDs,2) || nargin>=5, 'ID columns count mismatch: size(sourceIDs,2) does not match size(targetIDs,2)');
    %source IDs and source data must have the same number of rows:
      if(size(sourceIDs,1)~=size(sourceData,1))
        %Convenience syntax: If the row counts of sourceIDs and sourceData do not match, transpose to matching with warning, if possible:
          if(    size(sourceIDs,1)==size(sourceData,2) && size(sourceData,1)==1)
            warning('source records count mismatch: size(sourceData,1) did not match size(sourceIDs,1), but it matched after transposing sourceData');
            sourceData = sourceData';
              bTransposedSourceData = ~bTransposedSourceData;
          elseif(size(sourceIDs,1)==size(sourceData,2) && size(sourceIDs,1) ==1)
            warning('source records count mismatch: size(sourceData,1) did not match size(sourceIDs,1), but it matched after transposing sourceIDs');
            sourceIDs = sourceIDs';
              bTransposedSourceIDs = ~bTransposedSourceIDs;
          else %cannot autocorrect ambiguous cases with muticolumn IDs.
            error('source records count mismatch: size(sourceData,1) does not match size(sourceIDs,1)');
          end
      end
      assert(size(sourceIDs,1)==size(sourceData,1),'source records count mismatch: size(sourceData,1) does not match size(sourceIDs,1)');
      
    %type checks:
      if(isnumeric(sourceIDs) && ~isnumeric(targetIDs)) error('if sourceIDs is numeric, so must be targetIDs'); end
      if(iscell(sourceIDs) && all(cellfun(@isnumeric,sourceIDs(:))) && ~(iscell(targetIDs) && all(cellfun(@isnumeric,targetIDs(:))))) 
        if(~(isempty(sourceIDs) && isempty(sourceData)))
          error('if sourceIDs isa cellnumeric, so must be targetIDs'); 
        end
      end
      %if(~iscellstr(sourceIDs) && ~isnumeric(sourceIDs)) error('sourceIDs must be either cellstr or numeric'); end;
      if(iscellstr(sourceIDs) && ~iscellstr(targetIDs)) error('if sourceIDs isa cellstr, so must be targetIDs'); end
  %% Handle special case of empty targetIDs list => return empty targetData of appropriate size and type:
    if(isempty(targetIDs))
      if(isnumeric(sourceData))
        targetData = nan(0,size(sourceData,2));
      elseif(islogical(sourceData))
        targetData = false(0,size(sourceData,2));
      elseif(iscell(sourceData))
        targetData = cell(0,size(sourceData,2));
      elseif(isstruct(sourceData))
        targetData = repmat(struct([]), size(sourceData,2));
      else
        error('sourceData is of unsupported type');
      end    
      return;
    end
    %assert(nargin>=5 || isempty(sourceIDs) || size(sourceIDs,2)==size(targetIDs,2),'size(sourceIDs,2) does not match size(targetIDs,2)');
  %% connect the IDs (calculate the full cross join):
    %support for multicolumn cellstring IDs:
      %Note: numeric multicol IDs are supported by IFromTo natively.
      if(iscellstr(sourceIDs) && size(sourceIDs,2)>1)
        warning('multicolumn cellstring IDs like row 1 {%s} will be string-concatenated before computing the join', cellstring2separatedList(sourceIDs(1,:)));
        concatSourceIDs = cell(size(sourceIDs,1),1);
          for ci=1:size(sourceIDs,1)
            concatSourceIDs{ci} = cellstring2separatedList(sourceIDs(ci,:),'|'); %needs a separator to be unambiguous.
          end
        sourceIDs = concatSourceIDs;
        concatTargetIDs = cell(size(targetIDs,1),1);
          for ci=1:size(targetIDs,1)
            concatTargetIDs{ci} = cellstring2separatedList(targetIDs(ci,:),'|'); %needs a separator to be unambiguous.
          end
        targetIDs = concatTargetIDs;
      end
%Performance: Vercellung nicht mehr erforderlich.     %keys are expected by IFromTo to be cell column vectors: 
%       if(isnumeric(sourceIDs)) 
%         sourceIDs=num2cell(sourceIDs);
%         targetIDs=num2cell(targetIDs);
%       end
    if(nargin<5) fcnComparer = []; end; %use standard comparer inside of IFromTo
    if(nargin<6) customSortCols = ':'; end; %use standard comparer inside of IFromTo
    %[~,~,~,~,~,~,CITarget2Source,CISource2Target] = IFromTo(...
    [~,~,~,~,~,~,CITarget2Source] = IFromTo(... %Performance: do not request CISource2Target, as it is not needed here, but may take long to compute in IFromTo.
       sourceIDs, targetIDs...
      ,'rawJoinArraysOnly'...
      ,false,[],[]... %bAllowEmptyMatch, sourceAliasesFcn, targetAliasesFcn 
      ,fcnComparer...
      ,customSortCols ...
    );
  %unique matches only:
    if(nargin<4 || isempty(fcnAggregateSourceDataRows))
      BAvailable = cellfun(@(I)length(I)==1,CITarget2Source);
      ISource = [CITarget2Source{BAvailable}];
      ITarget = 1:length(CITarget2Source); ITarget(~BAvailable) = [];
      %create targetData:
        if(isnumeric(sourceData) || islogical(sourceData))
          targetData = nan(length(targetIDs),size(sourceData,2));
            targetData(ITarget,:) = sourceData(ISource,:);
          if(islogical(sourceData))
            if(any(isnan(targetData(:))))
              warning('Although sourceData is of class logical, targetData will be of class double, as it contained NaNs for targetIDs with no (or no unique) matches.');
            else
              targetData = logical(targetData);
            end
          end
        elseif(iscell(sourceData))
          targetData = cell(length(targetIDs),size(sourceData,2));
            targetData(ITarget,:) = sourceData(ISource,:);      
        elseif(isstruct(sourceData))
          targetData = repmat(cell2struct(cell(length(fieldnames(sourceData)),1), fieldnames(sourceData)), length(targetIDs), size(sourceData,2));
            targetData(ITarget,:) = sourceData(ISource,:);      
        else
          error('sourceData is of unsupported type');
        end    
      %Warning in case of unaggregates ambigious data:
        BAmbigious = cellfun(@(I)length(I)>1,CITarget2Source);
        if(any(BAmbigious))
          warning('data contains ambigious n:1 matches, which are not used, since no aggregation function was specified');
        end
      if(bTransposedTargetIDs)
        targetData = targetData'; %return in input format of targetIDs.
      end
      return;
  %n:1 matches with aggregation:
    else 
%       %preallocate targetData:
%         if(isnumeric(sourceData))
%           targetData = nan(length(targetIDs),size(sourceData,2));
%         elseif(iscell(sourceData))
%           targetData = cell(length(targetIDs),size(sourceData,2));
%         elseif(isstruct(sourceData))
%           targetData = repmat(cell2struct(cell(length(fieldnames(sourceData)),1), fieldnames(sourceData)), length(targetIDs), size(sourceData,2));
%         else
%           error('sourceData is of unsupported type');
%         end
      %aggregate source data and fill in:
        lastStatusOutputTime = datenummx(clock);
        for iTarget = 1:length(CITarget2Source)
          %status info:
            if((-lastStatusOutputTime+datenummx(clock))*24*60*60>7)
              disp([num2str(100*iTarget/length(CITarget2Source)),'% aggregated']);
              lastStatusOutputTime = datenummx(clock);
            end
          aggregatedRow = fcnAggregateSourceDataRows(sourceData(CITarget2Source{iTarget},:));
%             if(size(aggregatedRow,2)~=size(targetData,2))
%               warning(...
%                 [ 'The custom aggregation function for iTarget=%d did not result in a 1x%d row vector for the source data table of size %s.'...
%                  ,'\n=> the aggregation function changes the size of the output! (No future warning will follow, if all later elements are aggregated to compatible size.)'...
%                 ]...
%                ,iTarget,num2str(size(targetData,2)),mat2str(size(sourceData(CITarget2Source{iTarget},:)))...
%               );
%             end
          %infer targetData type from aggregatedRow type, if this is the first result:
%<-TODO: allow caller to specify the type, e.g. a cell array of jagged numeric arrays.
            if(iTarget==1)
              if(isnumeric(aggregatedRow))
                targetData = nan(length(targetIDs),size(aggregatedRow,2));
              elseif(islogical(aggregatedRow))
                targetData = false(length(targetIDs),size(aggregatedRow,2));
              elseif(iscell(aggregatedRow))
                targetData = cell(length(targetIDs),size(aggregatedRow,2));
              elseif(isstruct(aggregatedRow))
                targetData = repmat(cell2struct(cell(length(fieldnames(aggregatedRow)),1), fieldnames(aggregatedRow)), length(targetIDs), size(aggregatedRow,2));
              else
                error('aggregatedRow is of unsupported type; try wrapping with {} to produce a cell array as output.');
              end
            end
          targetData(iTarget,:) = aggregatedRow;
        end
        if(bTransposedTargetIDs)
          targetData = targetData'; %return in input format of targetIDs.
        end
        return;
    end
end
