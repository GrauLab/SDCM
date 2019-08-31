%ABSTRACT
% Inline replacement for subindexing on expressions evaluating to arrays.
%SYNTAX
% output = ilsub(A,I,subIndexType,...)
% I: a vector of indices or an evaluable index string like '1:end'. Can also
%    be a cell array of index vectors for multidim indexing (no eval support here).
%    Can also be a function handle of the form @(A)calcIndices(A), which allows
%    calculating indices based on the input array itself on the fly.
% subIndexType: default='()' with exceptions:
%               - for convenience, if ischar(I) and I(1)=='.', then default='.'.
%               - for convenience, if isstruct(A) and iscell(I), then default='*'.
%               optional='[{}]', 'vertcat({})', 'horzcat({})': different concatenation options.
%               optional '{}' (can return multiple outputs, but the caller must request them)
%               optional 'members'; returns A(ismember(A,I)). (keeping the order of A)
%               optional 'notMembers'; returns A(~ismember(A,I)). (keeping the order of A)
%               optional 'split'; I must be a numeric vector of split lengths, i.e. sum(I)=numel(A).
%                        Returns a cell array with length(output{i})=I(i) (splitting via linear sequential indexing)
%               <-convenience syntax: {} -> '{}' and [] -> '()'.
%               optional '.' or '*': The call is rerouted to mapSubPath for field or tree path subscripts (It becomes 
%                                    mapSubPath.subIndexingPath and any arguments after subIndexType are passed on.)
%               optional 'incomplete': supports NaN indices and returns a cell array (with empty cells for NaN indices)
%                but A with A(I) replaced by assignInplace. Use this to for exmaple to
%                make inline deletions like B(isnan(B))=[] via ilsub(B,isnan(B),[],[])
%                or left-assignments like X(isnan(X))=1:sum(isnan(X)) via ilsub(X,isnan(X),[],1:sum(isnan(X)))
%                (only (), {} and member sub indexing types are available when replacing values inplace).
%EXAMPLES
% concat and subindex in one line
%   ilsub([A1,A2],I)
% initialize multiple fields at once without for loop:
%   [structArray.newField] = ilsub(num2cell(nan(size(structArray))),':','{}')
% construct multiple-outputs function handle
%   fcn = @()ilsub({nargout,nargout},':',{})
%   [a,b] = fcn();
%NEAT EXAMPLE:
% [DB.measurements.BC.cellLines.CNs.extNorm.probes.geneID] = ilsub(joinAndFill(...
%   {knownHumanGenes.geneSymbol}',{knownHumanGenes.geneID}'...
%  ,{DB.measurements.BC.cellLines.CNs.extNorm.probes.HGNC_ID}'...
%  ,@(multiIDs)iif(true,{nanmin([NaN,multiIDs{:}])},inlineExec(@disp, iif(length(multiIDs)>1,[{'taking min of multi ID candidates: '},multiIDs{:}],[]) ))...
%  ,@(sourceIDsRow, targetIDsRow) ilsub({...
%       strcmpi(sourceIDsRow{1}, targetIDsRow{1})... %bEqual: we need to accept case-insensitive matches.
%      ,ilsub([ilsub( ... %bSourceRow_gte_TargetRow: a string is greater than the other, if its first non-equal letter is greater than the other at the same position.
%         sourceIDsRow{1}(1:min(end,length(targetIDsRow{1}))) >= targetIDsRow{1}(1:min(length(sourceIDsRow{1}),end))...
%        ,sourceIDsRow{1}(1:min(end,length(targetIDsRow{1}))) ~= targetIDsRow{1}(1:min(length(sourceIDsRow{1}),end))...
%       ),true],1)... 
%      ,isempty(sourceIDsRow{1})&&isempty(targetIDsRow{1})... %bIsEmptyMatch
%    },':',{}) ... %a working inline anonymous multi-output function for the signature: [bEqual, bSourceRow_gte_TargetRow, bIsEmptyMatch] = cmp(sourceIDsRow, targetIDsRow).
% ),':',{});
%AUTHOR
% (C) Michael Grau, 2012
% Library function for SDCM standalone deployment. Part of the Matlab inline 
% language toolbox. (Not to be redistributed separately.)

function varargout = ilsub(A,I,subIndexType,varargin)
  if(nargin<3) 
    if(ischar(I) && ~isempty(I) && I(1)=='.') subIndexType='.';
    elseif(isstruct(A) && iscell(I)) subIndexType='*';
    else subIndexType='()';
    end
  end
  if(iscell(subIndexType)&&isempty(subIndexType)) subIndexType='{}'; end %convenience syntax {}.
  if(isnumeric(subIndexType)&&isempty(subIndexType)) subIndexType='()'; end %convenience syntax [].
  %Support for field and treepath indexing by rerouting to mapSubPath:
    if(ischar(subIndexType) && (strcmp(subIndexType,'.')||strcmp(subIndexType,'*')))
      varargout = cell(1,max(1,nargout));
      [varargout{:}] = mapSubPath(A,I,varargin{:});
      return;
    end
  assert(ischar(subIndexType), 'subIndexType must be [], {} or a string as defined in the help sheet');
  bIWasFcnHandle = false;
  if(isa(I,'function_handle')&&nargin(I)==1) %support for reflection (calculating the indices based on A on the fly)
    bIWasFcnHandle = true;
    try
      I = I(A);
    catch ex
      error('Subindex input argument I was provided as function handle, but its execution with A as input argument errored: %s', ex.message);
    end
  end
  if(nargin<4) %retrieve values:
    switch(lower(subIndexType))
      case '()'
        if(ischar(I))
          try
            if(strcmp(I,':')) %performance for simple flattening (no eval needed)
              output = A(:);
            else
              output = eval(['A(',I,')']);
            end
          catch ex
            if(~bIWasFcnHandle)
              error('Subindex input argument I was provided as an eval string like ''1:end'', but this index evaluation errored: %s', ex.message);
            else
              error('Subindex input argument I was provided as function handle. Its execution lead to the string ''%s'', which then was treated as an eval subindex string like ''1:end''. This subindex evaluation errored: %s. (Usually a function handle returns a subindex vector instead of a string that has to be evalled again; check your input for I.)', I, ex.message);
            end
          end
        elseif(iscell(I))
          %if there are any chars in I, replace them by numeric index vectors:
            KtoEval = find(cellfun(@ischar,I));
            for k=KtoEval(:)'
              II = 1:size(A,k);
              I{k} = eval(['II(',I{k},')']);
            end
          output = A(I{:});
        else
          output = A(I);
        end
        varargout = {output};
      case 'members'
        output = A(ismember(A,I));
        varargout = {output};
      case 'notmembers'
        output = A(~ismember(A,I));
        varargout = {output};
      case 'split'
        groupLengths = I;
        startIndices = cumsum([1,groupLengths]);
        endIndices = startIndices(2:end)-1;
          assert(endIndices(end)==numel(A), 'if you use the split method, sum(I) must equal numel(A)');
        startIndices(end) = [];
        output = cellfun(@(iStart,iEnd)A(iStart:iEnd),num2cell(startIndices),num2cell(endIndices),'UniformOutput',false);
        varargout = {output};
      case '[{}]'
        if(ischar(I))
          try
            output = eval(['[A{',I,'}]']);
          catch ex
            if(~bIWasFcnHandle)
              error('Subindex input argument I was provided as an eval string like ''1:end'', but this index evaluation errored: %s', ex.message);
            else
              error('Subindex input argument I was provided as function handle. Its execution lead to the string ''%s'', which then was treated as an eval subindex string like ''1:end''. This subindex evaluation errored: %s. (Usually a function handle returns a subindex vector instead of a string that has to be evalled again; check your input for I.)', I, ex.message);
            end
          end
        elseif(iscell(I))
          output = [A{I{:}}];
        else
          output = [A{I}];
        end
        varargout = {output};
      case 'vertcat({})'
        if(ischar(I))
          try
            output = eval(['vertcat(A{',I,'})']);
          catch ex
            if(~bIWasFcnHandle)
              error('Subindex input argument I was provided as an eval string like ''1:end'', but this index evaluation errored: %s', ex.message);
            else
              error('Subindex input argument I was provided as function handle. Its execution lead to the string ''%s'', which then was treated as an eval subindex string like ''1:end''. This subindex evaluation errored: %s. (Usually a function handle returns a subindex vector instead of a string that has to be evalled again; check your input for I.)', I, ex.message);
            end
          end
        elseif(iscell(I))
          output = vertcat(A{I{:}});
        else
          output = vertcat(A{I});
        end
        varargout = {output};
      case 'horzcat({})'
        if(ischar(I))
          try
            output = eval(['horzcat(A{',I,'})']);
          catch ex
            if(~bIWasFcnHandle)
              error('Subindex input argument I was provided as an eval string like ''1:end'', but this index evaluation errored: %s', ex.message);
            else
              error('Subindex input argument I was provided as function handle. Its execution lead to the string ''%s'', which then was treated as an eval subindex string like ''1:end''. This subindex evaluation errored: %s. (Usually a function handle returns a subindex vector instead of a string that has to be evalled again; check your input for I.)', I, ex.message);
            end
          end
        elseif(iscell(I))
          output = horzcat(A{I{:}});
        else
          output = horzcat(A{I});
        end
        varargout = {output};
      case '{}' %support for eg. [structArray.newField] = ilsub(num2cell(nan(size(structArray))),':','{}')
        assert(iscell(A),'cannot {}-subindex into non-cell array');
        if(ischar(I))
          try
            varargout = eval(['A(',I,')']);
          catch ex
            if(~bIWasFcnHandle)
              error('Subindex input argument I was provided as an eval string like ''1:end'', but this index evaluation errored: %s', ex.message);
            else
              error('Subindex input argument I was provided as function handle. Its execution lead to the string ''%s'', which then was treated as an eval subindex string like ''1:end''. This subindex evaluation errored: %s. (Usually a function handle returns a subindex vector instead of a string that has to be evalled again; check your input for I.)', I, ex.message);
            end
          end
        elseif(iscell(I))
          varargout = A(I{:});
        else
          varargout = A(I);
        end
        if(nargout~=length(varargout))
          warning(['ilsub is used with {} indexing returning ',num2str(length(varargout)),' outputs, but the caller requested ',num2str(nargout),' output(s)!']);
        end
      case 'incomplete'
        output = cell(size(I));
        BIsNan = isnan(I);
        if(iscell(A))
          output(~BIsNan) = A(I(~BIsNan));
        else
          [output{~BIsNan}] = ilC2A( num2cell(A(I(~BIsNan))) );
        end
        varargout = {output};
      otherwise
        error('unknown sub index type "%s"; see help.',subIndexType);
    end
  else %use inplaceReplace to replace values:
    error('ilsub with inplaceReplace is deprecated; use illa instead');
    inplaceReplace = varargin{1};
    output = A;
    switch(lower(subIndexType))
      case '()'
        if(ischar(I)) %eval index
          try
            eval(['output(',I,')=inplaceReplace;']);
          catch ex
            if(~bIWasFcnHandle)
              error('Subindex input argument I was provided as an eval string like ''1:end'', but this index evaluation errored: %s', ex.message);
            else
              error('Subindex input argument I was provided as function handle. Its execution lead to the string ''%s'', which then was treated as an eval subindex string like ''1:end''. This subindex evaluation errored: %s. (Usually a function handle returns a subindex vector instead of a string that has to be evalled again; check your input for I.)', I, ex.message);
            end
          end
        elseif(iscell(I)) %multidim index
          output(I{:}) = inplaceReplace;
        else %standard numeric index vector
          output(I) = inplaceReplace;
        end
        varargout = {output};
      case 'members'
        output(ismember(A,I)) = inplaceReplace;
        varargout = {output};
      case '{}' %support for eg. [structArray.newField] = ilsub(num2cell(nan(size(structArray))),':','{}')
        assert(iscell(A),'cannot {}-subindex into non-cell array');
        if(ischar(I))
          try
            varargout = eval(['A(',I,')']);
          catch ex
            if(~bIWasFcnHandle)
              error('Subindex input argument I was provided as an eval string like ''1:end'', but this index evaluation errored: %s', ex.message);
            else
              error('Subindex input argument I was provided as function handle. Its execution lead to the string ''%s'', which then was treated as an eval subindex string like ''1:end''. This subindex evaluation errored: %s. (Usually a function handle returns a subindex vector instead of a string that has to be evalled again; check your input for I.)', I, ex.message);
            end
          end
        elseif(iscell(I)) %multidim index
          output(I{:}) = inplaceReplace;
          varargout = output;
        else %standard numeric index vector
          output(I) = inplaceReplace;
          varargout = output;
        end
        if(nargout~=length(varargout))
          warning(['ilsub is used with {} indexing returning ',num2str(length(varargout)),' outputs, but the caller requested ',num2str(nargout),' output(s)!']);
        end
      otherwise
        error('unspported sub index type "%s" for inline inplace replace command; see help.',subIndexType);
    end
  end
end
