% ABSTRACT
%   Returns a MATLAB structure outStructure composed of all fields from
%   inStructure1 and inStructure2, where the values of inStructure1 dominate, 
%   if there are fields with the same name.
%   If a field is of type structure and existent in both input structs,
%   MergeStruct is used recursively.
% SYNTAX:
%   [outStructure, unequalInformation] = mergeStructs(inStructure1,inStructure2)
% INPUT
%   inStructure1, inStructure2
%     MATLAB structures, for example CDOs.
% OUTPUT
%   outStructure
%     Matlab structure containing all fields from inStructure1 and all
%     fields from inStructure2, which do not cause a field name conflict.
%   unequalInformation
%     unequalInformation{1}.(inStructure1Fields{i}) = inStructure1.(inStructure1Fields{i}); %dominating value
%     unequalInformation{2}.(inStructure1Fields{i}) = outStructure.(inStructure1Fields{i}); %different value.
% AUTHOR
%   Michael Grau (C)2005.
%   (Included for SDCM standalone deployment; not to be redistributed separately.)

function [outStructure, unequalInformation] = mergeStructs(inStructure1, inStructure2, internal_treePathToParent, bDebug)
  %Check parameters:
    if(nargin<2)
      error('HelperMergeStructs parameter error: parameters 1 and 2 are not optional.');
    end
    if(nargin<3) internal_treePathToParent = {'root'}; end
    if(nargin<4) bDebug = false; end;
    if(~isstruct(inStructure1))
      error('HelperMergeStructs parameter error: parameters ''inStructure1'' must be a MATLAB structure');
    end
    if(~isstruct(inStructure2))
      error('HelperMergeStructs parameter error: parameters ''inStructure2'' must be a MATLAB structure');
    end
  %recurse on structure arrays:
    if(isstruct(inStructure1) && isstruct(inStructure2) && ~isscalar(inStructure1))
      if(length(size(inStructure1))~=length(size(inStructure2)) || any(size(inStructure1)~=size(inStructure2)))
        %Allow merging with an empty structure array:
          if(isempty(inStructure1))
            outStructure = inStructure2;
            unequalInformation = {};
            return;
          end
          if(isempty(inStructure2))
            outStructure = inStructure1;
            unequalInformation = {};
            return;
          end
        error('cannot merge structure arrays of different size at %s', cellstring2separatedList(internal_treePathToParent,''));
      end
      %Performance: in case of disjunct field names, merge directly:
        fieldNames1 = fieldnames(inStructure1);
        fieldNames2 = fieldnames(inStructure2);
        if(isempty(intersect(fieldNames1,fieldNames2)))
          SC1 = struct2cell(inStructure1);
          SC2 = struct2cell(inStructure2);
          outStructure = cell2struct([SC1;SC2], {fieldNames1{:},fieldNames2{:}}, 1);
          unequalInformation = {};
          return;
        end
      unequalInformation = cell(size(inStructure1));
      for i=1:numel(inStructure1)
        if(bDebug) disp(['-> recursing into struct array index i=',num2str(i),' at ', cellstring2separatedList([internal_treePathToParent, {sprintf('(%d)',i)}],'')]); end;
        if(i==1)
          [outStructure, unequalInformation{i}] = mergeStructs(inStructure1(i),inStructure2(i), [internal_treePathToParent, {sprintf('(%d)',i)}], bDebug);
          outStructure = repmat(outStructure, size(inStructure1)); %preallocate memory
        else
          [outStructure(i), unequalInformation{i}] = mergeStructs(inStructure1(i),inStructure2(i), [internal_treePathToParent, {sprintf('(%d)',i)}], bDebug);
        end
      end
      %outStructure = reshape(outStructure, size(inStructure1));
      return;
    end
  %Merging structures:
    outStructure = inStructure2;
    unequalInformation = {struct(),struct()}; 
    inStructure1Fields = fieldnames(inStructure1);
    for i=1:length(inStructure1Fields)           
        if( (isstruct(inStructure1.(inStructure1Fields{i})) || isempty(inStructure1.(inStructure1Fields{i}))) ... %Recursion, if structs in both or one isempty:
         && isfield(inStructure2,(inStructure1Fields{i})) ...
         && (isstruct(inStructure2.(inStructure1Fields{i})) || isempty(inStructure2.(inStructure1Fields{i}))) ...
         && (~isempty(inStructure1.(inStructure1Fields{i})) || ~isempty(inStructure2.(inStructure1Fields{i}))) ... %field at least one structure not empty
        )
          %Replace [] with struct():
            if(isempty(inStructure1.(inStructure1Fields{i})))
              inStructure1.(inStructure1Fields{i}) = repmat(struct(),size(inStructure2.(inStructure1Fields{i})));
            end
            if(isempty(inStructure2.(inStructure1Fields{i})))
              inStructure2.(inStructure1Fields{i}) = repmat(struct(),size(inStructure1.(inStructure1Fields{i})));
            end
          if(bDebug) disp(['-> recursing into struct array in field .',inStructure1Fields{i},' at ', cellstring2separatedList([internal_treePathToParent, {sprintf('.%s',inStructure1Fields{i})}],'')]); end;
          if(nargout>=2)
            [outStructure.(inStructure1Fields{i}), subUnequalInfo] = mergeStructs(inStructure1.(inStructure1Fields{i}),inStructure2.(inStructure1Fields{i}), [internal_treePathToParent, {sprintf('.%s',inStructure1Fields{i})}], bDebug);
              if(~isempty(subUnequalInfo))
                unequalInformation{1}.(inStructure1Fields{i}) = subUnequalInfo{1};
                unequalInformation{2}.(inStructure1Fields{i}) = subUnequalInfo{2};
              end
          else
            outStructure.(inStructure1Fields{i}) = mergeStructs(inStructure1.(inStructure1Fields{i}),inStructure2.(inStructure1Fields{i}), [internal_treePathToParent, {sprintf('.%s',inStructure1Fields{i})}], bDebug);
          end
        else %no recursion needed:
          %overwrite mode:
            if(nargout>=2 && isfield(outStructure,inStructure1Fields{i}) ...
              && ( any(size(outStructure.(inStructure1Fields{i})) ~= size(inStructure1.(inStructure1Fields{i}))) ...
                || ischar(outStructure.(inStructure1Fields{i})) && ~ischar(inStructure1.(inStructure1Fields{i})) ...
                || ~ischar(outStructure.(inStructure1Fields{i})) && ischar(inStructure1.(inStructure1Fields{i})) ...
                || ischar(outStructure.(inStructure1Fields{i})) && ischar(inStructure1.(inStructure1Fields{i})) && ~strcmp(outStructure.(inStructure1Fields{i}), inStructure1.(inStructure1Fields{i})) ...
                || isnumeric(outStructure.(inStructure1Fields{i})) && ~isnumeric(inStructure1.(inStructure1Fields{i})) ...
                || ~isnumeric(outStructure.(inStructure1Fields{i})) && isnumeric(inStructure1.(inStructure1Fields{i})) ...
                || isnumeric(outStructure.(inStructure1Fields{i})) && isnumeric(inStructure1.(inStructure1Fields{i})) && any(reshape(outStructure.(inStructure1Fields{i})~=inStructure1.(inStructure1Fields{i}),numel(outStructure.(inStructure1Fields{i})),1)) && ~(all(isnan(reshape(outStructure.(inStructure1Fields{i}),numel(outStructure.(inStructure1Fields{i})),1)))&&all(isnan(reshape(inStructure1.(inStructure1Fields{i}),numel(outStructure.(inStructure1Fields{i})),1)))) ...
            ))
              unequalInformation{1}.(inStructure1Fields{i}) = inStructure1.(inStructure1Fields{i}); %dominating value
              unequalInformation{2}.(inStructure1Fields{i}) = outStructure.(inStructure1Fields{i}); %different value.
              % warning(['unequal values (at least in field ',inStructure1Fields{i},'): ',13 ...
              % ,'  inStructure1=',evalc('disp(inStructure1)') ... 
              % ,'  inStructure2=',evalc('disp(inStructure2)') ... 
              % ]);
            elseif(nargout>=2 && ~isfield(outStructure,inStructure1Fields{i})) %field only exists in inStructure1
              unequalInformation{1}.(inStructure1Fields{i}) = 'missing in inStructure2, using value from inStructure1';
            end
          outStructure.(inStructure1Fields{i}) = inStructure1.(inStructure1Fields{i});
        end
    end
    if(nargout>=2)
      for fn=setdiff(fieldnames(inStructure2)',fieldnames(inStructure1)'); 
        fn=fn{1};
        unequalInformation{2}.(fn) = 'missing in inStructure1, using value from inStructure2';
      end
    end
    if(isempty(fieldnames(unequalInformation{1})) && isempty(fieldnames(unequalInformation{2})))
      unequalInformation = [];
    end
end
