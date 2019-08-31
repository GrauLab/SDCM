%ABSTRACT:
%  Inline-replacement for left-assignment.
%SYNTAX
%  result = illa(leftVar,leftIndex,rightSide)
%EXAMPLES
%  illa(1:10,2:3,nan)
%  illa(rand(100),{1:100,1:10},nan)
%  illa(rand(100),{1:100,1:10},@(previousValues)f(previousValues))
%  illa({1,2,3},2,{7}); %note that illa subindexes via () into leftVar, not via {}, so rightSide must be a cell array, if leftVar is a cell array.
%  Library function for SDCM standalone deployment. Part of the Matlab inline 
%  language toolbox. (Not be redistributed separately.)
%AUTHOR
%  (C) Michael Grau, 2012
%  Library function for SDCM standalone deployment. Part of the Matlab inline 
%  language toolbox. (Not to be redistributed separately.)

function result = illa(leftVar,leftIndex,rightSide)
  if(isa(leftIndex,'function_handle'))
    leftIndex = leftIndex(leftVar);
  end
  if(isnumeric(leftIndex) || islogical(leftIndex))
    if(isa(rightSide,'function_handle'))
      leftVar(leftIndex) = rightSide(leftVar(leftIndex));
    else
      leftVar(leftIndex) = rightSide;
    end
  elseif(iscell(leftIndex))
    assert(all(cellfun(@(c)isnumeric(c)||islogical(c),leftIndex)), 'index ranges of all dimensions must be either numeric or logical');
    if(isa(rightSide,'function_handle'))
      leftVar(leftIndex{:}) = rightSide(leftVar(leftIndex{:}));
    else
      leftVar(leftIndex{:}) = rightSide;
    end
  elseif(ischar(leftIndex))
    if(isa(rightSide,'function_handle'))
      newValues = rightSide(eval(['leftVar(',leftIndex,')']));
      eval(['leftVar(',leftIndex,') = newValues;']);
    else
      eval(['leftVar(',leftIndex,') = rightSide;']);
    end
  else
    error('unknown leftIndex syntax');
  end
  result = leftVar;
end
