%ABSTRACT
% Transforms a cell array into output arguments (as always the caller must request them).
% If C is numeric or logical, num2cell(C) will be processed.
%SYNTAX 
% varargout = ilC2A(C)
% varargout = ilC2A(C, nOutArgs)
%NOTE
% ilsub(cellArray,':',{}) does the same.
%AUTHOR
% (C) Michael Grau, 2012
% Library function for SDCM standalone deployment. Part of the Matlab inline 
% language toolbox. (Not to be redistributed separately.)

function varargout = ilC2A(C, nOutArgs)
  if(isnumeric(C)||islogical(C)) C= num2cell(C); end
  assert(iscell(C), 'C must be a cell array (or a numeric/logical array that will be converted via num2cell)');
  varargout = C;
  if(nargin>=2)
    if(nOutArgs>numel(C))
      error('ilC2A got a cell array of length %d, but was specified to return nOutArgs=%d output arguments', numel(C), nOutArgs);
    end
    C = C(1:nOutArgs);
  end
  if(nargout < numel(C))
    warning('ilC2A converted a cell array of %d elements to an equal amount of outputs, but the caller requested only %d outputs!', numel(C), nargout);
  elseif(nargout > numel(C))
    error('ilC2A converted a cell array of %d elements to an equal amount of outputs, but the caller requested %d, i.e. more outputs!', numel(C), nargout);
  end
end
