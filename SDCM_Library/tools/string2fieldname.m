%ABSTRACT
% Library function for SDCM. Convert a general text string to 
% a compatible Matlab fieldname or file system string.

function fn = string2fieldname(s,bDontTruncate, bDontLowercaseFirstLetter, bOnlyDisallowFilesystemSpecialChars)
  if(nargin<2)bDontTruncate=false; end
  if(nargin<3)bDontLowercaseFirstLetter=false; end
  if(nargin<4)bOnlyDisallowFilesystemSpecialChars=false; end
  
  %Support cellstrings via recursion:
    if(iscellstr(s))
      fn = cellfun(@string2fieldname, s, 'UniformOutput', false);
      return;
    end
  if(~bOnlyDisallowFilesystemSpecialChars)
    disallowedSigns = [' ',10,13,9,'.:,;#''+-*^~={[]}()/\&%$§"!@<>|',26,181,'±½'];
  else
    disallowedSigns = ['/\:*?"<>|',10,13,9,26,181];
  end
  fn = '';
  bUpper = false;
  for i=1:length(s)
    %disp(s(i)); disp(double(s(i)));
    if(~ismember(s(i),disallowedSigns) && double(s(i))<=127) %exclude explicitly disallowed signs and all Unicode signs
      if(bUpper)
        fn = [fn,upper(s(i))];
      else
        fn = [fn,s(i)];
      end
      bUpper = false;
    else
      switch(s(i)) %optional replacementStrings for some disallowed signs:
        case '>'; fn = [fn,iif(i<length(s)&&s(i+1)=='=', iif(length(fn)>1&&fn(end)==upper(fn(end)),'gte','Gte'), iif(length(fn)>1&&fn(end)==upper(fn(end)),'gt','Gt'))];
        case '<'; fn = [fn,iif(i<length(s)&&s(i+1)=='=', iif(length(fn)>1&&fn(end)==upper(fn(end)),'lte','Lte'), iif(length(fn)>1&&fn(end)==upper(fn(end)),'lt','Lt'))];
        case '/'; fn = [fn,'Per'];
        case '+'; fn = [fn,'Plus'];
        case '-'; fn = [fn,'Minus'];
        case '('; fn = [fn,'_'];
        case '['; fn = [fn,'_'];
        case '{'; fn = [fn,'_'];
        case ';'; fn = [fn,'_'];
        case '%'; fn = [fn,'Prc'];
        case '#'; fn = [fn,'Hash'];
        case 181; fn = [fn,'mu'];
        case '±'; fn = [fn,'Pm'];
        case '½'; fn = [fn,'Half'];
        case '|'; fn = [fn,'_'];
      end
      bUpper = true;
    end
  end
  %lower first sign, if it is a letter and the second and third are lower case letters:
    if(~bDontLowercaseFirstLetter)
      if(length(fn)>=3 && fn(1)>='A' && fn(1)<='Z' && fn(2)>='a' && fn(2)<='z' && fn(3)>='a' && fn(3)<='z')
        fn(1) = lower(fn(1));
      end
    end
  %if the fieldname does not start with a letter, prepend a _:
    if(isempty(fn) || (~bOnlyDisallowFilesystemSpecialChars && ~(fn(1)>='a'&&fn(1)<='z' || fn(1)>='A'&&fn(1)<='Z')))
      warning(['fieldname "',fn,'" does not start with a letter; prefixing a "n"']);
      fn = ['n',fn];
    end
  %length limit:
    if(~bDontTruncate && length(fn)>63)
      warning(['fieldname "',fn,'" is too long (must have 63 characters at the maximum).']);
    end
end
