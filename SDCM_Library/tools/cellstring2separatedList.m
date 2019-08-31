%ABSTRACT:
% Create list string from cellstring.
%SYNTAX: 
% S = cellstring2separatedList(C,separator,truncateLength4EachEntry,bSkipEmpty,bNumbered)
%AUTHOR
% (C) Michael Grau, 2012
% Library function for SDCM standalone deployment.

  function S = cellstring2separatedList(C,separator,truncateLength4EachEntry,bSkipEmpty,bNumbered)
    if(nargin<2) separator=', '; end;
    if(nargin<3 || isempty(truncateLength4EachEntry)) truncateLength4EachEntry = Inf; end
    if(nargin<4 || isempty(bSkipEmpty)) bSkipEmpty = true; end;
    if(nargin<5 || isempty(bNumbered)) bNumbered = false; end;
      bNumbering_onlyIfMoreThanOne = false;
      if(~islogical(bNumbered))
        if(ischar(bNumbered) && strcmp(bNumbered,'onlyIfMoreThanOne'))
          bNumbered = true;
          bNumbering_onlyIfMoreThanOne = true;
        else
          error('unexpected input parameter bNumbered; see help.');
        end
      end

    S = '';
    for i=1:numel(C)
      sEntry = C{i};
      if(bSkipEmpty && isempty(sEntry)) continue; end;
      if(length(sEntry)>truncateLength4EachEntry+2)
        sEntry = [sEntry(1:truncateLength4EachEntry),'..'];
      end
      if(~bNumbered || (bNumbered && bNumbering_onlyIfMoreThanOne && numel(C)<=1))
        if(i<numel(C)) %no separator after last token
          if(iscellstr(separator)) %exact number and types of separators provided (see splitters output of strsplit)
            S = [S, sEntry, separator{i}];
          else %same char seperator everywhere
            S = [S, sEntry, separator];
          end
        else
            S = [S, sEntry];
        end
      else
        if(i<numel(C)) %no separator after last token
          if(iscellstr(separator)) %exact number and types of separators provided (see splitters output of strsplit)
            S = [S, num2str(i),') ', sEntry, separator{i}];
          else %same char seperator everywhere
            S = [S, num2str(i),') ', sEntry, separator];
          end
        else
            S = [S, num2str(i),') ', sEntry];
        end
      end
    end
  end
