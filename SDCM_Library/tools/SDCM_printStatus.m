%ABSTRACT
% Library function for SDCM. Tool for status/log messages.
% Works like fprinf, but accepts a status level as first input
% parameter. Only if the status level is <= the one initially
% OutputLevelForThisSession' special syntax are displayed.

function s = SDCM_printStatus(varargin)
  %Configure status level for this session:
    persistent dissectSignal_statusOutputLevel;
      if(isempty(dissectSignal_statusOutputLevel) || ~isnumeric(dissectSignal_statusOutputLevel) ||~isscalar(dissectSignal_statusOutputLevel))
        dissectSignal_statusOutputLevel = Inf;
      end
      if(length(varargin)>=1 && ischar(varargin{1}) && strcmp(varargin{1},'INTERNAL_configureStatusOutputLevelForThisSession'))
        if(length(varargin)>=2)
          dissectSignal_statusOutputLevel = varargin{2};
        end
        s = dissectSignal_statusOutputLevel;
        return;
      end
    persistent dissectSignal_statusOutputLevel4checkpoints;
      if(isempty(dissectSignal_statusOutputLevel4checkpoints) || ~isnumeric(dissectSignal_statusOutputLevel4checkpoints) ||~isscalar(dissectSignal_statusOutputLevel4checkpoints))
        dissectSignal_statusOutputLevel4checkpoints = -1;
      end
      if(length(varargin)>=1 && ischar(varargin{1}) && strcmp(varargin{1},'INTERNAL_configureCheckpointLevelForThisSession'))
        if(length(varargin)>=2)
          dissectSignal_statusOutputLevel4checkpoints = varargin{2};
        end
        s = dissectSignal_statusOutputLevel4checkpoints;
        return;
      end
  %Display message, if its status priority is sufficiently high:
    if(isempty(varargin)) 
      return; 
    end
    if(isnumeric(varargin{1}) && isscalar(varargin{1}))
      nStatusLevel = varargin{1};
      if(nStatusLevel <= dissectSignal_statusOutputLevel)
        fprintf(varargin{2:end});
      end
      if(nargout>0)
        s = sprintf(varargin{2:end}); %cannot use this to save the fprintf above as the string might contain escape chars like \ or % in the text.
      end
    else
      nStatusLevel = 3; %default status level, if none is provided as first argument.
      if(nStatusLevel <= dissectSignal_statusOutputLevel)
        fprintf(varargin{1:end});
      end
      if(nargout>0)
        s = sprintf(varargin{1:end}); %cannot use this to save the fprintf above as the string might contain escape chars like \ or % in the text.
      end
    end
  %update console immediately after high-level status messages.
    %if(nStatusLevel==0) drawnow; end
    if(nStatusLevel<=1) drawnow; end 
  %support for interactive check points: inject a break point in the caller of this message, if the checkpoint figure has been closed:
    if(nStatusLevel<=dissectSignal_statusOutputLevel4checkpoints)
      interactiveCheckPoint();
    end
end
