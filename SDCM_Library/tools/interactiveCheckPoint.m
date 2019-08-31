%ABSTRACT
% Library function for SDCM. Inject a breakpoint in the stack 
% if the checkpoint figure has been closed (for development only).

function interactiveCheckPoint(bInitialize)
  global hStopFigure;
  if(nargin<1) bInitialize = false; end
  if(isempty(hStopFigure)||~ishandle(hStopFigure)||~strcmp(get(hStopFigure,'Name'),'(close me to break at the next checkpoint)'))
    clearInjectedBreakPoints(); %cleanup.
    if(~bInitialize)
      dbs = dbstack;
        if(isempty(dbs))
          error('cannot inject a debugger stop point into the base workspace');
        end
        dbs(ismember({dbs.name},{'interactiveCheckPoint','statusOutput'})) = [];
      bBreakPointInserted = false;
        for stackIndex=1:length(dbs)
          newStopPoint = dbs(stackIndex);
          s = dbstatus;
          initialChecksum = length(horzcat(s.line));
          for lineOffset=1:10
            newStopPoint.line = newStopPoint.line+lineOffset;
            newStopPoint.anonymous = false;
            newStopPoint.expression = {'length(''injected break point by interactiveCheckPoint'')>0'};
            newStopPoint.cond = '';
            newStopPoint.identifier = {};
            try
              dbstop(newStopPoint)
            catch ex
              %warning('Could not inject a stop point in %s at line %d; details: %s; retrying further down in the code or up in the stack...', newStopPoint.name, newStopPoint.line, ex.message);
              continue;
            end
            s = dbstatus;
            newChecksum = length(horzcat(s.line));
            bBreakPointInserted = newChecksum > initialChecksum;
            if(bBreakPointInserted)
              warning(sprintf('INTERACTIVE BREAK in %s at line %d: Paused execution, since the checkpoint figure was closed.\n - execute >>clearInjectedBreakPoints() to not pause again at the current position\n - use >>dbcont or F5 to continue execution\n - set inInfo.internal.bDevEnableInteractiveBreaks=false to disable the pause feature\n - to break again, close the checkpoint figure window during execution', newStopPoint.name, newStopPoint.line));
              break;
            end
          end
          if(bBreakPointInserted)
            break;
          else
            warning('Could not inject a break point in %s within 10 lines of line %d; trying further up in the stack...', newStopPoint.name, newStopPoint.line);
          end
        end
        if(~bBreakPointInserted)
          warning('Could not inject a break; the checkpoint figure will be reopened; try breaking again later.');
          %keyboard;
        end
    end
    hStopFigure = figure('Name','(close me to break at the next checkpoint)','Position',[68,0,355,1],'HandleVisibility','off');
    drawnow;
  end
end
