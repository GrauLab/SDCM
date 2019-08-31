%ABSTRACT
% Library function for SDCM. Clear interactive checkpoints (for development only).
  function clearInjectedBreakPoints()
    dbs = dbstatus;
    for i=1:length(dbs)
      J = find(strcmp(dbs(i).expression, 'length(''injected break point by interactiveCheckPoint'')>0'));
      for j=J(:)'
        dbclear('in',dbs(i).file,num2str(dbs(i).line(j)));
      end
    end
  end

