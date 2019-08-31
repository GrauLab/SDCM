%ABSTRACT
%Helper function that returns true, if the code runs on the worker and false
%if executed on the client.

function result = isWorker()
  result = ~isempty( getCurrentWorker() );
end
