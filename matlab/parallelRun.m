function results = parallelRun(cmds)
  results = true;
  numCPUs = str2num(getenv('NUMBER_OF_PROCESSORS'))/2; %hyperthreading is a lie
  if numCPUs < 1, numCPUs = 6; end
%  disp(sprintf('using %d CPUs',numCPUs));
  pids = [];
  maxI = numel(cmds);
  q=1;
  wb = waitbar(0,'0',"Name","Running child threads","createcancelbtn", "setappdata (gcbf,'interrupt', true)");
  for i=1:maxI
   if (numel(pids) - q) >= numCPUs % too many running processes, so wait
     waitpid(pids(q));
     q = q+1;
   endif
   if (! ishandle (wb))
     results = false; break;
   elseif (getappdata (wb, "interrupt"))
     close(wb);
     results = false; break;
   else
     waitbar(q/maxI , wb, sprintf('%d/%d',i,maxI));
   endif
   c = [regexprep(program_invocation_name(),'gui','cli') ' --no-gui --silent --eval "' regexprep(cmds{i},'\"','\\\"') '"'];
   disp(c)
   pids(i) = system(c,false,"async"); % system(cmds{7},false,"async")
  end
  if results
    while q <= numel(pids) % finish the queue
      waitpid(pids(q));
      if (! ishandle (wb))
        results = false; break;
      elseif (getappdata (wb, "interrupt"))
        close(wb);
        results = false; break;
      else
        waitbar(q/maxI, wb, sprintf('%d/%d',q,maxI ));
      end
      q = q+1;
    end
  end

  close(wb);
  end
