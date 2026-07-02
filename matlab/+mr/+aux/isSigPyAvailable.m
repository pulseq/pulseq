function [sigPyOK, pythonExe] = isSigPyAvailable()
%ISSIGPYAVAILABLE Test whether Python is available and has SigPy installed
%  If successfull stores a persistent state after the first call to accelerate 
%  subsequent checks; use clear GLOBALS to reset the persistent state.

persistent static_pythonExe;

if ~isempty (static_pythonExe)
    sigPyOK=true;
    pythonExe=static_pythonExe;
    return;
end

executables = {'python3', 'python', 'py'};
test_commands = {'which -a', 'where'};
exes = {};

if ispc
  null_out = 'nul';
  test_commands=test_commands(2:end); % there is no which on Windows anyway...
else
  null_out = '/dev/null';
end

for c=1:length(test_commands)
    for e=1:length(executables)
        cmd=sprintf('%s %s 2>%s',test_commands{c},executables{e},null_out);
        [status, output] = system(cmd);
        if status~=0, continue; end
        lines = regexp(output,'\n','split');
        for l=1:length(lines)
            exe=mr.aux.strstrip(lines{l});
            if isempty(exe), continue; end
            % test for Python with sigPy
            cmd=sprintf('"%s"  -c "import sigpy" 2>%s',exe,null_out);
            [status, output] = system(cmd);
            if status ~= 0, continue; end
            % found it!
            sigPyOK=true;
            pythonExe=exe;
            static_pythonExe=exe;
            return;
        end
    end
end

% did't find python or functional sigPy
sigPyOK = false;
pythonExe = '';

end
