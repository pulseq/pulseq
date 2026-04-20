function sigPyOK = isSigPyAvailable()
%ISSIGPYAVAILABLE Test whether Python and SigPy are available
    [status, ~] = system('python3 -c "import sigpy" 2>/dev/null');
    if status ~= 0
        [status, ~] = system('python -c "import sigpy" 2>/dev/null');
    end
    sigPyOK = (status == 0);
end


