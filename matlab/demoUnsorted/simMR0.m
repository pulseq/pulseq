function out=simMR0(seq,varargin)
% simulates the sequence seq using MR0 and Python
% to install MR0 run: pip3 install MRzeroCore
% on windows you may need to run: python -m pip MRzeroCore

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'simMR0';
     
    addParamValue(parser, 'pythonCmd', '', @(x)isstring(x)||ischar(x)); % optional Python executable path 
    addParamValue(parser, 'workPath', '', @(x)isstring(x)||ischar(x)); % optional path to the temporary directory for files stored by the function
    addParamValue(parser, 'baseName', 'mr0sim', @(x)isstring(x)||ischar(x)); % optional file name for the sequence and simulation result files in the temporary directory
    addParamValue(parser, 'noCleanUp', false, @(x)islogical(x)||isnumeric(x)); % leave sequence and simulated signal files in place
end
parse(parser, varargin{:});
opt = parser.Results;

if isempty(opt.workPath)
    % on unix-like systems we prefer to use memory-mapped temporary
    % storage (does this work on a Mac?)
    if exist('/dev/shm','dir')
        workPath='/dev/shm/mr0mat';
    else
        % otherwise fallback to Matlab's temp directory, which on windows
        % defaults to 'C:\Users\<name>\AppData\Local\Temp\' and on Linux to /tmp/
        workPath=[tempdir filesep 'mr0mat'];
    end
else
    workPath=opt.workPath;
end
seqpath=strrep([workPath filesep opt.baseName '.seq'], [filesep filesep], filesep);
sigpath=strrep([workPath filesep opt.baseName '.mat'], [filesep filesep], filesep);

% init the environment
if ~exist(workPath,'dir')
    mkdir(workPath);
end

if exist(seqpath,'file')
    delete(seqpath);
end
if exist(sigpath,'file')
    delete(sigpath);
end

% get phantom 
if ~exist([workPath filesep 'numerical_brain_cropped.mat'],'file')
    url='https://github.com/MRsources/MRzero-Core/raw/main/documentation/playground_mr0/numerical_brain_cropped.mat';
    % try wget (common in Lunix/Unix)
    [status,cmdout] = system(['cd ' workPath '&& wget ' url]);
    if 0~=status
        % try curl (usually included in Windows)
        [status,cmdout] = system(['cd ' workPath ' && curl -L -O ' url]); 
        if 0~=status
            warning("Failed to download numerical phantom, maybe your firewall is too cautious? \nThe simulation with MR0 will probably fail \nDownload the phantom from %s and store it under %s", url, workPath);

        end
    end
end

helper={
'import MRzeroCore as mr0'; 
'import scipy.io'
'dB0 = 0';
'sz = [128, 128]';
'# %% S4: SETUP SPIN SYSTEM/object on which we can run the MR sequence external.seq from above';
'# (i) load a phantom object from file';
'obj_p = mr0.VoxelGridPhantom.load_mat(''numerical_brain_cropped.mat'')';
'obj_p = obj_p.interpolate(sz[0], sz[1], 1)';
'# Manipulate loaded data';
'obj_p.T2dash[:] = 30e-3';
'obj_p.D *= 0';
'obj_p.B0 *= 1    # alter the B0 inhomogeneity';
'# Store PD and B0 for comparison';
'PD = obj_p.PD';
'B0 = obj_p.B0';
'# Manipulate loaded data';
'obj_p.B0+=dB0';
'obj_p.D*=0';
'#obj_p.plot()';
'# Convert Phantom into simulation data';
'obj_p=obj_p.build()';
'# %% SIMULATE  the external.seq file and add acquired signal to ADC plot';
'# Read in the sequence';
['seq0 = mr0.Sequence.import_file("' strrep(seqpath,'\','\\') '")'];
'#seq0.plot_kspace_trajectory()';
'# Simulate the sequence';
'graph=mr0.compute_graph(seq0, obj_p, 200, 1e-3)';
'signal=mr0.execute_graph(graph, seq0, obj_p, print_progress=False)';
['scipy.io.savemat("' strrep(sigpath,'\','\\') '", {''signal'':signal})']
};

% create the helper script
fid=fopen([workPath filesep 'helper.py'],'w');
for i=1:length(helper)
    fprintf(fid, '%s\n',helper{i});
end
fclose(fid);

% find/check python 
if ~isempty(opt.pythonCmd)
    [status, result]=system([opt.pythonCmd ' --version']);
    if status~=0
        error(['provided python executable ''' opt.pythonCmd ''' returns an error on the version check']);
    end
    python=opt.pythonCmd;
elseif ispc()
    % on Windows we rely on the PATH settings
    [status, result]=system('python --version');
    if status==0
        python='python';
    else
        [status, result]=system('py --version');
        if status~=0
            error('python executable not found, please check your system PATH settings');
        end
        python='py';
    end
else
    % this probably only works on linux and maybe also on mac
    [status, result]=system('which python3');
    if status==0
        python=mr.aux.strstrip(result);
    else
        [status, result]=system('which python');
        if status==0
            python=mr.aux.strstrip(result);
        else
            error('python executable not found');
        end
    end
end

seq.write_v141(seqpath);
%seq.write(seqpath); % as of May 2025, MR0 needs v141
[status,cmdout] = system(['cd ' workPath ' && ' python ' helper.py']);
data=load(sigpath);
out=data.signal;

if ~opt.noCleanUp
    delete(seqpath);
    delete(sigpath);
end
end
