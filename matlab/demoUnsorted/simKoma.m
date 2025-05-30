function out=simKoma(seq, varargin)
% simulates the sequence seq using KomaMRI.jl and Julia
% you need Julia be installed and available in the path
% also follow the official instructions for installing KomaMRI.jl
% additionally, we need the MAT package

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'simMR0';
     
    addParamValue(parser, 'workPath', '', @(x)isstring(x)||ischar(x)); % optional path to the temporary directory for files stored by the function
    addParamValue(parser, 'baseName', 'komasim', @(x)isstring(x)||ischar(x)); % optional file name for the sequence and simulation result files in the temporary directory
    addParamValue(parser, 'noCleanUp', false, @(x)islogical(x)||isnumeric(x)); % leave sequence and simulated signal files in place
end
parse(parser, varargin{:});
opt = parser.Results;

if isempty(opt.workPath)
    % on unix-like systems we prefer to use memory-mapped temporary
    % storage (does this work on a Mac?)
    if exist('/dev/shm','dir')
        workPath='/dev/shm/koma_mat';
    else
        % otherwise fallback to Matlab's temp directory, which on windows
        % defaults to 'C:\Users\<name>\AppData\Local\Temp\' and on Linux to /tmp/
        workPath=[tempdir filesep 'koma_mat'];
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

helper={
'using KomaMRI'; 
'using MAT';
'sys = Scanner()';
'obj = brain_phantom2D(; ss=1)'; % us=2
'sim_params = KomaMRICore.default_sim_params()';
['seq=read_seq("' strrep(seqpath,'\','\\') '")'];
'sim_params["return_type"] = "mat"';
'raw = simulate(obj, seq, sys; sim_params)';
['matwrite("' strrep(sigpath,'\','\\') '", Dict("raw" => raw))'];
};

% create the helper script
fid=fopen([workPath filesep 'helper.jl'],'w');
for i=1:length(helper)
    fprintf(fid, '%s\n',helper{i});
end
fclose(fid);

seq.write_v141(seqpath);
%seq.write(seqpath); % as of May 2025, koma.jl needs v141
if ispc()
    [status,cmdout] = system(['cd ' workPath ' && julia --threads auto helper.jl']);
else
    [status,cmdout] = system(['cd ' workPath '; LD_LIBRARY_PATH=; julia --threads auto helper.jl']);
end
data=load(sigpath);
out=data.raw;

if ~opt.noCleanUp
    delete(seqpath);
    delete(sigpath);
end
