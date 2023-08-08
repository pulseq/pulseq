%% Compile mex files from MATLAB source
% This script is designed to compile frequently used Pulseq functions into
% mex files for faster execution.
%
% The Pulseq toolbox has been designed such that these files are *not*
% necessary to use the toolbox, but if they are present a substantial
% performance improvement can be acheived.
%
% The script relies on the MATLAB Coder package, available in later
% versions of MATLAB.
%

if isempty(which('codegen'))
    error('codegen not found for this version of MATLAB. Try 2014b or higher');
end
currDir=pwd;
if ~strcmp(currDir(end-2:end),'+mr')
    error('Please run from directory: pulseq/matlab/+mr')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create configuration object of class 'coder.MexCodeConfig'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.SaturateOnIntegerOverflow = false;
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;
cfg.ExtrinsicCalls = false;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define argument types for entry-point 'compressShape_mat'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ARGS = cell(1,1);
ARGS{1} = cell(1,1);
ARGS{1}{1} = coder.typeof(0,[Inf  1],[1 0]);

codegen -config cfg -o compressShape_mex ./compressShape_mat.m -args ARGS{1}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define argument types for entry-point 'find_mat'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{2} = struct;
ARGS{1}{2}.array = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[1 Inf],[0 1]);
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]);

codegen -config cfg -o ./@EventLibrary/find_mex ./@EventLibrary/find_mat.m -args ARGS{1}

