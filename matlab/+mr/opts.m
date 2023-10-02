function out=optsOld(varargin)
%OPTS Set gradient limits of the MR system.
%   g=OPTS() Return the default amplitude and slew limits.
%
%   g=OPTS('maxGrad',30,'gradUnit','mT/m') Set the maximum gradient to
%   30mT/m.

persistent defaultUserOpts
persistent defaultStandardOpts
if isempty(defaultStandardOpts)
    defaultStandardOpts=struct(...
        'maxGrad',mr.convert(40,'mT/m'),...   % Default: 40 mT/m
        'maxSlew',mr.convert(170,'T/m/s'),...	% Default: 170 mT/m/ms
        'riseTime',[],...
        'rfDeadTime',0,...
        'rfRingdownTime',0,...
        'adcDeadTime',0,...
        'adcRasterTime',100e-9,...
        'rfRasterTime',1e-6,...
        'gradRasterTime',10e-6,...
        'blockDurationRaster',10e-6,...
        'gamma',42576000,...
        'B0',1.5...
    );
end
if ~isempty(defaultUserOpts)
    defaultOpts=defaultUserOpts;
else
    defaultOpts=defaultStandardOpts;
end

persistent parser
validGradUnits={'Hz/m','mT/m','rad/ms/mm'};
validSlewUnits={'Hz/m/s','mT/m/ms','T/m/s','rad/ms/mm/ms'};
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'opts';
    parser.addParamValue('gradUnit',validGradUnits{1},...
        @(x) any(validatestring(x,validGradUnits)));
    parser.addParamValue('slewUnit',validSlewUnits{1},...
        @(x) any(validatestring(x,validSlewUnits)));
    parser.addParamValue('maxGrad',[],@isnumeric);
    parser.addParamValue('maxSlew',[],@isnumeric);
    parser.addParamValue('riseTime',[],@isnumeric);
    parser.addParamValue('rfDeadTime',defaultOpts.rfDeadTime,@isnumeric);
    parser.addParamValue('rfRingdownTime',defaultOpts.rfRingdownTime,@isnumeric);
    parser.addParamValue('adcDeadTime',defaultOpts.adcDeadTime,@isnumeric);
    parser.addParamValue('adcRasterTime',defaultOpts.adcRasterTime,@isnumeric);
    parser.addParamValue('rfRasterTime',defaultOpts.rfRasterTime,@isnumeric);
    parser.addParamValue('gradRasterTime',defaultOpts.gradRasterTime,@isnumeric);
    parser.addParamValue('blockDurationRaster',defaultOpts.blockDurationRaster,@isnumeric);
    parser.addParamValue('gamma',defaultOpts.gamma,@isnumeric); % Hz/T
    parser.addParamValue('B0',defaultOpts.B0,@isnumeric); % T
    parser.addParamValue('setAsDefault',false,@islogical);
end
parse(parser,varargin{:});
opt = parser.Results;

if isempty(opt.maxGrad)
    maxGrad = defaultOpts.maxGrad;
else
    maxGrad = mr.convert(opt.maxGrad,opt.gradUnit,'Hz/m','gamma',opt.gamma);
end
if isempty(opt.maxSlew)
    maxSlew=defaultOpts.maxSlew;
else
    maxSlew = mr.convert(opt.maxSlew,opt.slewUnit,'Hz/m','gamma',opt.gamma);
end
if ~isempty(opt.riseTime)
    %maxSlew=[];
    maxSlew=maxGrad/opt.riseTime;
end

out.maxGrad = maxGrad;
out.maxSlew = maxSlew;
out.riseTime = opt.riseTime;
out.rfDeadTime = opt.rfDeadTime;
out.rfRingdownTime = opt.rfRingdownTime;
out.adcDeadTime = opt.adcDeadTime;
out.adcRasterTime = opt.adcRasterTime;
out.rfRasterTime = opt.rfRasterTime;
out.gradRasterTime = opt.gradRasterTime;
out.blockDurationRaster = opt.blockDurationRaster;
out.gamma=opt.gamma;
out.B0=opt.B0;

if opt.setAsDefault
    defaultUserOpts=out;
    parser=[];
end

end