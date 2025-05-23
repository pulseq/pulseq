function out=opts(varargin)
%OPTS Set gradient limits and other system properties of the MR system.
%   g=OPTS() Return the default amplitude and slew limits.
%
%   g=OPTS('maxGrad',30,'gradUnit','mT/m') Set the maximum gradient to
%   30mT/m.

persistent defaultUserOpts
persistent defaultStandardOpts
if isempty(defaultStandardOpts)
    defaultStandardOpts=struct(...
        'maxGrad',mr.convert(40,'mT/m'),...   % Default: 40 mT/m
        'maxSlew',mr.convert(170,'T/m/s'),... % Default: 170 mT/m/ms
        'maxB1',mr.convert(20,'uT'),...	      % Default: 20 uT
        'riseTime',[],...
        'rfDeadTime',0,...
        'rfRingdownTime',0,...
        'adcDeadTime',0,...
        'adcRasterTime',100e-9,...
        'rfRasterTime',1e-6,...
        'gradRasterTime',10e-6,...
        'blockDurationRaster',10e-6,...
        'adcSamplesLimit',0,... % 0 means no limit
        'rfSamplesLimit',0,... % 0 means no limit
        'adcSamplesDivisor',4,... % the number of which the adc.numSamples should be integer multiple 
        'gamma',42576000,...
        'B0',1.5...
    );
end

if ~isempty(defaultUserOpts)
    defaultOpts=defaultUserOpts;
else
    defaultOpts=defaultStandardOpts;
end

if isempty(varargin) % accelerate default constructor calls
    out=defaultOpts;
    return
end

persistent parser
validB1Units={'Hz','T','mT','uT'}; % todo: gauss?
validGradUnits={'Hz/m','mT/m','rad/ms/mm'};
validSlewUnits={'Hz/m/s','mT/m/ms','T/m/s','rad/ms/mm/ms'};
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'opts';
    parser.addParamValue('gradUnit',validGradUnits{1},...
        @(x) any(validatestring(x,validGradUnits)));
    parser.addParamValue('slewUnit',validSlewUnits{1},...
        @(x) any(validatestring(x,validSlewUnits)));
    parser.addParamValue('b1Unit',validB1Units{1},...
        @(x) any(validatestring(x,validSlewUnits)));
    parser.addParamValue('maxGrad',[],@isnumeric);
    parser.addParamValue('maxSlew',[],@isnumeric);
    parser.addParamValue('maxB1',[],@isnumeric);
    parser.addParamValue('riseTime',[],@isnumeric);
    parser.addParamValue('rfDeadTime',defaultOpts.rfDeadTime,@isnumeric);
    parser.addParamValue('rfRingdownTime',defaultOpts.rfRingdownTime,@isnumeric);
    parser.addParamValue('adcDeadTime',defaultOpts.adcDeadTime,@isnumeric);
    parser.addParamValue('adcRasterTime',defaultOpts.adcRasterTime,@isnumeric);
    parser.addParamValue('rfRasterTime',defaultOpts.rfRasterTime,@isnumeric);
    parser.addParamValue('gradRasterTime',defaultOpts.gradRasterTime,@isnumeric);
    parser.addParamValue('blockDurationRaster',defaultOpts.blockDurationRaster,@isnumeric);
    parser.addParamValue('adcSamplesLimit',defaultOpts.adcSamplesLimit,@isnumeric);
    parser.addParamValue('rfSamplesLimit',defaultOpts.rfSamplesLimit,@isnumeric);
    parser.addParamValue('adcSamplesDivisor',defaultOpts.adcSamplesDivisor,@isnumeric);
    parser.addParamValue('gamma',defaultOpts.gamma,@isnumeric); % Hz/T
    parser.addParamValue('B0',defaultOpts.B0,@isnumeric); % T
    parser.addParamValue('setAsDefault',false,@islogical);
end
parse(parser,varargin{:});
opt = parser.Results;

if isempty(opt.maxB1)
    maxB1 = defaultOpts.maxB1;
else
    maxB1 = mr.convert(opt.maxB1,opt.b1Unit,'Hz','gamma',opt.gamma);
end
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
out.maxB1 = maxB1;
out.riseTime = opt.riseTime;
out.rfDeadTime = opt.rfDeadTime;
out.rfRingdownTime = opt.rfRingdownTime;
out.adcDeadTime = opt.adcDeadTime;
out.adcRasterTime = opt.adcRasterTime;
out.rfRasterTime = opt.rfRasterTime;
out.gradRasterTime = opt.gradRasterTime;
out.blockDurationRaster = opt.blockDurationRaster;
out.adcSamplesLimit = opt.adcSamplesLimit;
out.rfSamplesLimit = opt.rfSamplesLimit;
out.adcSamplesDivisor = opt.adcSamplesDivisor;
out.gamma=opt.gamma;
out.B0=opt.B0;

if opt.setAsDefault
    defaultUserOpts=out;
    parser=[];
end

end