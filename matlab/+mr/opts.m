function out=opts(varargin)
%OPTS Set gradient limits of the MR system.
%   g=OPTS() Return the default amplitude and slew limits.
%
%   g=OPTS('MaxGradAmp',30,'AmpUnit','mT/m') Set the maximum gradient to
%   30mT/m.

persistent parser
validAmpUnits={'Hz/m','mT/m','rad/ms/mm'};
validSlewUnits={'Hz/m/s','mT/m/ms','T/m/s','rad/ms/mm/ms'};
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'opts';
    parser.addParamValue('ampUnit',validAmpUnits{1},...
        @(x) any(validatestring(x,validAmpUnits)));
    parser.addParamValue('slewUnit',validSlewUnits{1},...
        @(x) any(validatestring(x,validSlewUnits)));
    parser.addParamValue('maxGradAmp',[],@isnumeric);
    parser.addParamValue('maxGradSlew',[],@isnumeric);
end
parse(parser,varargin{:});
opt = parser.Results;

gamma = 42.576e6;      % Hz/T
if isempty(opt.maxGradAmp)
    maxAmp = 40*1e-3*gamma;     % Default: 40 mT/m
else
    
    switch opt.ampUnit
        case 'Hz/m'
            maxAmp = opt.maxGradAmp;
        case 'mT/m'
            maxAmp = opt.maxGradAmp*1e-3*gamma;
        case 'rad/ms/mm'
            maxAmp = opt.maxGradAmp*1e6/(2*pi);
    end
end
if isempty(opt.maxGradSlew)
    maxSlew=170*gamma;          % Default: 170 mT/m/ms
else
    switch opt.slewUnit
        case 'Hz/m/s'
            maxSlew = opt.maxGradSlew;
        case {'mT/m/ms','T/m/s'}
            maxSlew = opt.maxGradSlew*gamma;
        case 'rad/ms/mm/ms'
            maxSlew = opt.maxGradSlew*1e9/(2*pi);
    end
end
out.maxGradAmp = maxAmp;
out.maxGradSlew = maxSlew;

end