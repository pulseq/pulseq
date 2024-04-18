function out=convert(in,varargin)
%OPTS Set gradient limits of the MR system.
%   out=convert(in,from,to) Convert the numerical data in, given in
%   specified units 'from' to units specified in 'to'.
%
%   Valid unit strings are:
%    'Hz/m' 'mT/m' 'rad/ms/mm' 'Hz/m/s' 'mT/m/ms' 'T/m/s' 'rad/ms/mm/ms'

persistent parser
validGradUnits={'Hz/m','mT/m','rad/ms/mm'};
validSlewUnits={'Hz/m/s','mT/m/ms','T/m/s','rad/ms/mm/ms'};
validUnits=cat(2,validGradUnits,validSlewUnits);
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'convert';
    parser.addRequired('in',@isnumeric);
    parser.addRequired('fromUnit',...
        @(x) any(validatestring(x,validUnits)));
    parser.addOptional('toUnit',[],...
        @(x) any(validatestring(x,validUnits)));
    parser.addParamValue('gamma',42.576e6,@isnumeric); % Hz/T
end
parse(parser,in,varargin{:});
opt = parser.Results;

% Set default output unit if not given
if isempty(opt.toUnit)
    if ismember(opt.fromUnit,validGradUnits)
        opt.toUnit = validGradUnits{1};
    elseif ismember(opt.fromUnit,validSlewUnits)
        opt.toUnit = validSlewUnits{1};
    end
end

% Convert to standard units
switch opt.fromUnit
    % Grad units
    case 'Hz/m'
        standard = in;
    case 'mT/m'
        standard = in*1e-3*opt.gamma;
    case 'rad/ms/mm'
        standard = in*1e6/(2*pi);
    % Slew units
    case 'Hz/m/s'
        standard = in;
    case {'mT/m/ms','T/m/s'}
        standard = in*opt.gamma;
    case 'rad/ms/mm/ms'
        standard = in*1e9/(2*pi);
end

% Convert from standard units
switch opt.toUnit
    % Grad units
    case 'Hz/m'
        out = standard;
    case 'mT/m'
        out = 1e3*standard/opt.gamma;
    case 'rad/ms/mm'
        out = standard*2*pi*1e-6;
    % Slew units
    case 'Hz/m/s'
        out = standard;
    case {'mT/m/ms','T/m/s'}
        out = standard/opt.gamma;
    case 'rad/ms/mm/ms'
        out = standard*2*pi*1e-9;
end
