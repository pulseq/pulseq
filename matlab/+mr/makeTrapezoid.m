function grad=makeTrapezoid(channel,varargin)
%makeTrapezoid Create a trapezoid gradient event.
%   g=makeTrapezoid(channel, ...) Create trapezoid gradient on
%   the given channel.
%
%   g=makeTrapezoid(channel,lims,...) Create trapezoid with the specificed 
%   gradient limits (e.g. amplitude, slew).
%
%   g=makeTrapezoid(...,'Duration',d,'Area',a) Create a
%   trapezoid gradient with given duration (s) and total area (1/m)
%   including ramps.
%
%   g=makeTrapezoid(...,'FlatTime',d,'FlatArea',a) Create a
%   trapezoid gradient with given flat-top time and flat-top
%   area not including ramps.
%
%   g=makeTrapezoid(...,'Amplitude',a) Create a trapezoid gradient with
%   given amplitude (Hz/m).
%
%   See also  Sequence.addBlock  mr.opts

persistent parser

if isempty(parser)
    validChannels = {'x','y','z'};
    parser = inputParser;
    parser.FunctionName = 'makeTrapezoid';
    parser.addRequired('channel',...
        @(x) any(validatestring(x,validChannels)));
    parser.addOptional('system',mr.opts(),@isstruct);
    parser.addParamValue('duration',0,@(x)(isnumeric(x) && x>0));
    parser.addParamValue('area',0,@isnumeric);
    parser.addParamValue('flatTime',[],@isnumeric);
    parser.addParamValue('flatArea',[],@isnumeric);
    parser.addParamValue('amplitude',[],@isnumeric);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    parser.addParamValue('riseTime',0,@isnumeric);
    
end
parse(parser,channel,varargin{:});
opt = parser.Results;

maxSlew=opt.system.maxSlew;
riseTime=opt.system.riseTime;
maxGrad=opt.system.maxGrad;

if opt.maxGrad>0
    maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    maxSlew=opt.maxSlew;
end
if opt.riseTime>0
    riseTime=opt.riseTime;
end

if isempty(opt.area) && isempty(opt.flatArea) && isempty(opt.amplitude)
    error('makeTrapezoid:invalidArguments','Must supply either ''area'', ''flatArea'' or ''amplitude''');
end
if opt.flatTime>0
    if ~isempty(opt.amplitude)
        amplitude = opt.amplitude;
    else
        amplitude = opt.flatArea/opt.flatTime;
    end
    if isempty(riseTime)
        riseTime = abs(amplitude)/maxSlew;
        riseTime = ceil(riseTime/opt.system.gradRasterTime)*opt.system.gradRasterTime;
    end
    fallTime = riseTime;
    flatTime = opt.flatTime;
elseif opt.duration>0
    if ~isempty(opt.amplitude)
        amplitude = opt.amplitude;
    else
        if isempty(riseTime)
            dC = 1/abs(2*maxSlew) + 1/abs(2*maxSlew);
            possible = opt.duration^2 > 4*abs(opt.area)*dC;
            amplitude = ( opt.duration - sqrt(opt.duration^2 - 4*abs(opt.area)*dC) )/(2*dC);
        else
            amplitude = opt.area/(opt.duration-riseTime);
            possible = opt.duration>2*riseTime & abs(amplitude)<maxGrad;
        end    
        assert(possible,'Requested area is too large for this gradient.');
            
    end
    if isempty(riseTime)
        riseTime = ceil(amplitude/maxSlew/opt.system.gradRasterTime)*opt.system.gradRasterTime;
    end
    fallTime = riseTime;
    flatTime = opt.duration-riseTime-fallTime;
    if isempty(opt.amplitude)
        % Adjust amplitude (after rounding) to achieve given area
        amplitude = opt.area/(riseTime/2 + fallTime/2 + flatTime);
    end
else
    error('makeTrapezoid:invalidArguments','Must supply a duration');
end
assert(abs(amplitude)<=maxGrad,'makeTrapezoid:invalidAmplitude','Amplitude violation');

grad.type = 'trap';
grad.channel = opt.channel;
grad.amplitude = amplitude;
grad.riseTime = riseTime;
grad.flatTime = flatTime;
grad.fallTime = fallTime;
grad.area = amplitude*(flatTime + riseTime/2 + fallTime/2);
grad.flatArea = amplitude*flatTime;

end