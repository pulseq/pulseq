function grad=makeTrapezoid(channel,varargin)
%makeTrapezoid Create a trapezoid gradient event.
%   g=makeTrapezoid(channel, ...) Create trapezoid gradient on
%   the given channel.
%
%   g=makeTrapezoid(...,'Duration',d,'Area',a) Create a
%   trapezoid gradient with given duration and total area
%   including ramps.
%
%   g=makeTrapezoid(...,'FlatTime',d,'FlatArea',a) Create a
%   trapezoid gradient with given flat-top time and flat-top
%   area not including ramps.
%
%   See also  Sequence.addBlock

persistent parser
validChannels = {'x','y','z'};
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeTrapezoid';
    parser.addRequired('channel',...
        @(x) any(validatestring(x,validChannels)));
    parser.addParamValue('duration',0,@(x)(isnumeric(x) && x>0));
    parser.addParamValue('area',0,@isnumeric);
    parser.addParamValue('flatTime',[],@isnumeric);
    parser.addParamValue('flatArea',[],@isnumeric);
    defaultOpts=mr.opts();
    parser.addParamValue('maxGradAmp',defaultOpts.maxGradAmp,@isnumeric);
    parser.addParamValue('maxGradSlew',defaultOpts.maxGradSlew,@isnumeric);
end
parse(parser,channel,varargin{:});
opt = parser.Results;

maxSlew=opt.maxGradSlew;
maxAmp=opt.maxGradAmp; % TODO: use this when no duration is supplied


if isempty(opt.area) && isempty(opt.flatArea)
    error('makeTrapezoid:invalidArguments','Must supply either ''area'', ''flatArea'' or ''amplitude''');
end
if opt.flatTime>0
    amplitude = opt.flatArea/opt.flatTime;
    riseTime = amplitude/maxSlew;
    riseTime = ceil(riseTime/mr.Sequence.GradRasterTime)*mr.Sequence.GradRasterTime;
    fallTime = riseTime;
    flatTime = opt.flatTime;
elseif opt.duration>0
    dC = 1/abs(2*maxSlew) + 1/abs(2*maxSlew);
    possible = opt.duration^2 > 4*abs(opt.area)*dC;
    assert(possible,'Requested are is too large for this gradient.');
    amplitude = ( opt.duration - sqrt(opt.duration^2 - 4*abs(opt.area)*dC) )/(2*dC);
    
    riseTime = ceil(amplitude/maxSlew/mr.Sequence.GradRasterTime)*mr.Sequence.GradRasterTime;
    fallTime = riseTime;
    flatTime = opt.duration-riseTime-fallTime;
    % Adjust amplitude (after rounding) to achieve given area
    amplitude = opt.area/(riseTime/2 + fallTime/2 + flatTime);
else
    error('makeTrapezoid:invalidArguments','Must supply a duration');
end

grad.type = 'trap';
grad.channel = find(strcmp(opt.channel,validChannels));
grad.amplitude = amplitude;
grad.riseTime = riseTime;
grad.flatTime = flatTime;
grad.fallTime = fallTime;
grad.area = amplitude*(flatTime + riseTime/2 + fallTime/2);
grad.flatArea = amplitude*flatTime;

end