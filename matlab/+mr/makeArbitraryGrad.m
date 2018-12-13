function grad=makeArbitraryGrad(channel,varargin)
%makeArbitraryGrad Create an gradient event with arbitrary waveform.
%   g=makeArbitraryGrad(channel, waveform) Create gradient on
%   the given channel with the specified waveform.
%
%   g=makeArbitraryGrad(channel,waveform,lims) Ensure the waveform 
%   satisfies the gradient hardware constraints.
%
%   See also  Sequence.addBlock 

persistent parser

if isempty(parser)
    validChannels = {'x','y','z'};
    parser = inputParser;
    parser.FunctionName = 'makeArbitraryGrad';
    parser.addRequired('channel',...
        @(x) any(validatestring(x,validChannels)));
    parser.addRequired('waveform');
    parser.addOptional('system',mr.opts(),@isstruct);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    parser.addParamValue('delay',0,@isnumeric);
end
parse(parser,channel,varargin{:});
opt = parser.Results;

maxSlew=opt.system.maxSlew;
maxGrad=opt.system.maxGrad; % TODO: use this when no duration is supplied
if opt.maxGrad>0
    maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    maxSlew=opt.maxSlew;
end

g=opt.waveform;
slew=(g(2:end)-g(1:end-1))./opt.system.gradRasterTime;
if ~isempty(slew)
    assert(max(abs(slew))<maxSlew,'Slew rate violation (%.0f%%)',max(abs(slew))/maxSlew*100);
end
assert(max(abs(g))<maxGrad,'Gradient amplitude violation (%.0f%%)',max(abs(g))/maxGrad*100);


grad.type = 'grad';
grad.channel = opt.channel;
grad.waveform = g;
grad.delay = opt.delay;
grad.t = (0:length(g)-1)*opt.system.gradRasterTime;
grad.first = g(1); % MZ: eventually we should use extrapolation by 1/2 gradient rasters here
grad.last = g(end);

end