function trig = makeDigitalOutputPulse(channel, varargin)
%makeDigitalOutputPulse Create a digital output pulse event a.k.a. trigger.
%   trig_out=makeDigitalOutputPulse() Create an output trigger event on a  
%                            given channel with optional given delay and 
%                            duration.
%                            Possible channel values: 'osc0','osc1','ext1'
%
%   See also  Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeDigitalOutputPulse';
    
    addOptional(parser, 'delay', 0, @isnumeric);
    addOptional(parser, 'duration', 0, @isnumeric); % will replace with gradRadterTime below
    addOptional(parser, 'system', mr.opts(), @isstruct); 
end

if nargin<1
    error('makeDigitalOutputPulse:invalidArguments','Must supply a channel');
end

parse(parser, varargin{:});
opt = parser.Results;

channel_num=find(strcmp(channel,{'osc0','osc1','ext1'}));
assert(~isempty(channel_num) && channel_num>0,'makeDigitalOutputPulse:invalidChannel',...
    'Channel (%s) is invalid',channel);
trig.type = 'output';
trig.channel=channel;
trig.delay = opt.delay;
trig.duration = opt.duration;
if (trig.duration<=opt.system.gradRasterTime)
    trig.duration=opt.system.gradRasterTime;
end

end
