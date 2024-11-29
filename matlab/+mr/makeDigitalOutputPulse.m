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
    
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'duration', 0, @isnumeric); % will replace with gradRadterTime below
    addParamValue(parser, 'system', [], @isstruct); 
end

if nargin<1
    error('makeDigitalOutputPulse:invalidArguments','Must supply a channel');
end

parse(parser, varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

channel_num=find(strcmp(channel,{'osc0','osc1','ext1'}));
assert(~isempty(channel_num) && channel_num>0,'makeDigitalOutputPulse:invalidChannel',...
    'Channel (%s) is invalid',channel);
trig.type = 'output';
trig.channel=channel;
trig.delay = opt.delay;
trig.duration = opt.duration;
if (trig.duration<=system.gradRasterTime)
    trig.duration=system.gradRasterTime;
end

end
