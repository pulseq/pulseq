function trig = makeTrigger(channel, varargin)
%makeTrigger Create a trigger halt event.
%   trigger=makeTrigger() Create a trigger event for a synchronisation with
%                         an external signal from a given channel with an 
%                         optional given delay prio to the sync and
%                         duration after the sync.
%                         Possible channel values: 'physio1','physio2' (Siemens specific)
%
%   See also  Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeTrigger';
    
    addOptional(parser, 'delay', 0, @isnumeric);
    addOptional(parser, 'duration', 0, @isnumeric); % will replace with gradRadterTime below
    addOptional(parser, 'system', [], @isstruct); 
end

if nargin<1
    error('makeTrigger:invalidArguments','Must supply a channel');
end

parse(parser, varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

channel_num=find(strcmp(channel,{'physio1','physio2'}));
assert(~isempty(channel_num) && channel_num>0,'makeTrigger:invalidChannel',...
    'Channel (%s) is invalid',channel);
trig.type = 'trigger';
trig.channel=channel;
trig.delay = opt.delay;
trig.duration = opt.duration;
if (trig.duration<=system.gradRasterTime)
    trig.duration=system.gradRasterTime;
end

end
