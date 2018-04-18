function grad = addGradients(grads, varargin)
%addGradients Superposition of several gradients
%
%   [grad] = addGradients(grads, system) 
%   Returns the superposition of serveral gradients
%
%   See also  Sequence.addBlock  mr.opts  makeTrapezoid
%
%   Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>

persistent parser

if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'splitGradient';
	parser.addRequired('grads');
    parser.addOptional('system', mr.opts(), @isstruct);
    parser.addParamValue('maxGrad', 0, @isnumeric);
    parser.addParamValue('maxSlew', 0, @isnumeric);
end
parse(parser, grads, varargin{:});
opt = parser.Results;

maxSlew = opt.system.maxSlew;
maxGrad = opt.system.maxGrad;
if opt.maxGrad > 0
    maxGrad = opt.maxGrad;
end
if opt.maxSlew > 0
    maxSlew = opt.maxSlew;
end

% first gradient event defines channel:
channel = grads(1).channel;

% find out the general delay of all gradients
delays = [];
for ii = 1:length(grads)
    delays = [delays, grads(ii).delay];
end
common_delay = min(delays);


waveforms = {};
max_length = 0;
for ii = 1:length(grads)
    g = grads(ii);
    if strcmp(g.type, 'grad')
        waveforms{ii} = g.waveform;
    elseif strcmp(g.type, 'trap')
        times = [g.delay - common_delay ...
                 g.delay - common_delay + g.riseTime ...
                 g.delay - common_delay + g.riseTime + g.flatTime ...
                 g.delay - common_delay + g.riseTime + g.flatTime + g.fallTime];
        amplitudes = [0 g.amplitude g.amplitude 0];
        waveforms{ii} = mr.pts2waveform(times, amplitudes, ...
                                        opt.system.gradRasterTime);
    else
        error('Unknown gradient type.');
    end
    if g.delay - common_delay > 0
%         t_delay = 0:opt.system.gradRasterTime:g.delay-opt.system.gradRasterTime;
        t_delay = 0:opt.system.gradRasterTime:g.delay-common_delay-opt.system.gradRasterTime;
        waveforms{ii} = [t_delay waveforms{ii}];
    end
    num_points = length(waveforms{ii});
    if num_points > max_length
        max_length = num_points;
    end
end

w = zeros(max_length,1);
for ii = 1:length(grads)
    % SK: Matlab is so ridiculously cumbersome...
    wt = zeros(max_length, 1);
    wt(1:length(waveforms{ii})) = waveforms{ii};
    w = w + wt;
end
% figure;
% plot(w)
grad = mr.makeArbitraryGrad(channel, w, opt.system, ...
                            'maxSlew', maxSlew,...
                            'maxGrad', maxGrad,...
                            'delay', common_delay);
end