function grad = addGradients(grads, varargin)
%addGradients Superposition of several gradients
%
%   [grad] = addGradients(grads, system) 
%   Returns the superposition of serveral gradients
%   gradients have to be passed as a cell array e.g. {g1, g2, g3}
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

if ~iscell(grads)
    error('gradients have to be passed as cell array');
end

if length(grads)<2
    error('cannot add less then two gradients');
end

% first gradient event defines channel:
channel = grads{1}.channel;

% find out the general delay of all gradients and other statistics
delays = []; % TODO: preallocate instead of grow
firsts = [];
lasts = [];
durs=[];
is_trap=[];
is_arb=[];
for ii = 1:length(grads)
    if grads{ii}.channel~=channel
        error('cannot add gradients on different channels');
    end
    delays = [delays, grads{ii}.delay];
    firsts = [firsts, grads{ii}.first];
    lasts = [lasts, grads{ii}.last];
    durs = [durs, mr.calcDuration(grads{ii})];
    is_trap = [is_trap, strcmp(grads{ii}.type,'trap')];
    if is_trap(end)
        is_arb = [is_arb, false];
    else
        tt_rast=grads{ii}.tt/opt.system.gradRasterTime+0.5;
        is_arb = [is_arb, all(abs(tt_rast-(1:length(tt_rast)))<eps)];
    end
end
common_delay = min(delays);
total_duration = max(durs);

% check if we have a set of traps with the same timing
if all(is_trap)
    % now all fields are the same so we can convert cell to a normal array
    gradsa=cell2mat(grads);
    if 1==length(unique([gradsa.delay])) && ...
       1==length(unique([gradsa.riseTime])) && ...
       1==length(unique([gradsa.flatTime])) && ...
       1==length(unique([gradsa.fallTime])) 
        % TADA, all our gradients have the same timing, so we just add
        % the amplitudes!
        grad=gradsa(1);
        grad.amplitude = sum([gradsa.amplitude]);
        grad.area = sum([gradsa.area]);
        grad.flatArea = sum([gradsa.flatArea]);
        return;
    end    
end

% check if we only have arbitrary grads on irregular time samplings
% optionally mixed with trapezoids
if all(is_trap | ~is_arb)
    % we can do quite efficient calculations and keep the shapes still rather simple
    times=[];
    for ii = 1:length(grads)
        g=grads{ii};
        if strcmp(g.type,'trap')
            times = [times cumsum([g.delay g.riseTime g.flatTime g.fallTime])];
        else
            times = [times g.delay+g.tt];
        end
    end
    times=unique(times); % unique() also sorts the array
    %times=unique(round(times/opt.system.gradRasterTime)*opt.system.gradRasterTime); % rounding to raster would be too crude here
    dt=times(2:end)-times(1:end-1);
    ieps=find(dt<eps);
    if ~isempty(ieps)
        dtx=[times(1) dt];
        dtx(ieps)=dtx(ieps)+dtx(ieps+1); % this assumes that no more than two too similar values can occur
        dtx(ieps+1)=[];
        times=cumsum(dtx);
    end
    amplitudes=zeros(size(times));
    for ii = 1:length(grads)
        g=grads{ii};
        if strcmp(g.type,'trap')
            if g.flatTime>0 % trapezoid or triange
                g.tt=cumsum([0 g.riseTime g.flatTime g.fallTime]);
                g.waveform=[0 g.amplitude g.amplitude 0];
            else
                g.tt=cumsum([0 g.riseTime g.fallTime]);
                g.waveform=[0 g.amplitude 0];
            end
        end
        tt=g.delay+g.tt;
        % fix rounding for the first and last time points
        [tmin, imin]=min(abs(tt(1)-times));
        if tmin<eps
            tt(1)=times(imin);
        end
        [tmin, imin]=min(abs(tt(end)-times));
        if tmin<eps
            tt(end)=times(imin);
        end
        % give up the "ownership" of the first point of the shape if it starts at a non-zero value
        if abs(g.waveform(1))>eps && tt(1)>eps 
            tt(1)=tt(1)+eps;
        end
        amplitudes=amplitudes+interp1(tt,g.waveform,times,'linear',0);
    end
    grad=mr.makeExtendedTrapezoid(channel,'amplitudes',amplitudes,'times',times,'system',opt.system);
    return;
end

% OK, here we convert everything to a regularly-sampled waveform
waveforms = {};
max_length = 0;
for ii = 1:length(grads)
    g = grads{ii};
    if strcmp(g.type, 'grad')
        if is_arb(ii)
            waveforms{ii} = g.waveform;
        else
            waveforms{ii} = mr.pts2waveform(g.tt, g.waveform, opt.system.gradRasterTime);
        end
    elseif strcmp(g.type, 'trap')
        if (g.flatTime>0) % triangle or trapezoid
            times = [g.delay - common_delay ...
                     g.delay - common_delay + g.riseTime ...
                     g.delay - common_delay + g.riseTime + g.flatTime ...
                     g.delay - common_delay + g.riseTime + g.flatTime + g.fallTime];
            amplitudes = [0 g.amplitude g.amplitude 0];
        else
            times = [g.delay - common_delay ...
                     g.delay - common_delay + g.riseTime ...
                     g.delay - common_delay + g.riseTime + g.fallTime];
            amplitudes = [0 g.amplitude 0];
        end
        waveforms{ii} = mr.pts2waveform(times, amplitudes, opt.system.gradRasterTime);
    else
        error('Unknown gradient type.');
    end
    warning('addGradient(): potentially incorrect handling of delays... TODO: fixme!');
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

grad = mr.makeArbitraryGrad(channel, w, opt.system, ...
                            'maxSlew', maxSlew,...
                            'maxGrad', maxGrad,...
                            'delay', common_delay);
                       
% fix the first and the last values
% first is defined by the sum of firsts with the minimal delay (common_delay)
% last is defined by the sum of lasts with the maximum duration (total_duration)
grad.first=sum(firsts(delays==common_delay));
grad.last=sum(lasts(durs==total_duration));

end