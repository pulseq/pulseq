function [grads] = splitGradient(grad, varargin)
%SplitGradient Splits a trapezoidal Gradient into slew up, flat top and
%slew down.
%
%   [grads] = splitGradient(grad) 
%   Returns the individual gradient parts (slew up, flat top and slew down) 
%   as arbitrary gradient objects. The delays in the individual gradient
%   events are adapted such that addGradients(...) produces an gradient
%   equivalent to 'grad'.
%
%   See also  Sequence.addBlock  mr.opts  makeTrapezoid
%
%   Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>

persistent parser

if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'splitGradient';
	parser.addRequired('grad', @isstruct);
    parser.addOptional('system', mr.opts(), @isstruct);
end
parse(parser, grad, varargin{:});
opt = parser.Results;

gradRasterTime = opt.system.gradRasterTime;
mr.calcDuration(grad)
total_length = mr.calcDuration(grad);

if strcmp(grad.type, 'trap')
    ch = grad.channel;
    grad.riseTime = ceil(grad.riseTime/gradRasterTime)*gradRasterTime;
    grad.delay    = ceil(grad.delay   /gradRasterTime)*gradRasterTime;
    grad.flatTime = ceil(grad.flatTime/gradRasterTime)*gradRasterTime;
    grad.fallTime = ceil(grad.fallTime/gradRasterTime)*gradRasterTime;
    % ramp up
    times = [0, grad.riseTime];
    amplitudes = [0 grad.amplitude];
    rampup = mr.makeExtendedTrapezoid(ch, opt.system, 'times', times,...
                                      'amplitudes', amplitudes, ...
                                      'skip_check', true);
    rampup.delay = grad.delay;
%     rampup.t = rampup.t;


    % ramp down
    times = [0, grad.fallTime];
    amplitudes = [grad.amplitude 0];
    rampdown = mr.makeExtendedTrapezoid(ch, opt.system, 'times', times,...
                                        'amplitudes', amplitudes, ...
                                        'skip_check', true);
    rampdown.delay = total_length - grad.fallTime;
    rampdown.t = rampdown.t*gradRasterTime;
    
    % flattop
    flattop = struct;
    flattop.type = 'grad';
    flattop.channel = ch;
%     flattop.delay = (grad.delay + grad.riseTime + gradRasterTime); 
    flattop.delay = (grad.delay + grad.riseTime); 
    flattop.t = 0:gradRasterTime:(rampdown.delay-1*gradRasterTime-grad.delay-grad.riseTime);
    flattop.waveform = grad.amplitude*ones(size(flattop.t));
    flattop.first = grad.amplitude;
    flattop.last = grad.amplitude;

    grads = [rampup flattop rampdown];
elseif strcmp(grad.type, 'grad')
    error('Splitting of arbitrary gradients is not implemented yet.');
else
    error('Splitting of unsupported event.');
end

end