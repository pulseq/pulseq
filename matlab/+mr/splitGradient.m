function [grads] = splitGradient(grad, varargin)
%SplitGradient Splits a trapezoidal gradient into slew up, flat top and
%slew down.
%
%   [grads] = splitGradient(grad) 
%   Returns the individual gradient parts (slew up, flat top and slew down) 
%   as extended trapezoid gradient objects. The delays in the individual 
%   gradient events are adapted such that addGradients(...) produces an 
%   gradient equivalent to 'grad'.
%
%   See also  splitGradientAt makeExtendedTrapezoid makeTrapezoid
%             Sequence.addBlock  mr.opts
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
total_length = mr.calcDuration(grad);

if strcmp(grad.type, 'trap')
    ch = grad.channel;
    grad.delay    = round(grad.delay   /gradRasterTime)*gradRasterTime; % MZ: was ceil
    grad.riseTime = round(grad.riseTime/gradRasterTime)*gradRasterTime; % MZ: was ceil
    grad.flatTime = round(grad.flatTime/gradRasterTime)*gradRasterTime; % MZ: was ceil
    grad.fallTime = round(grad.fallTime/gradRasterTime)*gradRasterTime; % MZ: was ceil
    
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