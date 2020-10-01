function [grads] = splitGradientAt(grad, timepoint, varargin)
%SplitGradient Splits a trapezoidal gradient into two extended trapezoids
%(currently shaped gradients) defined by the cut line.
%
%   [grads] = splitGradient(grad) 
%   Returns the two gradient parts by cutting the original 'grad' at the 
%   'timepoint' . For the input type 'trapezoid' the results are trtyurned 
%   as extended trapezoids, for 'arb' as arbitrary gradient objects. The
%    delays in the individual gradient events are adapted such that
%   addGradients(...) produces an gradient equivalent to 'grad'.
%
%   See also  splitGradient makeExtendedTrapezoid makeTrapezoid
%             Sequence.addBlock  mr.opts  
%
%   Maxim Zaitsev <maxim.zaitsev@uniklinik-freiburg.de>
%   Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>

persistent parser

if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'splitGradientAt';
	parser.addRequired('grad', @isstruct);
    parser.addRequired('timepoint', @isnumeric);
    parser.addOptional('system', mr.opts(), @isstruct);
end
parse(parser, grad, timepoint, varargin{:});
opt = parser.Results;

gradRasterTime = opt.system.gradRasterTime;
    
% round the time point to the gradient raster;
timeindex = round(timepoint / gradRasterTime);
timepoint = timeindex * gradRasterTime;
timeindex = timeindex + 1; % convert to Matlab convention

if strcmp(grad.type, 'trap')    
    ch = grad.channel;
    grad.delay    = round(grad.delay   /gradRasterTime)*gradRasterTime; % MZ: was ceil
    grad.riseTime = round(grad.riseTime/gradRasterTime)*gradRasterTime; % MZ: was ceil
    grad.flatTime = round(grad.flatTime/gradRasterTime)*gradRasterTime; % MZ: was ceil
    grad.fallTime = round(grad.fallTime/gradRasterTime)*gradRasterTime; % MZ: was ceil
    
    % prepare the extended trapezoid structure
    if grad.flatTime == 0
        times      = [0 grad.riseTime  grad.riseTime+grad.fallTime];
        amplitudes = [0 grad.amplitude 0];
    else
        times      = [0 grad.riseTime  grad.riseTime+grad.flatTime grad.riseTime+grad.flatTime+grad.fallTime];
        amplitudes = [0 grad.amplitude grad.amplitude              0];
    end
    
    % if the cutline goes through the delay we need special treatment
    if timepoint < grad.delay
        times=[0 grad.delay+times];
        amplitudes = [0 amplitudes];
        grad.delay=0;
    end
    
    % sample at timepoint
    amp_tp=interp1(times, amplitudes, timepoint, 'linear'); % MZ: interp1() is not OK here for the corner situation TODO: fixme! (e.g. by restoring the corers as done in waveforms_and_times())
    % split the data
    times1 = [ times(times<timepoint) timepoint ];
    amplitudes1 = [ amplitudes(times<timepoint) amp_tp ];
    times2 = [ timepoint times(times>timepoint) ] - timepoint;
    amplitudes2 = [ amp_tp amplitudes(times>timepoint) ];
    
    % recareate gradients
    grad1 = mr.makeExtendedTrapezoid(ch, opt.system, 'times', times1,...
                                      'amplitudes', amplitudes1, ...
                                      'skip_check', true); 
    grad1.delay = grad.delay;
    grad2 = mr.makeExtendedTrapezoid(ch, opt.system, 'times', times2,...
                                      'amplitudes', amplitudes2, ...
                                      'skip_check', true); 
    grad2.delay = timepoint;

    grads = [grad1 grad2];
elseif strcmp(grad.type, 'grad')
    % arbitrary gradient -- the most trivial conversion
    % if timepoint is out of range we have nothing to do
    if timeindex == 1 || timeindex >= length(grad.t)
        grads = grad;
    else
        grad1=grad;
        grad2=grad;
        grad1.last=grad.waveform(timeindex);
        grad2.first=grad.waveform(timeindex);
        grad2.delay=grad.delay + grad.t(timeindex);
        grad1.t=grad.t(1:(timeindex-1));
        grad1.waveform=grad.waveform(1:(timeindex-1));
        grad2.t=grad.t(timeindex:end);
        grad2.waveform=grad.waveform(timeindex:end);
        
        grads = [grad1 grad2];
    end
else
    error('Splitting of unsupported event.');
end

end
