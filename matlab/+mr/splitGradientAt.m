function [varargout] = splitGradientAt(grad, timepoint, varargin)
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

ch = grad.channel;

if strcmp(grad.type, 'grad')
    % check if we have an arbitrary gradient or an exended trapezoid
    if abs(grad.tt(1)-0.5*gradRasterTime)<1e-10 && ... 
       all(abs(grad.tt(2:end)-grad.tt(1:end-1)-gradRasterTime)<1e-10)
        % arbitrary gradient -- the most trivial conversion
        % if timepoint is out of range we have nothing to do
        if timeindex == 1 || timeindex >= length(grad.tt)
            varargout{1} = grad;
        else
            grad1=grad;
            grad2=grad;
            grad1.last=0.5*(grad.waveform(timeindex-1)+grad.waveform(timeindex)); % FIXME: retrive the double-sampling point (e.g. the corner of the trapezoid)
            grad2.first=grad1.last;
            grad2.delay=grad.delay + timepoint;
            grad1.tt=grad.tt(1:(timeindex-1));
            grad1.waveform=grad.waveform(1:(timeindex-1));
            grad2.tt=grad.tt(timeindex:end) - timepoint;
            grad2.waveform=grad.waveform(timeindex:end);

            if nargout==1
                varargout{1} = [grad1 grad2];
            else
                varargout{1} = grad1;
                varargout{2} = grad2;
            end
        end
        return; % early return
    else
        % we have an extended trapezoid -- excellent choice!
        times      = grad.tt;
        amplitudes = grad.waveform;
    end
elseif strcmp(grad.type, 'trap')    
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
else
    error('Splitting of unsupported event.');
end

% if the cutline is behind the gradient there is no second gradient to create
if timepoint >= grad.delay+times(end)
    error('trying to place the splitting time point after the end of the gradient');
end

% now we have everything in the extended trapezoid structure

% if the cutline goes through the delay we need special treatment
if timepoint < grad.delay
    times=[0 grad.delay+times];
    amplitudes = [0 amplitudes];
    grad.delay=0;
else
    timepoint = timepoint - grad.delay;
end

% sample at timepoint
amp_tp=interp1(times, amplitudes, timepoint, 'linear'); % MZ: interp1() is not OK here for the corner situation TODO: fixme! (e.g. by restoring the corners as done in waveforms_and_times())
% split the data
teps=1e-10; % we need this because of the rounding problems
times1 = [ times(times<timepoint-teps) timepoint ];
amplitudes1 = [ amplitudes(times<timepoint-teps) amp_tp ];
times2 = [ timepoint times(times>timepoint+teps) ] - timepoint;
amplitudes2 = [ amp_tp amplitudes(times>timepoint+teps) ];

% recreate gradients
grad1 = mr.makeExtendedTrapezoid(ch, opt.system, 'times', times1,...
                                  'amplitudes', amplitudes1, ...
                                  'skip_check', true); 
grad1.delay = grad.delay;
grad2 = mr.makeExtendedTrapezoid(ch, opt.system, 'times', times2,...
                                  'amplitudes', amplitudes2, ...
                                  'skip_check', true); 
grad2.delay = timepoint + grad.delay;

%grads = [grad1 grad2];
if nargout==1
    varargout{1} = [grad1 grad2];
else
    varargout{1} = grad1;
    varargout{2} = grad2;
end

end
