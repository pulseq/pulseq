function [rf, gz, gzr, delay] = makeSincPulse(flip,varargin)
%makeSincPulse Create a slice selective since pulse.
%   rf=makeSincPulse(flip, 'Duration', dur) Create sinc pulse
%   with given flip angle (rad) and duration (s).
%
%   rf=makeSincPulse(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create sinc pulse with frequency offset (Hz) and phase offset (rad).
%
%   rf=makeSincPulse(..., 'ppmOffset')
%   Create arbitrary RF pulse with frequency offset specified in PPM (e.g.
%   actual frequency offset proportional to the true Larmor frequency); can
%   be combined with the 'freqOffset' specified in Hz.
%
%   [rf, gz]=makeSincPulse(...,'SliceThickness',st) Return the
%   slice select gradient corresponding to given slice thickness (m).
%
%   [rf, gz]=makeSincPulse(flip,lims,...) Create slice selection gradient 
%   with the specificed gradient limits (e.g. amplitude, slew).
%
%   [rf, gz, gzr]=makeSincPulse(flip,lims,...) Create slice selection and 
%   slice refocusing gradients with the specificed gradient limits 
%   (e.g. amplitude, slew) and taking into account 'centerpos' parameter
%
%   See also  Sequence.addBlock

validPulseUses = mr.getSupportedRfUse();

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeSincPulse';
    
    % RF params
    addRequired(parser, 'flipAngle', @isnumeric);
    addOptional(parser, 'system', [], @isstruct);
    addParamValue(parser, 'duration', 0, @isnumeric);
    addParamValue(parser, 'ppmOffset', 0, @isnumeric);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addParamValue(parser, 'timeBwProduct', 4, @isnumeric);
    addParamValue(parser, 'apodization', 0, @isnumeric);
    addParamValue(parser, 'centerpos', 0.5, @isnumeric);
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value
    % whether it is a refocusing pulse (for k-space calculation)
    addParamValue(parser, 'use', 'u', @(x) any(validatestring(x,validPulseUses)));
end
parse(parser, flip, varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

if opt.dwell==0
    opt.dwell=system.rfRasterTime;
end

if opt.duration<=0
    error('rf pulse duration must be positive');
end

BW = opt.timeBwProduct/opt.duration;
alpha = opt.apodization;
N = round(opt.duration/opt.dwell);
t = ((1:N)-0.5)'*opt.dwell;
tt = t - opt.duration*opt.centerpos;
window = (1.0-alpha+alpha*cos(2*pi*tt/opt.duration));
signal = window.*sinc(BW*tt);
flip = sum(signal)*opt.dwell*2*pi;
signal = signal*opt.flipAngle/flip;

rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.shape_dur=N*opt.dwell;
rf.ppmOffset = opt.ppmOffset;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = system.rfDeadTime;
rf.ringdownTime = system.rfRingdownTime;
rf.delay = opt.delay;
if ~isempty(opt.use)
    rf.use=opt.use;
end
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end
rf.center=mr.calcRfCenter(rf);

if nargout > 1
    assert(opt.sliceThickness > 0,'SliceThickness must be provided');
    if opt.maxGrad > 0
        system.maxGrad = opt.maxGrad;
    end
    if opt.maxSlew > 0
        system.maxSlew = opt.maxSlew;
    end
    
    amplitude = BW/opt.sliceThickness;
    area = amplitude*opt.duration;
    gz = mr.makeTrapezoid('z', system, 'flatTime', opt.duration, ...
                          'flatArea', area);
    gzr= mr.makeTrapezoid('z', system, 'Area', -area*(1-opt.centerpos)-0.5*(gz.area-area));
    if rf.delay > gz.riseTime
        gz.delay = ceil((rf.delay - gz.riseTime)/system.gradRasterTime)*system.gradRasterTime; % round-up to gradient raster
    end
    if rf.delay < (gz.riseTime+gz.delay)
        rf.delay = gz.riseTime+gz.delay; % these are on the grad raster already which is coarser 
    end
end

% v1.4 finally eliminates RF zerofilling
% if rf.ringdownTime > 0
%     tFill = (1:round(rf.ringdownTime/1e-6))*1e-6;  % Round to microsecond
%     rf.t = [rf.t rf.t(end)+tFill];
%     rf.signal = [rf.signal, zeros(size(tFill))];
% end
if rf.ringdownTime > 0 && nargout > 3
    delay=mr.makeDelay(mr.calcDuration(rf)+rf.ringdownTime);
end 

% RF amplitude check
rf_amplitude=max(abs(rf.signal));
if rf_amplitude>system.maxB1
    warning('WARNING: system maximum RF amplitude exceeded (%.01f%%)', rf_amplitude/system.maxB1*100);
end

function y = sinc(x)
    % sinc Calculate the sinc function:
    %   sinc(x) = sin(pi*x)/(pi*x)
    %
    % This is a useful helper function for those without the signal 
    % processing toolbox 
    
    i = find(x == 0);                                                              
    x(i) = 1;
    y = sin(pi*x)./(pi*x);                                                     
    y(i) = 1;   
end

end