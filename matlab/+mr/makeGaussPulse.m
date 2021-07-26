function [rf, gz, gzr, delay] = makeGaussPulse(flip,varargin)
%makeGaussPulse Create a [optionally slice selective] Gauss pulse.
%   rf=makeGaussPulse(flip, 'Duration', dur) Create Gauss pulse
%   with given flip angle (rad) and duration (s).
%
%   rf=makeGaussPulse(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create Gauss pulse with frequency offset (Hz) and phase offset (rad).
%
%   [rf, gz]=makeGaussPulse(...,'SliceThickness',st) Return the
%   slice select gradient corresponding to given slice thickness (m).
%
%   [rf, gz]=makeGaussPulse(flip,lims,...) Create slice selection gradient 
%   with the specificed gradient limits (e.g. amplitude, slew).
%
%   [rf, gz, gzr]=makeGaussPulse(flip,lims,...) Create slice selection and 
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
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'duration', 0, @isnumeric);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addOptional(parser, 'timeBwProduct', 3, @isnumeric);
    addOptional(parser, 'bandwidth', 0, @isnumeric);
    addParamValue(parser, 'apodization', 0, @isnumeric);
    addParamValue(parser, 'centerpos', 0.5, @isnumeric);
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
end
parse(parser, flip, varargin{:});
opt = parser.Results;

if opt.dwell==0
    opt.dwell=opt.system.rfRasterTime;
end

if opt.bandwidth == 0
    BW = opt.timeBwProduct/opt.duration;
else
    BW = opt.bandwidth;
end
alpha = opt.apodization;
N = round(opt.duration/opt.dwell);
t = ((1:N)-0.5)*opt.dwell;
tt = t - opt.duration*opt.centerpos;
window = (1.0-alpha+alpha*cos(2*pi*tt/opt.duration));
signal = window.*gauss(BW*tt);
flip = sum(signal)*opt.dwell*2*pi;
signal = signal*opt.flipAngle/flip;

rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.shape_dur=N*opt.dwell;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;
rf.delay = opt.delay;
if ~isempty(opt.use)
    rf.use=opt.use;
end
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end

if nargout > 1
    assert(opt.sliceThickness > 0,'SliceThickness must be provided');
    if opt.maxGrad > 0
        opt.system.maxGrad = opt.maxGrad;
    end
    if opt.maxSlew > 0
        opt.system.maxSlew = opt.maxSlew;
    end
    
    amplitude = BW/opt.sliceThickness;
    area = amplitude*opt.duration;
    gz = mr.makeTrapezoid('z', opt.system, 'flatTime', opt.duration, ...
                          'flatArea', area);
    gzr= mr.makeTrapezoid('z', opt.system, 'Area', -area*(1-opt.centerpos)-0.5*(gz.area-area));
    if rf.delay > gz.riseTime
        gz.delay = ceil((rf.delay - gz.riseTime)/opt.system.gradRasterTime)*opt.system.gradRasterTime; % round-up to gradient raster
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

function y = gauss(x)
    % gauss Calculate the Gaussian function:
    %   gauss(x) = exp(-pi*x^2)
    %
    % This is a useful helper function for those without the signal 
    % processing toolbox 
    
    y = exp(-pi*x.^2);
end

end