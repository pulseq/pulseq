function [rf, gz] = makeSincPulse(flip,varargin)
%makeSincPulse Create a slice selective since pulse.
%   rf=makeSincPulse(flip, 'Duration', dur) Create sinc pulse
%   with given flip angle (rad) and duration (s).
%
%   rf=makeSincPulse(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create sinc pulse with frequency offset (Hz) and phase offset (rad).
%
%   [rf, gz]=makeSincPulse(...,'SliceThickness',st) Return the
%   slice select gradient corresponding to given slice thickness (m).
%
%   [rf, gz]=makeSincPulse(flip,lims,...) Create slice selection gradient 
%   with the specificed gradient limits (e.g. amplitude, slew).
%
%   See also  Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeSincPulse';
    
    % RF params
    addRequired(parser,'flipAngle',@isnumeric);
    addOptional(parser,'system',mr.opts(),@isstruct);
    addParamValue(parser,'duration',0,@isnumeric);
    addParamValue(parser,'freqOffset',0,@isnumeric);
    addParamValue(parser,'phaseOffset',0,@isnumeric);
    addParamValue(parser,'timeBwProduct',4,@isnumeric);
    addParamValue(parser,'apodization',0,@isnumeric);
    % Slice params
    addParamValue(parser,'maxGrad',0,@isnumeric);
    addParamValue(parser,'maxSlew',0,@isnumeric);
    addParamValue(parser,'sliceThickness',0,@isnumeric);
end
parse(parser,flip,varargin{:});
opt = parser.Results;

BW=opt.timeBwProduct/opt.duration;
alpha=opt.apodization;
N=round(opt.duration/1e-6);
t = (1:N)*opt.system.rfRasterTime;
tt = t - opt.duration/2;
window = (1.0-alpha+alpha*cos(2*pi*tt/opt.duration));
signal = window.*sinc(BW*tt);
flip=sum(signal)*opt.system.rfRasterTime*2*pi;
signal=signal*opt.flipAngle/flip;

rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;

fillTime=0;
if nargout>1
    assert(opt.sliceThickness>0,'SliceThickness must be provided');
    if opt.maxGrad>0
        opt.system.maxGrad = opt.maxGrad;
    end
    if opt.maxSlew>0
        opt.system.maxSlew = opt.maxSlew;
    end
    
    amplitude = BW/opt.sliceThickness;
    area = amplitude*opt.duration;
    gz = mr.makeTrapezoid('z',opt.system,'flatTime',opt.duration,'flatArea',area);
    
    % Pad RF pulse with zeros during gradient ramp up
    fillTime=gz.riseTime;
    tFill = (1:round(fillTime/1e-6))*1e-6;   % Round to microsecond
    rf.t = [tFill rf.t+tFill(end) ];
    rf.signal=[zeros(size(tFill)), rf.signal];
end

% Add dead time to start of pulse, if required
if fillTime<rf.deadTime
    fillTime=rf.deadTime-fillTime;
    tFill = (1:round(fillTime/1e-6))*1e-6;   % Round to microsecond
    rf.t = [tFill rf.t+tFill(end) ];
    rf.signal=[zeros(size(tFill)), rf.signal];
end

if rf.ringdownTime>0
    tFill = (1:round(rf.ringdownTime/1e-6))*1e-6;  % Round to microsecond
    rf.t = [rf.t rf.t(end)+tFill];
    rf.signal = [rf.signal, zeros(size(tFill))];
end
    

function y=sinc(x)
    % sinc Calculate the sinc function:
    %   sinc(x) = sin(pi*x)/(pi*x)
    %
    % This is a useful helper function for those without the signal processing
    % toolbox 
    
    i=find(x==0);                                                              
    x(i)= 1;
    y = sin(pi*x)./(pi*x);                                                     
    y(i) = 1;   
end

end