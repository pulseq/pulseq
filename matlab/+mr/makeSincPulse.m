function [rf, gz] = makeSincPulse(flip,varargin)
%makeSincPulse Create a slice selective since pulse.
%   rf=makeSincPulse(flip, 'Duration', dur) Create sinc pulse
%   with given flip angle and duration.
%
%   rf=makeSincPulse(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create sinc pulse with frequency offset and phase offset.
%
%   [rf, gz]=makeSincPulse(...,'SliceThickness',st) Return the
%   corresponding slice select gradient suitable for MRI.
%
%   See also  Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeSincPulse';
    
    % RF params
    addRequired(parser,'flipAngle',@isnumeric);
    addParamValue(parser,'duration',0,@isnumeric);
    addParamValue(parser,'freqOffset',0,@isnumeric);
    addParamValue(parser,'phaseOffset',0,@isnumeric);
    addParamValue(parser,'timeBwProduct',4,@isnumeric);
    addParamValue(parser,'apodization',0,@isnumeric);
    % Slice params
    defaultOpts=mr.opts();
    addParamValue(parser,'maxGradAmp',defaultOpts.maxGradAmp,@isnumeric);
    addParamValue(parser,'maxGradSlew',defaultOpts.maxGradSlew,@isnumeric);
    addParamValue(parser,'sliceThickness',0,@isnumeric);
end
parse(parser,flip,varargin{:});
opt = parser.Results;

BW=opt.timeBwProduct/opt.duration;
alpha=opt.apodization;
N=round(opt.duration/1e-6);
t = (1:N)*mr.Sequence.RfRasterTime;
tt = t - opt.duration/2;
window = (1.0-alpha+alpha*cos(2*pi*tt/opt.duration));
signal = window.*sinc(BW*tt);
flip=sum(signal)*mr.Sequence.RfRasterTime*2*pi;
signal=signal*opt.flipAngle/flip;

rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;

if nargout>1
    assert(opt.sliceThickness>0,'SliceThickness must be provided');
    
    amplitude = BW/opt.sliceThickness;
    area = amplitude*opt.duration;
    gz = mr.makeTrapezoid('z','flatTime',opt.duration,'flatArea',area);
    
    tFill = (1:round(gz.riseTime/1e-6))*1e-6;   % Round to microsecond
    rf.t = [tFill rf.t+tFill(end) tFill+rf.t(end)+tFill(end)];
    rf.signal=[zeros(size(tFill)), rf.signal, zeros(size(tFill))];
end
end