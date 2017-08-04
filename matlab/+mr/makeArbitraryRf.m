function [rf, gz] = makeArbitraryRf(signal,flip,varargin)
%makeArbitraryRf Create an RF pulse with the given pulse shape.
%   rf=makeArbitraryRf(singal, flip) Create RF pulse with complex signal 
%   and given flip angle (in radians)
%
%   rf=makeArbitraryRf(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create block pulse with frequency offset and phase offset.
%
%   [rf, gz]=makeArbitraryRf(..., 'Bandwidth', bw, 'SliceThickness', st) 
%   Create RF pulse and corresponding slice select gradient. The bandwidth
%   of the pulse must be given for the specified shape.
%
%   See also  Sequence.makeSincPulse, Sequence.addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeArbitraryRf';
    
    % RF params
    addRequired(parser,'signal',@isnumeric);
    addRequired(parser,'flipAngle',@isnumeric);
    addOptional(parser,'system',mr.opts(),@isstruct);
    addParamValue(parser,'freqOffset',0,@isnumeric);
    addParamValue(parser,'phaseOffset',0,@isnumeric);
    addParamValue(parser,'timeBwProduct',0,@isnumeric);
    addParamValue(parser,'bandwidth',0,@isnumeric);
    % Slice params
    addParamValue(parser,'maxGrad',0,@isnumeric);
    addParamValue(parser,'maxSlew',0,@isnumeric);
    addParamValue(parser,'sliceThickness',0,@isnumeric);
end
parse(parser,signal,flip,varargin{:});
opt = parser.Results;

signal = signal./sum(signal.*opt.system.rfRasterTime)*flip/(2*pi);

N=length(signal);
duration = N*opt.system.rfRasterTime;
t = (1:N)*opt.system.rfRasterTime;

rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;

if nargout>1
    assert(opt.sliceThickness>0,'SliceThickness must be provided');
    assert(opt.bandwidth>0,'Bandwidth of pulse must be provided');
    if opt.maxGrad>0
        opt.system.maxGrad = opt.maxGrad;
    end
    if opt.maxSlew>0
        opt.system.maxSlew = opt.maxSlew;
    end
    
    BW=opt.bandwidth;
    if opt.timeBwProduct>0
        BW=opt.timeBwProduct/duration;
    end

    amplitude = BW/opt.sliceThickness;
    area = amplitude*opt.duration;
    gz = mr.makeTrapezoid('z',opt.system,'flatTime',opt.duration,'flatArea',area);
    
    tFill = (1:round(gz.riseTime/1e-6))*1e-6;   % Round to microsecond
    rf.t = [tFill rf.t+tFill(end) tFill+rf.t(end)+tFill(end)];
    rf.signal=[zeros(size(tFill)), rf.signal, zeros(size(tFill))];
end

if rf.ringdownTime>0
    tFill = (1:round(rf.ringdownTime/1e-6))*1e-6;  % Round to microsecond
    rf.t = [rf.t rf.t(end)+tFill];
    rf.signal = [rf.signal, zeros(size(tFill))];
end

end