function [rf, gz, gzr, delay] = makeAdiabaticPulse(type,varargin)
%makeAdiabaticPulse make an aiabatic inversion pulse
%     a wrapper to a python function(see below). See supported params below
%     in the 'parser' section. Currently it will probably only work on
%     Linux. On my system I could install the required Python library by
%     executing "pip3 install sigpy"
%     type must be one of {'hypsec','wurst'}
%     BE CAREFUL, some parameters only affect certain pulse types and are
%     ignored for other; e.g. bandwidth is ignored if type='hypsec'.
%
%     hypsec(n=512, beta=800, mu=4.9, dur=0.012)
%         Design a hyperbolic secant adiabatic pulse.
%         
%         mu * beta becomes the amplitude of the frequency sweep
%         
%         Args:
%             n (int): number of samples (should be a multiple of 4).
%             beta (float): AM waveform parameter.
%             mu (float): a constant, determines amplitude of frequency sweep.
%             dur (float): pulse time (s).
%         
%         Returns:
%             2-element tuple containing
%         
%             - **a** (*array*): AM waveform.
%             - **om** (*array*): FM waveform (radians/s).
%         
%         References:
%             Baum, J., Tycko, R. and Pines, A. (1985). 'Broadband and adiabatic
%             inversion of a two-level system by phase-modulated pulses'.
%             Phys. Rev. A., 32:3435-3447.
%     
%     wurst(n=512, n_fac=40, bw=40000.0, dur=0.002)
%         Design a WURST (wideband, uniform rate, smooth truncation) adiabatic
%          inversion pulse
%         
%         Args:
%             n (int): number of samples (should be a multiple of 4).
%             n_fac (int): power to exponentiate to within AM term. ~20 or greater is
%              typical.
%             bw (float): pulse bandwidth.
%             dur (float): pulse time (s).
%         
%         
%         Returns:
%             2-element tuple containing
%            - **a** (*array*): AM waveform.
%             - **om** (*array*): FM waveform (radians/s).
%         
%         References:
%             Kupce, E. and Freeman, R. (1995). 'Stretched Adiabatic Pulses for
%             Broadband Spin Inversion'.
%             J. Magn. Reson. Ser. A., 117:246-256.

            
validPulseTypes = {'hypsec','wurst'};
validPulseUses = mr.getSupportedRfUse();

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeAdiabaticPulse';
    
    % RF params
    addRequired(parser, 'type', @(x) any(validatestring(x,validPulseTypes)));
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'duration', 10e-3, @isnumeric);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addParamValue(parser, 'beta', 800, @isnumeric);
    addParamValue(parser, 'mu', 4.9, @isnumeric);
    addParamValue(parser, 'n_fac', 40, @isnumeric);
    addParamValue(parser, 'bandwidth', 40000, @isnumeric);
    addParamValue(parser, 'adiabaticity', 4, @isnumeric);
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
end

parse(parser, type, varargin{:});
opt = parser.Results;

if opt.dwell==0
    opt.dwell=opt.system.rfRasterTime;
end

% find python (probably only works on linux, maybe also mac)
[status, result]=system('which python');
if status==0
    python=strip(result);
else
    [status, result]=system('which python3');
    if status==0
        python=strip(result);
    else
        error('python executable not found');
    end
end

Nraw = round(opt.duration/opt.dwell+eps);
N = floor(Nraw/4)*4; % number of points must be divisible by four -- this is a requirement of the underlying library

switch type
    case 'hypsec'
        cmd=[python ' -c $''import sigpy.mri.rf\npulse=sigpy.mri.rf.hypsec(' ... % hypsec(n=512, beta=800, mu=4.9, dur=0.012)
            'n=' num2str(N) ',beta=' num2str(opt.beta) ',' ...
            'mu=' num2str(opt.mu) ',dur=' num2str(opt.duration) ...
            ')\nprint(*pulse[0])\nprint(*pulse[1])'''];
    case 'wurst'
        cmd=[python ' -c $''import sigpy.mri.rf\npulse=sigpy.mri.rf.wurst(' ... % wurst(n=512, n_fac=40, bw=40000.0, dur=0.002)
            'n=' num2str(N) ',n_fac=' num2str(opt.n_fac) ',' ...
            'bw=' num2str(opt.bandwidth) ',dur=' num2str(opt.duration) ...
            ')\nprint(*pulse[0])\nprint(*pulse[1])'''];
    otherwise
        error('unsupported adiabatic pulse type');
end

[status, result]=system(cmd);

if status~=0
    error('executing python command failed');
end

lines = regexp(result,'\n','split'); % the response from the python call contains some garbage 
am=str2num(lines{1});
fm=str2num(lines{2});
assert(length(am)==length(fm));

pm=cumsum(fm)*opt.dwell;

[dfm,ifm]=min(abs(fm)); % find the center of the pulse
% we will also use the ocasion to find the rate of change of the frequency
% at the center of the pulse
if dfm==0 
    pm0=pm(ifm);
    am0=am(ifm);
    roc_fm0=abs(fm(ifm+1)-fm(ifm-1))/2/opt.dwell;
else
    % we need to bracket the zero-crossing
    if fm(ifm)*fm(ifm+1) < 0
        b=1;
    else
        b=-1;
    end
    pm0=(pm(ifm)*fm(ifm+b)-pm(ifm+b)*fm(ifm))/(fm(ifm+b)-fm(ifm));
    am0=(am(ifm)*fm(ifm+b)-am(ifm+b)*fm(ifm))/(fm(ifm+b)-fm(ifm));
    roc_fm0=abs(fm(ifm)-fm(ifm+b))/opt.dwell;
end
pm=pm-pm0;
a=(roc_fm0*opt.adiabaticity)^0.5/2/pi/am0;

signal = a*am.*exp(1i*pm);

if (N~=Nraw)
    % we need to pad the signal vector
    Npad=Nraw-N;
    signal=[zeros(1,Npad-floor(Npad/2)) signal zeros(1,floor(Npad/2))];
    N=Nraw;
end

%BW = opt.timeBwProduct/opt.duration;
t = ((1:N)-0.5)*opt.dwell;
%flip = abs(sum(signal))*opt.dwell*2*pi;

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
else
    rf.use='inversion';
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
    
    switch type
        case 'hypsec'
            BW=mr.calcRfBandwidth(rf,0.1);
        case 'wurst'
            BW=opt.bandwidth;
        otherwise
            error('unsupported pulse type')
    end
    centerpos=mr.calcRfCenter(rf);
    
    amplitude = BW/opt.sliceThickness;
    area = amplitude*opt.duration;
    gz = mr.makeTrapezoid('z', opt.system, 'flatTime', opt.duration, ...
                          'flatArea', area);
    gzr= mr.makeTrapezoid('z', opt.system, 'Area', -area*(1-centerpos)-0.5*(gz.area-area));
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

end