function grad=makeTrapezoid(channel, varargin)
%makeTrapezoid Create a trapezoid gradient event.
%
%   PURPOSE
%     Build a trapezoidal gradient event struct (rise / flat-top / fall)
%     for a given logical channel. The returned struct is consumed by
%     mr.Sequence/addBlock to add the gradient to a sequence.
%
%   SIGNATURES
%     g = mr.makeTrapezoid(channel, ...)                       % uses mr.opts() defaults
%     g = mr.makeTrapezoid(channel, system, ...)               % system as 2nd positional arg
%     g = mr.makeTrapezoid(channel, ..., 'system', system)     % system as name/value
%     g = mr.makeTrapezoid(channel, ..., 'Duration', d, 'Area', a)
%     g = mr.makeTrapezoid(channel, ..., 'FlatTime', ft, 'FlatArea', fa)
%     g = mr.makeTrapezoid(channel, ..., 'FlatTime', ft, 'Amplitude', amp)
%     g = mr.makeTrapezoid(channel, ..., 'Area', a)            % shortest possible timing
%
%     Exactly one of 'Area', 'FlatArea', or 'Amplitude' must be supplied.
%     Timing is determined by 'FlatTime' if given, else 'Duration' if given,
%     else the shortest realizable timing for the requested 'Area'.
%     Parameter names are case-insensitive.
%
%   INPUTS
%     channel      [required]    char, 'x'|'y'|'z'
%     system       [optional]    struct from mr.opts; defaults to mr.opts() if omitted
%     'Duration'   [name/value]  double, total duration including ramps, seconds, >0
%     'Area'       [name/value]  double, total gradient area including ramps, 1/m
%     'FlatTime'   [name/value]  double, flat-top duration, seconds
%     'FlatArea'   [name/value]  double, flat-top-only area, 1/m
%     'Amplitude'  [name/value]  double, flat-top amplitude, Hz/m
%     'maxGrad'    [name/value]  double, override system.maxGrad, Hz/m
%     'maxSlew'    [name/value]  double, override system.maxSlew, Hz/m/s
%     'riseTime'   [name/value]  double, force rise time, seconds
%     'fallTime'   [name/value]  double, force fall time, seconds (requires riseTime)
%     'delay'      [name/value]  double, pre-event delay, seconds, default 0
%
%   OUTPUT
%     grad  struct with fields:
%       .type       char,    always 'trap' (includes degenerate triangle, flatTime=0)
%       .channel    char,    'x'|'y'|'z'
%       .amplitude  double,  flat-top amplitude, Hz/m
%       .riseTime   double,  ramp-up duration, seconds
%       .flatTime   double,  flat-top duration, seconds (0 for triangular)
%       .fallTime   double,  ramp-down duration, seconds
%       .area       double,  total area including ramps, 1/m
%                            (= amplitude * (flatTime + riseTime/2 + fallTime/2))
%       .flatArea   double,  flat-top-only area, 1/m  (= amplitude * flatTime)
%       .delay      double,  pre-event delay, seconds
%       .first      double,  gradient value at t=0,  Hz/m  (always 0 for trap)
%       .last       double,  gradient value at end,  Hz/m  (always 0 for trap)
%
%   ERRORS
%     makeTrapezoid:invalidArguments
%       - 'fallTime' specified without 'riseTime'.
%       - Not exactly one of 'Area' / 'FlatArea' / 'Amplitude' supplied.
%       - 'FlatTime' supplied without 'FlatArea' or 'Amplitude'.
%       - Neither 'Area' nor 'Duration' supplied.
%     makeTrapezoid:invalidDuration
%       - Requested area cannot be realized within the requested duration
%         under maxGrad/maxSlew. Error message reports the minimum
%         achievable duration in microseconds.
%     makeTrapezoid:invalidAmplitude
%       - Computed amplitude exceeds maxGrad.
%     Assertion failure
%       - With explicit riseTime+duration: duration < riseTime+fallTime,
%         or computed amplitude exceeds maxGrad.
%
%   NOTES
%     - All ramp/flat times are rounded up to system.gradRasterTime.
%     - Internal storage uses Hz/m and Hz/m/s regardless of the units
%       passed to mr.opts (mr.opts converts on input). Use mr.convert
%       if you need physical units (mT/m, T/m/s, etc.).
%     - Caches an inputParser in a persistent variable for performance;
%       no other global state.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     Nx = 256;  fov = 256e-3;  deltak = 1/fov;
%     % Readout gradient with fixed flat-top area
%     gx    = mr.makeTrapezoid('x', sys, 'FlatArea', Nx*deltak, 'FlatTime', 6.4e-3);
%     % Matching prephaser: half the (negative) area of the readout
%     gxPre = mr.makeTrapezoid('x', sys, 'Area', -gx.area/2, 'Duration', 1e-3);
%
%   SEE ALSO
%     mr.opts, mr.makeExtendedTrapezoid, mr.makeArbitraryGrad,
%     mr.calcDuration, mr.Sequence/addBlock

persistent parser

if isempty(parser)
    validChannels = {'x','y','z'};
    parser = mr.aux.InputParserCompat;
    parser.FunctionName = 'makeTrapezoid';
    parser.addRequired('channel',...
        @(x) any(validatestring(x,validChannels)));
    parser.addOptional('system',[],@isstruct);
    parser.addParamValue('duration',0,@(x)(isnumeric(x) && x>0));
    parser.addParamValue('area',[],@isnumeric);
    parser.addParamValue('flatTime',[],@isnumeric);
    parser.addParamValue('flatArea',[],@isnumeric);
    parser.addParamValue('amplitude',[],@isnumeric);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    parser.addParamValue('riseTime',0,@isnumeric);
    parser.addParamValue('fallTime',0,@isnumeric);
    parser.addParamValue('delay',0,@isnumeric);
end
parse(parser,channel,varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

maxSlew=system.maxSlew;
%riseTime=system.riseTime;
maxGrad=system.maxGrad;
fallTime = [];
riseTime = [];

if opt.maxGrad>0
    maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    maxSlew=opt.maxSlew;
end
if opt.riseTime>0
    riseTime=opt.riseTime;
end
if opt.fallTime>0
    if isempty(riseTime)
        error('makeTrapezoid:invalidArguments','Must always supply ''riseTime'' if ''fallTime'' is specified explicitly.');
    end
    fallTime=opt.fallTime;
end


if (isempty(opt.area)+isempty(opt.flatArea)+isempty(opt.amplitude))~=2
    error('makeTrapezoid:invalidArguments','Must supply either ''area'', ''flatArea'' or ''amplitude'', and only one of the three may be specified');
end
if ~isempty(opt.flatTime) % MZ was: opt.flatTime>0
    if ~isempty(opt.amplitude)
        amplitude = opt.amplitude;
    else
        if isempty(opt.flatArea)
            error('makeTrapezoid:invalidArguments','When ''flatTime'' is provided either ''flatArea'' or ''amplitude'' must be provided as well; you may consider providing ''duration'', ''area'' and optionally ramp times instead.');
        end
        amplitude = opt.flatArea/opt.flatTime;
    end
    if isempty(riseTime)
        riseTime = abs(amplitude)/maxSlew;
        riseTime = ceil(riseTime/system.gradRasterTime)*system.gradRasterTime;
        if riseTime==0
            riseTime=system.gradRasterTime;
        end
    end
    if isempty(fallTime)
        fallTime = riseTime;
    end
    flatTime = opt.flatTime;
elseif opt.duration>0
    if ~isempty(opt.amplitude)
        amplitude = opt.amplitude;
    else
        if isempty(riseTime)
            dC = 1/abs(2*maxSlew) + 1/abs(2*maxSlew);
            possible = opt.duration^2 > 4*abs(opt.area)*dC;
            if ~possible
                [~, t1, t2, t3]=calcShortestParamsForArea(opt.area,maxSlew,maxGrad,system.gradRasterTime);
                error('makeTrapezoid:invalidDuration',['Requested area is too large for this gradient. Minimum required duration for this area (accounting for the gradient raster time) is ' num2str((t1+t2+t3)*1e6) 'us']);
            end
            amplitude = ( opt.duration - sqrt(opt.duration^2 - 4*abs(opt.area)*dC) )/(2*dC);
        else
            if isempty(fallTime)
                fallTime = riseTime;
            end    
            amplitude = opt.area/(opt.duration-0.5*riseTime-0.5*fallTime);
            possible = opt.duration>=(riseTime+fallTime) & abs(amplitude)<maxGrad;
            assert(possible,['Requested area is too large for this gradient duration. Probably amplitude is violated (' num2str(round(abs(amplitude)/maxGrad*100)) '%)']);    
        end    
    end
    if isempty(riseTime)
        riseTime = ceil(abs(amplitude)/maxSlew/system.gradRasterTime)*system.gradRasterTime;
        if(riseTime==0)
            riseTime=system.gradRasterTime;
        end
    end
    if isempty(fallTime)
        fallTime = riseTime;
    end
    flatTime = opt.duration-riseTime-fallTime;
    if isempty(opt.amplitude)
        % Adjust amplitude (after rounding) to achieve given area
        amplitude = opt.area/(riseTime/2 + fallTime/2 + flatTime);
    end
else
    if isempty(opt.area)
        error('makeTrapezoid:invalidArguments','Must supply area or duration');
    else
        % call the local function to calculate the shortest timing
        [amplitude, riseTime, flatTime, fallTime]=calcShortestParamsForArea(opt.area,maxSlew,maxGrad,system.gradRasterTime);
    end
end
if abs(amplitude)>maxGrad
    if isempty(opt.area)
        error('makeTrapezoid:invalidAmplitude',['Amplitude violation (' num2str(round(abs(amplitude)/maxGrad*100)) '%%)']);
    else
        % this error can only be produced by the failed trapezoid
        % calculation with the specified area (leading to exceedingly high
        % amplitude), the triangular blip error should have occured around
        % line 103
        [~, t1, t2, t3]=calcShortestParamsForArea(opt.area,maxSlew,maxGrad,system.gradRasterTime);
        error('makeTrapezoid:invalidDuration',['Requested duration is too short for the area to be realized within system limits. Minimum duration for this trapezoid (accounting for the gradient raster time) is ' num2str((t1+t2+t3)*1e6) ' us']);
    end
end

grad.type = 'trap';
grad.channel = opt.channel;
grad.amplitude = amplitude;
grad.riseTime = riseTime;
grad.flatTime = flatTime;
grad.fallTime = fallTime;
grad.area = amplitude*(flatTime + riseTime/2 + fallTime/2);
grad.flatArea = amplitude*flatTime;
grad.delay = opt.delay;
grad.first = 0;
grad.last = 0;

end

function [amplitude, riseTime, flatTime, fallTime] = calcShortestParamsForArea(area,maxSlew,maxGrad,gradRasterTime)
    % find the shortest possible duration
    % first check if the area can be realized as a triangle
    % if not we calculate a trapezoid
    riseTime=ceil(sqrt(abs(area)/maxSlew)/gradRasterTime)*gradRasterTime;
    if riseTime < gradRasterTime % the "area" was probably 0 or almost 0 ...
        riseTime=gradRasterTime;
    end
    amplitude=area/riseTime;
    tEff=riseTime;
    if abs(amplitude)>maxGrad 
        tEff=ceil(abs(area)/maxGrad/gradRasterTime)*gradRasterTime;
        amplitude=area/tEff;
        riseTime=ceil(abs(amplitude)/maxSlew/gradRasterTime)*gradRasterTime;
        if(riseTime==0)
            riseTime=gradRasterTime;
        end
    end
    flatTime=tEff-riseTime;
    fallTime=riseTime;        
end
