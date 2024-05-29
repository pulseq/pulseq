function grad=makeTrapezoid(channel, varargin)
%makeTrapezoid Create a trapezoid gradient event.
%   g=makeTrapezoid(channel, ...) Create trapezoid gradient on
%   the given channel.
%
%   g=makeTrapezoid(channel,lims,...) Create trapezoid with the specificed 
%   gradient limits (e.g. amplitude, slew).
%
%   g=makeTrapezoid(...,'Duration',d,'Area',a) Create a
%   trapezoid gradient with given duration (s) and total area (1/m)
%   including ramps.
%
%   g=makeTrapezoid(...,'FlatTime',d,'FlatArea',a) Create a
%   trapezoid gradient with given flat-top time and flat-top
%   area not including ramps.
%
%   g=makeTrapezoid(...,'Amplitude',a) Create a trapezoid gradient with
%   given amplitude (Hz/m).
%
%   See also  Sequence.addBlock  mr.opts

persistent parser

if isempty(parser)
    validChannels = {'x','y','z'};
    parser = inputParser;
    parser.FunctionName = 'makeTrapezoid';
    parser.addRequired('channel',...
        @(x) any(validatestring(x,validChannels)));
    parser.addOptional('system',mr.opts(),@isstruct);
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

maxSlew=opt.system.maxSlew;
%riseTime=opt.system.riseTime;
maxGrad=opt.system.maxGrad;
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
        riseTime = ceil(riseTime/opt.system.gradRasterTime)*opt.system.gradRasterTime;
        if riseTime==0
            riseTime=opt.system.gradRasterTime;
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
                [~, t1, t2, t3]=calcShortestParamsForArea(opt.area,maxSlew,maxGrad,opt.system.gradRasterTime);
                error('makeTrapezoid:invalidDuration',['Requested area is too large for this gradient. Minimum required duration for this area (accounting for the gradient raster time) is ' num2str((t1+t2+t3)*1e6) 'us']);
            end
            amplitude = ( opt.duration - sqrt(opt.duration^2 - 4*abs(opt.area)*dC) )/(2*dC);
        else
            if isempty(fallTime)
                fallTime = riseTime;
            end    
            amplitude = opt.area/(opt.duration-0.5*riseTime-0.5*fallTime);
            possible = opt.duration>(riseTime+fallTime) & abs(amplitude)<maxGrad;
            assert(possible,['Requested area is too large for this gradient duration. Probably amplitude is violated (' num2str(round(abs(amplitude)/maxGrad*100)) '%)']);    
        end    
    end
    if isempty(riseTime)
        riseTime = ceil(abs(amplitude)/maxSlew/opt.system.gradRasterTime)*opt.system.gradRasterTime;
        if(riseTime==0)
            riseTime=opt.system.gradRasterTime;
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
        [amplitude, riseTime, flatTime, fallTime]=calcShortestParamsForArea(opt.area,maxSlew,maxGrad,opt.system.gradRasterTime);
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
        [~, t1, t2, t3]=calcShortestParamsForArea(opt.area,maxSlew,maxGrad,opt.system.gradRasterTime);
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
