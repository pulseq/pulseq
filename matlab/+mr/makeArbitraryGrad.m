function grad=makeArbitraryGrad(channel,varargin)
%makeArbitraryGrad Create an gradient event with arbitrary waveform.
%   g=makeArbitraryGrad(channel, waveform) Create gradient on
%   the given channel with the specified waveform.
%
%   g=makeArbitraryGrad(channel,waveform,lims) Ensure the waveform 
%   satisfies the gradient hardware constraints.
%
%   See also  Sequence.addBlock 

persistent parser

if isempty(parser)
    validChannels = {'x','y','z'};
    parser = inputParser;
    parser.FunctionName = 'makeArbitraryGrad';
    parser.addRequired('channel',...
        @(x) any(validatestring(x,validChannels)));
    parser.addRequired('waveform');    
    parser.addOptional('system', [], @isstruct);
    parser.addParamValue('oversampling',false,@islogical);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    parser.addParamValue('delay',0,@isnumeric);
    parser.addParamValue('first',NaN,@isnumeric);
    parser.addParamValue('last',NaN,@isnumeric);
end
parse(parser,channel,varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

maxSlew=system.maxSlew;
maxGrad=system.maxGrad; % TODO: use this when no duration is supplied
if opt.maxGrad>0
    maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    maxSlew=opt.maxSlew;
end

g=opt.waveform(:);

if isfinite(opt.first)
    first = opt.first;
else
    warning('it will be compulsory to provide the first point of the gradient shape in the future releases; finding the first by extrapolation for now...');
    if opt.oversampling
        first = 2*g(1)-g(2); % extrapolate by 1 gradient raster
    else
        first = (3*g(1)-g(2))*0.5; % extrapolate by 1/2 gradient of the raster
    end
end

if isfinite(opt.last)
    last = opt.last;
else
    warning('it will be compulsory to provide the last point of the gradient shape in the future releases; finding the last by extrapolation for now...');
    if opt.oversampling
        last = g(end)*2-g(end-1); % extrapolate by 1 gradient raster
    else
        last = (g(end)*3-g(end-1))*0.5; % extrapolate by 1/2 gradient of the raster
    end
end

if opt.oversampling
    slew=[(first-g(1)); (g(2:end)-g(1:end-1)); (last-g(end))]./system.gradRasterTime*2;
else
    slew=[(first-g(1))*2; (g(2:end)-g(1:end-1)); (g(end)-last)*2]./system.gradRasterTime;
end
if ~isempty(slew) && max(abs(slew))>maxSlew
    error('Slew rate violation (%.0f%%)',max(abs(slew))/maxSlew*100);
end
if max(abs(g))>maxGrad
    error('Gradient amplitude violation (%.0f%%)',max(abs(g))/maxGrad*100);
end

grad.type = 'grad';
grad.channel = opt.channel;
grad.waveform = g;
grad.delay = opt.delay;
grad.area=sum(grad.waveform)*system.gradRasterTime; % QC: Take gradient raster time into account. 20230719
% true timing and aux shape data
if opt.oversampling
    if (mod(length(g),2)~=1)
        error('when oversampling is active the gradient shape vector must contain an odd number of samples');
    end
    grad.tt = (1:length(g))'*0.5*system.gradRasterTime;
    grad.shape_dur = (length(g)+1)*0.5*system.gradRasterTime;    
else
    grad.tt = ((1:length(g))'-0.5)*system.gradRasterTime;
    grad.shape_dur = length(g)*system.gradRasterTime;
end
grad.first = first;
grad.last = opt.last;

end
