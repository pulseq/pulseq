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
    parser.addOptional('system', mr.opts(), @isstruct);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    parser.addParamValue('delay',0,@isnumeric);
    parser.addParamValue('first',NaN,@isnumeric);
    parser.addParamValue('last',NaN,@isnumeric);
end
parse(parser,channel,varargin{:});
opt = parser.Results;

maxSlew=opt.system.maxSlew;
maxGrad=opt.system.maxGrad; % TODO: use this when no duration is supplied
if opt.maxGrad>0
    maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    maxSlew=opt.maxSlew;
end

g=opt.waveform(:);
slew=(g(2:end)-g(1:end-1))./opt.system.gradRasterTime;
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
grad.area=sum(grad.waveform)*opt.system.gradRasterTime; % QC: Take gradient raster time into account. 20230719
% true timing and aux shape data
grad.tt = ((1:length(g))-0.5)*opt.system.gradRasterTime;
grad.shape_dur = length(g)*opt.system.gradRasterTime;

if isfinite(opt.first)
    grad.first = opt.first;
else
    grad.first = (3*g(1)-g(2))*0.5; % extrapolate by 1/2 gradient of the raster
end

if isfinite(opt.last)
    grad.last = opt.last;
else
    grad.last = (g(end)*3-g(end-1))*0.5; % extrapolate by 1/2 gradient of the raster
end

end
