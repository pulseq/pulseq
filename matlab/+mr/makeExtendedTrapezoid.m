function grad = makeExtendedTrapezoid(channel, varargin)
%makeExtendedTrapezoid Create an extended trapezoid gradient event.
%
%   g = makeExtendedTrapezoid(channel, lims, ...
%                             'times', times, ...
%                             'amplitudes', amplitudes) 
%   Create a gradient by specifying a set of points (amplitudes) at
%   specified time points(times) at a given channel with given system 
%   limits. This function returns an arbitrary gradient object. 
%
%   See also  Sequence.addBlock  mr.opts  makeTrapezoid
%
%   Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>

persistent parser

if isempty(parser)
    validChannels = {'x', 'y', 'z'};
    parser = inputParser;
    parser.FunctionName = 'makeExtendedTrapezoid';
    parser.addRequired('channel', ...
        @(x) any(validatestring(x, validChannels)));
    parser.addParamValue('times', 0, @isnumeric);
    parser.addParamValue('amplitudes', 0, @isnumeric);
    parser.addOptional('system', mr.opts(), @isstruct);
    parser.addParamValue('maxGrad', 0, @isnumeric);
    parser.addParamValue('maxSlew', 0, @isnumeric);
    parser.addParamValue('skip_check', false);
    parser.addParamValue('convert2arbitrary', false);
    
end
parse(parser,channel,varargin{:});
opt = parser.Results;

if any(size(opt.times) ~= size(opt.amplitudes))
    error('Times and amplitudes must have the same length.');
end

if all(opt.times == 0)
    error('At least one of the given times must be non-zero.');
end

if any(diff(opt.times)<=0)
    error('Times must be in ascending order and all times must be distinct.');
end

if abs(round(opt.times(end)/opt.system.gradRasterTime)*opt.system.gradRasterTime-opt.times(end))>eps
    error('The last time point must be on a gradient raster.');
end

%if all(opt.amplitudes == 0)
%    error('At least one of the given amplitudes must be non-zero.');
%end

if opt.skip_check == false && opt.times(1) > 0 && opt.amplitudes(1) ~= 0
    error('If first amplitude of a gradient is nonzero, it must connect to previous block!');
end

maxSlew = opt.system.maxSlew;
maxGrad = opt.system.maxGrad;
if opt.maxGrad > 0
    maxGrad = opt.maxGrad;
end
if opt.maxSlew > 0
    maxSlew = opt.maxSlew;
end

if (opt.convert2arbitrary)
    % represent the extended trapezoid on the regularly sampled time grid
    waveform = mr.pts2waveform(opt.times, opt.amplitudes, opt.system.gradRasterTime);
    grad = mr.makeArbitraryGrad(channel, waveform, opt.system, ...
                                'maxSlew', maxSlew,...
                                'maxGrad', maxGrad,...
                                'delay', opt.times(1));
else
    % keep the original possibly irregular sampling
    if any(abs(round(opt.times/opt.system.gradRasterTime)*opt.system.gradRasterTime-opt.times)>eps)
        error('All time points must be on a gradient raster or "convert2arbitrary" option must be used.');
    end
    %
    grad.type = 'grad';
    grad.channel = opt.channel;
    grad.waveform = opt.amplitudes;
    grad.delay = round(opt.times(1)/opt.system.gradRasterTime)*opt.system.gradRasterTime;
    grad.tt = opt.times - grad.delay;
    grad.shape_dur = round(opt.times(end)/opt.system.gradRasterTime)*opt.system.gradRasterTime;
end

% MZ: although makeArbitraryGrad sets the .first and .last for extended 
% trapezoids we can do it better
grad.first=opt.amplitudes(1);
grad.last=opt.amplitudes(end);

end



% figure; plot(waveform)
% waveform = waveform(1:end-2);
% opt.times = round(opt.times/opt.system.gradRasterTime)*opt.system.gradRasterTime; % round onto grid
% times_diff = diff(opt.times);
% amplitudes_diff = diff(opt.amplitudes);
% waveform = [];
% for ii = 1:length(opt.times)-1
%     % SK: there are no new points after the end, therefore we dont need to
%     % handle the overlap situation.
%     if ii == length(opt.times)-1
%         crop = 0;
%     else
%         crop = opt.system.gradRasterTime;
%     end
%     y = amplitudes_diff(ii)/times_diff(ii)*...
%         (0:opt.system.gradRasterTime:(opt.times(ii+1)-opt.times(ii)-crop))...
%         + opt.amplitudes(ii);
%     waveform = [waveform y(1:end)];
% end
