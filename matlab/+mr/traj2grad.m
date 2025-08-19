function [g sr]=traj2grad(k,varargin)
%traj2grad Convert a k-space trajectory to gradient waveform.
%   g=traj2grad(k) Convert k into gradient waveform g using finite
%   differences. The trajectory is in units of 1/m. The k-space points are
%   assumed to be sampled on the raster edges.
%   The size of k = [nChannel nTime].
%
%   g=traj2grad(k,'RasterTime',T) Calculate gradient waveforms assuming
%   the given raster time.
%
%   See also  Sequence.makeArbitraryGrad

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'traj2grad';
    parser.addRequired('k',@isnumeric);
    parser.addParamValue('first',[],@isnumeric);
    parser.addParamValue('firstGradStepHalfRaster',true,@islogical);
    parser.addParamValue('conservativeSlewEstimate',false,@islogical);    
    parser.addParamValue('system',[],@isstruct);
    parser.addParamValue('RasterTime',[],@isnumeric);
end
parse(parser,k,varargin{:});
opt = parser.Results;
if isempty(opt.system)
    opt.system=mr.opts();
end
if isempty(opt.RasterTime)
    opt.RasterTime=opt.system.gradRasterTime;
end
if isempty(opt.first)
    opt.first=zeros(size(k,1),1); % QC: if the first gradient point is not given, set it to zero. 2025.01.03
end

% Compute finite difference for gradients in Hz/m
%g=([k(:,2:end)-k(:,1:end-1) zeros(size(k,1),1)])/opt.RasterTime; % MZ: with zero-padding
g=[(k(:,2:end)-k(:,1:end-1))/opt.RasterTime]; % MZ: no zero-padding!

% Compute the slew rate (time derivative of the gradient)
sr0=(g-[opt.first g(:,1:end-1)])/opt.RasterTime;
if opt.firstGradStepHalfRaster
    sr0(:,1)=sr0(:,1)*2; % account for the half-step in the beginning of the shape
end

% now we think how to post-process the results 
% gradient is now sampled between the k-points (on raster cell centers)
% whilst the slew rate is between the gradient points, except of the first
% point, which relies on the opt.first value (and may be a bit off anyway,
% but this is the best estimate that we have)
sr=zeros(size(sr0));
sr(:,1)=sr0(:,1);
if (opt.conservativeSlewEstimate)
    if opt.firstGradStepHalfRaster
        sr(:,2)=sr0(:,2);
        sr(:,3:end)=max_abs(sr0(:,2:end-1),sr0(:,3:end));
    else
        sr(:,2:end)=max_abs(sr0(:,1:end-1),sr0(:,2:end));
    end
else
    if opt.firstGradStepHalfRaster
        sr(:,2)=sr0(:,2);
        sr(:,3:end)=0.5*(sr0(:,2:end-1)+sr0(:,3:end));
    else
        sr(:,2:end)=0.5*(sr0(:,1:end-1)+sr0(:,2:end));
    end
end

end

function out=max_abs(in1, in2)
    if size(in2)~=size(in2)
        error('arrays of incompatible sizes');
    end
    abs1gtoe=abs(in1)>=abs(in2);
    out=in1.*abs1gtoe+in2.*(~abs1gtoe);
end
