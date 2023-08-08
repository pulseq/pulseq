function [g sr]=traj2grad(k,varargin)
%traj2grad Convert a k-space trajectory to gradient waveform.
%   g=traj2grad(k) Convert k into gradient waveform g using finite
%   differences. The trajectory is in units of 1/m. 
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
    system=mr.opts();
    parser.addParamValue('RasterTime',system.gradRasterTime,@isnumeric);
end
parse(parser,k,varargin{:});
opt = parser.Results;

% Compute finite difference for gradients in Hz/m
%g=([k(:,2:end)-k(:,1:end-1) zeros(size(k,1),1)])/opt.RasterTime; % MZ: with zero-padding
g=(k(:,2:end)-k(:,1:end-1))/opt.RasterTime; % MZ: no zero-padding!

% Compute the slew rate (time derivative of the gradient)
sr0=(g(:,2:end)-g(:,1:end-1))/opt.RasterTime;

% now we think how to post-process the results 
% gradient is now sampled between the k-points
% whilst the slew rate is between the gradient points
sr=zeros(size(sr0,1),size(sr0,2)+1);
sr(:,1)=sr0(:,1);
sr(:,2:end-1)=0.5*(sr0(:,1:end-1)+sr0(:,2:end));
sr(:,end)=sr0(:,end);

end