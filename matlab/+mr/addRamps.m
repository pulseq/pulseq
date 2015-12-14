function varargout=addRamps(k,varargin)
%addRamps Add segment to kspace trajectory to ramp to and from the given 
%   trajectory.
%   kout=addRamps(k) Add a segment to k so kout travels from 0 to k(1) and
%   a segment so kout goes from k(end) back to 0 without violating the
%   gradient and slew constraints. 
%
%   [kx,ky,...]=addRamps({kx,ky,...}) Add segments of the same length
%   for each trajectory in the cell array.
%
%   [...,rf]=addRamps(...,'rf',rf) Add a segment of zeros over the ramp
%   times to an RF shape.
%
%   See also  Sequence.makeArbitraryGrad

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'addRamps';
    parser.addRequired('k',@(x)(isnumeric(x)||iscell(x)));
    parser.addOptional('system',mr.opts(),@isstruct);
    parser.addParamValue('rf',[],@isnumeric);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    
end
parse(parser,k,varargin{:});
opt = parser.Results;

if opt.maxGrad>0
    opt.system.maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    opt.system.maxSlew=opt.maxSlew;
end

if iscell(opt.k)
    k=cell2mat(opt.k(:));
else
    k=opt.k;
end

nChannels=size(k,1);
k=[k; zeros(3-nChannels,size(k,2))];    % Pad out with zeros if needed

[kUp, ok1]   = mr.calcRamp(zeros(3,2),k(:,1:2),opt.system);
[kDown, ok2] = mr.calcRamp(k(:,end-1:end),zeros(3,2),opt.system);
assert(ok1 & ok2,'Failed to calculate gradient ramps');

kUp = [zeros(3,2), kUp];            % Add start and end points to ramps
kDown = [kDown, zeros(3,1)];

k = [kUp, k, kDown];                % Add ramps to trajectory

if isnumeric(opt.k)
    varargout{1} = k(1:nChannels,:);
else
    for i=1:nChannels
        varargout{i}=k(i,:);
    end
end
if ~isempty(opt.rf)
    varargout{end+1}=[zeros(1,size(kUp,2)*10), opt.rf, zeros(1,size(kDown,2)*10)];
end


end