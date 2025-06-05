function [varargout] = rotate3D(rotation, varargin)
%rotate3D rotate all objects (gradients) in the block by a rotation matrix
%
%   [...] = rotate(rotMat, obj <, obj> ...);
%
%   Rotates the corresponding gradinet object(s) by the provided rotation
%   represented either as a 3x3 rotation matrix or a unit quaternion with
%   the scalar component at the first position; non-gradient objects are
%   not affected.   
%
%   Optional parameter list may include the keyword 'system' followed by a
%   system limits struct. The system can only be provided in the beginning
%   or at the ent of the list of optional parameters. 
%
%   Returns either a cell-array of objects if one return parameter is
%   provided or an explicit list of objects if multiple parameters are
%   given. Can be used directly as a parameter of seq.addBlock().
%
%   See also  mr.rotate, Sequence.addBlock
%
%   Maxim Zaitsev <maxim.zaitsev@uniklinik-freiburg.de>

% parse out the optional parameter 'system', which can only be at the beginning
% or in the end of the optional parameters
system=[];
if ischar(varargin{1}) && strcmp(lower(varargin{1}),'system')
    if ~isstruct(varargin{2}) || ~isfield(varargin{2},'gradRasterTime')
        error('Error parsing input parameters, keyword ''system'' is not followed by a valid system struct');
    end
    system=varargin{2};
    varargin=varargin(3:end);
elseif length(varargin)>1 && ischar(varargin{end-1}) && strcmp(lower(varargin{end-1}),'system')
    if ~isstruct(varargin{end}) || ~isfield(varargin{end},'gradRasterTime')
        error('Error parsing input parameters, keyword ''system'' is not followed by a valid system struct');
    end
    system=varargin{end};
    varargin=varargin(1:end-2);
end

% detect rotation matrix or quaternion formats
if size(rotation)==[3 3]
    rotMat=rotation;
elseif length(rotation)==4
    rotMat=mr.aux.quat.toRotMat(rotation);
else
    error('The parameter ''rotation'' must either bi a 3x3 matrix or a quaternion');
end

% make this function accept ready-made blocks
if isstruct(varargin{1}) && isfield(varargin{1}, 'rf') 
    varargin=mr.block2events(varargin);    
end
% we need this to allow for nested mr.rotate() calls
if 1==length(varargin) && iscell(varargin{1})
    va=varargin{1}; 
else
    va=varargin;
end

% first create indexes of the objects to be bypassed or rotated
ibypass=[];
grads3_in=cell(1,3);
axes={'x', 'y', 'z'};

for i=1:length(va)
    event = va{i};
    if isempty(event)
        continue;
    end
    if isnumeric(event) || ...
            ((~strcmp(event.type,'grad') && ...
             ~strcmp(event.type,'trap')))
        ibypass=[ibypass i];
    else
        iAxis=find(strcmp(event.channel,axes));
        if ~isempty(grads3_in{iAxis})
            error('More than one gradient on the same axis %s provided', event.channel);
        end
        grads3_in{iAxis}=event;
    end
end

max_mag=0; % measure of the relevant amplitude
for i=1:3
    if ~isempty(grads3_in{i})
        max_mag=max(max_mag, getGradAbsMag(grads3_in{i}));
    end
end
fthresh=1e-6;
thresh=fthresh*max_mag;

grads_out={};
for j=1:3
    grad_out_curr=[];
    for i=1:3
        if isempty(grads3_in{i}) || ...
                abs(rotMat(j,i))<fthresh
            continue;
        end
        g=mr.scaleGrad(grads3_in{i},rotMat(j,i));
        g.channel=axes{j};
        if isempty(grad_out_curr)
            grad_out_curr=g;
        else
            if isempty(system)
                grad_out_curr=mr.addGradients({grad_out_curr,g});
            else
                grad_out_curr=mr.addGradients({grad_out_curr,g},system);
            end
        end
    end
    % only output non-zero-amplitude gradients
    if ~isempty(grad_out_curr) && getGradAbsMag(grad_out_curr) >= thresh
        grads_out{end+1}=grad_out_curr;
    end
end

% export
bypass=va(ibypass);
out={bypass{:},grads_out{:}}; 

nout = nargout;
varargout = cell(1,nout);
if nout==1 
    varargout{1}=out;
else
    nr=min(nout,length(out));
    if nout<length(out) 
        warning('insufficient number of return parameters, some rotated gradient components might go lost');
    end
    for k=1:nr
        varargout{k}=out{k};
    end
end

end


function [out] = getGradAbsMag(grad)
    if strcmp(grad.type,'trap')
        out=abs(grad.amplitude);
    else
        out=max(abs(grad.waveform));
    end
end
