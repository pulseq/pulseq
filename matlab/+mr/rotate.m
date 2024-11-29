function [varargout] = rotate(raxis, angle, varargin)
%align set alignment of the objects in the block
%
%   [...] = rotate(axis, angle, obj <, obj> ...);
%
%   Rotates the corresponding gradinet object(s) about the given axis by 
%   the specified amount. Gradients parallel to the rotation axis and 
%   non-gradient objects are not affected. 
%   Possible rotation axes are 'x', 'y' or 'z'.
%
%   Optional parameter list may include the keyword 'system' followed by a
%   system limits struct. The system can only be provided in the beginning or 
%   at the ent of the list of optional parameters.
%
%   Returns either a cell-array of objects if one return parameter is
%   provided or an explicit list of objects if multiple parameters are
%   given. Can be used directly as a parameter of seq.addBlock().
%
%   See also  mr.rotate3D, Sequence.addBlock
%
%   Maxim Zaitsev <maxim.zaitsev@uniklinik-freiburg.de>

axes={'x', 'y', 'z'};

% cycle through the objects and rotate gradients non-parallel to the
% given rotation axis. Rotated gradients assigned to the same axis are then
% added together.

% first create indexes of the objects to be bypassed or rotated
irotate1=[];
irotate2=[];
ibypass=[];
axes2rot=axes(~strcmp(axes,raxis));
if length(axes2rot)~=2
    error('incorrect axis specification');
end
if strcmp('y',raxis)
    axes2rot=axes2rot(end:-1:1); % we need to reverse the list to preserve the correct handiness of the rotation matrix
end

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

for i=1:length(va)
    par = va{i};
    if isempty(par)
        continue;
    end
    if isnumeric(par) || ...
            ((~strcmp(par.type,'grad') && ...
             ~strcmp(par.type,'trap')) || ...
            strcmp(par.channel, raxis))%['g' axis]
        ibypass=[ibypass i];
    else
        if strcmp(par.channel, axes2rot(1)) %['g' axes2rot(1)]
            irotate1=[irotate1 i];
        else
            if (strcmp(par.channel, axes2rot(2))) %['g' axes2rot(2)]
                irotate2=[irotate2 i];
            else
                ibypass=[ibypass i]; % should never happen
            end
        end
    end
end

% now every gradient to be rotated generates two new gradients, one on the
% original axis and one on the other from the axes2rot list

rotated1=cell(1,length(irotate1)+length(irotate2));
rotated2=cell(1,length(irotate1)+length(irotate2));
max_mag=0; % measure of the relevant amplitude
for i=1:length(irotate1)
    g=va{irotate1(i)};
    max_mag=max(max_mag, getGradAbsMag(g));
    rotated1{i}=mr.scaleGrad(g,cos(angle));
    g=mr.scaleGrad(g,sin(angle));
    g.channel=axes2rot{2};
    rotated2{i}=g;
end
o=length(irotate1);
for i=1:length(irotate2)
    g=va{irotate2(i)};
    max_mag=max(max_mag, getGradAbsMag(g));
    rotated2{i+o}=mr.scaleGrad(g,cos(angle));
    g=mr.scaleGrad(g,-sin(angle));
    g.channel=axes2rot{1};
    rotated1{i+o}=g;
end

% eliminate zero-amplitude gradients
thresh=1e-6*max_mag;
for i=length(rotated1):-1:1
    if getGradAbsMag(rotated1{i})<thresh
        rotated1(i)=[];
    end
end
for i=length(rotated2):-1:1
    if getGradAbsMag(rotated2{i})<thresh
        rotated2(i)=[];
    end
end

g=cell(1,2);
% now we add gradients on the corresponding axis together
if (length(rotated1)>1)
    if isempty(system)
        g{1}=mr.addGradients(rotated1);
    else
        g{1}=mr.addGradients(rotated1,system);
    end
else
    if (~isempty(rotated1))
        g{1}=rotated1{1};
    end
end

if (length(rotated2)>1)
    if isempty(system)
        g{2}=mr.addGradients(rotated2);
    else
        g{2}=mr.addGradients(rotated2,system);
    end
else
    if (~isempty(rotated2))
        g{2}=rotated2{1};
    end
end

% eliminate zero-amplitude gradients
for i=length(g):-1:1
    if isempty(g{i}) || getGradAbsMag(g{i})<thresh
        g(i)=[];
    end
end

% export
bypass=va(ibypass);
out={bypass{:},g{:}}; 

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
