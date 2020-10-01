function [varargout] = align(varargin)
%align set alignment of the objects in the block
%
%   align(align_spec, obj <, obj> <, align_spec, obj> ...);
%
%   sets delays of the objects within the block to achieve the desired
%   alignment. 
%   All previously configured delays within objects are taken into account 
%   during calculating of the block duration but then reset according to 
%   the selected alignment.
%   Possible values for align_spec are 'left', 'center', 'right'
%   WARNING: 'center' may break graient raster alignment
%
%   See also  Sequence.addBlock
%
%   Maxim Zaitsev <maxim.zaitsev@uniklinik-freiburg.de>

alignment_options={'left', 'center', 'right'};
% parse parameters
if ~ischar(varargin{1})
    error('first parameter must be a string');
end
curr_align=find(strcmp(varargin{1},alignment_options));
iobjects=[];
alignments=[];
for i=2:length(varargin)
    if isempty(curr_align)
        error('invalid alignment spec');
    end
    if ischar(varargin{i})
        curr_align=find(strcmp(varargin{i},alignment_options));
        continue;
    end
    iobjects=[iobjects i];
    alignments=[alignments curr_align];
end

objects={varargin{iobjects}};

dur=mr.calcDuration(objects);

% set new delays
for i=1:length(objects)
    switch alignments(i)
        case 1
            objects{i}.delay=0;
        case 2
            objects{i}.delay=(dur - mr.calcDuration(objects{i}))/2; % FIXME check how to handle the existing delay
        case 3
            objects{i}.delay=dur - mr.calcDuration(objects{i}) + objects{i}.delay;
    end
end

varargout=objects;

end