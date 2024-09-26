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
%   When a numerical parameter is passed amongst the events it is
%   interpreted as a predefined duration of the block. If the duration of
%   any of the events exceeds the desired duration an error will be thrown.
%   The predified block duration is optionally returned as the last
%   parameter in the output list
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
required_duration=[];

for i=2:length(varargin)
    if isempty(curr_align)
        error('invalid alignment spec');
    end
    if ischar(varargin{i})
        curr_align=find(strcmp(varargin{i},alignment_options));
        continue;
    end
    if isnumeric(varargin{i})
        if ~isempty(required_duration)
            error('More than one numeric parameter given to align()');
        end
        required_duration=varargin{i};
        continue;
    end
    iobjects=[iobjects i];
    alignments=[alignments curr_align];
end

objects={varargin{iobjects}};

dur=mr.calcDuration(objects);
if ~isempty(required_duration)
    if dur-required_duration>eps
        error('Required block duration is %g s but actuall block duration is %g s', required_duration, dur);
    end
    dur=required_duration;
end

% set new delays
for i=1:length(objects)
    switch alignments(i)
        case 1
            objects{i}.delay=0;
        case 2
            objects{i}.delay=(dur - mr.calcDuration(objects{i}))/2; % FIXME check how to handle the existing delay
        case 3
            ev=objects{i};
            ev_dur=mr.calcDuration(ev);
            %if isfield(ev,'ringdownTime')
            %    ev_dur=ev_dur+ev.ringdownTime;
            %end
            objects{i}.delay=dur - ev_dur + objects{i}.delay;
            if objects{i}.delay < 0
                error('aligh() attempts to set a negative delay, probably some RF pulses ignore rfRingdownTime');
            end
    end
end

if nargout==length(objects)
    varargout=objects;
elseif ~isempty(required_duration) && nargout==length(objects)+1
    varargout=objects;
    varargout{end+1}=required_duration;
elseif nargout==1    
    varargout={objects};
    if ~isempty(required_duration)
        varargout{end+1}=required_duration;
    end
elseif nargout<length(objects)
    warning('not all objects can be assigned to the output argiments; we recommend using ~ to discard output arguments explicitly.');
    nout=min(nargout,length(objects));
    varargout=objects(1:nout);
else
    error(['the number of output arguments (' num2str(nargout) ') exceeds the number of sequence objects (' num2str(length(objects)) ') passed to the function.']);
end

end