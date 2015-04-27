function read(obj,filename)
%READ Load sequence from file.
%   READ(seqObj, filename) Read the given filename and load sequence
%   data into sequence object.
%
%   Examples:
%   Load the sequence defined in gre.txt in my_sequences directory
%
%       read(seqObj,'my_sequences/gre.txt')
%
% See also  write

fid=fopen(filename);
while true
    section = skipComments(fid);
    if section==-1
        break
    end
    
    switch section
        case '[DEFINITIONS]'
            obj.definitions = readDefinitions(fid);
        case '[BLOCKS]'
            obj.blockEvents = readBlocks(fid);
        case '[RF]'
            obj.rfLibrary = readEvents(fid,1);
        case '[GRAD]'
            obj.gradLibrary = readEvents(fid,1,'grad',obj.gradLibrary);
        case '[TRAP]'
            obj.gradLibrary = readEvents(fid,[1 1e-6 1e-6 1e-6],'trap',obj.gradLibrary);
        case '[ADC]'
            obj.adcLibrary = readEvents(fid,[1 1e-9 1e-6 1 1]);
        case '[DELAYS]'
            obj.delayLibrary = readEvents(fid,1e-6);
        case '[SHAPES]'
            obj.shapeLibrary = readShapes(fid);
    end
end
fclose(fid);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Helper functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function def = readDefinitions(fid)
        %readDefinitions Read the [DEFINITIONS] section of a sequence file.
        %   defs=readDefinitions(fid) Read user definitions from file
        %   identifier of an open MR sequence file and return a map of
        %   key/value entries.
        
        def = containers.Map;
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1)=='#')
            tok=textscan(line,'%s');
            def(tok{1}{1})=str2double(tok{1}(2:end));
            line=fgetl(fid);
        end
    end

    function eventTable = readBlocks(fid)
        %readBlocks Read the [BLOCKS] section of a sequence file.
        %   library=readBlocks(fid) Read blocks from file identifier of an
        %   open MR sequence file and return the event table.
        
        eventTable=[];
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1)=='#')
            blockEvents=sscanf(line,'%f')';
            eventTable=[eventTable; blockEvents(2:end)];
            line=fgetl(fid);
        end
    end

    function eventLibrary = readEvents(fid,scale,type,eventLibrary)
        %readEvents Read an event section of a sequence file.
        %   library=readEvents(fid) Read event data from file identifier of
        %   an open MR sequence file and return a library of events.
        %
        %   library=readEvents(fid,scale) Read event data and scale
        %   elements according to column vector scale.
        %
        %   library=readEvents(fid,scale,type) Attach the type string to
        %   elements of the library.
        %
        %   library=readEvents(...,library) Append new events to the given
        %   library.
        if nargin<2
            scale=1;
        end
        if nargin<4
            eventLibrary=containers.Map('KeyType','double','ValueType','any');
        end
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1)=='#')
            data=sscanf(line,'%f')';
            id=data(1);
            event.data=scale.*data(2:end);
            if nargin>2
                event.type=type;
            end
            eventLibrary(id)=event;
            line=fgetl(fid);
        end
    end

    function shapeLibrary = readShapes(fid)
        %readShapes Read the [SHAPES] section of a sequence file.
        %   library=readShapes(fid) Read shapes from file identifier of an
        %   open MR sequence file and return a library of shapes.

        shapeLibrary=containers.Map('KeyType','double','ValueType','any');
        line = skipComments(fid);
        while ~(~ischar(line) || isempty(line) || ~strcmp(line(1:8),'shape_id'))
            tok=textscan(line,'%s');
            id = str2double(tok{1}(2));
            line = skipComments(fid);
            tok=textscan(line,'%s');
            num_samples = str2double(tok{1}(2));
            data=[];
            line = skipComments(fid);   % first sample
            while ischar(line) && ~(isempty(line) || line(1)=='#')
                data=[data sscanf(line,'%f')];
                line=fgetl(fid);
            end
            line = skipComments(fid);
            shape.num_samples = num_samples;
            shape.data = data;
            shapeLibrary(id) = shape;
        end
    end

    function nextLine = skipComments(fid)
        %skipComments Read lines of skipping blank lines and comments.
        %   line=skipComments(fid) Read lines from valid file identifer and
        %   return the next non-comment line.
        
        line = fgetl(fid);
        while ischar(line) && (isempty(line) || line(1)=='#')
            line = fgetl(fid);
        end
        if ischar(line)
            nextLine=line;
        else
            nextLine=-1;
        end
    end
end