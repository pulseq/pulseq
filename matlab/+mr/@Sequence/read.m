function read(obj,filename,varargin)
%READ Load sequence from file.
%   READ(seqObj, filename, ...) Read the given filename and load sequence
%   data into sequence object.
%
%   optional parwameter 'detectRFuse' can be given to let the function
%   infer the currently missing flags concerning the intended use of the RF
%   pulses (excitation, refocusing, etc). These are important for the
%   k-space trajectory calculation
%
%   Examples:
%   Load the sequence defined in gre.seq in my_sequences directory
%
%       read(seqObj,'my_sequences/gre.seq')
%
% See also  write

detectRFuse=false;
if ~isempty(varargin) && ~isempty(strfind(varargin{:},'detectRFuse'))
    detectRFuse=true;
end

fid = fopen(filename);

% Clear sequence data
%obj.blockEvents = [];
obj.blockEvents = {};
obj.definitions = containers.Map();
obj.gradLibrary = mr.EventLibrary();
obj.shapeLibrary = mr.EventLibrary();
obj.rfLibrary = mr.EventLibrary();
obj.adcLibrary = mr.EventLibrary();
obj.delayLibrary = mr.EventLibrary();

% Load data from file
while true
    section = skipComments(fid);
    if section == -1
        break
    end
    
    switch section
        case '[DEFINITIONS]'
            obj.definitions = readDefinitions(fid);
        case '[VERSION]'
%             obj.definitions = readDefinitions(fid); % MZ: I guess this should be OK
                                                      % SK: Nope. Overwrites definitions.
            % SK: This works:
            [version_major, ...
             version_minor, ...
             version_revision] = readVersion(fid);
            assert(version_major == obj.version_major, ...
                'Unsupported version_major %d', version_major)
            assert(version_minor == obj.version_minor, ...
                'Unsupported version_minor %d', version_minor)
            assert(version_revision == obj.version_revision, ...
                'Unsupported version_revision %d', version_revision)
            obj.version_major = version_major;
            obj.version_minor = version_minor;
            obj.version_revision = version_revision;
        case '[BLOCKS]'
            obj.blockEvents = readBlocks(fid);
        case '[RF]'
            obj.rfLibrary = readEvents(fid, [1 1 1 1e-6 1 1]);
        case '[GRADIENTS]'
            obj.gradLibrary = readEvents(fid, [1 1 1e-6], 'g' ,obj.gradLibrary);
        case '[TRAP]'
            obj.gradLibrary = readEvents(fid, [1 1e-6 1e-6 1e-6 1e-6], 't', obj.gradLibrary);
        case '[ADC]'
            obj.adcLibrary = readEvents(fid, [1 1e-9 1e-6 1 1]);
        case '[DELAYS]'
            obj.delayLibrary = readEvents(fid, 1e-6);
        case '[SHAPES]'
            obj.shapeLibrary = readShapes(fid);
        otherwise
            error('Unknown section code: %s', section);
    end
end
fclose(fid);

gradChannels={'gx','gy','gz'};

if detectRFuse
    % find the RF pulses, list flip angles
    % and work around the current (rev 1.2.0) Pulseq file format limitation
    % that the RF pulse use is not stored in the file
    for k=obj.rfLibrary.keys
        libData=obj.rfLibrary.data(k).array;
        rf=obj.rfFromLibData(libData);
        flipAngleDeg=abs(sum(rf.signal))*rf.t(1)*360; %we use rfex.t(1) in place of opt.system.rfRasterTime
        % fix libData
        if length(libData) < 9
            if flipAngleDeg <= 90.0
                libData(9) = 0; % or 1 ?
            else
                libData(9) = 2; % or 1 ?
            end
            obj.rfLibrary.data(k).array=libData;
        end
    end
end

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
        while ischar(line) && ~(isempty(line) || line(1) == '#')
            tok = textscan(line, '%s');
            def(tok{1}{1}) = str2double(tok{1}(2:end));
            line = fgetl(fid);
        end
    end    

    function [major, minor, revision] = readVersion(fid)
        %readVersion Read the [VERSION] section of a sequence file.
        %   defs=readVersion(fid) Read Pulseq version from file
        %   identifier of an open MR sequence file and return it
        
        major = [];
        minor = [];
        revision = [];
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1)=='#')
            tok = textscan(line,'%s');
            switch tok{1}{1}
                case 'major'
                    major = str2double(tok{1}(2:end));
                case 'minor'
                    minor = str2double(tok{1}(2:end));
                case 'revision'
                    revision = str2double(tok{1}(2:end));
            end
            line = fgetl(fid);
        end
    end

    function eventTable = readBlocks(fid)
        %readBlocks Read the [BLOCKS] section of a sequence file.
        %   library=readBlocks(fid) Read blocks from file identifier of an
        %   open MR sequence file and return the event table.
        
        eventTable = {};
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1) == '#')
            blockEvents = sscanf(line, '%f')';
            %eventTable = [eventTable; blockEvents(2:end)];
            eventTable{blockEvents(1)} = blockEvents(2:end);
            line = fgetl(fid);
        end
    end

    function eventLibrary = readEvents(fid, scale, type, eventLibrary)
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
        if nargin < 2
            scale = 1;
        end
        if nargin < 4
            eventLibrary = mr.EventLibrary();
        end
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1) == '#')
            data = sscanf(line,'%f')';
            id = data(1);
            data = scale.*data(2:end);
            if nargin < 3
                eventLibrary.insert(id, data);
            else
                eventLibrary.insert(id, data, type);
            end
            
            line=fgetl(fid);
        end
    end

    function shapeLibrary = readShapes(fid)
        %readShapes Read the [SHAPES] section of a sequence file.
        %   library=readShapes(fid) Read shapes from file identifier of an
        %   open MR sequence file and return a library of shapes.

        shapeLibrary = mr.EventLibrary();
        line = skipComments(fid);
        while ~(~ischar(line) || isempty(line) || ~strcmp(line(1:8), 'shape_id'))
            tok = textscan(line, '%s');
            id = str2double(tok{1}(2));
            line = skipComments(fid);
            tok = textscan(line, '%s');
            num_samples = str2double(tok{1}(2));
            data = [];
            line = skipComments(fid);   % first sample
            while ischar(line) && ~(isempty(line) || line(1) == '#')
                data = [data sscanf(line, '%f')];
                line = fgetl(fid);
            end
            line = skipComments(fid);
            data = [num_samples data];
            shapeLibrary.insert(id, data);
        end
    end

    function nextLine = skipComments(fid)
        %skipComments Read lines of skipping blank lines and comments.
        %   line=skipComments(fid) Read lines from valid file identifer and
        %   return the next non-comment line.
        
        line = fgetl(fid);
        while ischar(line) && (isempty(line) || line(1) == '#')
            line = fgetl(fid);
        end
        if ischar(line)
            nextLine = line;
        else
            nextLine = -1;
        end
    end
end