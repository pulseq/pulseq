function readBinary(obj,filename)
%READBINARY Load sequence from binary file.
%   READBINARY(seqObj, filename) Read the given filename and load sequence
%   data stored in the binary version of the Pulseq open file format. The
%   binary format is described in the specficiation available at 
%   http://pulseq.github.io
%
%   Examples:
%   Load the sequence defined in gre.bin in sequences directory
%
%       readBinary(seqObj,'sequences/gre.bin')
%
% See also  writeBinary

binaryCodes = obj.getBinaryCodes();
fid=fopen(filename);
magicNum = fread(fid,8,'uchar');
assert(all(magicNum==binaryCodes.fileHeader(:)),'Not a binary file');
version_major = fread(fid,1,'int64');
version_minor = fread(fid,1,'int64');
version_revision = fread(fid,1,'int64');
assert(version_major==binaryCodes.version_major,'Unsupported version_major %d', version_major)
assert(version_minor==binaryCodes.version_minor,'Unsupported version_minor %d', version_minor)
assert(version_revision==binaryCodes.version_revision,'Unsupported version_revision %d', version_revision)

% set version
obj.version_major = version_major;
obj.version_minor = version_minor;
obj.version_revision = version_revision;

% Clear sequence data
obj.blockEvents={};
obj.definitions=containers.Map();
obj.gradLibrary=mr.EventLibrary();
obj.shapeLibrary=mr.EventLibrary();
obj.rfLibrary=mr.EventLibrary();
obj.adcLibrary=mr.EventLibrary();
obj.delayLibrary=mr.EventLibrary();

% Load data from file
while true
    section = int64(fread(fid,1,'int64'));
    if isempty(section)
        break
    end
    
    switch section
        case binaryCodes.section.definitions
            obj.definitions = readDefinitions(fid);
        case binaryCodes.section.blocks
            obj.blockEvents = readBlocks(fid);
        case binaryCodes.section.rf
            format = {'float64','int32','int32','float64','float64'};
            obj.rfLibrary = readEvents(fid,format,1);
        case binaryCodes.section.gradients
            format = {'float64','int32'};
            obj.gradLibrary = readEvents(fid,format,1,'g',obj.gradLibrary);
        case binaryCodes.section.trapezoids
            format = {'float64','int64','int64','int64'};
            obj.gradLibrary = readEvents(fid,format,[1 1e-6 1e-6 1e-6],'t',obj.gradLibrary);
        case binaryCodes.section.adc
            format = {'int64','int64','int64','float64','float64'};
            obj.adcLibrary = readEvents(fid,format,[1 1e-9 1e-6 1 1]);
        case binaryCodes.section.delays
            format = {'int64'};
            obj.delayLibrary = readEvents(fid,format,1e-6);
        case binaryCodes.section.shapes
            obj.shapeLibrary = readShapes(fid);
        otherwise
            error('Unknown section code: %s',dec2hex(section));
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
        
        def = containers.Map();
        numDefs = fread(fid,1,'int64');
        for i=1:numDefs
            c = fread(fid,1,'char');
            key=[];
            while c~=0 && length(key)<255
                key=[key c];
                c = fread(fid,1,'char');
            end
            key = char(key);
            type = fread(fid,1,'int8');
            count = fread(fid,1,'int8');
            switch type
                case 1
                    type = 'int64';
                case 2
                    type = 'float64';
                otherwise
                    error('Unknown definition type');
            end
            values = double(fread(fid,count,type));
            def(key) = values;
        end
    end

    function eventTable = readBlocks(fid)
        %readBlocks Read the [BLOCKS] section of a sequence file.
        %   library=readBlocks(fid) Read blocks from file identifier of an
        %   open MR sequence file and return the event table.
        
        numBlocks = double(fread(fid,1,'int64'));
        eventIds = double(fread(fid,6*numBlocks,'int32'));
        eventTable_tmp = reshape(eventIds,6,numBlocks)';
        eventTable = cell(1, numBlocks);
        for ii = 1:numBlocks
            eventTable{ii} = eventTable_tmp(ii,:);
        end
    end

    function eventLibrary = readEvents(fid,format,scale,type,eventLibrary)
        %readEvents Read an event section of a sequence file.
        %   library=readEvents(fid,fmt,scale) Read event data of given
        %   format and scale elements according to column vector scale.
        %
        %   library=readEvents(fid,fmt,scale,type) Attach the type string
        %   to elements of the library.
        %
        %   library=readEvents(...,library) Append new events to the given
        %   library.
        
        if nargin<3
            scale=1;
        end
        if nargin<5
            eventLibrary=mr.EventLibrary();
        end
        numEvents = fread(fid,1,'int64');
        for i=1:numEvents
            id = double(fread(fid,1,'int32'));
            data = zeros(1,length(format));
            for j=1:length(format)
                data(j) = double(fread(fid,1,format{j}));
            end
            data=scale.*data;
            if nargin<4
                eventLibrary.insert(id,data);
            else
                eventLibrary.insert(id,data,type);
            end
        end
    end

    function shapeLibrary = readShapes(fid)
        %readShapes Read the [SHAPES] section of a sequence file.
        %   library=readShapes(fid) Read shapes from file identifier of an
        %   open MR sequence file and return a library of shapes.
        
        shapeLibrary=mr.EventLibrary();
        numShapes = fread(fid,1,'int64');
        for i=1:numShapes
            id = double(fread(fid,1,'int32'));
            numUncompressed = double(fread(fid,1,'int64'));
            numCompressed = double(fread(fid,1,'int64'));
            data = double(fread(fid,numCompressed,'float64'))';
            shapeData = [numUncompressed data];
            shapeLibrary.insert(id,shapeData);
        end
    end

end