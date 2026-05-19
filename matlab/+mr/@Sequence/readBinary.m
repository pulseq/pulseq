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
obj.blockDurations=[];
obj.definitions=containers.Map();
obj.gradLibrary=mr.EventLibrary();
obj.shapeLibrary=mr.EventLibrary();
obj.rfLibrary=mr.EventLibrary();
obj.adcLibrary=mr.EventLibrary();
obj.trigLibrary=mr.EventLibrary();
obj.labelsetLibrary=mr.EventLibrary();
obj.labelincLibrary=mr.EventLibrary();
obj.extensionLibrary=mr.EventLibrary();
obj.rfShimLibrary=mr.EventLibrary();
obj.softDelayLibrary=mr.EventLibrary();
obj.softDelayHints1=containers.Map();
obj.softDelayHints2={};
obj.rotationLibrary=mr.EventLibrary();
obj.extensionStringIDs={};
obj.extensionNumericIDs=[];
obj.signatureType='';
obj.signatureFile='';
obj.signatureValue='';

% Load data from file
while true
    section = int64(fread(fid,1,'int64'));
    if isempty(section)
        break
    end

    switch section
        case binaryCodes.section.definitions
            obj.definitions = readDefinitions(fid);
            v=obj.getDefinition('GradientRasterTime');
            if ~isempty(v), obj.gradRasterTime=v; end
            v=obj.getDefinition('RadiofrequencyRasterTime');
            if ~isempty(v), obj.rfRasterTime=v; end
            v=obj.getDefinition('AdcRasterTime');
            if ~isempty(v), obj.adcRasterTime=v; end
            v=obj.getDefinition('BlockDurationRaster');
            if ~isempty(v), obj.blockDurationRaster=v; end

        case binaryCodes.section.blocks
            [obj.blockEvents, obj.blockDurations] = readBlocks(fid, obj.blockDurationRaster);

        case binaryCodes.section.rf
            % array: [amp mag_id phase_id time_shape_id center delay freqPPM phasePPM freq phase]
            % type (use) stored separately as a char
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id    = double(fread(fid,1,'int32'));
                amp   = double(fread(fid,1,'float64'));
                ids   = double(fread(fid,3,'int32'))';        % mag_id, phase_id, time_shape_id
                ctr   = double(fread(fid,1,'int64')) * 1e-12; % center (ps -> s)
                dly   = double(fread(fid,1,'int64')) * 1e-12; % delay  (ps -> s)
                fpp   = double(fread(fid,4,'float64'))';      % freqPPM, phasePPM, freq, phase
                use   = char(fread(fid,1,'char'));
                obj.rfLibrary.insert(id, [amp ids(1) ids(2) ids(3) ctr dly fpp(1) fpp(2) fpp(3) fpp(4)], use);
            end

        case binaryCodes.section.gradients
            % array: [amp first last amp_shape_id time_shape_id delay]
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id   = double(fread(fid,1,'int32'));
                afl  = double(fread(fid,3,'float64'))';       % amp, first, last
                ids2 = double(fread(fid,2,'int32'))';         % amp_shape_id, time_shape_id
                dly  = double(fread(fid,1,'int64')) * 1e-12;  % delay (ps -> s)
                obj.gradLibrary.insert(id, [afl ids2 dly], 'g');
            end

        case binaryCodes.section.trapezoids
            % array: [amp rise flat fall delay]
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id   = double(fread(fid,1,'int32'));
                amp  = double(fread(fid,1,'float64'));
                t4   = double(fread(fid,4,'int64'))' * 1e-12;  % rise,flat,fall,delay (ps->s)
                obj.gradLibrary.insert(id, [amp t4], 't');
            end

        case binaryCodes.section.adc
            % array: [num dwell delay freqPPM phasePPM freq phase phase_id]
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id     = double(fread(fid,1,'int32'));
                num    = double(fread(fid,1,'int64'));
                dwell  = double(fread(fid,1,'int64')) * 1e-12;  % ps -> s
                delay  = double(fread(fid,1,'int64')) * 1e-12;  % ps -> s
                f4     = double(fread(fid,4,'float64'))';       % freqPPM, phasePPM, freq, phase
                phid   = double(fread(fid,1,'int32'));
                obj.adcLibrary.insert(id, [num dwell delay f4(1) f4(2) f4(3) f4(4) phid]);
            end

        case binaryCodes.section.delays
            readLegacyDelays(fid);

        case binaryCodes.section.shapes
            obj.shapeLibrary = readShapes(fid);

        case binaryCodes.section.extensions
            % array per entry: [type ref next_id]
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id   = double(fread(fid,1,'int32'));
                data = double(fread(fid,3,'int32'))';
                obj.extensionLibrary.insert(id, data);
            end

        case binaryCodes.section.triggers
            % type id(i32) then events: id(i32) type(i32) channel(i32) delay(i32,us) duration(i32,us)
            ext_id = double(fread(fid,1,'int32'));
            obj.setExtensionStringAndID('TRIGGERS', ext_id);
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id   = double(fread(fid,1,'int32'));
                tc   = double(fread(fid,2,'int32'))';          % type, channel
                dd   = double(fread(fid,2,'int64'))' * 1e-12;  % delay, duration (ps->s)
                obj.trigLibrary.insert(id, [tc dd]);
            end

        case binaryCodes.section.labelset
            ext_id = double(fread(fid,1,'int32'));
            obj.setExtensionStringAndID('LABELSET', ext_id);
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id   = double(fread(fid,1,'int32'));
                data = double(fread(fid,2,'int32'))';         % value, label_index
                obj.labelsetLibrary.insert(id, data);
            end

        case binaryCodes.section.labelinc
            ext_id = double(fread(fid,1,'int32'));
            obj.setExtensionStringAndID('LABELINC', ext_id);
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id   = double(fread(fid,1,'int32'));
                data = double(fread(fid,2,'int32'))';         % value, label_index
                obj.labelincLibrary.insert(id, data);
            end

        case binaryCodes.section.softdelays
            ext_id = double(fread(fid,1,'int32'));
            obj.setExtensionStringAndID('DELAYS', ext_id);
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id      = double(fread(fid,1,'int32'));
                num     = double(fread(fid,1,'int32'));
                offset  = double(fread(fid,1,'int64')) * 1e-12;   % ps -> s
                factor  = double(fread(fid,1,'float64'));
                hlen    = double(fread(fid,1,'int32'));
                hint    = char(fread(fid,hlen,'char')');
                % register hint string and get its index
                if obj.softDelayHints1.isKey(hint)
                    hint_idx = obj.softDelayHints1(hint);
                else
                    hint_idx = length(obj.softDelayHints2) + 1;
                    obj.softDelayHints1(hint) = hint_idx;
                    obj.softDelayHints2{hint_idx} = hint;
                end
                obj.softDelayLibrary.insert(id, [num offset factor hint_idx]);
            end

        case binaryCodes.section.rfshims
            ext_id = double(fread(fid,1,'int32'));
            obj.setExtensionStringAndID('RF_SHIMS', ext_id);
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id       = double(fread(fid,1,'int32'));
                num_chan  = double(fread(fid,1,'int32'));
                chan_data = double(fread(fid,2*num_chan,'float64'))';
                obj.rfShimLibrary.insert(id, chan_data);
            end

        case binaryCodes.section.rotations
            ext_id = double(fread(fid,1,'int32'));
            obj.setExtensionStringAndID('ROTATIONS', ext_id);
            numEvents = double(fread(fid,1,'int64'));
            for i=1:numEvents
                id   = double(fread(fid,1,'int32'));
                quat = double(fread(fid,4,'float64'))';
                obj.rotationLibrary.insert(id, mr.aux.quat.normalize(quat));
            end

        case binaryCodes.section.signature
            type_len = double(fread(fid,1,'int32'));
            sig_type = char(fread(fid,type_len,'char')');
            hash_len = double(fread(fid,1,'int32'));
            hash_raw = uint8(fread(fid,hash_len,'uint8'));
            fread(fid,1,'int64'); % original file length prior to signature append

            obj.signatureType = sig_type;
            obj.signatureFile = 'bin';
            if isempty(hash_raw)
                obj.signatureValue = '';
            else
                obj.signatureValue = lower(reshape(dec2hex(hash_raw,2)',1,[]));
            end

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
        for iDef=1:numDefs
            c = fread(fid,1,'char');
            key=[];
            while c~=0 && length(key)<255
                key=[key c];
                c = fread(fid,1,'char');
            end
            key = char(key);
            count = fread(fid,1,'int8');
            type  = char(fread(fid,1,'char'));
            switch type
                case 'f'
                    values = double(fread(fid,count,'float64'));
                case 'i'
                    values = int32(fread(fid,count,'int32'));
                case 'c'
                    values = char(fread(fid,count,'char')');
                    if ~isempty(values) && values(end)==0
                        values=values(1:(end-1));
                    end
                otherwise
                    error('Unknown definition type: %s', type);
            end
            def(key) = values;
        end
    end

    function [eventTable, blockDurations] = readBlocks(fid, blockDurationRaster)
        %readBlocks Read the [BLOCKS] section of a binary sequence file.
        %   Each block stores: duration(int64) + 6 event IDs (int32):
        %   rf, gx, gy, gz, adc, ext
        %   Returns eventTable (cell array) with 7 elements per block
        %   [0 rf gx gy gz adc ext] and blockDurations in seconds.

        numBlocks = double(fread(fid,1,'int64'));
        eventTable = cell(1, numBlocks);
        blockDurations = zeros(1, numBlocks);
        for ii = 1:numBlocks
            dur_raster   = double(fread(fid,1,'int64'));
            event_ids    = double(fread(fid,6,'int32'))';
            blockDurations(ii) = dur_raster * blockDurationRaster;
            eventTable{ii} = [0 event_ids];  % prepend 0 for legacy delay placeholder
        end
    end

    function shapeLibrary = readShapes(fid)
        %readShapes Read the [SHAPES] section of a binary sequence file.

        shapeLibrary=mr.EventLibrary();
        numShapes = fread(fid,1,'int64');
        for iShape=1:numShapes
            id = double(fread(fid,1,'int32'));
            numUncompressed = double(fread(fid,1,'int64'));
            numCompressed   = double(fread(fid,1,'int64'));
            data = double(fread(fid,numCompressed,'float32'))';
            shapeData = [numUncompressed data];
            shapeLibrary.insert(id,shapeData);
        end
    end

    function readLegacyDelays(fid)
        % readLegacyDelays Read and ignore legacy [DELAYS] binary section.
        % Delay events are no longer represented via obj.delayLibrary.
        numEvents = double(fread(fid,1,'int64'));
        for ii=1:numEvents
            fread(fid,1,'int32'); % id
            fread(fid,1,'int64'); % delay value
        end
    end

end