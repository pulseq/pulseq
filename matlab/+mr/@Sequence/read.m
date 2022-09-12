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

if fid<0
    error('filed to open file ''%s''', filename);
end

% Clear sequence data
%obj.blockEvents = [];
obj.blockEvents = {};
obj.definitions = containers.Map();
obj.gradLibrary = mr.EventLibrary();
obj.shapeLibrary = mr.EventLibrary();
obj.rfLibrary = mr.EventLibrary();
obj.adcLibrary = mr.EventLibrary();
%obj.delayLibrary = mr.EventLibrary();
obj.trigLibrary = mr.EventLibrary();
obj.labelsetLibrary = mr.EventLibrary();
obj.labelincLibrary = mr.EventLibrary();
obj.extensionStringIDs={};
obj.extensionNumericIDs=[];

version_combined=0;

% Load data from file
while true
    section = skipComments(fid);
    if section == -1
        break
    end
    
    switch section
        case '[DEFINITIONS]'
            obj.definitions = readDefinitions(fid);
            v=obj.getDefinition('GradientRasterTime');
            if ~isempty(v)
                obj.gradRasterTime=v;
            end
            v=obj.getDefinition('RadiofrequencyRasterTime');
            if ~isempty(v)
                obj.rfRasterTime=v;
            end
            v=obj.getDefinition('AdcRasterTime');
            if ~isempty(v)
                obj.adcRasterTime=v;
            end
            v=obj.getDefinition('BlockDurationRaster');
            if ~isempty(v)
                obj.blockDurationRaster=v;
            end
        case '[SIGNATURE]'
            tmpSignDefs = readDefinitions(fid);
            if isKey(tmpSignDefs,'Type')
                obj.signatureType=tmpSignDefs('Type');
            end
            if isKey(tmpSignDefs,'Hash')
                obj.signatureValue=tmpSignDefs('Hash');
                obj.signatureFile='Text'; % we are reading a text file, so much is known for sure
            end
        case '[VERSION]'
            [version_major, ...
             version_minor, ...
             version_revision] = readVersion(fid);
            assert(version_major == obj.version_major, ...
                    'Unsupported version_major %d', version_major)
            %
            version_combined=1000000*version_major+1000*version_minor+version_revision;
            %
            if version_combined < 1002000 
                error('Unsupported version %d.%d.%d, only file format revision 1.2.0 and above are supported', version_major, version_minor, version_revision);
            end
            if version_combined < 1003001 
                warning('Loading older Pulseq format file (version %d.%d.%d) some code may function not as expected', version_major, version_minor, version_revision);
            end
        case '[BLOCKS]'
            if ~exist('version_major')
                error('Pulseq file MUST include [VERSION] section prior to [BLOCKS] section');
            end
            [obj.blockEvents,obj.blockDurations,delayInd_tmp] = readBlocks(fid, obj.blockDurationRaster, version_combined);
        case '[RF]'
            if version_combined >= 1004000 
                obj.rfLibrary = readEvents(fid, [1 1 1 1 1e-6 1 1]); % this is 1.4.x format 
            else
                obj.rfLibrary = readEvents(fid, [1 1 1 1e-6 1 1]); % this is 1.3.x and below  
                % we will have to scan through the library later after all the shapes have been loaded
            end
        case '[GRADIENTS]'
            if version_combined >= 1004000 
                obj.gradLibrary = readEvents(fid, [1 1 1 1e-6], 'g' ,obj.gradLibrary); % this is 1.4.x format 
            else
                obj.gradLibrary = readEvents(fid, [1 1 1e-6], 'g' ,obj.gradLibrary); % this is 1.3.x and below 
            end
        case '[TRAP]'
            obj.gradLibrary = readEvents(fid, [1 1e-6 1e-6 1e-6 1e-6], 't', obj.gradLibrary);
        case '[ADC]'
            obj.adcLibrary = readEvents(fid, [1 1e-9 1e-6 1 1]);
        case '[DELAYS]'
            if version_combined >= 1004000 
                error('Pulseq file revision 1.4.0 and above MUST NOT contain [DELAYS] section');
            end
            tmp_delayLibrary = readEvents(fid, 1e-6);
        case '[SHAPES]'
            obj.shapeLibrary = readShapes(fid, (version_major==1 && version_minor<4));
        case '[EXTENSIONS]'
            obj.extensionLibrary = readEvents(fid);
        otherwise
            if     strncmp('extension TRIGGERS', section, 18) 
                id=str2num(section(19:end));
                obj.setExtensionStringAndID('TRIGGERS',id);
                obj.trigLibrary = readEvents(fid, [1 1 1e-6 1e-6]);
            elseif strncmp('extension LABELSET', section, 18) 
                id=str2num(section(19:end));
                obj.setExtensionStringAndID('LABELSET',id);
                obj.labelsetLibrary = readAndParseEvents(fid,@str2num,@(s)find(strcmp(mr.getSupportedLabels,s)));
            elseif strncmp('extension LABELINC', section, 18) 
                id=str2num(section(19:end));
                obj.setExtensionStringAndID('LABELINC',id);
                obj.labelincLibrary = readAndParseEvents(fid,@str2num,@(s)find(strcmp(mr.getSupportedLabels,s)));
            else
                error('Unknown section code: %s', section);
            end
    end
end
fclose(fid);

%
if version_combined < 1002000 
    error('Unsupported version %07d, only file format revision 1.2.0 (1002000) and above are supported', version_combined);
end            
% fix blocks, gradients and RF objects imported from older versions
if version_combined < 1004000  % MZ: FIXME: the code below does not update key-to-id mapping, so the libraries are not fully functional... these are only partially updated by the shaped gradients patching below (first/last detection)
    % scan through the RF objects
    for i=1:length(obj.rfLibrary.data) 
        % % need to (partially) decode the magnitude shape to find out the pulse duration
        %magSamples = obj.shapeLibrary.data(obj.rfLibrary.data(i).array(2)).array(1);
        % % create time shape
        %timeShape = mr.compressShape((1:magSamples)-0.5); % time shape is stored in units of RF raster
        %data = [timeShape.num_samples timeShape.data];
        %timeID = obj.shapeLibrary.find_or_insert(data);
        obj.rfLibrary.data(i).array = [obj.rfLibrary.data(i).array(1:3), 0, obj.rfLibrary.data(i).array(4:end)];
        obj.rfLibrary.lengths(i) = obj.rfLibrary.lengths(i) + 1;
    end
    
    % scan through the gradient objects and update 't'-s (trapezoids) und 'g'-s (free-shape gradients)
    for i=1:length(obj.gradLibrary.data)
        if obj.gradLibrary.type(i)=='t'
            if obj.gradLibrary.data(i).array(2)==0
                if abs(obj.gradLibrary.data(i).array(1))==0 && obj.gradLibrary.data(i).array(3) > 0
                    obj.gradLibrary.data(i).array(3)=obj.gradLibrary.data(i).array(3)-obj.gradRasterTime;
                    obj.gradLibrary.data(i).array(2)=obj.gradRasterTime;
                end
            end
            if obj.gradLibrary.data(i).array(4)==0
                if abs(obj.gradLibrary.data(i).array(1))==0 && obj.gradLibrary.data(i).array(3) > 0
                    obj.gradLibrary.data(i).array(3)=obj.gradLibrary.data(i).array(3)-obj.gradRasterTime;
                    obj.gradLibrary.data(i).array(4)=obj.gradRasterTime;
                end
            end
        end
        if obj.gradLibrary.type(i)=='g'
            % % need to (partially) decode the shape to find out the duration
            %nSamples = obj.shapeLibrary.data(obj.gradLibrary.data(i).array(2)).array(1);
            % % create time shape
            %timeShape = mr.compressShape((1:nSamples)-0.5); % time shape is stored in units of grad raster
            %data = [timeShape.num_samples timeShape.data];
            %timeID = obj.shapeLibrary.find_or_insert(data);
            obj.gradLibrary.data(i).array = [obj.gradLibrary.data(i).array(1:2), 0, obj.gradLibrary.data(i).array(3:end)];
            obj.gradLibrary.lengths(i) = obj.gradLibrary.lengths(i) + 1;
        end
    end
    
    % for versions prior to 1.4.0 blockDurations have not been initialized
    obj.blockDurations=zeros(1,length(obj.blockEvents));
    % scan trhough blocks and calculate durations
    for iB = 1:length(obj.blockEvents)
        b=obj.getBlock(iB);
        if delayInd_tmp(iB) > 0
            b.delay.type = 'delay';
            b.delay.delay = tmp_delayLibrary.data(delayInd_tmp(iB)).array;
        end
        obj.blockDurations(iB)=mr.calcDuration(b);
    end
end

gradChannels={'gx','gy','gz'};
gradPrevLast=zeros(1,length(gradChannels));
for iB = 1:length(obj.blockEvents)
    b=obj.getBlock(iB);
    block_duration=obj.blockDurations(iB);
    %obj.blockDurations(iB)=block_duration;
    % we also need to keep track of the event IDs because some Pulseq files written by external software may contain repeated entries so searching by content will fail 
    eventIDs=obj.blockEvents{iB};
    % update the objects by filling in the fields not contained in the
    % pulseq file
    for j=1:length(gradChannels)
        grad=b.(gradChannels{j});
        if isempty(grad)
            gradPrevLast(j)=0;
            continue;
        end
        if strcmp(grad.type,'grad')
            if grad.delay>0 
                gradPrevLast(j)=0;
            end
            if isfield(grad,'first')
                continue;
            end
            grad.first = gradPrevLast(j);
            % is this an extended trapezoid?
            if grad.time_id~=0
                grad.last=grad.waveform(end);
                grad_duration=grad.delay+grad.tt(end);
            else
                % restore samples on the edges of the gradient raster intervals
                % for that we need the first sample
                odd_step1=[grad.first 2*grad.waveform'];
                odd_step2=odd_step1.*(mod(1:length(odd_step1),2)*2-1);
                waveform_odd_rest=(cumsum(odd_step2).*(mod(1:length(odd_step2),2)*2-1))';
%                 delta_odd=waveform_odd_rest(2:end-1)-0.5*(grad.waveform(1:end-1)+grad.waveform(2:end));
%                 delta_odd_signed=delta_odd.*(mod(1:length(delta_odd),2)*2-1)';
%                 delta_odd_signed_flt=medfilt1(delta_odd_signed,49);
%                 delta_odd_flt=delta_odd_signed_flt.*(mod(1:length(delta_odd),2)*2-1)';
%                 waveform_odd_rest1=waveform_odd_rest-[0; delta_odd_flt; 0];
%                 waveform_odd_rest1(end)=2*grad.waveform(end)-waveform_odd_rest1(end-1); %restore the final sample based on the recurrent relation 
%                 waveform_odd_rest0=waveform_odd_rest;
%                 waveform_odd_rest=waveform_odd_rest1;
                grad.last = waveform_odd_rest(end);
                grad_duration=grad.delay+length(grad.waveform)*obj.gradRasterTime;
            end
            % bookkeeping
            gradPrevLast(j) = grad.last;
            if grad_duration+eps<block_duration
                gradPrevLast(j)=0;
            end
            %b.(gradChannels{j})=grad;
            % update library object
% this does not work s we don't know how the amplitude was defined
%             amplitude = max(abs(grad.waveform));
%             if amplitude>0
%                 [~,~,fnz]=find(grad.waveform,1); % find the first non-zero value and make it positive
%                 amplitude=amplitude*sign(fnz);
%             end
            % need to recover the amplidute from the library data directly...
            id=eventIDs(j+2);
            amplitude=obj.gradLibrary.data(id).array(1);
            %
            if version_combined>=1004000
                old_data = [amplitude grad.shape_id grad.time_id grad.delay];
            else
                old_data = [amplitude grad.shape_id grad.delay];
            end
            new_data = [amplitude grad.shape_id grad.time_id grad.delay grad.first grad.last];
            update_data(obj.gradLibrary, id, old_data, new_data,'g');
        else
            gradPrevLast(j)=0;
        end
    end
    %% copy updated objects back into the event library
    %obj.setBlock(iB,b);
    
    
%for iB=1:size(obj.blockEvents,1)
%     % update the objects by filling in the fields not contained in the
%     % pulseq file
%     for j=1:length(gradChannels)
%         grad=b.(gradChannels{j});
%         if isempty(grad)
%             continue;
%         end
%         if strcmp(grad.type,'grad')
%             grad.first = grad.waveform(1); % MZ: eventually we should use extrapolation by 1/2 gradient rasters here
%             grad.last = grad.waveform(end);
%             b.(gradChannels{j})=grad;
%         end;
%     end
%     % copy updated objects back into the event library
%     obj.setBlock(iB,b);
end

if detectRFuse
    % find the RF pulses, list flip angles
    % and work around the current (rev 1.2.0) Pulseq file format limitation
    % that the RF pulse use is not stored in the file
    for k=obj.rfLibrary.keys
        libData=obj.rfLibrary.data(k).array;
        rf=obj.rfFromLibData(libData);
        %flipAngleDeg=abs(sum(rf.signal))*rf.t(1)*360; %we use rfex.t(1) in place of opt.system.rfRasterTime
        flipAngleDeg=abs(sum(rf.signal(1:end-1).*(rf.t(2:end)-rf.t(1:end-1))))*360;
        offresonance_ppm=1e6*rf.freqOffset/obj.sys.B0/obj.sys.gamma;
        % fix library %%%% if length(obj.rfLibrary.type)>=eventInd(2)
        if flipAngleDeg < 90.01 % we add 0.01 degree to account for rounding errors which we've experienced for very short RF pulses
            obj.rfLibrary.type(k) = 'e';
        else
            if rf.shape_dur > 6e-3 && offresonance_ppm >= -3.5 && offresonance_ppm <= -3.4 % approx -3.45 ppm
                obj.rfLibrary.type(k) = 's'; % saturation (fat-sat)
            else
                obj.rfLibrary.type(k) = 'r';
            end
        end
%         % fix libData
%         if length(libData) < 9
%             if flipAngleDeg < 90.01 % we add 0.01 degree to account for rounding errors which we've experienced for very short RF pulses
%                 libData(9) = 0; % or 1 ?
%             else
%                 libData(9) = 2; % or 1 ?
%             end
%             obj.rfLibrary.data(k).array=libData;
%         end
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
        %line = strip(fgetl(fid));
        line = strip(skipComments(fid));
        while ischar(line) && ~(isempty(line) || line(1) == '#')
            tok = textscan(line, '%s');
            def(tok{1}{1}) = str2double(tok{1}(2:end));
            if ~all(isfinite(def(tok{1}{1})))
                def(tok{1}{1}) = strtrim(line((length(tok{1}{1})+2):end));
            end
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

    function [eventTable,blockDurations,delayIDs_tmp] = readBlocks(fid, blockDurationRaster, version_combined)
        %readBlocks Read the [BLOCKS] section of a sequence file.
        %   library=readBlocks(fid) Read blocks from file identifier of an
        %   open MR sequence file and return the event table.
        
        eventTable = {};
        blockDurations = [];
        delayIDs_tmp = [];
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1) == '#')
            blockEvents = sscanf(line, '%f')';
            %eventTable = [eventTable; blockEvents(2:end)];
            if version_combined<=1002001
                eventTable{blockEvents(1)} = [0 blockEvents(3:end) 0];
            else
                eventTable{blockEvents(1)} = [0 blockEvents(3:end)];
            end
            if version_combined>=1004000
                blockDurations(blockEvents(1)) = blockEvents(2)*blockDurationRaster;
            else
                delayIDs_tmp(blockEvents(1)) = blockEvents(2);
            end
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

    function eventLibrary = readAndParseEvents(fid, varargin)
        %readAndParseEvents Read an event section of a sequence file.
        %   library=readAndParseEvents(fid) Read event data from file 
        %   identifier of an open MR sequence file and return a library of 
        %   events.
        %
        %   library=readAndParseEvents(fid,parser1,parser2,...) Read event  
        %   data and convert the elements using to the provided parser. 
        %   Default parser is str2num()
        %
        eventLibrary = mr.EventLibrary();
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1) == '#')
            datas=regexp(line, '(\s+)','split');
            data=zeros(1,length(datas)-1);
            id = str2num(datas{1});
            for i=2:length(datas)
                if i>nargin
                    data(i-1) = str2num(datas{i});
                else
                    data(i-1) = varargin{i-1}(datas{i});
                end
            end
            eventLibrary.insert(id, data);

            line=fgetl(fid);
        end
    end

    function shapeLibrary = readShapes(fid, forceConvertUncompressed)
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
                %data = [data single(sscanf(line, '%f'))]; % C-code uses single precision and we had problems already due to the rounding during reading in of the shapes...
                line = fgetl(fid);
            end
            
            line = skipComments(fid, true); % MZ: second parameter to prevent readShapes from reading into the next section (long-standing bug)

            % check if conversion is needed: in v1.4.x we use length(data)==num_samples 
            % as a marker for the uncompressed (stored) data. In older versions this condition could occur by chance 
            if forceConvertUncompressed && length(data)==num_samples
                shape.data=data;
                shape.num_samples=num_samples;
                shape = mr.compressShape(mr.decompressShape(shape,true));
                data = [shape.num_samples shape.data];
            else
                data = [num_samples data];
            end
            shapeLibrary.insert(id, data);
        end
    end

    function nextLine = skipComments(fid, stopBeforeSection)
        %skipComments Read lines of skipping blank lines and comments.
        %   line=skipComments(fid) Read lines from valid file identifer and
        %   return the next non-comment line.
        
        if (nargin<2)
            stopBeforeSection=false;
        end
        
        tmpPos=ftell(fid);
        line = fgetl(fid);
        while ischar(line) && (isempty(line) || line(1) == '#')
            tmpPos=ftell(fid);
            line = fgetl(fid);
        end
        if ischar(line)
            if stopBeforeSection && line(1)=='['
                fseek(fid,tmpPos,-1); % restore the file position
                nextLine = ''; % feasible (non-error) dummy return
            else
                nextLine = line;
            end
        else
            nextLine = -1;
        end
    end
end
