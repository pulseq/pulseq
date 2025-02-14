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
old_definitions = obj.definitions;
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
requiredDefs=struct('GradientRasterTime',false,'RadiofrequencyRasterTime',false,'AdcRasterTime',false,'BlockDurationRaster',false);

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
                requiredDefs.GradientRasterTime=true;
            end
            v=obj.getDefinition('RadiofrequencyRasterTime');
            if ~isempty(v)
                obj.rfRasterTime=v;
                requiredDefs.RadiofrequencyRasterTime=true;
            end
            v=obj.getDefinition('AdcRasterTime');
            if ~isempty(v)
                obj.adcRasterTime=v;
                requiredDefs.AdcRasterTime=true;
            end
            v=obj.getDefinition('BlockDurationRaster');
            if ~isempty(v)
                obj.blockDurationRaster=v;
                requiredDefs.BlockDurationRaster=true;
            end
            if version_combined >= 1004000 && ~all(struct2array(requiredDefs))
                fn=fieldnames(requiredDefs);
                fn=fn(struct2array(requiredDefs)==0);
                error(['Required definitions ' sprintf('%s ',fn{:}) 'are missing in the file']);
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
            if version_combined >= 1005000 && detectRFuse
                warning('Option ''detectRFuse'' is not supported for file format version 1.5.0 and above');
                detectRFuse=false;
            end
        case '[BLOCKS]'
            if ~exist('version_major')
                error('Pulseq file MUST include [VERSION] section prior to [BLOCKS] section');
            end
            [obj.blockEvents,obj.blockDurations,delayInd_tmp] = readBlocks(fid, obj.blockDurationRaster, version_combined);
        case '[RF]'
            if version_combined >= 1005000 
                obj.rfLibrary = readEvents(fid, [1 1 1 1 1e-6 1e-6 1 1 1 1 NaN]); % this is 1.5.x format 
            elseif version_combined >= 1004000 
                obj.rfLibrary = readEvents(fid, [1 1 1 1 1e-6 1 1]); % this is 1.4.x format 
                % we fix it below
            else
                obj.rfLibrary = readEvents(fid, [1 1 1 1e-6 1 1]); % this is 1.3.x and below  
                % we will have to scan through the library later after all the shapes have been loaded
            end
        case '[GRADIENTS]'
            if version_combined >= 1005000 
                obj.gradLibrary = readEvents(fid, [1 1 1 1 1 1e-6], 'g' ,obj.gradLibrary); % this is 1.5.x format 
            elseif version_combined >= 1004000 
                obj.gradLibrary = readEvents(fid, [1 1 1 1e-6], 'g' ,obj.gradLibrary); % this is 1.4.x format 
            else
                obj.gradLibrary = readEvents(fid, [1 1 1e-6], 'g' ,obj.gradLibrary); % this is 1.3.x and below 
            end
        case '[TRAP]'
            obj.gradLibrary = readEvents(fid, [1 1e-6 1e-6 1e-6 1e-6], 't', obj.gradLibrary);
        case '[ADC]'
            if version_combined >= 1005000 
                obj.adcLibrary = readEvents(fid, [1 1e-9 1e-6 1 1 1 1 1]); % this is 1.5.x format 
            else
                obj.adcLibrary=readEvents(fid, [1 1e-9 1e-6 1 1]); % this is 1.4.x and older format
                % for now we don't have the phase vector in the ADC library
                %obj.adcLibrary.data = [obj.adcLibrary.data(:,1:3) 0 obj.adcLibrary.data(:,4:5)]; % import from the older format
            end
        case '[DELAYS]'
            if version_combined >= 1004000 
                error('Pulseq file revision 1.4.0 and above MUST NOT contain the [DELAYS] section');
            end
            tmp_delayLibrary = readEvents(fid, 1e-6);
        case '[SHAPES]'
            obj.shapeLibrary = readShapes(fid, (version_major==1 && version_minor<4));
        case '[EXTENSIONS]'
            obj.extensionLibrary = readEvents(fid);
        otherwise
            if strncmp('extension', section, 9)
                extension=section(11:end);                    
                if strncmp('TRIGGERS', extension, 8) 
                    id=str2num(extension(9:end));
                    obj.setExtensionStringAndID('TRIGGERS',id);
                    obj.trigLibrary = readEvents(fid, [1 1 1e-6 1e-6]);
                elseif strncmp('LABELSET', extension, 8) 
                    id=str2num(extension(9:end));
                    obj.setExtensionStringAndID('LABELSET',id);
                    obj.labelsetLibrary = readAndParseEvents(fid,@str2num,@(s)find(strcmp(mr.getSupportedLabels,s)));
                elseif strncmp('LABELINC', extension, 8) 
                    id=str2num(extension(9:end));
                    obj.setExtensionStringAndID('LABELINC',id);
                    obj.labelincLibrary = readAndParseEvents(fid,@str2num,@(s)find(strcmp(mr.getSupportedLabels,s)));
                elseif strncmp('DELAYS', extension, 6) 
                    id=str2num(extension(7:end));
                    obj.setExtensionStringAndID('DELAYS',id);
                    obj.softDelayLibrary = readAndParseEvents(fid,@str2num,@(s) 1e-6*str2num(s),@str2num,@(s) parseSoftDelayHint(s, obj));
                else
                    warning('Ignoring unknown extension, input string: %s', extension);
                    exts=regexp(extension, '(\s+)','split');                    
                    obj.setExtensionStringAndID(exts{1}, str2num(exts{2}));
                    skipSection(fid);
                end
            else
                error('Unknown section code: %s', section);
            end
    end
end
fclose(fid);

% fix sequence data imported from older verisons
if version_combined < 1002000 
    error('Unsupported version %07d, only file format revision 1.2.0 (1002000) and above are supported', version_combined);
end

% a special case for ADCs as the format for them has only been updated once (in v1.5.0)
% we have to do it first because seq.getBlock is used in the next version porting code section (version_combined < 1004000)
if version_combined < 1005000 
    % scan though the ADCs and add empty phase shape IDs
    for i=1:length(obj.adcLibrary.data)
        obj.adcLibrary.update_data(...
            obj.adcLibrary.keys(i), ...
            obj.adcLibrary.data(i).array, ...
            [obj.adcLibrary.data(i).array(1:3) 0 0 obj.adcLibrary.data(i).array(4:5) 0]); % add empty freqPPM, phasePPM and phase_id fields
    end
end

% fix blocks, gradients and RF objects imported from older versions (< v1.4.0)
if version_combined < 1004000  
    % fix definitions which are be missing in older files
    if ~obj.definitions.isKey('GradientRasterTime') 
        obj.setDefinition('GradientRasterTime', obj.gradRasterTime);
    end
    if ~obj.definitions.isKey('RadiofrequencyRasterTime') 
        obj.setDefinition('RadiofrequencyRasterTime', obj.rfRasterTime);
    end
    if ~obj.definitions.isKey('AdcRasterTime') 
        obj.setDefinition('AdcRasterTime', obj.adcRasterTime);
    end
    if ~obj.definitions.isKey('BlockDurationRaster') 
        obj.setDefinition('BlockDurationRaster', obj.blockDurationRaster);
    end
    
    % scan through the RF objects
    obj.rfLibrary.type(obj.rfLibrary.keys) = 'u'; % undefined for now, we'll attempt the type detection later (see below)
    for i=1:length(obj.rfLibrary.data)
        % % need to (partially) decode the magnitude shape to find out the pulse duration
        %magSamples = obj.shapeLibrary.data(obj.rfLibrary.data(i).array(2)).array(1);
        % % create time shape
        %timeShape = mr.compressShape((1:magSamples)-0.5); % time shape is stored in units of RF raster
        %data = [timeShape.num_samples timeShape.data];
        %timeID = obj.shapeLibrary.find_or_insert(data);
        rf=rmfield(obj.rfFromLibData([obj.rfLibrary.data(i).array(1:3) 0 0 obj.rfLibrary.data(i).array(4) 0 0 obj.rfLibrary.data(i).array(5:6)],'u'),'center');
        center=mr.calcRfCenter(rf);
        obj.rfLibrary.update_data(...
            obj.rfLibrary.keys(i), ...
            obj.rfLibrary.data(i).array, ...
            [obj.rfLibrary.data(i).array(1:3) 0 center obj.rfLibrary.data(i).array(4) 0 0 obj.rfLibrary.data(i).array(5:6)]); % 0 between (4) and (5:6) are the freqPPM and phasePPM
    end
    
    % scan through the gradient objects and update 't'-s (trapezoids) und 'g'-s (free-shape gradients)
    for i=1:length(obj.gradLibrary.data)
        if obj.gradLibrary.type(i)=='t' % we need to fix some trapezoids, namely ones having zero amplitude and zero ramp times
            if obj.gradLibrary.data(i).array(2)==0
                if abs(obj.gradLibrary.data(i).array(1))==0 && obj.gradLibrary.data(i).array(3) > 0
                    obj.gradLibrary.update_data(...
                        obj.gradLibrary.keys(i), ...
                        obj.gradLibrary.data(i).array, ...
                        [obj.gradLibrary.data(i).array(1) obj.gradRasterTime obj.gradLibrary.data(i).array(3)-obj.gradRasterTime obj.gradLibrary.data(i).array(4:5)],...
                        obj.gradLibrary.type(i));
                end
            end
            if obj.gradLibrary.data(i).array(4)==0
                if abs(obj.gradLibrary.data(i).array(1))==0 && obj.gradLibrary.data(i).array(3) > 0
                    obj.gradLibrary.update_data(...
                        obj.gradLibrary.keys(i), ...
                        obj.gradLibrary.data(i).array, ...
                        [obj.gradLibrary.data(i).array(1:2) obj.gradLibrary.data(i).array(3)-obj.gradRasterTime obj.gradRasterTime obj.gradLibrary.data(i).array(5)],...
                        obj.gradLibrary.type(i));
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
            obj.gradLibrary.update_data(...
                obj.gradLibrary.keys(i), ...
                obj.gradLibrary.data(i).array, ...
                [obj.gradLibrary.data(i).array(1) NaN NaN obj.gradLibrary.data(i).array(2) 0 obj.gradLibrary.data(i).array(3)], ... % we use NaNs to label the non-initialized first/last fields. These will be restored in the code below
                'g');
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
elseif version_combined < 1005000 
    % port from v1.4.x : RF, ADC and GRAD objects need to be updated
    % this needs to be done on the level of the libraries, because getBlock will fail
    
    % scan though the RFs and add center, freqPPM, phasePPM and use fields
    obj.rfLibrary.type(obj.rfLibrary.keys) = 'u'; % undefined for now, we'll attemp the type detection later (see below)
    for i=1:length(obj.rfLibrary.data)
        % use goes into the type field, and this is done separately
        rf=rmfield(obj.rfFromLibData([obj.rfLibrary.data(i).array(1:4) 0 obj.rfLibrary.data(i).array(5) 0 0 obj.rfLibrary.data(i).array(6:7)],'u'),'center');
        center=mr.calcRfCenter(rf);
        obj.rfLibrary.update_data(...
            obj.rfLibrary.keys(i), ...
            obj.rfLibrary.data(i).array, ...
            [obj.rfLibrary.data(i).array(1:4) center obj.rfLibrary.data(i).array(5) 0 0 obj.rfLibrary.data(i).array(6:7)]); % 0 between (5) and (6:7) are the freqPPM and phasePPM 
    end
    % scan through the gradient objects and update 'g'-s (free-shape gradients)
    for i=1:length(obj.gradLibrary.data)
        if obj.gradLibrary.type(i)=='g'
            obj.gradLibrary.update_data(...
                obj.gradLibrary.keys(i), ...
                obj.gradLibrary.data(i).array, ...
                [obj.gradLibrary.data(i).array(1) NaN NaN obj.gradLibrary.data(i).array(2:4)], ... % we use NaNs to label the non-initialized first/last fields. These will be restored in the code below
                'g');
        end
    end
end


% another run through for all older versions
if version_combined < 1005000 
    gradChannels={'gx','gy','gz'};
    gradPrevLast=zeros(1,length(gradChannels));
    for iB = 1:length(obj.blockEvents)
        b=obj.getBlock(iB);
        block_duration=obj.blockDurations(iB);
        %obj.blockDurations(iB)=block_duration;
        % we also need to keep track of the event IDs because some Pulseq files written by external software may contain repeated entries so searching by content will fail 
        eventIDs=obj.blockEvents{iB};
        processedGradIDs=zeros(1,length(gradChannels));
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
                if isfield(grad,'first') && isfinite(grad.first)
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
                if j>1 && any(processedGradIDs(1:j)==id)
                    continue; % avoid repeated updates if the same gradient is applied on differen gradient axes
                end
                processedGradIDs(j)=id;
                amplitude=obj.gradLibrary.data(id).array(1);
                %
                old_data = [amplitude NaN NaN grad.shape_id grad.time_id grad.delay];
                new_data = [amplitude grad.first grad.last grad.shape_id grad.time_id grad.delay];
                update_data(obj.gradLibrary, id, old_data, new_data,'g');
            else
                gradPrevLast(j)=0;
            end
        end
        %% copy updated objects back into the event library
        %obj.setBlock(iB,b);
    end
    
    
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
        rf=obj.rfFromLibData(libData,'u');
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
        %line = localStrip(fgetl(fid));
        line = localStrip(skipComments(fid));
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

    function str = format_helper(scale)
        if isfinite(scale)
            str='%f ';
        else
            str='%c ';
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
            format='%f';
            type_idx=[];
            data_mask=[];
        else
            % new in v1.5.0 : generate format string; NaN labels character param(s)
            format=['%f ' cell2mat(arrayfun(@format_helper,scale,'UniformOutput',false))];
            format(end)=[]; % matlab is so incredibly ugly!
            data_mask=isfinite(scale);
            type_idx=find(~data_mask);
            if length(type_idx)>2
                error('Only one type field (marked as NaN) can be provided');
            end
            if isempty(type_idx)
                data_mask=[];
            end
        end
        if nargin < 4
            eventLibrary = mr.EventLibrary();
        end        
        %
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1) == '#')
            data = sscanf(line,format)';
            id = data(1);
            if ~isempty(type_idx)
                type=char(data(type_idx+1)); % need +1 because of the eventID in the first position
            end
            data = scale.*data(2:end);
            if nargin < 3 && isempty(type_idx)
                if isempty(data_mask)
                    eventLibrary.insert(id, data);
                else
                    eventLibrary.insert(id, data(data_mask));
                end
            else
                if isempty(data_mask)
                    eventLibrary.insert(id, data, type);
                else
                    eventLibrary.insert(id, data(data_mask), type);
                end
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

    function skipSection(fid)
        %skipSection Read an event section of a sequence file without 
        %   interpreting it.
        %
        line = fgetl(fid);
        while ischar(line) && ~(isempty(line) || line(1) == '#')
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

    function id=parseSoftDelayHint(s, seq) 
        try
            id=seq.softDelayHints1(s);
        catch
            id=seq.softDelayHints1.length()+1;
            seq.softDelayHints1(s)=id;
            seq.softDelayHints2{id}=s;
        end
    end
    
    % for compatibility with Octave which has no strip()
    function s=localStrip(s)
        a=1;
        b=length(s);
        while a<=b && isspace(s(a)), a=a+1; end
        while a<=b && isspace(s(b)), b=b-1; end
        s=s(a:b);
    end
end
