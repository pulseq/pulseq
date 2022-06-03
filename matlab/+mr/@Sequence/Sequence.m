classdef Sequence < handle
    % Sequence   Generate sequences and read/write sequence files.
    % This class defines properties and methods to define a complete
    % MR sequence including RF pulses, gradients, ADC events, etc.
    %
    % The class provides an implementation of the open MR sequence format
    % defined by the Pulseq project.
    %   See http://pulseq.github.io/
    %
    % Sequence Properties:
    %    definitions - A list of custom definitions
    %
    % Sequence Methods:
    %    read - Load sequence from open MR sequence format
    %    write - Write sequence to open MR sequence format
    %
    % Sequence Static Methods:
    %    makeTrapezoid - Create a trapezoid gradient structure
    %
    % Examples:
    %
    % To read a sequence from file:
    %     read(seqObj,'my_sequences/gre.seq');
    %
    % To plot a sequence:
    %     plot(seqObj)
    %
    % See also   demoRead.m, demoWrite.m
    % Examples defining an MRI sequence and reading/writing files
    %
    % Kelvin Layton <kelvin.layton@uniklinik-freiburg.de>
    % Maxim Zaitsev <maxim.zaitsev@uniklinik-freiburg.de>

    % Private properties
    %
    properties(GetAccess = public, SetAccess = private)
        version_major;
        version_minor;
        version_revision;
        rfRasterTime;        % RF raster time (system dependent)
        gradRasterTime;      % Gradient raster time (system dependent)
        adcRasterTime;       % minimum unit/increment of the ADC dwell time (system dependent) 
        blockDurationRaster; % unit/increment of the block duration (system dependent)         
        definitions       % Optional sequence definitions
        
        blockEvents;      % Event table (references to events)
        blockDurations;   % Cache of block durations
        rfLibrary;        % Library of RF events
        gradLibrary;      % Library of gradient events
        adcLibrary;       % Library of ADC readouts
        trigLibrary;      % Library of trigger events ( referenced from the extentions library )
        labelsetLibrary;  % Library of Label(set) events ( reference from the extensions library )
        labelincLibrary;  % Library of Label(inc) events ( reference from the extensions library )
        extensionLibrary; % Library of extension events. Extension events form single-linked zero-terminated lists
        shapeLibrary;     % Library of compressed shapes
        extensionStringIDs;  % string IDs of the used extensions (cell array)
        extensionNumericIDs; % numeric IDs of the used extensions (numeric array)
        
        signatureType; % type of the hashing function used, currently 'md5'
        signatureFile; % which data were hashed, currently 'text' or 'bin' (used file format of the save function)
        signatureValue; % the hash of the exported Pulse sequence
        
        sys;
    end
    
    methods
        
        function obj = Sequence(varargin)
            obj.version_major = 1;
            obj.version_minor = 4; % version minor 3 will now support control events (8th column in the event table) mv4 supports/expects timing vectors for arbitrary grads
            obj.version_revision = 0;
            obj.definitions = containers.Map();
            obj.gradLibrary = mr.EventLibrary();
            obj.shapeLibrary = mr.EventLibrary();
            obj.rfLibrary = mr.EventLibrary();
            obj.adcLibrary = mr.EventLibrary();
            obj.trigLibrary = mr.EventLibrary();
            obj.labelsetLibrary = mr.EventLibrary();
            obj.labelincLibrary = mr.EventLibrary();
            obj.extensionLibrary = mr.EventLibrary();
            obj.extensionStringIDs={};
            obj.extensionNumericIDs=[];
            obj.blockEvents = {};
            
            if nargin<1
                sys=mr.opts();
            else
                sys=varargin{1};
            end
            obj.sys = sys;
            obj.rfRasterTime = sys.rfRasterTime;
            obj.gradRasterTime = sys.gradRasterTime;
            obj.adcRasterTime = sys.adcRasterTime;
            obj.blockDurationRaster = sys.blockDurationRaster;
            obj.setDefinition('GradientRasterTime', obj.gradRasterTime);
            obj.setDefinition('RadiofrequencyRasterTime', obj.rfRasterTime);
            obj.setDefinition('AdcRasterTime', obj.adcRasterTime);
            obj.setDefinition('BlockDurationRaster', obj.blockDurationRaster);
            obj.signatureType=''; 
            obj.signatureFile=''; 
            obj.signatureValue='';
        end
        
        
        % See read.m
        read(obj,filename,varargin)
        
        % See write.m
        write(obj,filename,create_signature)
        
        % See write.m
        write_file(obj,filename)
        
        % See readBinary.m
        readBinary(obj,filename);
        
        % See writeBinary.m
        writeBinary(obj,filename);
        
        
        % See calcPNS.m
        [ok, pns_norm, pns_comp, t_axis]=calcPNS(obj,hardware,doPlots)
        
        % See testReport.m
        %testReport(obj);
        
        function [duration, numBlocks, eventCount]=duration(obj)
            % duration() 
            %     Returns the total duration of the sequence
            %     optionally returns the total count of events
            %
            
            % Loop over blocks and gather statistics
            numBlocks = length(obj.blockEvents);
            eventCount=zeros(size(obj.blockEvents{1}));
            duration=0;
            for iB=1:numBlocks
                %b=obj.getBlock(iB);
                eventCount = eventCount + (obj.blockEvents{iB}>0);
                duration=duration+obj.blockDurations(iB);%mr.calcDuration(b);
            end
        end
        
        function [is_ok, errorReport]=checkTiming(obj)
            % checkTiming() 
            %     Checks timing of all blocks and objects in the sequence 
            %     optionally returns the detailed error log as cell array
            %     of strings. This function also modifies the sequence
            %     object by adding the field "TotalDuration" to sequence
            %     definitions
            %
            
            % Loop over blocks and gather statistics
            numBlocks = length(obj.blockEvents);
            is_ok=true;
            errorReport={};
            totalDuration=0;
            for iB=1:numBlocks
                b=obj.getBlock(iB);
                % assemble cell array of events
                %ev={b.rf, b.gx, b.gy, b.gz, b.adc, b.delay, b.ext}; 
                %ind=~cellfun(@isempty,ev);
                % the above does not work for ext because it may be
                % missing from some blocks and may have multiple entries in
                % others. 
                ind=~structfun(@isempty,b);
                fn=fieldnames(b);
                ev=cellfun(@(f) b.(f), fn(ind), 'UniformOutput', false);
                [res, rep, dur] = mr.checkTiming(obj.sys,ev{:}); %ev{ind});
                
                is_ok = (is_ok && res); 
                
                % check the stored total block duration
                if abs(dur-obj.blockDurations(iB))>eps
                    rep = [rep 'inconsistency between the stored block duration and the duration of the block content'];
                    is_ok = false;
                    dur=obj.blockDurations(iB);
                end
                
                % check RF dead times
                if ~isempty(b.rf)
                    if b.rf.delay-b.rf.deadTime < -eps
                        rep = [rep 'delay of ' num2str(b.rf.delay*1e6) 'us is smaller than the RF dead time ' num2str(b.rf.deadTime*1e6) 'us'];
                        is_ok = false;
                    end
                    if b.rf.delay+b.rf.t(end)+b.rf.ringdownTime-dur > eps
                        rep = [rep 'time between the end of the RF pulse at ' num2str((b.rf.delay+b.rf.t(end))*1e6) ' and the end of the block at ' num2str(dur*1e6) 'us is shorter than rfRingdownTime'];
                        is_ok = false;
                    end
                end
                
                % check ADC dead times
                if ~isempty(b.adc) 
                    if b.adc.delay-obj.sys.adcDeadTime < -eps
                        rep = [rep ' adc.delay<system.adcDeadTime'];
                        is_ok=false;
                    end
                    if b.adc.delay+b.adc.numSamples*b.adc.dwell+obj.sys.adcDeadTime-dur > eps
                        rep = [rep ' adc: system.adcDeadTime (post-adc) violation'];
                        is_ok=false;
                    end
                end

                % update report
                if ~isempty(rep)
                    errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' ' rep '\n' ] };
                end
                %
                totalDuration = totalDuration+dur;
            end
            
            % check whether all gradients in the last block are ramped down properly
            if ~isempty(ev) && isstruct(ev)
                for en=1:length(ev)
                    if length(ev{en})==1 && strcmp(ev{en}.type,'grad') % length(ev{en})==1 excludes arrays of extensions 
                        if ev{en}.last~=0 % must be > sys.slewRate*sys.gradRasterTime
                            errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' gradients do not ramp to 0 at the end of the sequence\n' ] };
                        end
                    end
                end
            end
            
            obj.setDefinition('TotalDuration', totalDuration);%sprintf('%.9g', totalDuration));
        end
        
        function value=getDefinition(obj,key)
            %getDefinition Return the values of custom definition.
            %   val=getDefinitions(seqObj,key) Return value of the
            %   definition specified by the key.
            
            %   These definitions can be added manually or read from the
            %   header of a sequence file defined in the sequence header.
            %   An empty array is return if the key is not defined.
            %
            %   See also setDefinition
            if isKey(obj.definitions,key)
                value = obj.definitions(key);
            else
                value = [];
            end
        end
        
        function setDefinition(seqObj,key,val)
            %setDefinition Modify a custom definition of the sequence.
            %   setDefinition(seqObj,def,val) Set the user definition 'key'
            %   to value 'val'. If the definition does not exist it will be
            %   created.
            %
            %   See also getDefinition
            if strcmp(key,'FOV')
                % issue a warning if FOV is too large e.g. is in mm
                if max(val)>1
                    warning('WARNING: definition FOV uses values exceeding 1m. New Pulseq interpreters expect values in units of meters!\n');
                end
            end
            seqObj.definitions(key)=val;
        end
        
        function addBlock(obj,varargin)
            %addBlock Add a new block to the sequence.
            %   addBlock(obj, blockStruct) Adds a sequence block with
            %   provided as a block struture
            %
            %   addBlock(obj, e1, e2, ...) Adds a block with multiple
            %   events e1, e2, etc.
            %
            %   See also  setBlock, makeAdc, makeTrapezoid, makeSincPulse
            %setBlock(obj,size(obj.blockEvents,1)+1,varargin{:});
            setBlock(obj,length(obj.blockEvents)+1,varargin{:});            
        end
                
        function modGradAxis(obj,axis,modifier)
            %modGradAxis Invert or scale all gradinents along the corresponding
            %   axis/channel. The function acts on all gradient objects 
            %   already added to the sequence object
            %
            channelNum = find(strcmp(axis, ...
                                                 {'x', 'y', 'z'}));
            otherChans = find(~strcmp(axis, ...
                                                 {'x', 'y', 'z'}));
            % go through all event table entries and list gradient
            % objects in the library
            %paren = @(x, varargin) x(varargin{:}); % anonymous function to access the array on the fly
            paren2 = @(x, varargin) x(:,varargin{:}); % anonymous function to access the array on the fly
            %allGradEvents = paren(vertcat(obj.blockEvents{:}),:,3:5);
            allGradEvents = paren2(vertcat(obj.blockEvents{:}),3:5);
            
            selectedEvents=unique(allGradEvents(:,channelNum));
            selectedEvents=selectedEvents(0~=selectedEvents); % elliminate 0
            otherEvents=unique(allGradEvents(:,otherChans));
            assert(isempty(intersect(selectedEvents,otherEvents)),'ERROR: the same gradient event is used on multiple axes, this is not yet supported by modGradAxis()');
                                    
            for i = 1:length(selectedEvents)                
                %type = obj.gradLibrary.type(i);
                %libData = obj.gradLibrary.data(i).array;
                %if strcmp(grad.type,'grad')
                %    amplitude = libData(1);
                %else
                %    %grad.amplitude = libData(1);
                %end
                % 
                % based on the above we just patch the first element of the
                % gradient library data entries
                obj.gradLibrary.data(selectedEvents(i)).array(1)=modifier*obj.gradLibrary.data(selectedEvents(i)).array(1);
                if obj.gradLibrary.type(selectedEvents(i))=='g' && obj.gradLibrary.lengths(selectedEvents(i))==5
                    % need to update .first and .last fields
                    obj.gradLibrary.data(selectedEvents(i)).array(4)=modifier*obj.gradLibrary.data(selectedEvents(i)).array(4);
                    obj.gradLibrary.data(selectedEvents(i)).array(5)=modifier*obj.gradLibrary.data(selectedEvents(i)).array(5);
                end
            end
        end
        
        function flipGradAxis(obj, axis)
            %flipGradAxis Invert all gradinents along the corresponding
            %   axis/channel. The function acts on all gradient objects 
            %   already added to the sequence object
            %
            modGradAxis(obj,axis,-1);
        end
        
        function rf = rfFromLibData(obj, libData, use)                
            rf.type = 'rf';

            amplitude = libData(1);
            magShape = libData(2);
            phaseShape = libData(3);
            shapeData = obj.shapeLibrary.data(magShape).array;
            compressed.num_samples = shapeData(1);
            compressed.data = shapeData(2:end);
            mag = mr.decompressShape(compressed);
            shapeData = obj.shapeLibrary.data(phaseShape).array;
            compressed.num_samples = shapeData(1);
            compressed.data = shapeData(2:end);
            phase = mr.decompressShape(compressed);
            rf.signal = amplitude*mag.*exp(1j*2*pi*phase);
            timeShape = libData(4);
            if timeShape>0
                shapeData = obj.shapeLibrary.data(timeShape).array;
                compressed.num_samples = shapeData(1);
                compressed.data = shapeData(2:end);
                rf.t = mr.decompressShape(compressed)*obj.rfRasterTime;
                rf.shape_dur=ceil((rf.t(end)-eps)/obj.rfRasterTime)*obj.rfRasterTime;
            else
                % generate default time raster on the fly
                rf.t = ((1:length(rf.signal))-0.5)'*obj.rfRasterTime;
                rf.shape_dur=length(rf.signal)*obj.rfRasterTime;
            end

            rf.delay = libData(5);
            rf.freqOffset = libData(6);
            rf.phaseOffset = libData(7);
            
            rf.deadTime = obj.sys.rfDeadTime;
            rf.ringdownTime = obj.sys.rfRingdownTime;

%             % SK: Is this a hack? (MZ: see below)
%             if length(libData) < 8
%                 libData(8) = 0;
%             end
%             rf.deadTime = libData(9);
%             % SK: Using the same hack here
%             if length(libData) < 9
%                 libData(9) = 0;
%             end
%             rf.ringdownTime = libData(9);

            if nargin>2
                switch use
                    case 'e'
                        rf.use='excitation';
                    case 'r'
                        rf.use='refocusing';
                    case 'i'
                        rf.use='inversion';
                    case 's'
                        rf.use='saturation';
                    case 'p'
                        rf.use='preparation';
                    otherwise
                        rf.use='undefined';
                end
            else
                rf.use='undefined';
            end
        end
        
        function [id shapeIDs]=registerRfEvent(obj, event)
            % registerRfEvent Add the event to the libraries (object,
            % shapes, etc and retur the event's ID. This I can be stored in
            % the object to accelerate addBlock()
            
            mag = abs(event.signal);
            amplitude = max(mag);
            mag = mag / amplitude;
            phase = angle(event.signal);
            phase(phase < 0) = phase(phase < 0) + 2*pi;
            phase = phase / (2*pi);
            may_exist=true;
            
            if isfield(event,'shapeIDs')
                shapeIDs=event.shapeIDs;
            else
                shapeIDs=[0 0 0];

                magShape = mr.compressShape(mag(:));
                data = [magShape.num_samples magShape.data];
                [shapeIDs(1),found] = obj.shapeLibrary.find_or_insert(data);
                may_exist=may_exist & found;
                
                phaseShape = mr.compressShape(phase);
                data = [phaseShape.num_samples phaseShape.data];
                [shapeIDs(2),found] = obj.shapeLibrary.find_or_insert(data);
                may_exist=may_exist & found;
                
                timeShape = mr.compressShape(event.t/obj.rfRasterTime); % time shape is stored in units of RF raster
                if length(timeShape.data)==4 && all(timeShape.data == [0.5 1 1 timeShape.num_samples-3]) 
                    shapeIDs(3)=0;
                else
                    data = [timeShape.num_samples timeShape.data];
                    [shapeIDs(3),found] = obj.shapeLibrary.find_or_insert(data);
                    may_exist=may_exist & found;
                end
            end

            use = 'u';
            if isfield(event,'use')
                switch event.use
                    case {'excitation','refocusing','inversion','saturation','preparation'}
                        use = event.use(1);
                    otherwise
                        use = 'u'; % undefined
                end
            end

            data = [amplitude shapeIDs(1) shapeIDs(2) shapeIDs(3) ...
                    event.delay event.freqOffset event.phaseOffset ];%...
                    %event.deadTime event.ringdownTime];
            if may_exist
                id = obj.rfLibrary.find_or_insert(data,use);
            else
                id = obj.rfLibrary.insert(0,data,use);
            end
        end
        
        function [id,shapeIDs]=registerGradEvent(obj, event)
            may_exist=true;
            switch event.type
                case 'grad'
                    amplitude = max(abs(event.waveform));
                    if amplitude>0
                        [~,~,fnz]=find(event.waveform,1); % find the first non-zero value and make it positive
                        amplitude=amplitude*sign(fnz);
                    end
                    if isfield(event,'shapeIDs')
                        shapeIDs=event.shapeIDs;
                    else
                        shapeIDs=[0 0];
                        % fill the shape IDs
                        if amplitude~=0
                            g = event.waveform./amplitude;
                        else
                            g = event.waveform;
                        end
                        c_shape = mr.compressShape(g);
                        s_data = [c_shape.num_samples c_shape.data];
                        [shapeIDs(1),found] = obj.shapeLibrary.find_or_insert(s_data);
                        may_exist=may_exist & found;
                        c_time = mr.compressShape(event.tt/obj.gradRasterTime);
                        if ~(length(c_time.data)==4 && all(c_time.data == [0.5 1 1 c_time.num_samples-3])) 
                            t_data = [c_time.num_samples c_time.data];
                            [shapeIDs(2),found] = obj.shapeLibrary.find_or_insert(t_data);
                            may_exist=may_exist & found;
                        end
                    end
                    data = [amplitude shapeIDs event.delay event.first event.last];
                case 'trap'
                    data = [event.amplitude event.riseTime ...
                            event.flatTime event.fallTime ...
                            event.delay];
                otherwise
                    error('unknown grdient type passed to registerGradEvent()');
            end
            if may_exist
                id = obj.gradLibrary.find_or_insert(data,event.type(1));
            else
                id = obj.gradLibrary.insert(0,data,event.type(1));
            end
        end
        
        function id=registerAdcEvent(obj, event)
            data = [event.numSamples event.dwell max(event.delay,event.deadTime) ... % MZ: replaced event.delay+event.deadTime with a max(...) because we allow for overlap of the delay and the dead time
                event.freqOffset event.phaseOffset event.deadTime];
            id = obj.adcLibrary.find_or_insert(data);
        end
        
        function id=registerControlEvent(obj, event)
            event_type=find(strcmp(event.type,{'output','trigger'}));
            if (event_type==1)
                event_channel=find(strcmp(event.channel,{'osc0','osc1','ext1'})); % trigger codes supported by the Siemens interpreter as of May 2019
            elseif (event_type==2)
                event_channel=find(strcmp(event.channel,{'physio1','physio2'})); % trigger codes supported by the Siemens interpreter as of June 2019
            else
                error('unsupported control event type');
            end
            data = [event_type event_channel event.delay event.duration];
            id = obj.trigLibrary.find_or_insert(data);
        end
        
        function id=registerLabelEvent(obj, event)
            label_id=find(strcmp(event.label,mr.getSupportedLabels()));
            data=[event.value label_id];
            switch event.type
                case 'labelset'
                    id = obj.labelsetLibrary.find_or_insert(data);
                case 'labelinc'
                    id = obj.labelincLibrary.find_or_insert(data);
                otherwise
                    error('unknown label type passed to registerLabelEvent()');                    
            end
        end
        
        %TODO: Replacing blocks in the middle of sequence can cause unused
        %events in the libraries. These can be detected and pruned.
        function setBlock(obj, index, varargin)
            %setBlock Replace or add sequence block.
            %   setBlock(obj, index, bStruct) Replace block at index with new
            %   block provided as block structure.
            %
            %   setBlock(obj, index, e1, e2, ...) Create a new block from
            %   events and store at position given by index.
            %
            %   The block or events are provided in uncompressed form and
            %   will be stored in the compressed, non-redundant internal
            %   libraries.
            %
            %   See also  getBlock, addBlock
            
            %block_duration = mr.calcDuration(varargin);% don't seem to be needed
            
            % Convert block structure to cell array of events
            varargin=mr.block2events(varargin);    
            
            obj.blockEvents{index}=zeros(1,7);
            duration = 0;
            
            check_g = {}; % cell-array containing a structure, each with the index and pairs of gradients/times
            extensions = [];
            
            % Loop over events adding to library if necessary and creating
            % block event structure.
            for i = 1:length(varargin)
                event = varargin{i};
                switch event.type
                    case 'rf'
                        if isfield(event,'id')
                            id=event.id;
                        else
                            id = obj.registerRfEvent(event);
                        end
                        obj.blockEvents{index}(2) = id;
                        duration = max(duration, event.shape_dur + event.delay + event.ringdownTime);
                    case 'grad'
                        channelNum = find(strcmp(event.channel, ...
                                                 {'x', 'y', 'z'}));
                        idx = 2 + channelNum;
                                        
                        grad_start = event.delay + floor(event.tt(1)/obj.gradRasterTime+1e-10)*obj.gradRasterTime;
                        grad_duration = event.delay + ceil(event.tt(end)/obj.gradRasterTime-1e-10)*obj.gradRasterTime;
                        
                        check_g{channelNum}.idx = idx;
                        check_g{channelNum}.start = [grad_start, event.first];
                        check_g{channelNum}.stop  = [grad_duration, event.last]; 
                        

                        if isfield(event,'id')
                            id=event.id;
                        else
                            id = obj.registerGradEvent(event);
                        end
                        obj.blockEvents{index}(idx) = id;
                        duration = max(duration, grad_duration);

                    case 'trap'
                        channelNum = find(strcmp(event.channel,{'x','y','z'}));
                        
                        idx = 2 + channelNum;
                        
                        check_g{channelNum}.idx = idx;
                        check_g{channelNum}.start = [0, 0];
                        check_g{channelNum}.stop  = [event.delay + ...
                                                     event.riseTime + ...
                                                     event.fallTime + ...
                                                     event.flatTime, 0];
                        
                        if isfield(event,'id')
                            id=event.id;
                        else
                            id=obj.registerGradEvent(event);
                        end
                        obj.blockEvents{index}(idx)=id;
                        duration=max(duration,event.delay+event.riseTime+event.flatTime+event.fallTime);

                    case 'adc'
                        if isfield(event,'id')
                            id=event.id;
                        else
                            id=obj.registerAdcEvent(event);
                        end
                        obj.blockEvents{index}(6)=id;
                        duration=max(duration,event.delay+event.numSamples*event.dwell+event.deadTime);
                    case 'delay' 
                        %if isfield(event,'id')
                        %    id=event.id;
                        %else
                        %    id = obj.registerDelayEvent(event);
                        %end
                        %obj.blockEvents{index}(1)=id;
                        % delay is not a true event any more so we account
                        % for the duration but do not add anything
                        duration=max(duration,event.delay);
                    case {'output','trigger'} 
                        if isfield(event,'id')
                            id=event.id;
                        else
                            id=obj.registerControlEvent(event);
                        end
                        %obj.blockEvents{index}(7)=id; % now we just
                        % collect the list of extension objects and we will
                        % add it to the event table later
                        % ext=struct('type', 1, 'ref', id);
                        ext=struct('type', obj.getExtensionTypeID('TRIGGERS'), 'ref', id);
                        extensions=[extensions ext];
                        duration=max(duration,event.delay+event.duration);
                    case {'labelset','labelinc'}
                        if isfield(event,'id')
                            id=event.id;
                        else
                            id=obj.registerLabelEvent(event);
                        end
% %                         label_id=find(strcmp(event.label,mr.getSupportedLabels()));
% %                         data=[event.value label_id];
% %                         [id,found] = obj.labelsetLibrary.find(data);
% %                         if ~found
% %                             obj.labelsetLibrary.insert(id,data);
% %                         end
                        
                        % collect the list of extension objects and we will
                        % add it to the event table later
                        %ext=struct('type', 2, 'ref', id);
                        ext=struct('type', obj.getExtensionTypeID(upper(event.type)), 'ref', id);
                        extensions=[extensions ext];
                end
            end
            
            if ~isempty(extensions)
                % add extensions now... but it's tricky actually
                % we need to check whether the exactly the same list if
                % extensions already exists, otherwise we have to create a
                % new one... ooops, we have a potential problem with the 
                % key mapping then... The trick is that we rely on the
                % sorting of the extension IDs and then we can always find
                % the last one in the list by setting the reference to the
                % next to 0 and then proceed with the otehr elements.
                [~,I]=sort([extensions(:).ref]);
                extensions=extensions(I);
                all_found=true;
                id=0;
                for i=1:length(extensions)
                    data=[extensions(i).type extensions(i).ref id];
                    [id,found] = obj.extensionLibrary.find(data);
                    all_found = all_found && found;
                    if ~found
                        break;
                    end 
                end
                if ~all_found
                    % add the list
                    id=0;
                    for i=1:length(extensions)
                        data=[extensions(i).type extensions(i).ref id];
                        [id,found] = obj.extensionLibrary.find(data);
                        if ~found
                            obj.extensionLibrary.insert(id,data);
                        end 
                    end
                end
                % now we add the ID
                obj.blockEvents{index}(7)=id;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% PERFORM GRADIENT CHECKS                                 %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % check if connection to the previous block is correct
%             check_g
            for cg_temp = check_g
                cg=cg_temp{1}; % cg_temp is still a cell-array with a single element here...
                if isempty(cg) 
                    continue
                end
                % check the start 
                %if event.delay ~= 0 && event.first ~= 0 
                %if cg.start(1) ~= 0 && cg.start(2) ~= 0 
                %    error('No delay allowed for gradients which start with a non-zero amplitude.');
                %end
        
                if abs(cg.start(2)) > obj.sys.maxSlew * obj.sys.gradRasterTime % MZ: we only need the following check if the current gradient starts at non-0
                    if cg.start(1) ~= 0
                        error('No delay allowed for gradients which start with a non-zero amplitude.');
                    end
                    if index > 1
                        prev_id = obj.blockEvents{index-1}(cg.idx);
                        if prev_id ~= 0
                            prev_lib = obj.gradLibrary.get(prev_id);
                            prev_dat = prev_lib.data;
                            prev_type = prev_lib.type;
                            if prev_type == 't'
                                error('Two consecutive gradients need to have the same amplitude at the connection point');
                            elseif prev_type == 'g'
                                last = prev_dat(6); % '6' means last
                                if abs(last - cg.start(2)) > obj.sys.maxSlew * obj.sys.gradRasterTime
                                    error('Two consecutive gradients need to have the same amplitude at the connection point');
                                end
                            end
                        end
                    else                   
                        error('First gradient in the the first block has to start at 0.');
                    end
                end
                
                % Check if gradients, which do not end at 0, are as long as the
                % block itself.
                %assert(abs(duration-block_duration)<eps); % TODO: if this never fails we should remove mr.calcDuration at the beginning
                if cg.stop(2) > obj.sys.maxSlew * obj.sys.gradRasterTime && abs(cg.stop(1)-duration) > 1e-7
                    error('A gradient that doesn''t end at zero needs to be aligned to the block boundary.');
                end
            end
       
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% GRADIENT CHECKS DONE                                    %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            %assert(abs(duration-block_duration)<eps); % TODO: if this never fails we should remove mr.calcDuration at the beginning
            obj.blockDurations(index)=duration;
        end
        
        function block = getBlock(obj, index)
            %getBlock Return a block of the sequence.
            %   b=getBlock(obj, index) Return the block specified by the
            %   index.
            %
            %   The block is created from the sequence data with all
            %   events and shapes decompressed.
            %
            %   See also  setBlock, addBlock
            
            block=struct('blockDuration', 0, 'rf', {}, 'gx', {}, 'gy', {}, 'gz', {}, 'adc', {} );

            block(1).rf = [];
            eventInd = obj.blockEvents{index};
            
            if eventInd(7) > 0
                % we have extensions -- triggers, labels, etc
                % we will eventually isolate this into a separate function
                nextExtID=eventInd(7);
                while nextExtID~=0
                    extData = obj.extensionLibrary.data(nextExtID).array;
                    % format: extType, extID, nextExtID
                    extTypeStr=obj.getExtensionTypeString(extData(1));
                    switch extTypeStr
                        case 'TRIGGERS'
                            trigger_types={'output','trigger'};
                            data = obj.trigLibrary.data(extData(2)).array;
                            trig.type = trigger_types{data(1)};
                            if (data(1)==1)
                                trigger_channels={'osc0','osc1','ext1'};
                                trig.channel=trigger_channels{data(2)};
                            elseif (data(1)==2)
                                trigger_channels={'physio1','physio2'}; 
                                trig.channel=trigger_channels{data(2)};;
                            else
                                error('unsupported trigger event type');
                            end
                            trig.delay = data(3);
                            trig.duration = data(4);
                            % allow for multiple triggers per block
                            if(isfield(block, 'trig'))
                                block.trig(length(block.trig)+1) = trig;
                            else
                                block.trig=trig;
                            end
                        case {'LABELSET','LABELINC'} 
                            label.type=lower(extTypeStr);
                            supported_labels=mr.getSupportedLabels();
                            if strcmp(extTypeStr,'LABELSET')
                                data = obj.labelsetLibrary.data(extData(2)).array;
                            else
                                data = obj.labelincLibrary.data(extData(2)).array;
                            end
                            label.label=supported_labels{data(2)};
                            label.value=data(1);
                            % allow for multiple labels per block
                            if(isfield(block, 'label'))
                                block.label(length(block.label)+1) = label;
                            else
                                block.label=label;
                            end
                        otherwise
                            error('unknown extension ID %d', extData(1));
                    end
                    % now update nextExtID
                    nextExtID=extData(3);
                end
            end
            if eventInd(2) > 0 
                if length(obj.rfLibrary.type)>=eventInd(2)
                    block.rf = obj.rfFromLibData(obj.rfLibrary.data(eventInd(2)).array,obj.rfLibrary.type(eventInd(2)));
                else
                    block.rf = obj.rfFromLibData(obj.rfLibrary.data(eventInd(2)).array); % undefined type/use
                end
            end
            gradChannels = {'gx', 'gy', 'gz'};
            for i = 1:length(gradChannels)
                if eventInd(2+i) > 0
                    type = obj.gradLibrary.type(eventInd(2+i));
                    libData = obj.gradLibrary.data(eventInd(2+i)).array;
                    if type == 't'
                        grad.type = 'trap';
                    else
                        grad.type = 'grad';
                    end
                    grad.channel = gradChannels{i}(2);
                    if strcmp(grad.type,'grad')
                        amplitude = libData(1);
                        shapeId = libData(2);
                        timeId = libData(3);
                        delay = libData(4);
                        shapeData = obj.shapeLibrary.data(shapeId).array;
                        compressed.num_samples = shapeData(1);
                        compressed.data = shapeData(2:end);
                        try
                            g = mr.decompressShape(compressed);
                        catch
                            fprintf('  mr.decompressShape() failed for shapeId %d\n', shapeId);
                            error('mr.decompressShape() failed for shapeId %d', shapeId);
                        end
                        grad.waveform = amplitude*g;
                        % SK: This looks like a bug to me.
%                         grad.t = (1:length(g))'*obj.gradRasterTime;
                        if (timeId==0)
                            grad.tt = ((1:length(g))-0.5)'*obj.gradRasterTime; % TODO: evetually we may remove these true-times
                            t_end=length(g)*obj.gradRasterTime;
                            %grad.t = (0:length(g)-1)'*obj.gradRasterTime;
                        else
                            tShapeData = obj.shapeLibrary.data(timeId).array;
                            compressed.num_samples = tShapeData(1);
                            compressed.data = tShapeData(2:end);
                            try
                                grad.tt = mr.decompressShape(compressed)*obj.gradRasterTime;
                            catch
                                fprintf('  mr.decompressShape() failed for shapeId %d\n', shapeId);
                                error('mr.decompressShape() failed for shapeId %d', shapeId);
                            end
                            assert(length(grad.waveform) == length(grad.tt));
                            t_end=grad.tt(end);                            
                        end
                        grad.shape_id=shapeId; % needed for the second pass of read()
                        grad.time_id=timeId; % needed for the second pass of read()
                        grad.delay = delay;
                        grad.shape_dur = t_end;
                        if length(libData)>5
                            grad.first = libData(5);
                            grad.last = libData(6);
                        else
%			    assert(false); % this should never happen, as now we recover the correct first/last values during reading
%                             % for the data read from a file we need to
%                             % infer the missing fields here
%                             grad.first = grad.waveform(1); % MZ: eventually we should use extrapolation by 1/2 gradient rasters here
%                             grad.last = grad.waveform(end);
%                             % true / extrapolated values
%                             grad.tfirst = (3*g(1)-g(2))*0.5; % extrapolate by 1/2 gradient rasters 
%                             grad.tlast = (g(end)*3-g(end-1))*0.5; % extrapolate by 1/2 gradient rasters
                        end
                    else
                        grad.amplitude = libData(1);
                        grad.riseTime = libData(2);
                        grad.flatTime = libData(3);
                        grad.fallTime = libData(4);
                        grad.delay = libData(5);                        
                        grad.area = grad.amplitude*(grad.flatTime + ...
                                                    grad.riseTime/2 + ...
                                                    grad.fallTime/2);
                        grad.flatArea = grad.amplitude*grad.flatTime;
                    end
                    
                    block.(gradChannels{i}) = grad;
                end
            end
            if eventInd(6) > 0
                libData = obj.adcLibrary.data(eventInd(6)).array;
                if length(libData) < 6
                    libData(end+1) = 0;
                end
                adc = cell2struct(num2cell(libData), ...
                                  {'numSamples', 'dwell', 'delay', ...
                                   'freqOffset', 'phaseOffset', ...
                                   'deadTime'}, 2);
                adc.type = 'adc';
                block.adc = adc;
            end
            block.blockDuration=obj.blockDurations(index);
%             % now that delays in v1.4 and later Pulseq revisions are not
%             % stored, we need to see whether we need a delay to explain the
%             % current block duration
%             if (mr.calcDuration(block)~=obj.blockDurations(index))
%                 tmpDelay.type = 'delay';
%                 tmpDelay.delay = obj.blockDurations(index);
%                 block.delay = tmpDelay;
%             end
        end

        function [ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = calculateKspaceUnfunc(obj, varargin)
            % calculate the k-space trajectory of the entire pulse sequence
            %   optional parameter 'trajectory_delay' sets the compensation
            %   factor to align ADC and gradients in the reconstruction
            %   Return values: ktraj_adc, ktraj, t_excitation, t_refocusing
        
            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'calculateKspace';
                parser.addParamValue('trajectory_delay',0,@(x)(isnumeric(x)));
            end
            parse(parser,varargin{:});
            opt = parser.Results;
            
            if any(abs(opt.trajectory_delay)>100e-6)
                warning('trajectory delay of (%s) us is suspiciously high',num2str(opt.trajectory_delay*1e6));
            end
          
            % initialise the counters and accumulator objects
            c_excitation=0;
            c_refocusing=0;
            c_adcSamples=0;
            % loop throught the blocks to prepare preallocations
            for iB=1:length(obj.blockEvents)
                block = obj.getBlock(iB);
                if ~isempty(block.rf)
                    if (~isfield(block.rf,'use') || strcmp(block.rf.use,'excitation') || strcmp(block.rf.use,'undefined'))
                        c_excitation=c_excitation+1;
                    elseif strcmp(block.rf.use,'refocusing')
                        c_refocusing=c_refocusing+1;
                    end
                end
                if ~isempty(block.adc)
                    c_adcSamples=c_adcSamples+block.adc.numSamples;
                end
            end
            
            %
            t_excitation=zeros(c_excitation,1);
            t_refocusing=zeros(c_refocusing,1);
            ktime=zeros(c_adcSamples,1);
            current_dur=0;
            c_excitation=1;
            c_refocusing=1;
            kcouter=1;
            traj_recon_delay=opt.trajectory_delay;  
            
            % go through the blocks and collect RF and ADC timing data
            for iB=1:length(obj.blockEvents)
                block = obj.getBlock(iB);
                if ~isempty(block.rf)
                    rf=block.rf;
                    t=rf.delay+mr.calcRfCenter(rf);
                    if (~isfield(block.rf,'use') || strcmp(block.rf.use,'excitation') || strcmp(block.rf.use,'undefined'))
                        t_excitation(c_excitation) = current_dur+t;
                        c_excitation=c_excitation+1;
                    elseif strcmp(block.rf.use,'refocusing')
                        t_refocusing(c_refocusing) = current_dur+t;
                        c_refocusing=c_refocusing+1;
                    end
                end
                if ~isempty(block.adc)
                    ktime(kcouter:(kcouter-1+block.adc.numSamples)) = ((0:(block.adc.numSamples-1))+0.5)... % according to the information from Klaus Scheffler and indirectly from Siemens this is the present convention (the samples are shifted by 0.5 dwell)
                        *block.adc.dwell + block.adc.delay + current_dur + traj_recon_delay;
                    kcouter=kcouter+block.adc.numSamples;
                end
                current_dur=current_dur+obj.blockDurations(iB);%mr.calcDuration(block);
            end
            
            % now calculate the actual k-space trajectory based on the
            % gradient waveforms
            gw=obj.gradient_waveforms();
            i_excitation=round(t_excitation/obj.gradRasterTime);
            i_refocusing=round(t_refocusing/obj.gradRasterTime);
%             ii_next_excitation=min(length(i_excitation),1);
%             ii_next_refocusing=min(length(i_refocusing),1);
%             ktraj=zeros(size(gw));
%             k=[0;0;0];
%             % TODO: replace this plain stupid loop with a segment-wise
%             % integration (with segments defined by the RF pulses)
%             for i=1:size(gw,2)
%                 k=k+gw(:,i)*obj.gradRasterTime;
%                 ktraj(:,i)=k;
%                 %if find(i_excitation==i,1)
%                 if ii_next_excitation>0 && i_excitation(ii_next_excitation)==i
%                     k=0;
%                     ktraj(:,i)=NaN; % we use NaN-s to mark the excitation point, they interrupt the plots
%                     ii_next_excitation = min(length(i_excitation),ii_next_excitation+1);
%                 end
%                 %if find(i_refocusing==i,1)
%                 if ii_next_refocusing>0 && i_refocusing(ii_next_refocusing)==i
%                     k=-k;
%                     ii_next_refocusing = min(length(i_refocusing),ii_next_refocusing+1);
%                 end
%             end
            i_periods=sort([1; i_excitation+1; i_refocusing+1; size(gw,2)+1]); % we need thise +1 for compatibility with the above code which prooved to be correct
            ii_next_excitation=min(length(i_excitation),1);
            ii_next_refocusing=min(length(i_refocusing),1);
            ktraj=zeros(size(gw));
            k=[0;0;0];
            for i=1:(length(i_periods)-1)                
                %k=k+gw(:,i)*obj.gradRasterTime;
                i_period_end=(i_periods(i+1)-1);
                % here we use a trick to add current k value to the cumsum()
                k_period=cumsum([k,gw(:,i_periods(i):i_period_end)*obj.gradRasterTime],2);
                ktraj(:,i_periods(i):i_period_end)=k_period(:,2:end); % remove the first 'dummy' sample (see the trick above)
                k=k_period(:,end);
                if ii_next_excitation>0 && i_excitation(ii_next_excitation)==i_period_end
                    k(:)=0;
                    ktraj(:,i_period_end)=NaN; % we use NaN-s to mark the excitation point, they interrupt the plots
                    ii_next_excitation = min(length(i_excitation),ii_next_excitation+1);
                end
                if ii_next_refocusing>0 && i_refocusing(ii_next_refocusing)==i_period_end
                    k=-k;
                    ii_next_refocusing = min(length(i_refocusing),ii_next_refocusing+1);
                end
            end

            % now calculate the k-space positions at the ADC time points
            % sample the k-space positions at the ADC time points
            ktraj_adc=interp1((1:(size(ktraj,2)))*obj.gradRasterTime, ktraj', ktime)';
            t_adc=ktime; % we now also return the sampling time points
        end
        
        function f = plot(obj, varargin)
            %plot Plot the sequence in a new figure.
            %   plot(seqObj) Plot the sequence
            %
            %   plot(...,'TimeRange',[start stop]) Plot the sequence
            %   between the times specified by start and stop.
            %
            %   plot(...,'TimeDisp',unit) Display time in:
            %   's', 'ms' or 'us'.
            %
            %   plot(...,'Label','LIN,REP') Plot label values for ADC events:
            %   in this example for LIN and REP labels; other valid labes are 
            %   accepted as a comma-separated list.
            %
            %   plot(...,'showBlocks',1) Plot grid and tick labels at the
            %   block boundaries. Accebts a numeric or a boolean parameter.
            %
            %   f=plot(...) Return the new figure handle.
            %
            validTimeUnits = {'s','ms','us'};
            validLabel = mr.getSupportedLabels();
            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'plot';
                parser.addParamValue('showBlocks',false,@(x)(isnumeric(x) || islogical(x)));
                parser.addParamValue('timeRange',[0 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('timeDisp',validTimeUnits{1},...
                    @(x) any(validatestring(x,validTimeUnits)));
                parser.addOptional('label',[]);%,@(x)(isstr(x)));%@(x) any(validatestring(x,validLabel))
            end
            parse(parser,varargin{:});
            opt = parser.Results;
            
            fig=figure;
            if nargout>0
                f=fig;
            end
            ax=zeros(1,6);
            for i=1:6
                ax(i)=subplot(3,2,i);
            end
            ax=ax([1 3 5 2 4 6]);   % Re-order axes
            arrayfun(@(x)hold(x,'on'),ax);
            arrayfun(@(x)grid(x,'on'),ax);
            labels={'ADC/labels','RF mag (Hz)','RF/ADC ph (rad)','Gx (kHz/m)','Gy (kHz/m)','Gz (kHz/m)'};
            arrayfun(@(x)ylabel(ax(x),labels{x}),1:6);
            
            tFactorList = [1 1e3 1e6];
            tFactor = tFactorList(strcmp(opt.timeDisp,validTimeUnits));
            xlabel(ax(3),['t (' opt.timeDisp ')']);
            xlabel(ax(6),['t (' opt.timeDisp ')']);
            
            t0=0;
            label_defined=false;
            label_indexes_2plot=[];
            label_legend_2plot=[];
            for i=1:length(validLabel)
                label_store.(validLabel{i})=0;
                if ~isempty(strfind(upper(opt.label),validLabel{i}))
                    label_indexes_2plot=[label_indexes_2plot i];
                    label_legend_2plot=[label_legend_2plot; validLabel{i}];
                end
            end 
            if ~isempty(label_indexes_2plot)
                label_colors_2plot=parula(length(label_indexes_2plot)+1);
                label_colors_2plot=label_colors_2plot(1:end-1,:);
            end
            
            % block timings
            blockEdges=[0 cumsum(obj.blockDurations)];
            blockEdgesInRange=blockEdges(logical((blockEdges>=opt.timeRange(1)).*(blockEdges<=opt.timeRange(2))));
            if (opt.timeDisp=='us')
                for i=1:6
                    xax=get(ax(i),'XAxis');
                    xax.ExponentMode='manual';
                    xax.Exponent=0;
                end
            end
            if opt.showBlocks
                % show block edges in plots
                for i=1:6
                    xax=get(ax(i),'XAxis');
                    xax.TickValues=tFactor.*blockEdgesInRange;
                    set(ax(i),'XTickLabelRotation',90);
                    %xax.MinorTickValues=tFactor.*blockEdgesInRange;
                    %set(ax(i),'XMinorTick', 'on');
                    %set(ax(i),'XMinorGrid', 'on');
                    %set(ax(i),'GridColor',0.8*[1 1 1]); 
                    %set(ax(i),'MinorGridColor',0.6*[1 1 1]); 
                    %set(ax(i),'MinorGridLineStyle','-'); 
                end
            end            
            %for iB=1:size(obj.blockEvents,1)
            for iB=1:length(obj.blockEvents)
                block = obj.getBlock(iB);
                isValid = t0+obj.blockDurations(iB)>opt.timeRange(1) && t0<=opt.timeRange(2);
                if isValid
                    if isfield(block,'label') %current labels, works on the curent or next adc
                        for i=1:length(block.label)
                            if strcmp(block.label(i).type,'labelinc')
                                label_store.(block.label(i).label)=...
                                    label_store.(block.label(i).label)+block.label(i).value;
                            else
                                label_store.(block.label(i).label)=block.label(i).value;
                            end
                        end
                        label_defined=true;
                    end                  
                    if ~isempty(block.adc)
                        adc=block.adc;
                        t=adc.delay + ((0:adc.numSamples-1)+0.5)*adc.dwell; % according to the imformation from Klaus Scheffler and indirectly from Siemens this is the present convention (the samples are shifted by 0.5 dwell)
                        plot(tFactor*(t0+t),zeros(size(t)),'rx','Parent',ax(1));
                        plot(tFactor*(t0+t), angle(exp(1i*adc.phaseOffset).*exp(1i*2*pi*t*adc.freqOffset)),'b.','MarkerSize',1,'Parent',ax(3)); % plot ADC phase
                        if label_defined && ~isempty(label_indexes_2plot)
                            set(ax(1),'ColorOrder',label_colors_2plot);
                            label_store_cell=struct2cell(label_store);
                            lbl_vals=[label_store_cell{label_indexes_2plot}];
                            t=t0+adc.delay + (adc.numSamples-1)/2*adc.dwell;
                            p=plot(tFactor*t,lbl_vals,'.','markersize',5,'Parent',ax(1));
                            if ~isempty(label_legend_2plot)
                                legend(ax(1),p,label_legend_2plot,'location','Northwest','AutoUpdate','off');
                                label_legend_2plot=[];
                            end
                        end
                    end
                    if ~isempty(block.rf)
                        rf=block.rf;
                        [tc,ic]=mr.calcRfCenter(rf);
                        t=rf.t;
                        s=rf.signal;
                        if abs(s(1))~=0 % fix strangely looking phase / amplitude in the beginning
                            s=[0; s];
                            t=[t(1); t];
                            ic=ic+1;
                        end
                        if abs(s(end))~=0 % fix strangely looking phase / amplitude in the beginning
                            s=[s; 0];
                            t=[t; t(end)];
                        end
                        plot(tFactor*(t0+t+rf.delay),  abs(s),'Parent',ax(2));
                        plot(tFactor*(t0+t+rf.delay),  angle(s    *exp(1i*rf.phaseOffset).*exp(1i*2*pi*t    *rf.freqOffset)),...
                             tFactor*(t0+tc+rf.delay), angle(s(ic)*exp(1i*rf.phaseOffset).*exp(1i*2*pi*t(ic)*rf.freqOffset)),'xb',...
                             'Parent',ax(3));
                    end
                    gradChannels={'gx','gy','gz'};
                    for j=1:length(gradChannels)
                        grad=block.(gradChannels{j});
                        if ~isempty(grad)
                            if strcmp(grad.type,'grad')
                                % we extend the shape by adding the first 
                                % and the last points in an effort of 
                                % making the display a bit less confusing...
                                %t=grad.delay + [0; grad.t + (grad.t(2)-grad.t(1))/2; grad.t(end) + grad.t(2)-grad.t(1)];
                                t= grad.delay+[0; grad.tt; grad.shape_dur];
                                waveform=1e-3* [grad.first; grad.waveform; grad.last];
                            else
                                t=cumsum([0 grad.delay grad.riseTime grad.flatTime grad.fallTime]);
                                waveform=1e-3*grad.amplitude*[0 0 1 1 0];
                            end
                            plot(tFactor*(t0+t),waveform,'Parent',ax(3+j));
                        end
                    end                
                end
                t0=t0+obj.blockDurations(iB);%mr.calcDuration(block);
            end
            
            % Set axis limits and zoom properties
            dispRange = tFactor*[opt.timeRange(1) min(opt.timeRange(2),t0)];
            arrayfun(@(x)xlim(x,dispRange),ax);
            linkaxes(ax(:),'x')
            h = zoom(fig);
            setAxesZoomMotion(h,ax(1),'horizontal');
        end
        
        function grad_waveforms=gradient_waveforms1(obj) % currently disfunctional (feature_ExtTrap)
            % gradient_waveforms()
            %   Decompress the entire gradient waveform
            %   Returns an array of gradient_axes x timepoints
            %   gradient_axes is typically 3.
            %
             
            [duration, numBlocks, ~]=obj.duration();
            
            wave_length = ceil(duration / obj.gradRasterTime);
            grad_channels=3;
            grad_waveforms=zeros(grad_channels, wave_length);
            gradChannels={'gx','gy','gz'};
            
            t0=0;
            t0_n=0;
            for iB=1:numBlocks
                block = obj.getBlock(iB);
                for j=1:length(gradChannels)
                    grad=block.(gradChannels{j});
                    if ~isempty(block.(gradChannels{j}))
                        if strcmp(grad.type,'grad')
                            %nt_start=round((grad.delay+grad.t(1))/obj.gradRasterTime);
                            nt_start=round((grad.delay)/obj.gradRasterTime);
                            waveform=grad.waveform;
                        else
                            nt_start=round(grad.delay/obj.gradRasterTime);
                            if (abs(grad.flatTime)>eps) % interp1 gets confused by triangular gradients
                                t=cumsum([0 grad.riseTime grad.flatTime grad.fallTime]);
                                trapform=grad.amplitude*[0 1 1 0];
                            else
                                t=cumsum([0 grad.riseTime grad.fallTime]);
                                trapform=grad.amplitude*[0 1 0];
                            end
                            %
                            tn=floor(t(end)/obj.gradRasterTime);
                            
                            % it turns out that we need an additional zero-
                            % padding at the end otherwise interp1() 
                            % generates NaNs at the end of the shape
                            t=[t t(end)+obj.gradRasterTime];
                            trapform=[trapform 0];
                            
                            %fprintf('%g : %g | ', [t*1e6 ;trapform]);
                            %fprintf('\n');
                            
                            if abs(grad.amplitude)>eps 
                                % MZ: for consistency we change it to the
                                % corresponding mr. function
                                %waveform=interp1(t,trapform,obj.gradRasterTime*(0:tn),'linear');
                                %waveform=interp1(t,trapform,'linear','pp');
                                waveform=mr.pts2waveform(t,trapform,obj.gradRasterTime);
                            else
                                waveform=zeros(1,tn+1);
                            end
                        end
                        if numel(waveform)~=sum(isfinite(waveform(:)))
                            fprintf('Warning: not all elements of the generated waveform are finite!\n');
                        end
                        %plot(tFactor*(t0+t),waveform,'Parent',ax(3+j));
                        grad_waveforms(j,(t0_n+1+nt_start):(t0_n+nt_start+length(waveform)))=waveform;
                    end
                end                

                t0=t0+obj.blockDurations(iB);%mr.calcDuration(block);
                t0_n=round(t0/obj.gradRasterTime);
            end
        end
               
        function [wave_data, tfp_excitation, tfp_refocusing, t_adc, fp_adc]=waveforms_and_times(obj, appendRF)
            % waveforms_and_times()
            %   Decompress the entire gradient waveform
            %   Returns gradient wave forms as a cell array with
            %   gradient_axes (typically 3) dimension; each cell contains
            %   time points and the correspndig gradient amplitude values.
            %   Additional return values are time points of excitations,
            %   refocusings and ADC sampling points.
            %   If the optionl parameter 'appendRF' is set to true the RF 
            %   wave shapes are appended after the gradients
            %   Optional output parameters: tfp_excitation contains time
            %   moments, frequency and phase offsets of the excitation RF
            %   pulses (similar for tfp_refocusing); t_adc contains times 
            %   of all ADC sample points; fp_adc contains frequency and
            %   phase offsets of each ADC object (not sample).
            %   TODO: return RF frequency offsets and RF waveforms and t_preparing (once its available)
            
            if nargin < 2
                appendRF=false;
            end
             
            grad_channels=3;
            gradChannels={'gx','gy','gz'}; % FIXME: this is not OK for matrix gradient systems
            
            t0=0;
            t0_n=0;
            
            numBlocks=length(obj.blockEvents);

            % collect the shape pieces into a cell array
            if appendRF
                shape_channels=length(gradChannels)+1; % the last "channel" is RF
            else
                shape_channels=length(gradChannels);
            end
            shape_pieces=cell(shape_channels,numBlocks); 
            % also collect RF and ADC timing data
            % t_excitation t_refocusing t_adc
            tfp_excitation=[];
            tfp_refocusing=[];
            t_adc=[];
            fp_adc=[];
            %block_durations=zeros(1,numBlocks);
            curr_dur=0;
            out_len=zeros(1,shape_channels); % the last "channel" is RF
            for iB=1:numBlocks
                block = obj.getBlock(iB);
                for j=1:length(gradChannels)
                    grad=block.(gradChannels{j});
                    if ~isempty(block.(gradChannels{j}))
                        if strcmp(grad.type,'grad')
                            % check if we have an extended trapezoid or an arbitrary gradient on a regular raster
                            tt_rast=grad.tt/obj.gradRasterTime+0.5;
                            if all(abs(tt_rast-(1:length(tt_rast)))<eps)
                                % arbitrary gradient
                                %
                                % restore & recompress shape: if we had a
                                % trapezoid converted to shape we have to find
                                % the "corners" and we can eliminate internal
                                % samples on the straight segments
                                % but first we have to restore samples on the
                                % edges of the gradient raster intervals
                                % for that we need the first sample
                                max_abs=max(abs(grad.waveform));
                                odd_step1=[grad.first 2*grad.waveform'];
                                odd_step2=odd_step1.*(mod(1:length(odd_step1),2)*2-1);
                                waveform_odd_rest=(cumsum(odd_step2).*(mod(1:length(odd_step2),2)*2-1))';
                                waveform_odd_interp=[grad.first; 0.5*(grad.waveform(1:end-1)+grad.waveform(2:end)); grad.last];
                                assert(abs(waveform_odd_rest(end)-grad.last)<=2e-5*max_abs,['last restored point of shaped gradient differs too much from the recorded last, deviation: ' num2str(abs(waveform_odd_rest(end)-grad.last)) 'Hz/m (' num2str(abs(waveform_odd_rest(end)-grad.last)/max_abs*100) '%). Block number: ' num2str(iB)]); % what's the reasonable threshold? 
                                %figure; plot([0,10e-6+grad.t'],waveform_odd_rest-waveform_odd_interp);
                                waveform_odd_mask=abs(waveform_odd_rest-waveform_odd_interp)<=eps+2e-5*max_abs; % threshold ???
                                waveform_odd=waveform_odd_interp.*waveform_odd_mask+waveform_odd_rest.*(1-waveform_odd_mask);

                                % combine odd & even
                                comb=[ 0 grad.waveform' ; waveform_odd' ];
                                waveform_os=comb(2:end)';

                                tt_odd=(0:(length(waveform_odd_rest)-1))*obj.gradRasterTime;
                                tt_os=(0:(length(waveform_os)-1))*obj.gradRasterTime*0.5;

                                waveform_even_reint=0.5*(waveform_odd_rest(1:end-1)+waveform_odd_rest(2:end));

                                maskChanges = abs([1; diff(waveform_os,2); 1])>1e-8;   % TRUE if values change
                                waveform_chg = waveform_os(maskChanges);                     % Elements without repetitions
                                tt_chg=tt_os(maskChanges);
                                %figure;plot(grad.tt,grad.waveform);hold on; plot(tt_chg,waveform_chg); plot(tt_chg,waveform_chg,'o');
                                tgc=[tt_chg; waveform_chg'];
                                out_len(j)=size(tgc,2);
                                shape_pieces{j,iB}=curr_dur+grad.delay+tgc;
                            else
                                % extended trapezoid (the easy case!)
                                out_len(j)=out_len(j)+length(grad.tt);
                                shape_pieces{j,iB}=[curr_dur+grad.delay+grad.tt'; grad.waveform'];
                            end
                        else
                            if (abs(grad.flatTime)>eps) % interp1 gets confused by triangular gradients (repeating sample)
                                out_len(j)=out_len(j)+4;
                                shape_pieces{j,iB}=[
                                    curr_dur+grad.delay+cumsum([0 grad.riseTime grad.flatTime grad.fallTime]);...
                                    grad.amplitude*[0 1 1 0]];
                            else
                                if (abs(grad.riseTime)>eps && abs(grad.fallTime)>eps) % we skip 'empty' gradients
                                    out_len(j)=out_len(j)+3;
                                    shape_pieces{j,iB}=[
                                        curr_dur+grad.delay+cumsum([0 grad.riseTime grad.fallTime]);...
                                        grad.amplitude*[0 1 0]];
                                else
                                    if abs(grad.amplitude)>eps
                                        warning('''empty'' gradient with non-zero magnitude detected in block %d',iB);
                                    end
                                end
                            end
                        end
                    end
                end
                if ~isempty(block.rf)
                    rf=block.rf;
                    tc=mr.calcRfCenter(rf);
                    t=rf.delay+tc;
                    if (~isfield(block.rf,'use') || strcmp(block.rf.use,'excitation') || strcmp(block.rf.use,'undefined'))
                        tfp_excitation(:,end+1) = [curr_dur+t; block.rf.freqOffset; block.rf.phaseOffset+block.rf.freqOffset*tc];
                    elseif strcmp(block.rf.use,'refocusing')
                        tfp_refocusing(:,end+1) = [curr_dur+t; block.rf.freqOffset; block.rf.phaseOffset+block.rf.freqOffset*tc];
                    end
                    if appendRF
                        pre=[];
                        post=[];
                        if abs(rf.signal(1))>0
                            pre=[curr_dur+rf.delay+rf.t(1)-eps;0];
                        end
                        if abs(rf.signal(end))>0
                            post=[curr_dur+rf.delay+rf.t(end)+eps;0];
                        end
%                         pre=[curr_dur+rf.delay+rf.t(1)-eps;NaN];
%                         post=[curr_dur+rf.delay+rf.t(end)+eps;NaN];
                        out_len(end)=out_len(j)+length(rf.t)+size(pre,2)+size(post,2);
                        shape_pieces{end,iB}=[pre [curr_dur+rf.delay+rf.t'; (rf.signal.*exp(1i*(rf.phaseOffset+2*pi*rf.freqOffset*rf.t)))'] post];
                    end
                end
                if ~isempty(block.adc)
                    ta=block.adc.dwell*((0:(block.adc.numSamples-1))+0.5); % according to the information from Klaus Scheffler and indirectly from Siemens this is the present convention (the samples are shifted by 0.5 dwell) % according to the information from Klaus Scheffler and indirectly from Siemens this is the present convention (the samples are shifted by 0.5 dwell)
                    t_adc((end+1):(end+block.adc.numSamples)) = ta + block.adc.delay + curr_dur;
                    fp_adc(:,(end+1):(end+block.adc.numSamples)) = [block.adc.freqOffset*ones(1,block.adc.numSamples); block.adc.phaseOffset+block.adc.freqOffset*ta];
                end
                curr_dur=curr_dur+obj.blockDurations(iB);%mr.calcDuration(block);
            end
            
            % collect wave data
            wave_data=cell(1,shape_channels);
            for j=1:shape_channels
                wave_data{j}=zeros(2,out_len(j));
            end
            wave_cnt=zeros(1,shape_channels);
            curr_dur=0;
            for iB=1:numBlocks
                for j=1:shape_channels
                    if ~isempty(shape_pieces{j,iB})
                        wave_data_local=shape_pieces{j,iB};
                        len=size(wave_data_local,2);
                        if wave_cnt(j)~=0 && wave_data{j}(1,wave_cnt(j))+obj.gradRasterTime < wave_data_local(1,1)
                            if  wave_data{j}(2,wave_cnt(j))~=0
                                warning('waveforms_and_times(): forcing ramp-down from a non-zero gradient sample on axis %d at t=%d us \ncheck your sequence, some calculations are probably wrong.', j, round(1e6*wave_data{j}(1,wave_cnt(j))));
                                wave_data{j}(:,wave_cnt(j)+1)=[wave_data{j}(1,wave_cnt(j))+obj.gradRasterTime/2; 0]; % this is likely to cause memory reallocations
                                wave_cnt(j)=wave_cnt(j)+1;
                            end
                            if wave_data_local(2,1)~=0
                                warning('waveforms_and_times(): forcing ramp-up to a non-zero gradient sample on axis %d at t=%d us \ncheck your sequence, some calculations are probably wrong.', j, round(1e6*wave_data_local(1,1)));
                                wave_data_local=[[wave_data_local(1,1)-obj.gradRasterTime/2; 0] wave_data_local]; % this is likely to cause memory reallocations also later on
                                len=len+1;
                            end
                        end
                        if wave_cnt(j)==0 || wave_data{j}(1,wave_cnt(j))~=wave_data_local(1,1)
                            wave_data{j}(:,wave_cnt(j)+(1:len))=wave_data_local;
                            wave_cnt(j)=wave_cnt(j)+len;
                        else
                            wave_data{j}(:,wave_cnt(j)+(1:len-1))=wave_data_local(:,2:end);
                            wave_cnt(j)=wave_cnt(j)+len-1;
                        end
                        if wave_cnt(j)~=length(unique(wave_data{j}(1,1:wave_cnt(j))))
                            warning('Warning: not all elements of the generated time vector are unique!\n');
                        end
                    end
                end
            end
            
            % trim the output data
            for j=1:shape_channels
                if wave_cnt(j)<size(wave_data{j},2)
                    wave_data{j}(:,(wave_cnt(j)+1):end)=[];
                end
            end
            
%             % convert wave data to piecewise polynomials              
%             wave_pp=cell(1,length(gradChannels));
%             for j=1:length(gradChannels)
%                 if (wave_cnt(j)<=0)
%                     continue;
%                 end
%                 if ~all(isfinite(wave_data{j}(:)))
%                    fprintf('Warning: not all elements of the generated waveform are finite!\n');
%                 end
%                 wave_pp{j} = interp1(wave_data{j}(1,1:wave_cnt(j)),wave_data{j}(2,1:wave_cnt(j)),'linear','pp');
%             end            
        end
        
        function [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = calculateKspacePP(obj, varargin)
            % calculate the k-space trajectory of the entire pulse sequence
	        %   using piecewise-polynomial gradient wave representation 
	        %   which is much faster for simple shapes and large delays
            %   optional parameter 'trajectory_delay' sets the compensation
            %   factor to align ADC and gradients in the reconstruction
            %   optional parameter 'gradient_offset' allows to simulate 
            %   background gradients or verifz spin-echo conditions
            %   Return values: ktraj_adc, t_adc, ktraj, t_ktraj,
            %   t_excitation, t_refocusing, slicepos
        
            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'calculateKspacePP';
                parser.addParamValue('trajectory_delay',0,@(x)(isnumeric(x)));
                parser.addParamValue('gradient_offset',0,@(x)(isnumeric(x)));
            end
            parse(parser,varargin{:});
            opt = parser.Results;
            
            if any(abs(opt.trajectory_delay)>100e-6)
                warning('trajectory delay of (%s) us is suspiciously high',num2str(opt.trajectory_delay*1e6));
            end
                      
            total_duration=sum(obj.blockDurations);
            
            [gw_data, tfp_excitation, tfp_refocusing, t_adc]=obj.waveforms_and_times();
            
            ng=length(gw_data);
            % gradient delay handling
            if length(opt.trajectory_delay)==1
                gradient_delays(1:ng)=opt.trajectory_delay;
            else
                assert(length(opt.trajectory_delay)==ng); % we need to have the same number of gradient channels
                gradient_delays=opt.trajectory_delay;
            end
            % gradient offset handling
            if length(opt.gradient_offset)==1
                gradient_offset(1:ng)=opt.gradient_offset;
            else
                assert(length(opt.gradient_offset)==ng); % we need to have the same number of gradient channels
                gradient_offset=opt.gradient_offset;
            end
                        
            % convert wave data to piecewise polynomials              
            gw_pp=cell(1,ng);
            for j=1:ng
                wave_cnt=size(gw_data{j},2);
                if wave_cnt==0 
                    if abs(gradient_offset(j))<=eps % gradient offset support, part 1
                        continue;
                    else
                        gw=[0,total_duration; 0, 0];
                    end
                else
                    gw=gw_data{j};
                end
                % now gw contains the wave form for the current axis
                if abs(gradient_delays(j))>eps 
                    gw(1,:)=gw(1,:)-gradient_delays(j); % (anisotropic) gradient delay support
                end
                if ~all(isfinite(gw(:)))
                   fprintf('Warning: not all elements of the generated waveform are finite!\n');
                end
                teps=1e-12; % eps is too small and may go lost due to rounding errors (e.g. total_duration+eps==total_duration)
                if gw(1,1)>0 && gw(1,end) < total_duration
                    gw=[ [-teps gw(1,1)-teps;0 0] gw [gw(1,end)+teps total_duration+teps;0 0] ]; % we need these "eps" terms to avoid integration errors over extended periods of time
                elseif gw(1,1)>0 
                    gw=[ [-teps gw(1,1)-teps;0 0] gw ]; % we need these "eps" terms to avoid integration errors over extended periods of time
                elseif gw(1,end) < total_duration 
                    gw=[ gw [gw(1,end)+teps total_duration+teps;0 0] ]; % we need these "eps" terms to avoid integration errors over extended periods of time
                end
                %
                if abs(gradient_offset(j))>eps 
                    gw(2,:)=gw(2,:)+gradient_offset(j); % gradient offset support, part 2
                end
                gw_pp{j} = interp1(gw(1,:),gw(2,:),'linear','pp');
            end
            
            % calculate slice positions. for now we entirely rely on the
            % excitation -- ignoring complicated interleaved refocused sequences
            if ~isempty(tfp_excitation)
                slicepos=zeros(length(gw_data),size(tfp_excitation,2)); % position in x,y,z;
                for j=1:length(gw_data)
                    if isempty(gw_pp{j})
                        slicepos(j,:)=NaN;
                    else
                        slicepos(j,:)=tfp_excitation(2,:)./ppval(gw_pp{j},tfp_excitation(1,:));
                    end
                end
                slicepos(~isfinite(slicepos))=0; % reset undefined to 0 (or is NaN better?)
                t_slicepos=tfp_excitation(1,:);
            else
                slicepos=[];
                t_slicepos=[];
            end
            
            %t_adc = t_adc + opt.trajectory_delay;
            % this was wrong because it dod not shift RF events (which are
            % intrinsically well-synchronized) and did not allow for
            % anisotropic delays for different gradient axes. For the new
            % implementation see "gradient_delays" vector above 
            
            % integrate waveforms as PPs to produce gadient moments
            gm_pp=cell(1,ng);
            tc = {};
            for i=1:ng
                if isempty(gw_pp{i})
                    continue;
                end
                gm_pp{i}=fnint(gw_pp{i});
                tc{end+1}=gm_pp{i}.breaks;
                % "sample" ramps for display purposes otherwise piecewise-linear diplay (plot) fails (looks stupid)
                ii=find(abs(gm_pp{i}.coefs(:,1))>eps);
                if ~isempty(ii)
                    tca=cell(1,length(ii));
                    for j=1:length(ii)
                        tca{j}=(floor(gm_pp{i}.breaks(ii(j))/obj.gradRasterTime):ceil(gm_pp{i}.breaks(ii(j)+1)/obj.gradRasterTime))*obj.gradRasterTime;
                    end
                    tc{end+1}=[tca{:}];
                end
            end
            %t = unique([tc{:}, 0, t_excitation-obj.gradRasterTime, t_excitation, t_refocusing, t_adc]);
            % we round to 100ns, otherwise unique() fails...
            
            if isempty(tfp_excitation)
                t_excitation=[];
            else
                t_excitation=tfp_excitation(1,:);
            end
            if isempty(tfp_refocusing)
                t_refocusing=[];
            else
                t_refocusing=tfp_refocusing(1,:);
            end
            
            tacc=1e-10; % temporal accuracy
            taccinv=1/tacc;
            t_ktraj = tacc*unique(round(taccinv*[tc{:}, 0, t_excitation-2*obj.rfRasterTime, t_excitation-obj.rfRasterTime, t_excitation, t_refocusing-obj.rfRasterTime, t_refocusing, t_adc, total_duration]));
            % % the "proper" matlab's function ismember() is slow and returns a bool array, but builtin('_ismemberhelper'...) is not compatible accross versions (known not to work on the windows 2021a version)
            %[~,i_excitation]=builtin('_ismemberhelper',tacc*round(taccinv*t_excitation),t_ktraj);
            %[~,i_refocusing]=builtin('_ismemberhelper',tacc*round(taccinv*t_refocusing),t_ktraj);
            %[~,i_adc]=builtin('_ismemberhelper',tacc*round(taccinv*t_adc),t_ktraj);
            % this is nother undocumented solution, see https://undocumentedmatlab.com/blog_old/ismembc-undocumented-helper-function
            i_excitation=ismembc2(tacc*round(taccinv*t_excitation),t_ktraj);
            i_refocusing=ismembc2(tacc*round(taccinv*t_refocusing),t_ktraj);
            i_adc=ismembc2(tacc*round(taccinv*t_adc),t_ktraj);
            %
            i_periods=unique([1, i_excitation, i_refocusing, length(t_ktraj)]);
            if ~isempty(i_excitation)
                ii_next_excitation=1;
            else
                ii_next_excitation=0;
            end
            if ~isempty(i_refocusing)
                ii_next_refocusing=1;
            else
                ii_next_refocusing=0;
            end            
            ktraj=zeros(3, length(t_ktraj));
            for i=1:ng
                if isempty(gw_pp{i})
                    continue;
                end
                %[~,it]=builtin('_ismemberhelper',[gm_pp{i}.breaks(1),gm_pp{i}.breaks(end)],t);
                %ktraj(i,it(1):it(2))=ppval(gm_pp{i},t(it(1):it(2)));
                it=find(t_ktraj>=tacc*round(taccinv*gm_pp{i}.breaks(1)) & t_ktraj<=tacc*round(taccinv*gm_pp{i}.breaks(end)));
                ktraj(i,it)=ppval(gm_pp{i},t_ktraj(it));
                if t_ktraj(it(end))<t_ktraj(end)
                    ktraj(i,(it(end)+1):end)=ktraj(i,it(end));
                end
            end
            % convert gradient moments to k-space
            dk=-ktraj(:,1);%[0;0;0];
            for i=1:(length(i_periods)-1)
                i_period=i_periods(i);
                i_period_end=i_periods(i+1);
                if ii_next_excitation>0 && i_excitation(ii_next_excitation)==i_period
                    if abs(t_ktraj(i_period)-t_excitation(ii_next_excitation))>tacc
                        warning('abs(t_ktraj(i_period)-t_excitation(ii_next_excitation))<%g failed for ii_next_excitation=%d error=%g', tacc, ii_next_excitation, t_ktraj(i_period)-t_excitation(ii_next_excitation));
                    end
                    dk=-ktraj(:,i_period);
                    if (i_period>1)
                        ktraj(:,i_period-1)=NaN; % we use NaN-s to mark the excitation point, they interrupt the plots
                    end
                    ii_next_excitation = min(length(i_excitation),ii_next_excitation+1);
                elseif ii_next_refocusing>0 && i_refocusing(ii_next_refocusing)==i_period
                    %dk=-ktraj(:,i_period);
                    dk=-2*ktraj(:,i_period)-dk;
                    %if (i_period>1)
                    %    ktraj(:,i_period-1)=NaN; % we use NaN-s to mark the excitation point, they interrupt the plots
                    %end
                    ii_next_refocusing = min(length(i_refocusing),ii_next_refocusing+1);
                end
                
                ktraj(:,i_period:(i_period_end-1))=ktraj(:,i_period:(i_period_end-1))+dk;
            end
            ktraj(:,i_period_end)=ktraj(:,i_period_end)+dk;
            ktraj_adc=ktraj(:,i_adc);
%             % add first and last slicepos
%             if t_slicepos(1)>0
%                 t_slicepos=[0 t_slicepos];
%                 slicepos=[ zeros(length(gw_data),1) slicepos ];
%             end
%             if t_slicepos(end)<t_ktraj(end)
%                 t_slicepos=[t_slicepos t_ktraj(end)];
%                 slicepos=[ slicepos slicepos(:,end)];
%             end
        end
        
        function sound_data=sound(obj)
            %sound()
            %   "play out" the sequence through the system speaker
            %
            
            gw_data=obj.waveforms_and_times();
            total_duration=sum(obj.blockDurations);
            
            sample_rate=44100; %Hz
            dwell_time=1/sample_rate;
            sound_length=floor(total_duration/dwell_time)+1;
            
            sound_data(2,sound_length)=0; %preallocate
            %sound_data(1,:)=interp1((0:(grad_wavelen-1))*obj.gradRasterTime,grad_waveforms(1,:)+0.5*grad_waveforms(3,:),(0:(sound_length-1))*dwell_time);
            sound_data(1,:)=interp1(gw_data{1}(1,:),gw_data{1}(2,:),(0:(sound_length-1))*dwell_time,'linear',0) + ...
                            interp1(gw_data{3}(1,:),gw_data{3}(2,:),(0:(sound_length-1))*dwell_time,'linear',0)*0.5;
            %sound_data(2,:)=interp1((0:(grad_wavelen-1))*obj.gradRasterTime,grad_waveforms(2,:)+0.5*grad_waveforms(3,:),(0:(sound_length-1))*dwell_time);
            sound_data(1,:)=interp1(gw_data{2}(1,:),gw_data{2}(2,:),(0:(sound_length-1))*dwell_time,'linear',0) + ...
                            interp1(gw_data{3}(1,:),gw_data{3}(2,:),(0:(sound_length-1))*dwell_time,'linear',0)*0.5;
            
            % filter like we did it in the gradient music project
            %b = fir1(40, 10000/sample_rate);
            %sound_data = filter(b, 1, sound_data,[],2);
            % use Gaussian convolution instead to supress ringing
            gw=gausswin(round(sample_rate/6000)*2+1);
            gw=gw/sum(gw(:));
            sound_data(1,:) = conv(sound_data(1,:), gw, 'same');
            sound_data(2,:) = conv(sound_data(2,:), gw, 'same');
            
            sound_data_max=max(sound_data(:));            
            sound_data = 0.95 * sound_data / sound_data_max;
            
            % info
            fprintf('playing out the sequence waveform, duration %.1gs\n', sound_length*dwell_time);
            
            % play out the sound
            % we have to zero-pad the weveform due to the limitations of
            % matlab-to-sound interface
            sound([zeros(2,sample_rate/2) sound_data zeros(2,sample_rate/2)], sample_rate); 
        end
        
        function ok=install(seq,dest,name)
            %install Install sequence on RANGE system.
            %   install(seq) Install sequence by copying files to Siemens
            %   host and RANGE controller
            %
            %   install(seq,'siemens') Install Siemens Numaris4 files.
            %   install(seq,'siemensNX') Install Siemens NumarisX files.
            %
            
            if nargin<3
                name='external';
            end

            if ispc 
                % windows
                ping_command='ping -w 1000 -n 1';
            elseif isunix || ismac
                % unix-like
                ping_command='ping -q -n -W1 -c1';
            end

            ok = true;
            if any(strcmpi(dest, {'both', 'siemens', 'siemensNX'}))
                seq.write('external.seq.tmp')
                % we assume we are in the internal network but ICE computer
                % on VB and VE has different IPs
                if ~strcmpi(dest, 'siemensNX')
                    ice_ips={'192.168.2.3', '192.168.2.2'};
                    ice_ip=[];
                    for i=1:length(ice_ips)
                        [status, ~] = system([ping_command ' ' ice_ips{i}]);
                        if status == 0
                            ice_ip=ice_ips{i};
                            break;
                        end                
                    end
                else
                    ice_ip='192.168.2.2';
                end
                if isempty(ice_ip)
                    error('Scanner not found, sequence install failed.')
                end
                pulseq_seq_path='/opt/medcom/MriCustomer/seq/pulseq';
                if strcmpi(dest,'siemensNX')
                    pulseq_seq_path='/opt/medcom/MriCustomer/CustomerSeq/pulseq';
                end
                [status, ~] = system(['scp -oBatchMode=yes -oStrictHostKeyChecking=no external.seq.tmp root@' ice_ip ':' pulseq_seq_path '/external_tmp.seq']);
                ok = ok & status == 0;
                if ok
                    system(['ssh -oBatchMode=yes -oStrictHostKeyChecking=no root@' ice_ip ' "chmod a+rw ' pulseq_seq_path '/external_tmp.seq"']);
                    system(['ssh -oBatchMode=yes -oStrictHostKeyChecking=no root@' ice_ip ' "rm -f ' pulseq_seq_path '/' name '.seq"']);
                    system(['ssh -oBatchMode=yes -oStrictHostKeyChecking=no root@' ice_ip ' "mv ' pulseq_seq_path '/external_tmp.seq ' pulseq_seq_path '/' name '.seq"']);
                end
            end
            
            
            if strcmp(dest, 'dropbox')
                seq.write('external.seq.tmp')
                [status, ~] = system(['cp external.seq.tmp ~/Dropbox/inbox/maxim/' datestr(datetime, 'yyyymmddHHMMSS') '.seq']);
                ok = ok & status == 0;
	        end
            
            if ok
                fprintf('Sequence installed as %s.seq\n',name)
            else
                error('Sequence install failed.')
            end
        end
                
        function codes = getBinaryCodes(obj)
            %getBinaryCodes Return binary codes for section headers in
            %   in a binary sequence file.
            %
            %   See also  writeBinary

            codes.fileHeader = [1 'pulseq' 2];
            codes.version_major = int64(obj.version_major);
            codes.version_minor = int64(obj.version_minor);
            codes.version_revision = int64(obj.version_revision);
            prefix = bitshift(int64(hex2dec('FFFFFFFF')), 32);
            codes.section.definitions = bitor(prefix, int64(1));
            codes.section.blocks      = bitor(prefix, int64(2));
            codes.section.rf          = bitor(prefix, int64(3));
            codes.section.gradients   = bitor(prefix, int64(4));
            codes.section.trapezoids  = bitor(prefix, int64(5));
            codes.section.adc         = bitor(prefix, int64(6));
            codes.section.delays      = bitor(prefix, int64(7));
            codes.section.shapes      = bitor(prefix, int64(8));
        end
        
        function id = getExtensionTypeID(obj, str) 
            % get numeric ID for the given string extention ID
            % will automatically create a new ID if unknown 
            num=find(strcmp(obj.extensionStringIDs,str));
            if isempty(num)
                if isempty(obj.extensionNumericIDs)
                    id=1;
                else
                    id=1+max(obj.extensionNumericIDs);
                end
                obj.extensionNumericIDs(1+length(obj.extensionNumericIDs))=id;
                obj.extensionStringIDs{1+length(obj.extensionStringIDs)}=str;
                assert(length(obj.extensionNumericIDs)==length(obj.extensionStringIDs));
            else
                id=obj.extensionNumericIDs(num);
            end
        end
        
        function str = getExtensionTypeString(obj, id)
            % get numeric ID for the given string extention ID
            % may fail 
            num=find(obj.extensionNumericIDs==id);
            if isempty(num)
                error(['Extension for the given ID ' num2str(id) ' is unknown']);
            end
            str=obj.extensionStringIDs{num};
        end
        
        function setExtensionStringAndID(obj, str, id) 
            % set numeric ID for the given string extention ID
            % may fail if not unique
            if any(strcmp(obj.extensionStringIDs,str)) || any(obj.extensionNumericIDs==id) 
                error('Numeric or String ID is not unique');
            end
            obj.extensionNumericIDs(1+length(obj.extensionNumericIDs))=id;
            obj.extensionStringIDs{1+length(obj.extensionStringIDs)}=str;
            assert(length(obj.extensionNumericIDs)==length(obj.extensionStringIDs))
        end
    end
end % classdef
