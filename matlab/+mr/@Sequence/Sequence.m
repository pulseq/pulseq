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
        softDelayLibrary; % Library of 'soft delay' extension events.
        softDelayHints1;  % Map of string hints that are the part of the 'soft delay' extension objects
        softDelayHints2;  % cell array of string hints that are the part of the 'soft delay' extension objects
        extensionStringIDs;  % string IDs of the used extensions (cell array)
        extensionNumericIDs; % numeric IDs of the used extensions (numeric array)

        gradCheckData;    % struct caching date used for checking of extended gradients cthat cross block boundaries
        
        signatureType; % type of the hashing function used, currently 'md5'
        signatureFile; % which data were hashed, currently 'text' or 'bin' (used file format of the save function)
        signatureValue; % the hash of the exported Pulse sequence

        rfID2NameMap;   % optional names of objects in the plot
        adcID2NameMap;  % optional names of objects in the plot
        gradID2NameMap; % optional names of objects in the plot
        
        sys;
    end
    
    methods
        
        function obj = Sequence(varargin)
            [obj.version_major, obj.version_minor, obj.version_revision] = mr.aux.version();
            % version minor 3 will now support control events (8th column in the event table) mv4 supports/expects timing vectors for arbitrary grads
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
            obj.softDelayLibrary = mr.EventLibrary();
            obj.softDelayHints1 = containers.Map();
            obj.softDelayHints2 = {};
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
            obj.rfID2NameMap = containers.Map('KeyType', 'int32', 'ValueType', 'char'); 
            obj.adcID2NameMap = containers.Map('KeyType', 'int32', 'ValueType', 'char'); 
            obj.gradID2NameMap = containers.Map('KeyType', 'int32', 'ValueType', 'char'); 
            
            obj.gradCheckData=struct('validForBlockNum',0,'lastGradVals', [0 0 0]);

        end
        
        
        % See read.m
        read(obj,filename,varargin)
        
        % See write.m
        write(obj,filename,create_signature)
        
        % See write_v141.m
        write_v141(obj,filename,create_signature)
        
        % See write.m
        write_file(obj,filename)
        
        % See readBinary.m
        readBinary(obj,filename);
        
        % See writeBinary.m
        writeBinary(obj,filename);
        
        
        % See calcPNS.m
        [ok, pns_norm, pns_comp, t_axis]=calcPNS(obj,hardware,doPlots,calcCNS)

        %see calcMomentsBtensor.m
        [B, m1, m2, m3] = calcMomentsBtensor(obj, calcB, calcm1, calcm2, Ndummy, calcm3)
        
        % See testReport.m
        [ report ] = testReport( obj, varargin )
        
        function [duration, numBlocks, eventCount]=duration(obj)
            % duration() 
            %     Returns the total duration of the sequence
            %     optionally returns the total count of events
            %
            
            % Loop over blocks and gather statistics
            numBlocks = length(obj.blockEvents);
            if numBlocks>0 && nargout>2
                eventCount=zeros(size(obj.blockEvents{1}));
            end
            duration=0;
            for iB=1:numBlocks
                if nargout>2
                    eventCount = eventCount + (obj.blockEvents{iB}>0);
                end
                duration=duration+obj.blockDurations(iB);
            end
        end
        
        function [is_ok, errorReport]=checkTiming(obj)
            % checkTiming() 
            %     Checks timing (and some other parameters) of all blocks 
            %     and objects in the sequence optionally returns a detailed
            %     error log as cell array of strings. This function also 
            %     modifies the sequence object by adding the field 
            %     "TotalDuration" to sequence definitions
            %
            
            % Loop over blocks and gather statistics
            numBlocks = length(obj.blockEvents);
            is_ok=true;
            errorReport={};
            totalDuration=0;
            gradBook=struct();
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
                
                % check the stored block duration
                if abs(dur-obj.blockDurations(iB))>eps
                    rep = [rep ' inconsistency between the stored block duration and the duration of the block content'];
                    is_ok = false;
                    dur=obj.blockDurations(iB);
                end
                
                % check that block duration is aligned to the blockDurationRaster
                bd=obj.blockDurations(iB)/obj.blockDurationRaster;
                bdr=round(bd);
                if abs(bdr-bd)>=1e-6
                    rep = [rep ' block duration is not aligned to the blockDurationRaster'];
                    is_ok = false;
                end
                
                % check RF dead times
                if ~isempty(b.rf)
                    if b.rf.delay-b.rf.deadTime < -eps
                        rep = [rep ' delay of ' num2str(b.rf.delay*1e6) 'us is smaller than the RF dead time ' num2str(b.rf.deadTime*1e6) 'us'];
                        is_ok = false;
                    end
                    if b.rf.delay+b.rf.t(end)+b.rf.ringdownTime-dur > eps
                        rep = [rep ' time between the end of the RF pulse at ' num2str((b.rf.delay+b.rf.t(end))*1e6) ' and the end of the block at ' num2str(dur*1e6) 'us is shorter than rfRingdownTime'];
                        is_ok = false;
                    end
                end
                
                % check ADC dead times, dwell times and numbers of samples
                if ~isempty(b.adc) 
                    if b.adc.delay-obj.sys.adcDeadTime < -eps
                        rep = [rep ' adc.delay<system.adcDeadTime'];
                        is_ok=false;
                    end
                    if b.adc.delay+b.adc.numSamples*b.adc.dwell+obj.sys.adcDeadTime-dur > eps
                        rep = [rep ' adc: system.adcDeadTime (post-adc) violation'];
                        is_ok=false;
                    end
                    if abs(b.adc.dwell/obj.sys.adcRasterTime-round(b.adc.dwell/obj.sys.adcRasterTime)) > 1e-10 % the check against eps was too strict 
                        rep = [rep ' adc: dwell time is not an integer multiple of sys.adcRasterTime'];
                        is_ok=false;
                    end
                    if abs(b.adc.numSamples/obj.sys.adcSamplesDivisor-round(b.adc.numSamples/obj.sys.adcSamplesDivisor)) > eps
                        rep = [rep ' adc: numSamples is not an integer multiple of sys.adcSamplesDivisor'];
                        is_ok=false;
                    end
                end

                % check shaped gradients that may potentially end/start at non-zero values
                gradBookCurr=struct();
                if ~isempty(ev) && iscell(ev)
                    for en=1:length(ev)
                        if length(ev{en})==1 && isstruct(ev{en}) && strcmp(ev{en}.type,'grad') % length(ev{en})==1 excludes arrays of extensions 
                            g=ev{en};
                            if g.first~=0 
                                if g.delay~=0
                                    errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' ' g.channel ' gradient starts at a non-zero value but defines a delay\n' ] };
                                end
                                if ~isfield(gradBook, g.channel) || gradBook.(g.channel)~=g.first
                                    errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' ' g.channel ' gradient''s start value ' num2str(g.first) ' differs from the previous block end value\n' ] };
                                else
                                    gradBook.(g.channel)=0; % reset as properly consumed
                                end
                            end
                            if g.last~=0 
                                if g.delay+g.shape_dur ~= dur
                                    errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' ' g.channel ' gradient ends at a non-zero value but does not last until the end of the block\n' ] };
                                end
                                gradBookCurr.(g.channel)=g.last; % update bookkeeping
                            end
                        end
                    end
                end

                % check soft delays
                if isfield(b, 'softDelay') && ~isempty(b.softDelay)
                    if ~exist('softDelayState','var')
                        softDelayState={};
                    end
                    if b.softDelay.factor==0
                        errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' soft delay ' b.softDelay.hint '/' num2str(b.softDelay.num) ' has factor parameter of 0 which is invalid\n' ] };
                        is_ok=false;
                    end
                    % calculate the default delay value based on the current block duration
                    def_del=(obj.blockDurations(iB)-b.softDelay.offset)*b.softDelay.factor;
                    if (b.softDelay.num>=0)
                        % remember or check for consistency
                        if length(softDelayState)<b.softDelay.num+1 || isempty(softDelayState{b.softDelay.num+1})
                            softDelayState{b.softDelay.num+1}=struct('def',def_del,'hint',b.softDelay.hint, 'blk', iB);
                        else
                            if abs(def_del-softDelayState{b.softDelay.num+1}.def)>1e-7 % what is the reasonable threshold?
                                errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' soft delay ' b.softDelay.hint '/' num2str(b.softDelay.num) ': default duration derived from this block (' num2str(def_del*1e6) 'us) is inconsistent with the previous default (' num2str(softDelayState{b.softDelay.num+1}.def*1e6) 'us) that was derived from block ' num2str(softDelayState{b.softDelay.num+1}.blk) '\n' ] };
                                is_ok=false;
                            end
                            if ~strcmp(b.softDelay.hint, softDelayState{b.softDelay.num+1}.hint) 
                                errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' soft delay ' b.softDelay.hint '/' num2str(b.softDelay.num) ': soft delays with the same numeric ID are expected to share the same text hint but previous hint recorded in block ' num2str(softDelayState{b.softDelay.num+1}.blk) ' is ' softDelayState{b.softDelay.num+1}.hint '\n' ] };
                                is_ok=false;
                            end
                        end
                    else
                        errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' contains a soft delay ' b.softDelay.hint ' with an invalid numeric ID' num2str(b.softDelay.num) '\n' ] };
                        is_ok=false;
                    end
                end

                % check whether all gradient bookkeeping values have been properly consumed
                if dur~=0 
                    if any(0~=struct2array(gradBook))
                        errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' some gradients in the previous non-empty block are ending at non-zero values but are not continued here\n' ] };
                    end
                    gradBook=gradBookCurr;
                end

                % update report
                if ~isempty(rep)
                    errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' ' rep '\n' ] };
                end
                %
                totalDuration = totalDuration+dur;
            end
            
            % check whether all gradients in the last block are ramped down properly
            if ~isempty(ev) && iscell(ev)
                for en=1:length(ev)
                    if length(ev{en})==1 && isstruct(ev{en}) && strcmp(ev{en}.type,'grad') % length(ev{en})==1 excludes arrays of extensions 
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
            %   addBlock(obj, duration, e1, e2, ...) Create a new block
            %   with the given predefined duration populated with events
            %   e1, e2, etc. If the duration of any of the events exceeds
            %   the desired duration an error will be thrown.
            %
            %   See also  setBlock, makeAdc, makeTrapezoid, makeSincPulse
            %setBlock(obj,size(obj.blockEvents,1)+1,varargin{:});
            setBlock(obj,length(obj.blockEvents)+1,varargin{:});            
        end
        
        function iB=findBlockByTime(obj,t)
            if nargin<3
                nonEmpty=true;
            end
            iB=find(cumsum(obj.blockDurations)>t,1); 
            if iB>length(obj.blockDurations);
                iB=[]; % or length(obj.blockDurations)
            end
            assert(obj.blockDurations(iB)>0);
            %if nonEmpty && ~isempty(iB)
            %    iB=find(obj.blockDurations(1:iB)>0,1,'last');
            %end
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
                    obj.gradLibrary.data(selectedEvents(i)).array(2)=modifier*obj.gradLibrary.data(selectedEvents(i)).array(2);  % change in v150
                    obj.gradLibrary.data(selectedEvents(i)).array(3)=modifier*obj.gradLibrary.data(selectedEvents(i)).array(3);  % change in v150
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

            rf.center = libData(5); % new in v150
            rf.delay = libData(6); % changed in v150
            rf.freqPPM = libData(7); % changed in v150
            rf.phasePPM = libData(8); % new changed v150
            rf.freqOffset = libData(9); % changed in v150
            rf.phaseOffset = libData(10); % new changed v150
            
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

            if nargin<=2
                error('Parameter ''use'' is not optional since v1.5.0');
            end
            %TODO: fixme : use map built from mr.getSupportedRfUse();
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
        end
        
        function [id shapeIDs]=registerRfEvent(obj, event)
            % registerRfEvent : Add the event to the libraries (object,
            % shapes, etc and return the event's ID. This ID should be 
            % stored in the object to accelerate addBlock()
            
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

            if isfield(event,'use')
                % todo: fixme: use map from getSupportedRfUse
                switch event.use
                    case {'excitation','refocusing','inversion','saturation','preparation','other'}
                        use = event.use(1);
                    otherwise
                        if strcmp(event.use,'u')
                            event.use='undefined'; % make it little more user-friendly
                        end
                        warning('Unknown or undefined RF pulse parameter ''use''=%s, which is not optional since v1.5.0',event.use);
                        use = 'u'; % undefined
                end
            else
                error('Parameter ''use'' is not optional since v1.5.0');
            end

            data = [amplitude shapeIDs(1) shapeIDs(2) shapeIDs(3) ...
                    event.center event.delay event.freqPPM event.phasePPM event.freqOffset event.phaseOffset ];%...
                    %event.deadTime event.ringdownTime];
            if may_exist
                id = obj.rfLibrary.find_or_insert(data,use);
            else
                id = obj.rfLibrary.insert(0,data,use);
            end

            if isfield(event,'name')
                obj.rfID2NameMap(id) = event.name; 
            end
        end
        
        function [id,shapeIDs]=registerGradEvent(obj, event)
            % registerGradEvent : Add the event to the libraries (object,
            % shapes, etc and return the event's ID. This ID should be 
            % stored in the object to accelerate addBlock()
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
                        if (length(c_time.data)==4 && all(c_time.data == [0.5 1 1 c_time.num_samples-3])) 
                            % conventional grad on standard raster: shapeID
                            % is readily 0, nothing needs to be done
                            %shapeIDs(2)=0;
                        elseif (length(c_time.data)==3 && all(c_time.data == [0.5 0.5 c_time.num_samples-2])) 
                            % grad on a half-raster (oversampling): shapeID
                            % needs to be set to -1 as a flag for oversampling
                            shapeIDs(2)=-1;
                        else
                            t_data = [c_time.num_samples c_time.data];
                            [shapeIDs(2),found] = obj.shapeLibrary.find_or_insert(t_data);
                            may_exist=may_exist & found;
                        end
                    end
                    data = [amplitude event.first event.last shapeIDs event.delay];
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

            if isfield(event,'name')
                obj.gradID2NameMap(id) = event.name; 
            end
        end
        
        function [id,shapeID]=registerAdcEvent(obj, event)
            % registerAdcEvent : Add the event to the libraries (object,
            % shapes, etc and return the event's ID. This ID should be 
            % stored in the object to accelerate addBlock()

            surely_new=false;
            if isempty(event.phaseModulation)
                shapeID=0;
            else
                if isfield(event,'shapeID')
                    shapeID=event.shapeID;
                else
                    phaseShape = mr.compressShape(event.phaseModulation(:));
                    data = [phaseShape.num_samples phaseShape.data];
                    [shapeID,shape_found] = obj.shapeLibrary.find_or_insert(data);
                    if ~shape_found
                        surely_new=true;
                    end
                end
            end

            data = [event.numSamples event.dwell max(event.delay,event.deadTime), ... % MZ: replaced event.delay+event.deadTime with a max(...) because we allow for overlap of the delay and the dead time
                event.freqPPM event.phasePPM event.freqOffset event.phaseOffset shapeID]; % event.deadTime];
            if surely_new
                id = obj.adcLibrary.insert(0,data);
            else
                id = obj.adcLibrary.find_or_insert(data);
            end

            if isfield(event,'name')
                obj.adcID2NameMap(id) = event.name;
            end
        end
        
        function id=registerControlEvent(obj, event)
            % registerControlEvent : Add the event to the libraries (object,
            % shapes, etc and return the event's ID. This ID should be 
            % stored in the object to accelerate addBlock()
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
            % registerLabelEvent : Add the event to the libraries (object,
            % shapes, etc and return the event's ID. This ID should be 
            % stored in the object to accelerate addBlock()
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
        
        function id=registerSoftDelayEvent(obj, event)
            % registerDeleyEvent : Add the event to the libraries (object,
            % shapes, etc and return the event's ID. This ID should be 
            % stored in the object to accelerate addBlock()
            try
                hintID=obj.softDelayHints1(event.hint);
            catch
                hintID=obj.softDelayHints1.length()+1;
                obj.softDelayHints1(event.hint)=hintID;
                obj.softDelayHints2{hintID}=event.hint;
            end
            data = [event.num event.offset event.factor hintID];
            id = obj.softDelayLibrary.find_or_insert(data);
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
            %   setBlock(obj, index, duration, e1, e2, ...) Create a new 
            %   block with the given predefined duration populated with
            %   events e1, e2, etc. and store at position given by index.
            %   If the duration of any of the events exceeds the desired
            %   duration an error will be thrown.
            %
            %   The block or events are provided in uncompressed form and
            %   will be stored in the compressed, non-redundant internal
            %   libraries.
            %
            %   See also  getBlock, addBlock
            
            % Convert block structure to cell array of events
            varargin=mr.block2events(varargin);

            obj.blockEvents{index}=zeros(1,7);
            duration = 0;
            
            check_g = cell(1,3); % cell-array containing a structure, each with the index and pairs of gradients/times
            extensions = [];
            required_duration=[];
            
            % Loop over events adding to library if necessary and creating
            % block event structure.
            for i = 1:length(varargin)
                event = varargin{i};
                if isstruct(event)
                    switch event(1).type % we accept multiple extensions and one of the possibilities is an array of extensions
                        case 'rf'
                            if isfield(event,'id')
                                obj.blockEvents{index}(2)=event.id;
                            else
                                obj.blockEvents{index}(2) = obj.registerRfEvent(event);
                            end
                            duration = max(duration, event.shape_dur + event.delay + event.ringdownTime);
                        case 'grad'
                            channelNum = find(strcmp(event.channel, ...
                                                     {'x', 'y', 'z'}));

                            idx = 2 + channelNum;                            
                            grad_duration = event.delay + ceil(event.tt(end)/obj.gradRasterTime-1e-10)*obj.gradRasterTime;

                            grad_start = event.delay + floor(event.tt(1)/obj.gradRasterTime+1e-10)*obj.gradRasterTime;                                
                            
                            %check_g{channelNum}.idx = idx;
                            check_g{channelNum}.start = [grad_start, event.first];
                            check_g{channelNum}.stop  = [grad_duration, event.last]; 
                            
    
                            if isfield(event,'id')
                                obj.blockEvents{index}(idx) = event.id;
                            else
                                obj.blockEvents{index}(idx) = obj.registerGradEvent(event);
                            end
                            duration = max(duration, grad_duration);
    
                        case 'trap'
                            channelNum = find(strcmp(event.channel,{'x','y','z'}));
                            

                            idx = 2 + channelNum;

                            % MZ: the checks as implemented below only make sense for non-trapezoid gradients, commented out the coe below
                            % %check_g{channelNum}.idx = idx;
                            % check_g{channelNum}.start = [0, 0];
                            % check_g{channelNum}.stop  = [event.delay + ...
                            %                              event.riseTime + ...
                            %                              event.fallTime + ...
                            %                              event.flatTime, 0];
                            
                            if isfield(event,'id')
                                obj.blockEvents{index}(idx) = event.id;
                            else
                                obj.blockEvents{index}(idx) = obj.registerGradEvent(event);
                            end
                            duration=max(duration,event.delay+event.riseTime+event.flatTime+event.fallTime);
    
                        case 'adc'
                            if isfield(event,'id')
                                obj.blockEvents{index}(6) = event.id;
                            else
                                obj.blockEvents{index}(6) = obj.registerAdcEvent(event);
                            end
                            duration=max(duration,event.delay+event.numSamples*event.dwell+event.deadTime); % adcDeadTime is added after the sampling period (mr.makeADC also adds a delay before the actual sampling if it was shorter)
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
                            for e=event % allow multiple extensions as an array
                                if isfield(e,'id')
                                    id=e.id;
                                else
                                    id=obj.registerControlEvent(e);
                                end
                                %obj.blockEvents{index}(7)=id; % now we just
                                % collect the list of extension objects and we will
                                % add it to the event table later
                                % ext=struct('type', 1, 'ref', id);
                                ext=struct('type', obj.getExtensionTypeID('TRIGGERS'), 'ref', id);
                                extensions=[extensions ext];
                                duration=max(duration,e.delay+e.duration);
                            end
                        case {'labelset','labelinc'}
                            for e=event % allow multiple extensions as an array
                                if isfield(e,'id')
                                    id=e.id;
                                else
                                    id=obj.registerLabelEvent(e);
                                end
        % %                         label_id=find(strcmp(e.label,mr.getSupportedLabels()));
        % %                         data=[e.value label_id];
        % %                         [id,found] = obj.labelsetLibrary.find(data);
        % %                         if ~found
        % %                             obj.labelsetLibrary.insert(id,data);
        % %                         end
                                
                                % collect the list of extension objects and we will
                                % add it to the event table later
                                %ext=struct('type', 2, 'ref', id);
                                ext=struct('type', obj.getExtensionTypeID(upper(e.type)), 'ref', id);
                                extensions=[extensions ext];
                            end
                        case 'softDelay'
                            if isfield(event,'id')
                                id=event.id;
                            else
                                id=obj.registerSoftDelayEvent(event);
                            end
                            % collect the list of extension objects and
                            % add it to the event table later
                            % ext=struct('type', 1, 'ref', id);
                            ext=struct('type', obj.getExtensionTypeID('DELAYS'), 'ref', id);
                            extensions=[extensions ext];
                        otherwise
                            error('Attempting to add an unknown event to the block.');
                    end
                else
                    if isnumeric(event)
                        % interpret the single numeric parameter as a
                        % requested duration, but throw an error if
                        % multiple numbers are provided
                        if isempty(required_duration)
                            required_duration=event;
                        else
                            error('More than one numeric parameter given to setBlock()');
                        end
                    end
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
                % sanity checks for the softDelay
                nSoftDelays=sum([extensions(:).type]==obj.getExtensionTypeID('DELAYS'));
                if nSoftDelays
                    if nSoftDelays>1
                        error('Only one ''softDelay'' extension event can be added per block');
                    end
                    if duration==0 && isempty(required_duration)
                        error('Soft delay extension can only be used in conjunstion with blocks of non-zero duration'); % otherwise the gradient checks get tedious
                    end
                    if any(obj.blockEvents{index}(2:6)~=0)
                        error('Soft delay extension can only be used in empty blocks (blocks containing no conventional events such as RF, adc or gradients).')
                    end
                end
                % now we add the ID
                obj.blockEvents{index}(7)=id;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% PERFORM GRADIENT CHECKS                                 %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % see if we have the valid preceding gradient value data 
            %gradCheckData % struct('validForBlockNum',0,'lastGradVals', [0 0 0]);
            if duration>0
                if index > 1 && obj.gradCheckData.validForBlockNum ~= index-1
                    % need to update gradCheckData
                    obj.gradCheckData.validForBlockNum = index-1;
                    obj.gradCheckData.lastGradVals(:)=0;
                    [~,prev_nonempty_block]=find(obj.blockDurations(1:(index-1))>0, 1, 'last');
                    if ~isempty(prev_nonempty_block)
                        for i= 1:length(obj.gradCheckData.lastGradVals) % TODO: MZ: check this with external gradient channels !!!
                            prev_id = obj.blockEvents{prev_nonempty_block}(i+2); % careful! direct eventLib access
                            if prev_id ~= 0
                                prev_lib = obj.gradLibrary.get(prev_id); % MZ: for performance reasons we access the gradient library directly. I know, this is not elegant 
                                prev_dat = prev_lib.data;
                                prev_type = prev_lib.type;
                                if prev_type == 'g'
                                    obj.gradCheckData.lastGradVals(i) = prev_dat(3);  % change in v150 - now it is '3' % '6' means last; MZ: I know, this is a real hack...
                                end
                            end
                        end
                    end
                end
                % check if connection to the previous block is correct using check_g
                for i= 1:3 %length(check_g) % TODO: MZ: check this with external gradient channels !!!
                    cg=check_g{i}; % cg_temp is still a cell-array with a single element here...
                    % connection to the previous block in case of extended or shaped gradients
                    if isempty(cg)
                        if abs(obj.gradCheckData.lastGradVals(i)) > obj.sys.maxSlew * obj.sys.gradRasterTime
                            error('Error in block %d on gradient axis %d: previous block ended with non-zero amplitude but the current block has no compatible gradient.', index, i);
                        end
                        % update the gradCheckData.lastGradVals(i)
                        obj.gradCheckData.lastGradVals(i)=0;
                        continue; 
                    end
                    
                    % check the start 
                    if abs(cg.start(2)) > obj.sys.maxSlew * obj.sys.gradRasterTime % MZ: we only need the following check if the current gradient starts at non-0
                        if cg.start(1) ~= 0
                            error('Error in block %d: No delay allowed for gradients which start with a non-zero amplitude', index);
                        end
                        if index > 1
                            if abs(obj.gradCheckData.lastGradVals(i) - cg.start(2)) > obj.sys.maxSlew * obj.sys.gradRasterTime
                                error('Error in block %d on gradient axis %d: Two consecutive gradients need to have the same amplitude at the connection point', index, i);
                            end
                        else                   
                            error('First gradient in the the first block has to start at 0.');
                        end
                    end
                    
                    % Check if gradients, which do not end at 0, are as long as the block itself.
                    if cg.stop(2) > obj.sys.maxSlew * obj.sys.gradRasterTime && abs(cg.stop(1)-duration) > 1e-7
                        error('Error in block %d: A gradient that doesn''t end at zero needs to be aligned to the block boundary', index);
                    end
    
                    % update the gradCheckData.lastGradVals(i)
                    obj.gradCheckData.lastGradVals(i)=cg.stop(2);
                end
            end
            % finish updating gradCheckData (if current block duration is 0 we simply update the validity indicator)
            obj.gradCheckData.validForBlockNum = index;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% GRADIENT CHECKS DONE                                    %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            if ~isempty(required_duration)
                if duration-required_duration>eps
                    error('Required block duration is %g s but actuall block duration is %g s', required_duration, duration);
                end
                duration=required_duration;
            end
                
            obj.blockDurations(index)=duration;
        end

        function raw_block = getRawBlockContentIDs(obj, index)
            %getRawBlockContentIDs Return a block content of the sequence.
            %   b=getRawBlockContentIDs(obj, index) Return the block 
            %   content IDs specified by the index of the block.
            %
            %   No block events are created, only the IDs of the objects 
            %   are returned.
            %
            %   See also  getBlock, setBlock, addBlock
            
            raw_block=struct('blockDuration', 0, 'rf', [], 'gx', [], 'gy', [], 'gz', [], 'adc', [], 'ext', [] );

            eventInd = obj.blockEvents{index};
            
            if eventInd(7) > 0
                % we have extensions -- triggers, labels, etc
                % first find how many and preallocate raw_block.ext array
                nextExtID=eventInd(7);
                cExt=0;
                while nextExtID~=0
                    cExt=cExt+1;
                    % now update nextExtID
                    nextExtID=obj.extensionLibrary.data(nextExtID).array(3);
                end
                raw_block.ext=zeros(2,cExt);
                % now scan through the list again and extract the extension data
                nextExtID=eventInd(7);
                cExt=0;
                while nextExtID~=0
                    cExt=cExt+1;
                    extData = obj.extensionLibrary.data(nextExtID).array;
                    % format: extType, extID, nextExtID
                    raw_block.ext(:,cExt)=extData(1:2);                    
                    % now update nextExtID
                    nextExtID=extData(3);
                end
            end
            if eventInd(2) > 0 
                raw_block.rf=eventInd(2);
            end
            gradChannels = {'gx', 'gy', 'gz'};
            for i = 1:length(gradChannels)
                if eventInd(2+i) > 0
                    raw_block.(gradChannels{i})=eventInd(2+i);
                end
            end
            if eventInd(6) > 0
                raw_block.adc = eventInd(6);
            end
        end        
        
        function block = getBlock(obj, index, addIDs)
            %getBlock Return a block of the sequence.
            %   b=getBlock(obj, index) Return the block specified by the
            %   index.
            %
            %   The block is created from the sequence data with all
            %   events and shapes decompressed.
            %
            %   See also  setBlock, addBlock
            
            if nargin < 3
                addIDs=false;
            end

            block=struct('blockDuration', 0, 'rf', [], 'gx', [], 'gy', [], 'gz', [], 'adc', [] );

            %block(1).rf = [];
            %eventInd = obj.blockEvents{index};
            raw_block = obj.getRawBlockContentIDs(index);
            
            if ~isempty(raw_block.ext)
                % we have extensions -- triggers, labels, etc
                % ext field format: extType, extID

                % unpack trigger(s)
                trig_ext=raw_block.ext(2,raw_block.ext(1,:)==obj.getExtensionTypeID('TRIGGERS'));
                if ~isempty(trig_ext)
                    trigger_types={'output','trigger'};
                    for i=length(trig_ext):-1:1 % backwards for preallocation
                        data = obj.trigLibrary.data(trig_ext(i)).array;
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
                        if addIDs
                            trig.id=trig_ext(i);
                        end
                        % we allow for multiple triggers per block
                        block.trig(i)=trig;
                    end
                end
                % unpack labels
                lid_set=obj.getExtensionTypeID('LABELSET');
                lid_inc=obj.getExtensionTypeID('LABELINC');
                supported_labels=mr.getSupportedLabels();
                label_ext=raw_block.ext(:,raw_block.ext(1,:)==lid_set | raw_block.ext(1,:)==lid_inc);
                if ~isempty(label_ext)
                    for i=size(label_ext,2):-1:1 % backwards for preallocation
                        if label_ext(1,i)==lid_set
                            label.type='labelset';
                            data = obj.labelsetLibrary.data(label_ext(2,i)).array;
                        else
                            label.type='labelinc';
                            data = obj.labelincLibrary.data(label_ext(2,i)).array;
                        end
                        label.label=supported_labels{data(2)};
                        label.value=data(1);
                        if addIDs
                            label.id=label_ext(2,i);
                        end
                        block.label(i) = label;
                    end
                end
                % unpack delay
                delay_ext=raw_block.ext(:,raw_block.ext(1,:)==obj.getExtensionTypeID('DELAYS'));
                if ~isempty(delay_ext)
                    if size(delay_ext,2)>1
                        error('Only one soft delay extension object per block is allowed');
                    end
                    data = obj.softDelayLibrary.data(delay_ext(2,1)).array;
                    if addIDs
                        block.softDelay=struct('type','softDelay','num',data(1),'offset',data(2),'factor',data(3),'hint',obj.softDelayHints2{data(4)},'id',delay_ext(2,1));
                    else
                        block.softDelay=struct('type','softDelay','num',data(1),'offset',data(2),'factor',data(3),'hint',obj.softDelayHints2{data(4)});
                    end
                end
                if length(trig_ext)+size(label_ext,2)~=size(raw_block.ext,2)
                    for i=1:size(raw_block.ext,2)
                        if raw_block.ext(1,i)~=obj.getExtensionTypeID('TRIGGERS') && ...
                           raw_block.ext(1,i)~=obj.getExtensionTypeID('LABELSET') && ...
                           raw_block.ext(1,i)~=obj.getExtensionTypeID('LABELSET') && ...
                           raw_block.ext(1,i)~=obj.getExtensionTypeID('DELAYS')
                            warning('unknown extension ID %d', raw_block.ext(1,i));
                        end
                    end
                end
            end
            % finished extensions

            if ~isempty(raw_block.rf) 
                if length(obj.rfLibrary.type)>=raw_block.rf
                    block.rf = obj.rfFromLibData(obj.rfLibrary.data(raw_block.rf).array,obj.rfLibrary.type(raw_block.rf));
                else
                    block.rf = obj.rfFromLibData(obj.rfLibrary.data(raw_block.rf).array); % undefined type/use
                end
                if addIDs
                    block.rf.id=raw_block.rf;
                    % this is a bit of a hack because we have to access the rfLibrary directly
                    block.rf.shapeIDs=obj.rfLibrary.data(raw_block.rf).array(2:4); % ampl_shape_id phase_shape_id time_shape_id
                end
            end
            gradChannels = {'gx', 'gy', 'gz'};
            for i = 1:length(gradChannels)
                gid=raw_block.(gradChannels{i});
                if ~isempty(gid)
                    type = obj.gradLibrary.type(gid);
                    libData = obj.gradLibrary.data(gid).array;
                    if type == 't'
                        grad.type = 'trap';
                    else
                        grad.type = 'grad';
                    end
                    grad.channel = gradChannels{i}(2);
                    if strcmp(grad.type,'grad')
                        amplitude = libData(1);
                        shapeId = libData(4); % change in v150
                        timeId = libData(5);  % change in v150
                        delay = libData(6);   % change in v150
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
                        elseif (timeId==-1)
                            % gradient with oversampling by a factor of 2
                            grad.tt = ((1:length(g)))'/2*obj.gradRasterTime; 
                            assert(length(grad.tt)==length(grad.waveform));
                            assert(mod(length(g),2)==1);
                            t_end=(length(g)+1)/2*obj.gradRasterTime;
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
                        grad.first = libData(2); % change in v150 - we always have first/last now
                        grad.last = libData(3);  % change in v150 - we always have first/last now                       
                        if addIDs
                            grad.shapeIDs = [shapeId timeId];
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
                    if addIDs
                        grad.id=gid;
                    end
                    block.(gradChannels{i}) = grad;
                end
            end
            if ~isempty(raw_block.adc) 
                libData = obj.adcLibrary.data(raw_block.adc).array;
                shapeIdPhaseModulation=libData(end);
                if shapeIdPhaseModulation                    
                    shapeData = obj.shapeLibrary.data(shapeIdPhaseModulation).array;
                    compressed.num_samples = shapeData(1);
                    compressed.data = shapeData(2:end);
                    try
                        phaseShape = mr.decompressShape(compressed);
                    catch
                        error('mr.decompressShape() failed for shapeId %d', shapeIdPhaseModulation);
                    end
                else
                    phaseShape=0; % wee need a 0 trick here because [] did not work
                end
                adc = cell2struct(num2cell([libData(1:end-1) 0 obj.sys.adcDeadTime]), ...
                                  {'numSamples', 'dwell', 'delay', ...
                                   'freqPPM', 'phasePPM', 'freqOffset', 'phaseOffset', ...
                                   'phaseModulation','deadTime'}, 2);
                if shapeIdPhaseModulation
                    adc.phaseModulation=phaseShape;
                else
                    adc.phaseModulation=[]; % replace 0 with an empty array
                end
                adc.type = 'adc';
                if addIDs
                    adc.id=raw_block.adc;
                end
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
        
        function labels = evalLabels(obj, varargin)
            %Evaluate Label values of the entire sequence or its part
            %   evalLabels(seqObj) Returns the label values at the end of 
            %   the sequence. Return value of the function is the structure
            %   'labels' with fields named after the labels used in the
            %   sequence. Only the fiels corresponding to the lables
            %   actually used are created.
            %
            %   evalLabels(...,'blockRange',[first last]) Evaluate label
            %   values starting from the first specified block to the last
            %   one. 
            %
            %   evalLabels(...,'init',labels_struct) Evaluate labels
            %   assuming the initial values from 'labels_struct'. Useful if
            %   evaluating labes block-by-block.
            %
            %   evalLabels(...,'evolution',<flag>) Evaluate labels and
            %   return the evolution depending on the provided <flag>,
            %   which is one of 'none','adc','label','blocks': 'blocks'
            %   means return label values for all blocks; 'adc' - only for
            %   blocks containig ADC objects; 'label' - only for blocks
            %   where labels are manipulated.
            %
            validEvolutionValues = {'none','adc','label','blocks'};
            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'evalLabels';
                parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('init',struct([]),@(x)(isempty(x) || isstruct(x)));
                parser.addParamValue('evolution','none',@(x)any(validatestring(x,validEvolutionValues)));
            end
            parse(parser,varargin{:});
            opt = parser.Results;

            if isempty(opt.init)
                labels=struct();
            else
                labels=opt.init;
            end

            if ~strcmp(opt.evolution,'none')
                label_evol={};
            end

            if ~isfinite(opt.blockRange(2))
                opt.blockRange(2)=length(obj.blockEvents);
            end

            for iB=opt.blockRange(1):opt.blockRange(2)
                block = obj.getBlock(iB);
                if isfield(block,'label') %current block has labels
                    for i=1:length(block.label)
                        if strcmp(block.label(i).type,'labelinc')
                            if ~isfield(labels,block.label(i).label)
                                labels.(block.label(i).label)=0;
                            end
                            labels.(block.label(i).label)=...
                                labels.(block.label(i).label)+block.label(i).value;
                        else
                            labels.(block.label(i).label)=block.label(i).value;
                        end
                    end
                    if strcmp(opt.evolution,'label')
                        label_evol{end+1}=labels;
                    end
                end
                if strcmp(opt.evolution,'blocks') || ...
                   (strcmp(opt.evolution,'adc') && ~isempty(block.adc))
                    label_evol{end+1}=labels;
                end
            end
            n=length(label_evol);
            if n>1 && ~strcmp(opt.evolution,'none')
                % convert cell array of structures to a structure of arrays
                % %l = cell2mat(label_evol); 
                % l = [label_evol{:}]; % step1: convert to array of structures
                % f = fields(label_store);
                % a = cell(2,length(f)); % step2: prepare argumet 
                % for i=1:length(f)
                %     a{1,i}=f{i};
                %     a{2,i}=[l.(f{i})];
                % end
                % label_store=struct(a{:}); % step3: create the final structure
                f = fields(labels);
                for i=1:length(f)
                    labels.(f{i})=zeros(1,n);
                    for j=1:n
                        if isfield(label_evol{j},f{i})
                            labels.(f{i})(j)=label_evol{j}.(f{i});
                        end
                    end
                end
            end
        end
        
        function sp = plot(obj, varargin)
            %plot Plot the sequence in a new figure.
            %   plot(seqObj) Plot the sequence
            %
            %   plot(...,'timeRange',[start stop]) Plot the sequence
            %   between the times specified by start and stop.
            %
            %   plot(...,'blockRange',[first last]) Plot the sequence
            %   starting from the first specified block to the last one.
            %
            %   plot(...,'timeDisp',unit) Display time in:
            %   's', 'ms' or 'us'.
            %
            %   plot(...,'label','LIN,REP') Plot label values for ADC events:
            %   in this example for LIN and REP labels; other valid labes are 
            %   accepted as a comma-separated list.
            %
            %   plot(...,'showBlocks',1) Plot grid and tick labels at the
            %   block boundaries. Accepts a numeric or a boolean parameter.
            %
            %   plot(...,'stacked',1) Rearrange the plots such they are vertically
            %   stacked and share the same x-axis. Accepts a numeric or a boolean
            %   parameter.
            %
            %   plot(...,'showGuides',1) How dynamic hairline guides that follow 
            %   the data cursor to help verifying event alignment. Accepts a 
            %   numeric or a boolean parameter.
            % 
            %   f=plot(...) Return the new figure handle.
            %

            if nargout == 1
                sp = mr.aux.SeqPlot(obj, varargin{:});
            else
                mr.aux.SeqPlot(obj, varargin{:});
            end

        end

        function sp = paperPlot(obj, varargin)
            %paperPlot Plot the sequence in a stzle similar to that used in 
            %          scientific papers.
            %   paperPlot(seqObj) Plot the sequence
            %
            %   paperPlot(...,'blockRange',[first last]) Plot the sequence
            %   starting from the first specified block to the last one.
            %
            %   paperPlot(...,'lineWidth', w) Plot the sequence
            %   using the specified line width.
            %
            %   paperPlot(...,'axesColor', w) Plot the sequence
            %   using the specified color for the horisontal axes.
            %
            %   paperPlot(...,'rfColor', w) Plot the sequence
            %   using the specified color for the RF and ADC events.
            %
            %   paperPlot(...,'gxColor', w) Plot the sequence
            %   using the specified color for the X gradients.
            %
            %   paperPlot(...,'gyColor', w) Plot the sequence
            %   using the specified color for the Y gradients.
            %
            %   paperPlot(...,'gzColor', w) Plot the sequence
            %   using the specified color for the Z gradients.
            %
            %   paperPlot(...,'rfPlot', <'abs', 'real', 'imag'>) Plot the
            %   RF pulses as the magnitude or real or imaginary part.
            %
            %   Color parameters can be provided as a common color names,
            %   e.g. 'red', 'blue', 'black', character strings starting
            %   from '#' followed by a hexadecimal RGB values ranging from
            %   00 to ff or as an 1x3 vector of doubles ranging from 0 to 1
            %   containing RGB values. 
            %
            %   f=paperPlot(...) Return the new figure handle.
            %
            
            function c=my_validatecolor(c)
                try 
                    c=validatecolor(c); 
                catch 
                    c=[]; 
                end
            end
            validRfPlotValues = {'abs','real','imag'};
            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'paperPlot';
                parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('lineWidth',1.2,@(x)(isnumeric(x)));
                parser.addParamValue('axesColor',[0.5 0.5 0.5],@(x)~isempty(validatecolor(x)));
                parser.addParamValue('rfColor','black',@(x)~isempty(my_validatecolor(x)));
                parser.addParamValue('gxColor','blue',@(x)~isempty(my_validatecolor(x)));
                parser.addParamValue('gyColor','red',@(x)~isempty(my_validatecolor(x)));
                parser.addParamValue('gzColor',[0 0.5 0.3],@(x)~isempty(my_validatecolor(x)));            
                parser.addParamValue('rfPlot','abs',@(x)any(validatestring(x,validPlotRfValues)));
            end
            parse(parser,varargin{:});
            opt = parser.Results;
            
            if mr.aux.isOctave()
              warning('Function paperPlot() does not (yet) work on Octave.');
              return;
            end

            lw=opt.lineWidth;
            axes_clr=opt.axesColor;

            blockRange=opt.blockRange;
            if ~isfinite(blockRange(2))
                blockRange(2)=length(obj.blockDurations);
            end
            
            [wave_data,~,~,t_adc]=obj.waveforms_and_times(true,blockRange); % also export RF
            
            gwm=max(abs([wave_data{1:3}]'));
            rfm=max(abs([wave_data{4}]'));
            gwm(1)=max(gwm(1),t_adc(end));
            
            % remove horizontal lines with 0s. we detect 0 0 and insert a NaN in between
            for i=1:4
                j=size(wave_data{i},2);
                while j>1
                    if wave_data{i}(2,j)==0 && wave_data{i}(2,j-1)==0 
                        wave_data{i}(:,j:end+1)=[ [0.5*(wave_data{i}(1,j-1)+wave_data{i}(1,j));NaN] wave_data{i}(:,j:end)];
                    end
                    j=j-1;
                end
            end
            
            f=figure; 
            %f=colordef(f,'white'); %Set color scheme
            
            f.Color='w'; %Set background color of figure window
            
            t = tiledlayout(4,1,'TileSpacing','none');
            ax=[];
            
            nexttile
            % plot the 'axis'
            plot([-0.01*gwm(1),1.01*gwm(1)],[0 0],'Color',axes_clr,'LineWidth',lw/5); hold on; 
            % plot the RF waveform
            %wave_data{4}(2,wave_data{4}(2,:)==0)=NaN; % hide 0s
            switch opt.rfPlot
                case 'real'
                    plot(wave_data{4}(1,:), real(wave_data{4}(2,:)),'Color',opt.rfColor,'LineWidth',lw); 
                case 'imag'
                    plot(wave_data{4}(1,:), imag(wave_data{4}(2,:)),'Color',opt.rfColor,'LineWidth',lw); 
                otherwise
                    plot(wave_data{4}(1,:), abs(wave_data{4}(2,:)),'Color',opt.rfColor,'LineWidth',lw); 
            end
            
            % plot ADCs
            t_adc_x3=repmat(t_adc,[3 1]);
            y_adc_x3=repmat([0; rfm(2)/5; NaN],[1 length(t_adc)]);
            plot(t_adc_x3(:),y_adc_x3(:),'Color',opt.rfColor,'LineWidth',lw/4);
            
            xlim([-0.03*gwm(1),1.03*gwm(1)]);
            ylim([-1.03*rfm(2),1.03*rfm(2)]);
            set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
            set(get(gca, 'XAxis'), 'Visible', 'off');
            set(get(gca, 'YAxis'), 'Visible', 'off');
            ax(end+1)=gca;
            
            nexttile
            % plot the 'axis'
            plot([-0.01*gwm(1),1.01*gwm(1)],[0 0],'Color',axes_clr,'LineWidth',lw/5); hold on; 
            % plot the entire gradient waveforms
            plot(wave_data{3}(1,:), wave_data{3}(2,:),'Color',opt.gzColor,'LineWidth',lw); 
            
            xlim([-0.03*gwm(1),1.03*gwm(1)]);
            ylim([-1.03*gwm(2),1.03*gwm(2)]);
            set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
            set(get(gca, 'XAxis'), 'Visible', 'off');
            set(get(gca, 'YAxis'), 'Visible', 'off');
            ax(end+1)=gca;
            
            nexttile
            % plot the 'axis'
            plot([-0.01*gwm(1),1.01*gwm(1)],[0 0],'Color',axes_clr,'LineWidth',lw/5); hold on; 
            % plot the entire gradient waveforms
            plot(wave_data{2}(1,:), wave_data{2}(2,:),'Color',opt.gyColor,'LineWidth',lw);
            
            xlim([-0.03*gwm(1),1.03*gwm(1)]);
            ylim([-1.03*gwm(2),1.03*gwm(2)]);
            set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
            set(get(gca, 'XAxis'), 'Visible', 'off');
            set(get(gca, 'YAxis'), 'Visible', 'off');
            ax(end+1)=gca;
            
            nexttile
            % plot the 'axis'
            plot([-0.01*gwm(1),1.01*gwm(1)],[0 0],'Color',axes_clr,'LineWidth',lw/5); hold on; 
            % plot the entire gradient waveforms
            plot(wave_data{1}(1,:), wave_data{1}(2,:),'Color',opt.gxColor,'LineWidth',lw);
            
            xlim([-0.03*gwm(1),1.03*gwm(1)]);
            ylim([-1.03*gwm(2),1.03*gwm(2)]);
            set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
            set(get(gca, 'XAxis'), 'Visible', 'off');
            set(get(gca, 'YAxis'), 'Visible', 'off');
            ax(end+1)=gca;
            
            % link zooming on the time axis
            linkaxes(ax(:),'x')
            
            if nargout == 1
                sp = f;
            end
        end
                       
        function [wave_data, tfp_excitation, tfp_refocusing, t_adc, fp_adc]=waveforms_and_times(obj, appendRF, blockRange)
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
            
            if nargin < 3
                blockRange=[1, length(obj.blockEvents)];
            else
                if length(blockRange)~=2
                    error('parameter ''blockRange'' must contain exactly two numbers: first and last blocks from the range to consider');
                end
            end
            
            if nargin < 2
                appendRF=false;
            end
             
            grad_channels=3;
            gradChannels={'gx','gy','gz'}; % FIXME: this is not OK for matrix gradient systems
            
            t0=0;
            t0_n=0;
            
            numBlocks=blockRange(2)-blockRange(1)+1;

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
            iP=0;
            out_len=zeros(1,shape_channels); % the last "channel" is RF
            for iBc=blockRange(1):blockRange(2)
                block = obj.getBlock(iBc);
                iP=iP+1;
                for j=1:length(gradChannels)
                    grad=block.(gradChannels{j});
                    if ~isempty(block.(gradChannels{j}))
                        if strcmp(grad.type,'grad')
                            % check if we have an extended trapezoid or an arbitrary gradient 
                            % on a regular raster. Arbitrary gradient on a pure centers raster 
                            % (shifted by 0.5) needs special processing
                            tt_rast=grad.tt/obj.gradRasterTime;
                            if all(abs(tt_rast-((1:length(tt_rast))-0.5)')<1e-6)
                                % arbitrary gradient on a centers raster (no oversampling)
                                % restore shape: if we had a trapezoid converted to shape we 
                                % have to find the "corners" and we can eliminate internal
                                % samples on the straight segments but first we have to 
                                % restore samples on the edges of the gradient raster 
                                % intervals - for that we need the first sample
                                
                                [tt_chg, waveform_chg] = mr.restoreAdditionalShapeSamples(grad.tt,grad.waveform,grad.first,grad.last,obj.gradRasterTime,iBc);
                                
                                out_len(j)=out_len(j)+length(tt_chg);
                                shape_pieces{j,iP}=[curr_dur+grad.delay+tt_chg; waveform_chg];%curr_dur+grad.delay+tgc;
                            else
                                % extended trapezoid or sampled gradient with oversampling (the easy case!)
                                % the only caveat is that we need to add the first and last poins to the 
                                % shape in case of the oversampled rasterized gradient
                                if abs(tt_rast(1)-0.5)<1e-6 % rasterized gradient's first sample is always on half-raster, extended trapezoid is always on a raster edje
                                    out_len(j)=out_len(j)+length(grad.tt)+2;
                                    shape_pieces{j,iP}=[curr_dur+grad.delay+[0 grad.tt' grad.shape_dur]; [grad.first grad.waveform' grad.last]];
                                else
                                    out_len(j)=out_len(j)+length(grad.tt);
                                    shape_pieces{j,iP}=[curr_dur+grad.delay+grad.tt'; grad.waveform'];
                                end
                            end
                        else
                            if (abs(grad.flatTime)>eps) % interp1 gets confused by triangular gradients (repeating sample)
                                out_len(j)=out_len(j)+4;
                                shape_pieces{j,iP}=[
                                    curr_dur+grad.delay+cumsum([0 grad.riseTime grad.flatTime grad.fallTime]);...
                                    grad.amplitude*[0 1 1 0]];
                            else
                                if (abs(grad.riseTime)>eps && abs(grad.fallTime)>eps) % we skip 'empty' gradients
                                    out_len(j)=out_len(j)+3;
                                    shape_pieces{j,iP}=[
                                        curr_dur+grad.delay+cumsum([0 grad.riseTime grad.fallTime]);...
                                        grad.amplitude*[0 1 0]];
                                else
                                    if abs(grad.amplitude)>eps
                                        warning('''empty'' gradient with non-zero magnitude detected in block %d',iBc);
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
                    full_freqOffset=rf.freqOffset+rf.freqPPM*1e-6*obj.sys.gamma*obj.sys.B0;
                    full_phaseOffset=rf.phaseOffset+rf.phasePPM*1e-6*obj.sys.gamma*obj.sys.B0;
                    if (~isfield(rf,'use') || strcmp(rf.use,'excitation') || strcmp(rf.use,'undefined'))
                        tfp_excitation(:,end+1) = [curr_dur+t; full_freqOffset; full_phaseOffset+2*pi*full_freqOffset*tc];
                    elseif strcmp(rf.use,'refocusing')
                        tfp_refocusing(:,end+1) = [curr_dur+t; full_freqOffset; full_phaseOffset+2*pi*full_freqOffset*tc];
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
                        shape_pieces{end,iP}=[pre [curr_dur+rf.delay+rf.t'; (rf.signal.*exp(1i*(full_phaseOffset+2*pi*full_freqOffset*rf.t)))'] post];
                    end
                end
                if ~isempty(block.adc)
                    ta=block.adc.dwell*((0:(block.adc.numSamples-1))+0.5); % according to the information from Klaus Scheffler and indirectly from Siemens this is the present convention (the samples are shifted by 0.5 dwell) % according to the information from Klaus Scheffler and indirectly from Siemens this is the present convention (the samples are shifted by 0.5 dwell)
                    t_adc((end+1):(end+block.adc.numSamples)) = ta + block.adc.delay + curr_dur;
                    full_freqOffset=block.adc.freqOffset+block.adc.freqPPM*1e-6*obj.sys.gamma*obj.sys.B0;
                    full_phaseOffset=block.adc.phaseOffset+block.adc.phasePPM*1e-6*obj.sys.gamma*obj.sys.B0;
                    if isempty(block.adc.phaseModulation)
                        block.adc.phaseModulation=0;
                    end                        
                    fp_adc(:,(end+1):(end+block.adc.numSamples)) = [full_freqOffset*ones(1,block.adc.numSamples); full_phaseOffset+block.adc.phaseModulation+full_freqOffset*ta];
                end
                curr_dur=curr_dur+obj.blockDurations(iBc);%mr.calcDuration(block);
            end
            
            % collect wave data
            wave_data=cell(1,shape_channels);
            for j=1:shape_channels
                wave_data{j}=zeros(2,out_len(j));
            end
            wave_cnt=zeros(1,shape_channels);
            curr_dur=0;
            for iP=1:numBlocks
                for j=1:shape_channels
                    if ~isempty(shape_pieces{j,iP})
                        wave_data_local=shape_pieces{j,iP};
                        len=size(wave_data_local,2);
                        if wave_cnt(j)~=0 && wave_data{j}(1,wave_cnt(j))+obj.gradRasterTime < wave_data_local(1,1)
                            if  wave_data{j}(2,wave_cnt(j))~=0
                                if abs(wave_data{j}(2,wave_cnt(j)))>1e-6 % todo: real physical tolarance
                                    warning('waveforms_and_times(): forcing ramp-down from a non-zero gradient sample on axis %d at t=%d us \ncheck your sequence, some calculations are possibly wrong. If using mr.makeArbitraryGrad() consider using explicit values for ''first'' and ''last'' and setting them correctly.', j, round(1e6*wave_data{j}(1,wave_cnt(j))));
                                    wave_data{j}(:,wave_cnt(j)+1)=[wave_data{j}(1,wave_cnt(j))+obj.gradRasterTime/2; 0]; % this is likely to cause memory reallocations
                                    wave_cnt(j)=wave_cnt(j)+1;
                                else
                                    % we are within the tolorance, just set it to 0 quietly
                                    wave_data{j}(2,wave_cnt(j))=0.0;
                                end
                            end
                            if wave_data_local(2,1)~=0
                                if abs(wave_data_local(2,1))>1e-6 % todo: real physical tolarance
                                    warning('waveforms_and_times(): forcing ramp-up to a non-zero gradient sample on axis %d at t=%d us \ncheck your sequence, some calculations are probably wrong.  If using mr.makeArbitraryGrad() consider using explicit values for ''first'' and ''last'' and setting them correctly.', j, round(1e6*wave_data_local(1,1)));
                                    wave_data_local=[[wave_data_local(1,1)-obj.gradRasterTime/2; 0] wave_data_local]; % this is likely to cause memory reallocations also later on
                                    len=len+1;
                                else
                                    % we are wihin the tolorance, just set it to 0 quaietly
                                    wave_data_local(2,1)=0.0;
                                end
                            end
                        end
                        if wave_cnt(j)==0 || wave_data{j}(1,wave_cnt(j))<wave_data_local(1,1)
                            wave_data{j}(:,wave_cnt(j)+(1:len))=wave_data_local;
                            wave_cnt(j)=wave_cnt(j)+len;
                        else
                            if (wave_data_local(1,1)<wave_data{j}(1,wave_cnt(j))-1e-9) % TODO consistent time tolerance
                                warning('Warning: looks like rounding errors for some elements exceed the acceptable tolerance!\n');
                            end
                            [~,d]=find(wave_data_local(1,:)>wave_data{j}(1,wave_cnt(j)),1);
                            wave_data{j}(:,wave_cnt(j)+(1:(len-d+1)))=wave_data_local(:,d:end);
                            wave_cnt(j)=wave_cnt(j)+len-d+1;
                        end
                    end
                end
            end
            for j=1:shape_channels
                if any(diff(wave_data{j}(1,1:wave_cnt(j)))<=0.0) %&& ... % quick pre-check whether the time vector is monotonously increasing to avoid too often unique() calls 
                    %wave_cnt(j)~=length(unique(wave_data{j}(1,1:wave_cnt(j))))
                    warning('Warning: not all elements of the generated time vector are unique and sorted in accending order!\n');
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
        
        function [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos, gw_pp] = calculateKspacePP(obj, varargin)
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
                parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('externalWaveformsAndTimes',struct([]),@(x)(isstruct(x)));
            end
            parse(parser,varargin{:});
            opt = parser.Results;
            
            if any(abs(opt.trajectory_delay)>100e-6)
                warning('trajectory delay of (%s) us is suspiciously high',num2str(opt.trajectory_delay*1e6));
            end
            
            blockRange=opt.blockRange;
            if blockRange(1)<1
                blockRange(1)=1;
            end
            if ~isfinite(blockRange(2))
                blockRange(2)=length(obj.blockDurations);
            end            
                      
            total_duration=sum(obj.blockDurations(blockRange(1):blockRange(2)));
            
            if isempty(opt.externalWaveformsAndTimes)
                [gw_data, tfp_excitation, tfp_refocusing, t_adc]=obj.waveforms_and_times(false,blockRange);
            else
                gw_data=opt.externalWaveformsAndTimes.gw_data;
                tfp_excitation=opt.externalWaveformsAndTimes.tfp_excitation;
                tfp_refocusing=opt.externalWaveformsAndTimes.tfp_refocusing;
                t_adc=opt.externalWaveformsAndTimes.t_adc;
                % how do we verify that the total_duration is correct???
            end
            
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
                if mr.aux.isOctave()
                  gm_pp{i}=ppint(gw_pp{i});
                else
                  gm_pp{i}=fnint(gw_pp{i});
                end
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
            if mr.aux.isOctave()
              [~,i_excitation]=ismember(tacc*round(taccinv*t_excitation),t_ktraj);
              [~,i_refocusing]=ismember(tacc*round(taccinv*t_refocusing),t_ktraj);
              [~,i_adc]=ismember(tacc*round(taccinv*t_adc),t_ktraj);
            else
              i_excitation=ismembc2(tacc*round(taccinv*t_excitation),t_ktraj);
              i_refocusing=ismembc2(tacc*round(taccinv*t_refocusing),t_ktraj);
              i_adc=ismembc2(tacc*round(taccinv*t_adc),t_ktraj);
            end
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

        function [mean_pwr, peak_pwr, rf_rms, total_energy]=calcRfPower(obj, varargin)
            %calcRfPower : Calculate the relative** power of the RF pulse
            %   Returns the (relative) energy of the pulse expressed in the units of 
            %   RF amplitude squared multiplied by time, e.g. in Pulseq these are
            %   Hz * Hz * s = Hz. Sounds strange, but is true. Return
            %   parameter mean_pwr is closely related to the relative** SAR.
            %   Also returns peak power [Hz^2] and RMS B1 amplitude [Hz].
            %   Optional parameter 'blockRange' can be used to specify the
            %   part of the sequence for which the energy is calculated.
            %   Optional parameter 'windowDuration' can be used to specify
            %   the time window for the total_energy, mean_pwr and rf_rms
            %   calculation. The values returned in this case are the
            %   maximum values over all time windows. The time window is
            %   rounded up to a certain number of complete blocks.
            %   ** Note: the power and rf amplitude calculated by this function is 
            %   relative as it is calculated in units of Hz^2 or Hz. The rf amplitude 
            %   can be converted to T by dividing the resulting value by gamma. 
            %   Correspondingly, The power can be converted to mT^2*s by dividing
            %   the given value by gamma^2. Nonetheless, the absolute SAR is related to
            %   the electric field, so the further scaling coeficient is both tx-coil-
            %   dependent (e.g. depends on the coil design) and also subject-dependent 
            %   (e.g. depends on the reference voltage). 

            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'calcRfPower';
                parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('windowDuration',NaN,@(x)(isnumeric(x)));
            end
            parse(parser,varargin{:});
            opt = parser.Results;
            
            if ~isfinite(opt.blockRange(2))
                opt.blockRange(2)=length(obj.blockEvents);
            end
            
            dur=0;
            total_energy=0;
            peak_pwr=0;
            rf_ms=0;

            windowOn= isfinite(opt.windowDuration);
            if windowOn
                blockBookkeeping = zeros(2,opt.blockRange(2)-opt.blockRange(1)+1);
                currentWindowDur=0.0;
                currentWindowStartIdx=opt.blockRange(1);
                total_energy_max=0.0;
                rf_ms_max=0.0;
            end
            
            for iBc=opt.blockRange(1):opt.blockRange(2)
                block = obj.getBlock(iBc);
                dur=dur+obj.blockDurations(iBc);

                if ~isempty(block.rf)
                    rf=block.rf;
                    [e,pp,rms]=mr.calcRfPower(rf);
                    
                    total_energy=total_energy+e;
                    rf_ms=rf_ms+rms^2*rf.shape_dur;
                    peak_pwr=max(peak_pwr,pp);

                    if windowOn
                        blockBookkeeping(:,iBc-opt.blockRange(1)+1)=[e,rms^2*rf.shape_dur];
                        total_energy_max=max(total_energy_max,total_energy);
                        rf_ms_max=max(rf_ms_max,rf_ms);
                    end
                end
                % keep track of the window width and make it shorter if needed
                if windowOn
                    currentWindowDur=currentWindowDur+obj.blockDurations(iBc);
                    while currentWindowDur>opt.windowDuration
                        % remove the front side of the window from total_energy and rf_ms
                        total_energy=total_energy-blockBookkeeping(1,currentWindowStartIdx);
                        rf_ms=rf_ms-blockBookkeeping(2,currentWindowStartIdx);
                        currentWindowDur=currentWindowDur-obj.blockDurations(currentWindowStartIdx);
                        currentWindowStartIdx=currentWindowStartIdx+1;
                    end
                end
            end

            if windowOn
                total_energy=total_energy_max;
                mean_pwr=total_energy/opt.windowDuration; % we divide here by the nominal window duration, which may lead to a slight overestimation
                rf_rms=sqrt(rf_ms_max/opt.windowDuration); % we divide here by the nominal window duration, which may lead to a slight overestimation
            else
                mean_pwr=total_energy/dur;
                rf_rms=sqrt(rf_ms/dur);
            end
        end

        function applySoftDelay(obj, varargin) 
            % applies soft delays to the sequence by modifying the block
            % durations of the respective blocks. Input parameters are
            % pairs of soft delays and values, whereas the soft delay is
            % identified by its string hint and the value is the duration
            % in seconds. Not all soft delays defined in the sequence  need 
            % to be specified. Examples: 
            %
            %    seq.applySoftDelay('TE',40e-3); % set TE to 40ms
            %    seq.applySoftDelay('TE',50e-3,'TR',2); % set TE to 50ms and TR to 2 s
            %
            % See also: mr.makeSoftDelay()

            % parse arguments manually
            sdm_input=containers.Map('KeyType', 'char', 'ValueType', 'double');
            for i=1:2:length(varargin)
                if ~ischar(varargin{i})
                    error('Argument at the position %d must be a character string ID of the soft delay',i);
                end
                if i>length(varargin) || ~isnumeric(varargin{i+1})
                    error('Argument at the position %d must be the value of the soft delay ''%s''',i+1,varargin{i});
                end
                sdm_input(varargin{i})=varargin{i+1};
            end
            % go through all the blocks and update durations, at the same time 
            % checking the consistency of the soft delays 
            sdm_Str2numIDs=containers.Map('KeyType', 'char', 'ValueType', 'double');
            sdm_num2hint=containers.Map('KeyType', 'double', 'ValueType', 'char');
            sdm_warnings=containers.Map('KeyType', 'double', 'ValueType', 'logical');
            for iBc=1:length(obj.blockDurations)
                b = obj.getBlock(iBc);
                if isfield(b, 'softDelay') && ~isempty(b.softDelay)
                    % check the numeric ID consistency 
                    if ~sdm_Str2numIDs.isKey(b.softDelay.hint)
                        sdm_Str2numIDs(b.softDelay.hint)=b.softDelay.num;
                    else
                        if sdm_Str2numIDs(b.softDelay.hint)~=b.softDelay.num
                            error('Soft delay in block %d with numeric ID %d and string hint ''%s'' is inconsistent with the previous occurences of the same string hint', iBc, b.softDelay.num, b.softDelay.hint);
                        end
                    end
                    if ~sdm_num2hint.isKey(b.softDelay.num)
                        sdm_num2hint(b.softDelay.num)=b.softDelay.hint;
                    else
                        if ~strcmp(sdm_num2hint(b.softDelay.num),b.softDelay.hint)
                            error('Soft delay in block %d with numeric ID %d and string hint ''%s'' is inconsistent with the previous occurences of the same numeric ID', iBc, b.softDelay.num, b.softDelay.hint);
                        end
                    end
                    if sdm_input.isKey(b.softDelay.hint)
                        % calculate the new block duration 
                        new_dur_ru=(sdm_input(b.softDelay.hint)/b.softDelay.factor + b.softDelay.offset)/obj.sys.blockDurationRaster;
                        new_dur=round(new_dur_ru)*obj.sys.blockDurationRaster;
                        if abs(new_dur-new_dur_ru*obj.sys.blockDurationRaster)>0.5e-6 && ~sdm_warnings.isKey(b.softDelay.num)                            
                            warning('Block duration for block %d, soft delay ''%s'', had to be substantially rounded to become aligned to the raster time. This warning is only displayed for the first block where it occurs.', iBc, b.softDelay.hint);
                            sdm_warnings(b.softDelay.num)=true;
                        end
                        if new_dur<0 
                            error('Calculated new duration of the block %i, soft delay %s/%d is negative (%g s)', iBc, b.softDelay.hint, b.softDelay.num, new_dur);
                        end
                        obj.blockDurations(iBc)=new_dur;
                    end
                end
            end
            % now check if there are some input soft delays which haven't been found in the sequence
            all_input_hints=sdm_input.keys;
            for i=1:length(all_input_hints)
                if ~sdm_Str2numIDs.isKey(all_input_hints{i})
                    error('Specified soft delay ''%s'' does not exist in the sequence', all_input_hints{i});
                end
            end
        end

        function [easyStruct, errorReport, softDelayState] =getDefaultSoftDelayValues(obj) 
            % go through all the blocks checking the consistency of the soft delays 
            % the code below is copied from checkTiming; we should merge
            % the functionality eventually... TODO/FIXME
            errorReport={};
            softDelayState={};
            for iB=1:length(obj.blockDurations)
                b = obj.getBlock(iB);
                % check soft delays
                if isfield(b, 'softDelay') && ~isempty(b.softDelay)
                    if b.softDelay.factor==0
                        errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' soft delay ' b.softDelay.hint '/' num2str(b.softDelay.num) ' has factor parameter of 0 which is invalid\n' ] };
                        is_ok=false;
                    end                    
                    % calculate the default delay value based on the current block duration
                    def_del=(obj.blockDurations(iB)-b.softDelay.offset)*b.softDelay.factor;                    
                    if (b.softDelay.num>=0)
                        % remember or check for consistency
                        if length(softDelayState)<b.softDelay.num+1 || isempty(softDelayState{b.softDelay.num+1})
                            softDelayState{b.softDelay.num+1}=struct('def',def_del,'hint',b.softDelay.hint, 'blk', iB, 'min', 0.0, 'max', +Inf);
                        else
                            if abs(def_del-softDelayState{b.softDelay.num+1}.def)>1e-7 % what is the reasonable threshold?
                                errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' soft delay ' b.softDelay.hint '/' num2str(b.softDelay.num) ': default duration derived from this block (' num2str(def_del*1e6) 'us) is inconsistent with the previous default (' num2str(softDelayState{b.softDelay.num+1}.def*1e6) 'us) that was derived from block ' num2str(softDelayState{b.softDelay.num+1}.blk) '\n' ] };
                                is_ok=false;
                            end
                            if ~strcmp(b.softDelay.hint, softDelayState{b.softDelay.num+1}.hint) 
                                errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' soft delay ' b.softDelay.hint '/' num2str(b.softDelay.num) ': soft delays with the same numeric ID are expected to share the same text hint but previous hint recorded in block ' num2str(softDelayState{b.softDelay.num+1}.blk) ' is ' softDelayState{b.softDelay.num+1}.hint '\n' ] };
                                is_ok=false;
                            end
                        end
                        % calculate the delay value that would make the block duration of 0, which correponds to min/max
                        lim_del=(-b.softDelay.offset)*b.softDelay.factor;
                        if b.softDelay.factor>0
                            % lim_del corresponds to a minimum
                            if lim_del>softDelayState{b.softDelay.num+1}.min
                                softDelayState{b.softDelay.num+1}.min=lim_del;
                            end
                        else
                            % lim_del corresponds to a maximum
                            if lim_del<softDelayState{b.softDelay.num+1}.max
                                softDelayState{b.softDelay.num+1}.max=lim_del;
                            end
                        end
                    else
                        errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' contains a soft delay ' b.softDelay.hint ' with an invalid numeric ID' num2str(b.softDelay.num) '\n' ] };
                        is_ok=false;
                    end
                end
            end
            % re-package softDelayState into easyStruct
            easyStruct=struct;
            for i=1:length(softDelayState)
                if isempty(softDelayState{i}) 
                    warning('SoftDelay numeric ID %d is unused, we expect contiguous numbering of soft delays',i-1);
                    continue;
                end
                if isfield(easyStruct,softDelayState{i}.hint)
                    error('SoftDelay with numeric ID %d uses the rame hint ''%s'' as some previous SoftDelay',i-1,softDelayState{i}.hint);
                    continue;
                end
                easyStruct.(softDelayState{i}.hint)=softDelayState{i}.def;
            end
        end
        
        function soundData=sound(obj, varargin)
            %sound()
            %   "play out" the sequence through the system speaker
            %

            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'evalLabels';
                parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('channelWeights',[1 1 1],@(x)(isnumeric(x) && length(x)==3));
                parser.addParamValue('onlyProduceSoundData',false,@(x)(islogical(x)));
            end
            parse(parser,varargin{:});
            opt = parser.Results;
            
            if ~isfinite(opt.blockRange(2))
                opt.blockRange(2)=length(obj.blockEvents);
            end
            
            gw_data=obj.waveforms_and_times(false,opt.blockRange);
            total_duration=sum(obj.blockDurations);
            
            sample_rate=44100; %Hz
            dwell_time=1/sample_rate;
            sound_length=floor(total_duration/dwell_time)+1;
            
            soundData(2,sound_length)=0; %preallocate
            
            if ~isempty(gw_data{1})
                soundData(1,:)=interp1(gw_data{1}(1,:),gw_data{1}(2,:)*opt.channelWeights(1),(0:(sound_length-1))*dwell_time,'linear',0);
            end
            if ~isempty(gw_data{2})
                soundData(2,:)=interp1(gw_data{2}(1,:),gw_data{2}(2,:)*opt.channelWeights(2),(0:(sound_length-1))*dwell_time,'linear',0);
            end            
            if ~isempty(gw_data{3})
                tmp=interp1(gw_data{3}(1,:),0.5*gw_data{3}(2,:)*opt.channelWeights(3),(0:(sound_length-1))*dwell_time,'linear',0);
                soundData(1,:)=soundData(1,:)+tmp;
                soundData(2,:)=soundData(2,:)+tmp;
            end
            
            % filter like we did it in the gradient music project
            %b = fir1(40, 10000/sample_rate);
            %sound_data = filter(b, 1, sound_data,[],2);
            % use Gaussian convolution instead to supress ringing
            gw=gausswin(round(sample_rate/6000)*2+1);
            gw=gw/sum(gw(:));
            soundData(1,:) = conv(soundData(1,:), gw, 'same');
            soundData(2,:) = conv(soundData(2,:), gw, 'same');
            
            sound_data_max=max(abs(soundData(:))); 
            soundData = 0.95 * soundData / sound_data_max;
                        
            if ~opt.onlyProduceSoundData
                % info
                fprintf('playing out the sequence waveform, duration %.1gs\n', sound_length*dwell_time);            
                % play out the sound
                % we have to zero-pad the weveform due to the limitations of
                % matlab-to-sound interface            
                sound([zeros(2,sample_rate/2) soundData zeros(2,sample_rate/2)], sample_rate); 
            end
        end
        
        function ok=install(seq,param1,param2)
            %install Install sequence on RANGE system.
            %   install(seq) Install sequence by copying files to Siemens
            %   host and RANGE controller
            % 
            %   install(seq,'sequence_path_or_name') Auto-detect scanner
            %   environment and install the sequence under the given file
            %   name. If sub-directories are provided prior to the name
            %   they will be create automatically.
            %
            %   install(seq,'siemens') Install Siemens Numaris4 file as external.seq
            %   install(seq,'siemensNX') Install Siemens NumarisX file as external.seq
            %   install(seq,'siemens','sequence_path_or_name') Install
            %           Pulseq file assuming a Numaris4 Siemens system 
            %           under the given name and optinally path.
            %   install(seq,'siemens','sequence_path_or_name') Install
            %           Pulseq file assuming a NumarisX Siemens system 
            %           under the given name and optinally path.
            
            if ispc 
                % windows
                ping_command='ping -w 1000 -n 1';
            elseif isunix || ismac
                % unix-like
                ping_command='ping -q -n -W1 -c1';
            end

            % for compatibility with older versions
            if nargin==3
                name=param2;
                dest=param1;
            else
                switch param1
                    case {'siemens','siemensNX'...
                            }
                        name='external';
                        dest=param1;
                    otherwise
                        name=param1;
                        % auto-detect the scanner environment
                        cmd=[ping_command ' 192.168.2.2 && ssh -oBatchMode=yes -oStrictHostKeyChecking=no -oHostKeyAlgorithms=+ssh-rsa root@192.168.2.2 ls /opt/medcom/MriCustomer/CustomerSeq'];
                        [status, ~] = system(cmd);
                        if status == 0
                            fprintf('Siemens NumarisX environment detected\n');
                            dest='siemensNX';
                        else
                            fprintf('Assuming Siemens Numaris4 environment (not tested yet)\n');
                            dest='siemens';
                        end
                end                
            end

            ok = true;
            if any(strcmpi(dest, {'both', 'siemens', 'siemensNX'}))
                seq.write('external.seq.tmp');
                [filepath,filename,ext] = fileparts(name);
                filepath=strrep(filepath,'\','/'); % fix for windows users
                if contains(filepath,'../')
                    error('No relative path elements (like ..) are allowed.');
                end
                name=[filepath '/' filename]; % discard the extension
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
                [status, retmes] = system(['scp -oBatchMode=yes -oStrictHostKeyChecking=no -oHostKeyAlgorithms=+ssh-rsa external.seq.tmp root@' ice_ip ':' pulseq_seq_path '/external_tmp.seq']);
                ok = ok & status == 0;
                if ok
                    if ~isempty(filepath)
                        system(['ssh -oBatchMode=yes -oStrictHostKeyChecking=no -oHostKeyAlgorithms=+ssh-rsa root@' ice_ip ' "mkdir -p ' pulseq_seq_path '/' filepath '"']);
                    end
                    system(['ssh -oBatchMode=yes -oStrictHostKeyChecking=no -oHostKeyAlgorithms=+ssh-rsa root@' ice_ip ' "chmod a+rw ' pulseq_seq_path '/external_tmp.seq"']);
                    system(['ssh -oBatchMode=yes -oStrictHostKeyChecking=no -oHostKeyAlgorithms=+ssh-rsa root@' ice_ip ' "rm -f ' pulseq_seq_path '/' name '.seq"']);
                    system(['ssh -oBatchMode=yes -oStrictHostKeyChecking=no -oHostKeyAlgorithms=+ssh-rsa root@' ice_ip ' "mv ' pulseq_seq_path '/external_tmp.seq ' pulseq_seq_path '/' name '.seq"']);
                else
                    error(['Failed to copy the sequence file to the scanner, the returned error message is: ' strip(retmes)]);
                end
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
