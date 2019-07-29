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
        rfRasterTime;     % RF raster time (system dependent)
        gradRasterTime;   % Gradient raster time (system dependent)
        definitions       % Optional sequence definitions
        
        blockEvents;      % Event table (references to events)
        rfLibrary;        % Library of RF events
        gradLibrary;      % Library of gradient events
        adcLibrary;       % Library of ADC readouts
        delayLibrary;     % Library of delay events
        trigLibrary;      % Library of trigger events ( referenced from the extentions library )
        extensionLibrary; % Library of extension events. Extension events form single-linked zero-terminated lists
        shapeLibrary;     % Library of compressed shapes
        sys;
    end
    
    methods
        
        function obj = Sequence(varargin)
            obj.version_major = 1;
            obj.version_minor = 3; % version minor 3 will now support control events (8th column in the event table)
            obj.version_revision = 0;
            obj.definitions = containers.Map();
            obj.gradLibrary = mr.EventLibrary();
            obj.shapeLibrary = mr.EventLibrary();
            obj.rfLibrary = mr.EventLibrary();
            obj.adcLibrary = mr.EventLibrary();
            obj.delayLibrary = mr.EventLibrary();
            obj.trigLibrary = mr.EventLibrary();
            obj.extensionLibrary = mr.EventLibrary();
            obj.blockEvents = {};
            
            if nargin<1
                sys=mr.opts();
            else
                sys=varargin{1};
            end
            obj.sys = sys;
            obj.rfRasterTime = sys.rfRasterTime;
            obj.gradRasterTime = sys.gradRasterTime;
        end
        
        
        % See read.m
        read(obj,filename,varargin)
        
        % See write.m
        write(obj,filename)
        
        % See readBinary.m
        readBinary(obj,filename);
        
        % See writeBinary.m
        writeBinary(obj,filename);
        
        
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
                b=obj.getBlock(iB);
                eventCount = eventCount + (obj.blockEvents{iB}>0);
                duration=duration+mr.calcDuration(b);
            end
        end
        
        function [is_ok, errorReport]=checkTiming(obj)
            % checkTiming() 
            %     Checks timing of all blocks and objects in the sequence 
            %     optionally returns the detailed error log as cell array
            %     of strings. This function also modifies the sequence
            %     objectby adding the field "TotalDuration" to sequence
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
                if ~isempty(rep)
                    errorReport = { errorReport{:}, [ '   Block:' num2str(iB) ' ' rep '\n' ] };
                end
                totalDuration = totalDuration+dur;
            end
            
            obj.setDefinition('TotalDuration', sprintf('%.9g', totalDuration));
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
        
        function flip(obj,axis)
            %flip Invert all gradinents along the corresponding
            %   axis/channel. The function acts on all gradient objects 
            %   already added to the sequence object
            %
            channelNum = find(strcmp(axis, ...
                                                 {'x', 'y', 'z'}));
            otherChans = find(~strcmp(axis, ...
                                                 {'x', 'y', 'z'}));
            % go through all event table entries and list gradient
            % objects in the library
            paren = @(x, varargin) x(varargin{:}); % anonymous function to access the array on the fly
            allGradEvents = paren(vertcat(obj.blockEvents{:}),:,3:5);
            
            selectedEvents=unique(allGradEvents(:,channelNum));
            selectedEvents=selectedEvents(0~=selectedEvents); % elliminate 0
            otherEvents=unique(allGradEvents(:,otherChans));
            assert(isempty(intersect(selectedEvents,otherEvents)),'ERROR: the same gradient event is used on multiple axes, rhis is not yet supported by flip()');
                                    
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
                obj.gradLibrary.data(selectedEvents(i)).array(1)=-obj.gradLibrary.data(selectedEvents(i)).array(1);
            end
        end
        
        %TODO: Replacing blocks in the middle of sequence can cause unused
        %events in the libraries. These can be detected and pruned.
        function setBlock(obj, index, varargin)
            %setBlock Replace sequence block.
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
            
            block_duration = mr.calcDuration(varargin);
            
            % Convert block structure to cell array of events
            varargin=mr.block2events(varargin);    
            
            obj.blockEvents{index}=zeros(1,7);
            duration = 0;
            
            check_g = {};
            extensions = [];
            
            % Loop over events adding to library if necessary and creating
            % block event structure.
            for i = 1:length(varargin)
                event = varargin{i};
                switch event.type
                    case 'rf'
                        % TODO: Interpolate to 1us time grid using event.t
                        % if required.
                        
                        mag = abs(event.signal);
                        amplitude = max(mag);
                        mag = mag / amplitude;
                        phase = angle(event.signal);
                        phase(phase < 0) = phase(phase < 0) + 2*pi;
                        phase = phase / (2*pi);
                          
                        magShape = mr.compressShape(mag(:));
                        data = [magShape.num_samples magShape.data];
                        [magId,found] = obj.shapeLibrary.find(data);
                        if ~found
                            obj.shapeLibrary.insert(magId, data);
                        end
                        
                        phaseShape = mr.compressShape(phase);
                        data = [phaseShape.num_samples phaseShape.data];
                        [phaseId,found] = obj.shapeLibrary.find(data);
                        if ~found
                            obj.shapeLibrary.insert(phaseId, data);
                        end
                        
                        use = 0;
                        if isfield(event,'use')
                            switch event.use
                                case 'excitation'
                                    use = 1;
                                case 'refocusing'
                                    use = 2;
                                case 'inversion'
                                    use = 3;
                            end
                        end
                        
                        data = [amplitude magId phaseId event.delay ...
                                event.freqOffset event.phaseOffset ...
                                event.deadTime event.ringdownTime use];
                        [id, found] = obj.rfLibrary.find(data);
                        if ~found
                            obj.rfLibrary.insert(id, data);
                        end
                        
                        obj.blockEvents{index}(2) = id;
                        duration = max(duration, length(mag) * ...
                                   obj.rfRasterTime + ...
                                   event.deadTime + ...
                                   event.ringdownTime + event.delay);
                    case 'grad'
                        channelNum = find(strcmp(event.channel, ...
                                                 {'x', 'y', 'z'}));
                        idx = 2 + channelNum;
                                        
                        check_g{channelNum}.idx = idx;
                        check_g{channelNum}.start = [event.delay+min(event.t), event.first];
                        check_g{channelNum}.stop  = [event.delay+max(event.t)+obj.sys.gradRasterTime, event.last]; % MZ: we need to add this gradient raster time, otherwise the gradient appears to be one step too short
                        
                        amplitude = max(abs(event.waveform));
                        if amplitude>0
                            g = event.waveform./amplitude;
                        else
                            g = event.waveform;
                        end
                        shape = mr.compressShape(g);
                        data = [shape.num_samples shape.data];
                        [shapeId,found] = obj.shapeLibrary.find(data);
                        if ~found
                            obj.shapeLibrary.insert(shapeId,data);
                        end
                        data = [amplitude shapeId event.delay event.first event.last];
                        [id,found] = obj.gradLibrary.find(data);
                        if ~found
                            obj.gradLibrary.insert(id, data,'g');
                        end
                        obj.blockEvents{index}(idx) = id;
                        
                        grad_duration = event.delay + length(g)*obj.gradRasterTime; %MZ: was: (length(g)-1)
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
                        
                        data = [event.amplitude event.riseTime ...
                                event.flatTime event.fallTime ...
                                event.delay];
                        [id,found] = obj.gradLibrary.find(data);
                        if ~found
                            obj.gradLibrary.insert(id,data,'t');
                        end
                        obj.blockEvents{index}(idx)=id;
                        duration=max(duration,event.delay+event.riseTime+event.flatTime+event.fallTime);

                    case 'adc'
%                         data = [event.numSamples event.dwell event.delay ...
%                             event.freqOffset event.phaseOffset event.deadTime];
                        data = [event.numSamples event.dwell max(event.delay,event.deadTime) ... % MZ: replaced event.delay+event.deadTime with a max(...) because we allow for overlap of the delay and the dead time
                            event.freqOffset event.phaseOffset event.deadTime];
                        [id,found] = obj.adcLibrary.find(data);
                        if ~found
                            obj.adcLibrary.insert(id,data);
                        end
                        obj.blockEvents{index}(6)=id;
                        duration=max(duration,event.delay+event.numSamples*event.dwell+event.deadTime);
                    case 'delay'
                        data = [event.delay];
                        [id,found] = obj.delayLibrary.find(data);
                        if ~found
                            obj.delayLibrary.insert(id,data);
                        end
                        obj.blockEvents{index}(1)=id;
                        duration=max(duration,event.delay);
                    case {'output','trigger'} 
                        event_type=find(strcmp(event.type,{'output','trigger'}));
                        if (event_type==1)
                            event_channel=find(strcmp(event.channel,{'osc0','osc1','ext1'})); % trigger codes supported by the Siemens interpreter as of May 2019
                        elseif (event_type==2)
                            event_channel=find(strcmp(event.channel,{'physio1','physio2'})); % trigger codes supported by the Siemens interpreter as of June 2019
                        else
                            error('unsupported control event type');
                        end
                        data = [event_type event_channel event.delay event.duration];
                        [id,found] = obj.trigLibrary.find(data);
                        if ~found
                            obj.trigLibrary.insert(id,data);
                        end
                        %obj.blockEvents{index}(7)=id; % now we just
                        % collect the list of extension objects and we will
                        % add it to the event table later
                        ext=struct('type', 1, 'ref', id);
                        extensions=[extensions ext];
                        duration=max(duration,event.delay+event.duration);
                end
            end
            
            if ~isempty(extensions)
                % add extensions now... but it's tricky actually
                % we need to heck whether the exactly the same list if
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
                                last = prev_dat(5);
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
                if cg.stop(2) > obj.sys.maxSlew * obj.sys.gradRasterTime && abs(cg.stop(1)-block_duration) > 1e-7
                    error('A gradient that doesnt end at zero needs to be aligned to the block boundary.');
                end
            end
       
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% GRADIENT CHECKS DONE                                    %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
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
            
            block=struct('rf', {}, 'gx', {}, 'gy', {}, 'gz', {}, ...
                         'adc', {}, 'delay', {});
            block(1).rf = [];
            %eventInd = obj.blockEvents(index,:);
            eventInd = obj.blockEvents{index};
            
            if eventInd(1) > 0
                delay.type = 'delay';
                delay.delay = obj.delayLibrary.data(eventInd(1)).array;
                block.delay = delay;
            end
            if eventInd(7) > 0
                % we have extensions -- triggers for now, may be more to
                % come in the future, e.g. rotation matrices
                % we will eventually isolate this into a separate function
                nextExtID=eventInd(7);
                while nextExtID~=0 
                    extData = obj.extensionLibrary.data(nextExtID).array;
                    % format: extType, extID, nextExtID
                    if (extData(1)==1) % trigger
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
                        %block.trig = trig;
                        % generate extension-specific name
                        filedName=sprintf('trig%d', extData(2));
                        block.(filedName)=trig;
                    else
                        error('unknown extension ID %d', extData(1));
                    end
                    % now update nextExtID
                    nextExtID=extData(3);
                end
            end
            if eventInd(2) > 0
%                 rf.type = 'rf';
%                 libData = obj.rfLibrary.data(eventInd(2)).array;
%                 
%                 amplitude = libData(1);
%                 magShape = libData(2);
%                 phaseShape = libData(3);
%                 shapeData = obj.shapeLibrary.data(magShape).array;
%                 compressed.num_samples = shapeData(1);
%                 compressed.data = shapeData(2:end);
%                 mag = mr.decompressShape(compressed);
%                 shapeData = obj.shapeLibrary.data(phaseShape).array;
%                 compressed.num_samples = shapeData(1);
%                 compressed.data = shapeData(2:end);
%                 phase = mr.decompressShape(compressed);
%                 rf.signal = amplitude*mag.*exp(1j*2*pi*phase);
%                 rf.t = (1:length(mag))'*obj.rfRasterTime;
%                 
%                 rf.delay = libData(4);
%                 rf.freqOffset = libData(5);
%                 rf.phaseOffset = libData(6);
% 
%                 % SK: Is this a hack? (MZ: see below)
%                 if length(libData) < 7
%                     libData(7) = 0;
%                 end
%                 rf.deadTime = libData(7);
%                 % SK: Using the same hack here
%                 if length(libData) < 8
%                     libData(8) = 0;
%                 end
%                 rf.ringdownTime = libData(8);
%                 
%                 % MZ: I think this is needed for compatilbility with reading
%                 % (possibly older) seq-files
%                 if length(libData) < 9
%                     libData(9) = 0;
%                 end
%                 switch libData(9)
%                     case 1
%                         rf.use='excitation';
%                     case 2
%                         rf.use='refocusing';
%                     case 3
%                         rf.use='inversion';
%                 end
                
                block.rf = obj.rfFromLibData(obj.rfLibrary.data(eventInd(2)).array);
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
                        delay = libData(3);
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
                        grad.t = (0:length(g)-1)'*obj.gradRasterTime;
                        grad.delay = delay;
                        if length(libData)>4
                            grad.first = libData(4);
                            grad.last = libData(5);
                        else
                            % for the data read from a file we need to
                            % infer the missing fields here
                            grad.first = grad.waveform(1); % MZ: eventually we should use extrapolation by 1/2 gradient rasters here
                            grad.last = grad.waveform(end);
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
            
        end
        
        function rf = rfFromLibData(obj, libData)                
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
            rf.t = (1:length(mag))'*obj.rfRasterTime;

            rf.delay = libData(4);
            rf.freqOffset = libData(5);
            rf.phaseOffset = libData(6);

            % SK: Is this a hack? (MZ: see below)
            if length(libData) < 7
                libData(7) = 0;
            end
            rf.deadTime = libData(7);
            % SK: Using the same hack here
            if length(libData) < 8
                libData(8) = 0;
            end
            rf.ringdownTime = libData(8);

            % MZ: I think this is needed for compatilbility with reading
            % (possibly older) seq-files
            if length(libData) < 9
                libData(9) = 0;
            end
            switch libData(9)
                case 1
                    rf.use='excitation';
                case 2
                    rf.use='refocusing';
                case 3
                    rf.use='inversion';
            end
        end

        function [ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = calculateKspace(obj, varargin)
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
          
            % initialise the counters and accumulator objects
            c_excitation=0;
            c_refocusing=0;
            c_adcSamples=0;
            % loop throught the blocks to prepare preallocations
            for iB=1:length(obj.blockEvents)
                block = obj.getBlock(iB);
                if ~isempty(block.rf)
                    if (~isfield(block.rf,'use') || ~strcmp(block.rf.use,'refocusing'))
                        c_excitation=c_excitation+1;
                    else
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
                    if (~isfield(block.rf,'use') || ~strcmp(block.rf.use,'refocusing'))
                        t_excitation(c_excitation) = current_dur+t;
                        c_excitation=c_excitation+1;
                    else
                        t_refocusing(c_refocusing) = current_dur+t;
                        c_refocusing=c_refocusing+1;
                    end
                end
                if ~isempty(block.adc)
                    ktime(kcouter:(kcouter-1+block.adc.numSamples)) = (0:(block.adc.numSamples-1))*block.adc.dwell + block.adc.delay + current_dur + traj_recon_delay;
                    kcouter=kcouter+block.adc.numSamples;
                end
                current_dur=current_dur+mr.calcDuration(block);
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
            %   plot(...,'Type',type) Plot the sequence with gradients
            %   displayed according to type: 'Gradient' or 'Kspace'.
            %
            %   plot(...,'TimeRange',[start stop]) Plot the sequence
            %   between the times specified by start and stop.
            %
            %   plot(...,'TimeDisp',unit) Display time in:
            %   's', 'ms' or 'us'.
            %
            %   f=plot(...) Return the new figure handle.
            %
            validPlotTypes = {'Gradient','Kspace'};
            validTimeUnits = {'s','ms','us'};
            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'plot';
                parser.addParamValue('type',validPlotTypes{1},...
                    @(x) any(validatestring(x,validPlotTypes)));
                parser.addParamValue('timeRange',[0 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('timeDisp',validTimeUnits{1},...
                    @(x) any(validatestring(x,validTimeUnits)));
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
            labels={'ADC','RF mag (Hz)','RF ph (rad)','Gx (kHz/m)','Gy (kHz/m)','Gz (kHz/m)'};
            arrayfun(@(x)ylabel(ax(x),labels{x}),1:6);
            
            tFactorList = [1 1e3 1e6];
            tFactor = tFactorList(strcmp(opt.timeDisp,validTimeUnits));
            xlabel(ax(3),['t (' opt.timeDisp ')']);
            xlabel(ax(6),['t (' opt.timeDisp ')']);
            
            t0=0;
            %for iB=1:size(obj.blockEvents,1)
            for iB=1:length(obj.blockEvents)
                block = obj.getBlock(iB);
                isValid = t0>=opt.timeRange(1) && t0<=opt.timeRange(2);
                if isValid
                    if ~isempty(block.adc)
                        adc=block.adc;
                        t=adc.delay + (0:adc.numSamples-1)*adc.dwell;
                        plot(tFactor*(t0+t),zeros(size(t)),'rx','Parent',ax(1));
                    end
                    if ~isempty(block.rf)
                        rf=block.rf;
                        [tc,ic]=mr.calcRfCenter(rf);
                        t=rf.t + rf.delay;
                        tc=tc + rf.delay;
                        plot(tFactor*(t0+t),abs(rf.signal),'Parent',ax(2));
                        plot(tFactor*(t0+t), angle(rf.signal    *exp(1i*rf.phaseOffset).*exp(1i*2*pi*rf.t    *rf.freqOffset)),...
                             tFactor*(t0+tc),angle(rf.signal(ic)*exp(1i*rf.phaseOffset).*exp(1i*2*pi*rf.t(ic)*rf.freqOffset)),'xb',...
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
                                t=grad.delay + [0; grad.t + (grad.t(2)-grad.t(1))/2; grad.t(end) + grad.t(2)-grad.t(1)];
                                waveform=1e-3* [grad.first; grad.waveform; grad.last];
                            else
                                t=cumsum([0 grad.delay grad.riseTime grad.flatTime grad.fallTime]);
                                waveform=1e-3*grad.amplitude*[0 0 1 1 0];
                            end
                            plot(tFactor*(t0+t),waveform,'Parent',ax(3+j));
                        end
                    end                
                end
                t0=t0+mr.calcDuration(block);
            end
            
            % Set axis limits and zoom properties
            dispRange = tFactor*[opt.timeRange(1) min(opt.timeRange(2),t0)];
            arrayfun(@(x)xlim(x,dispRange),ax);
            linkaxes(ax(:),'x')
            h = zoom(fig);
            setAxesZoomMotion(h,ax(1),'horizontal');
        end
        
        function grad_waveforms=gradient_waveforms(obj)
            %gradient_waveforms()
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
                            nt_start=round((grad.delay+grad.t(1))/obj.gradRasterTime);
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
                                %waveform=interp1(t,trapform,obj.gradRasterTime*(0:tn));
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

                t0=t0+mr.calcDuration(block);
                t0_n=round(t0/obj.gradRasterTime);
            end
        end
        
        function sound_data=sound(obj)
            %sound()
            %   "play out" the sequence through the system speaker
            %
            
            grad_waveforms=obj.gradient_waveforms();
            grad_wavelen=size(grad_waveforms,2);
            
            sample_rate=44100; %Hz
            dwell_time=1/sample_rate;
            sound_length=floor((grad_wavelen-1)*obj.gradRasterTime/dwell_time)+1;
            
            sound_data(2,sound_length)=0; %preallocate
            sound_data(1,:)=interp1((0:(grad_wavelen-1))*obj.gradRasterTime,grad_waveforms(1,:)+0.5*grad_waveforms(3,:),(0:(sound_length-1))*dwell_time);
            sound_data(2,:)=interp1((0:(grad_wavelen-1))*obj.gradRasterTime,grad_waveforms(2,:)+0.5*grad_waveforms(3,:),(0:(sound_length-1))*dwell_time);
            
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
    end
end % classdef
