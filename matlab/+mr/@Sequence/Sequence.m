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

    % Private properties
    %
    properties(GetAccess = public, SetAccess = private)
        version_major;
        version_minor;
        version_revision;
        rfRasterTime;   % RF raster time (system dependent)
        gradRasterTime; % Gradient raster time (system dependent)
        definitions     % Optional sequence definitions
        
        blockEvents;    % Event table (references to events)
        rfLibrary;      % Library of RF events
        gradLibrary;    % Library of gradient events
        adcLibrary;     % Library of ADC readouts
        delayLibrary;   % Library of delay events
        shapeLibrary;   % Library of compressed shapes
    end
    
    methods
        
        function obj = Sequence(varargin)
            obj.version_major = 1;
            obj.version_minor = 1;
            obj.version_revision = 0;
            obj.definitions=containers.Map();
            obj.gradLibrary=mr.EventLibrary();
            obj.shapeLibrary=mr.EventLibrary();
            obj.rfLibrary=mr.EventLibrary();
            obj.adcLibrary=mr.EventLibrary();
            obj.delayLibrary=mr.EventLibrary();
            obj.blockEvents={};
            
            if nargin<1
                sys=mr.opts();
            else
                sys=varargin{1};
            end
            obj.rfRasterTime = sys.rfRasterTime;
            obj.gradRasterTime = sys.gradRasterTime;
        end
        
        
        % See read.m
        read(obj,filename)
        
        % See write.m
        write(obj,filename)
        
        % See readBinary.m
        readBinary(obj,filename);
        
        % See writeBinary.m
        writeBinary(obj,filename);
        
        
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
        
        %TODO: Replacing blocks in the middle of sequence can cause unused
        %events in the libraries. These can be detected and pruned.
        function setBlock(obj,index,varargin)
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
            
            % Convert block structure to cell array of events
            varargin=mr.block2events(varargin);
            %varargin_bak=varargin;
            %varargin(cellfun(@(C)isempty(C),varargin))=[]; % this costs time but does not seem to do anything
            %if (~isequal(varargin, varargin_bak)) then
            %    printf('it does something!');
            %    varargin_bak
            %    varargin
            %end                
            
            %obj.blockEvents(index,:) = zeros(1,6);
            obj.blockEvents{index}=zeros(1,6);
            duration = 0;
            
            % Loop over events adding to library if necessary and creating
            % block event structure.
            for i=1:length(varargin)
                event = varargin{i};
                switch event.type
                    case 'rf'
                        % TODO: Interpolate to 1us time grid using event.t
                        % if required.
                        
                        mag = abs(event.signal);
                        amplitude=max(mag);
                        mag = mag/amplitude;
                        phase = angle(event.signal);
                        phase(phase<0)=phase(phase<0)+2*pi;
                        phase=phase/(2*pi);
                        
                        magShape = mr.compressShape(mag(:));
                        data = [magShape.num_samples magShape.data];
                        [magId,found] = obj.shapeLibrary.find(data);
                        if ~found
                            obj.shapeLibrary.insert(magId,data);
                        end
                        
                        phaseShape = mr.compressShape(phase);
                        data = [phaseShape.num_samples phaseShape.data];
                        [phaseId,found] = obj.shapeLibrary.find(data);
                        if ~found
                            obj.shapeLibrary.insert(phaseId,data);
                        end
                        
                        data = [amplitude magId phaseId event.freqOffset event.phaseOffset event.deadTime event.ringdownTime];
                        [id,found] = obj.rfLibrary.find(data);
                        if ~found
                            obj.rfLibrary.insert(id,data);
                        end
                        
                        %obj.blockEvents(index,2)=id;
                        obj.blockEvents{index}(2)=id;
                        duration=max(duration,length(mag)*obj.rfRasterTime+event.deadTime+event.ringdownTime);
                    case 'grad'
                        channelNum = find(strcmp(event.channel,{'x','y','z'}));
                        amplitude = max(abs(event.waveform));
                        g = event.waveform./amplitude;
                        shape = mr.compressShape(g);
                        data = [shape.num_samples shape.data];
                        [shapeId,found] = obj.shapeLibrary.find(data);
                        if ~found
                            obj.shapeLibrary.insert(shapeId,data);
                        end
                        data = [amplitude shapeId];
                        [id,found] = obj.gradLibrary.find(data);
                        if ~found
                            obj.gradLibrary.insert(id,data,'g');
                        end
                        idx = 2+channelNum;
                        %obj.blockEvents(index,idx)=id;
                        obj.blockEvents{index}(idx)=id;
                        duration=max(duration,length(g)*obj.gradRasterTime);
                    case 'trap'
                        channelNum = find(strcmp(event.channel,{'x','y','z'}));
                        data = [event.amplitude event.riseTime event.flatTime event.fallTime];
                        [id,found] = obj.gradLibrary.find(data);
                        if ~found
                            obj.gradLibrary.insert(id,data,'t');
                        end
                        idx = 2+channelNum;
                        %obj.blockEvents(index,idx)=id;
                        obj.blockEvents{index}(idx)=id;
                        duration=max(duration,event.riseTime+event.flatTime+event.fallTime);
                    case 'adc'
                        data = [event.numSamples event.dwell event.delay ...
                            event.freqOffset event.phaseOffset event.deadTime];
                        [id,found] = obj.adcLibrary.find(data);
                        if ~found
                            obj.adcLibrary.insert(id,data);
                        end
                        %obj.blockEvents(index,6)=id;
                        obj.blockEvents{index}(6)=id;
                        duration=max(duration,event.delay+event.numSamples*event.dwell+event.deadTime);
                    case 'delay'
                        data = [event.delay];
                        [id,found] = obj.delayLibrary.find(data);
                        if ~found
                            obj.delayLibrary.insert(id,data);
                        end
                        %obj.blockEvents(index,1)=id;
                        obj.blockEvents{index}(1)=id;
                        duration=max(duration,event.delay);
                end
            end
        end
        
        function block = getBlock(obj,index)
            %getBlock Return a block of the sequence.
            %   b=getBlock(obj, index) Return the block specified by the
            %   index.
            %
            %   The block is created from the sequence data with all
            %   events and shapes decompressed.
            %
            %   See also  setBlock, addBlock
            
            block=struct('rf',{},'gx',{},'gy',{},'gz',{},'adc',{},'delay',{});
            block(1).rf=[];
            %eventInd = obj.blockEvents(index,:);
            eventInd = obj.blockEvents{index};
            
            if eventInd(1)>0
                delay.type = 'delay';
                delay.delay = obj.delayLibrary.data(eventInd(1)).array;
                block.delay = delay;
            end
            if eventInd(2)>0
                rf.type='rf';
                libData = obj.rfLibrary.data(eventInd(2)).array;
                
                amplitude = libData(1);
                magShape = libData(2);
                phaseShape = libData(3);
                shapeData = obj.shapeLibrary.data(magShape).array;
                compressed.num_samples = shapeData(1);
                compressed.data=shapeData(2:end);
                mag = mr.decompressShape(compressed);
                shapeData = obj.shapeLibrary.data(phaseShape).array;
                compressed.num_samples = shapeData(1);
                compressed.data=shapeData(2:end);
                phase = mr.decompressShape(compressed);
                rf.signal = amplitude*mag.*exp(1j*2*pi*phase);
                rf.t = (1:length(mag))'*obj.rfRasterTime;
                
                rf.freqOffset = libData(4);
                rf.phaseOffset = libData(5);
                % SK: Is this a hack?
                if length(libData)<6
                    libData(end+1)=0;
                end
                rf.deadTime = libData(6);
                % SK: Using the same hack here
                if length(libData) < 7
                    libData(end+1) = 0;
                end
                rf.ringdownTime = libData(7);
                
                block.rf = rf;
            end
            gradChannels = {'gx','gy','gz'};
            for i=1:length(gradChannels)
                if eventInd(2+i)>0
                    type = obj.gradLibrary.type(eventInd(2+i));
                    libData = obj.gradLibrary.data(eventInd(2+i)).array;
                    if type=='t'
                        grad.type = 'trap';
                    else
                        grad.type = 'grad';
                    end
                    grad.channel = gradChannels{i}(2);
                    if strcmp(grad.type,'grad')
                        amplitude = libData(1);
                        shapeId = libData(2);
                        shapeData = obj.shapeLibrary.data(shapeId).array;
                        compressed.num_samples = shapeData(1);
                        compressed.data=shapeData(2:end);
                        g = mr.decompressShape(compressed);
                        grad.waveform = amplitude*g;
                        grad.t = (1:length(g))'*obj.gradRasterTime;
                    else
                        grad.amplitude = libData(1);
                        grad.riseTime = libData(2);
                        grad.flatTime = libData(3);
                        grad.fallTime = libData(4);
                        grad.area = grad.amplitude*(grad.flatTime+grad.riseTime/2+grad.fallTime/2);
                        grad.flatArea = grad.amplitude*grad.flatTime;
                    end
                    
                    block.(gradChannels{i}) = grad;
                end
            end
            if eventInd(6)>0
                libData = obj.adcLibrary.data(eventInd(6)).array;
                if length(libData)<6
                    libData(end+1)=0;
                end
                adc = cell2struct(num2cell(libData),...
                    {'numSamples','dwell','delay','freqOffset','phaseOffset','deadTime'},2);
                adc.type='adc';
                block.adc = adc;
            end
            
        end
        
        
        function f=plot(obj,varargin)
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
                        t=rf.t;
                        plot(tFactor*(t0+t),abs(rf.signal),'Parent',ax(2));
                        plot(tFactor*(t0+t),angle(rf.signal),'Parent',ax(3));
                    end
                    gradChannels={'gx','gy','gz'};
                    for j=1:length(gradChannels)
                        grad=block.(gradChannels{j});
                        if ~isempty(block.(gradChannels{j}))
                            if strcmp(grad.type,'grad')
                                t=grad.t;
                                waveform=1e-3*grad.waveform;
                            else
                                t=cumsum([0 grad.riseTime grad.flatTime grad.fallTime]);
                                waveform=1e-3*grad.amplitude*[0 1 1 0];
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
                            nt_start=round(grad.t(1)/obj.gradRasterTime);
                            waveform=grad.waveform;
                        else
                            nt_start=0;
                            if (abs(grad.flatTime)>eps) % interp1 gets confused by triangular gradients
                                t=cumsum([0 grad.riseTime grad.flatTime grad.fallTime]);
                                trapform=grad.amplitude*[0 1 1 0];
                            else
                                t=cumsum([0 grad.riseTime grad.fallTime]);
                                trapform=grad.amplitude*[0 1 0];
                            end
                            tn=floor(t(end)/obj.gradRasterTime);
                            
                            % it turns out that we need an additional zero-
                            % padding at the endotherwise interp1() 
                            % generates NaNs at the end of the shape
                            t=[t t(end)+obj.gradRasterTime];
                            trapform=[trapform 0];
                            
                            %fprintf('%g : %g | ', [t*1e6 ;trapform]);
                            %fprintf('\n');
                            
                            if abs(grad.amplitude)>eps 
                                waveform=interp1(t,trapform,obj.gradRasterTime*(0:tn));
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
        
                
        function codes=getBinaryCodes(obj)
            %getBinaryCodes Return binary codes for section headers in
            %   in a binary sequence file.
            %
            %   See also  writeBinary

            codes.fileHeader = [1 'pulseq' 2];
            codes.version_major = int64(obj.version_major);
            codes.version_minor = int64(obj.version_minor);
            codes.version_revision = int64(obj.version_revision);
            prefix = bitshift(int64(hex2dec('FFFFFFFF')),32);
            codes.section.definitions = bitor(prefix,int64(1));
            codes.section.blocks      = bitor(prefix,int64(2));
            codes.section.rf          = bitor(prefix,int64(3));
            codes.section.gradients   = bitor(prefix,int64(4));
            codes.section.trapezoids  = bitor(prefix,int64(5));
            codes.section.adc         = bitor(prefix,int64(6));
            codes.section.delays      = bitor(prefix,int64(7));
            codes.section.shapes      = bitor(prefix,int64(8));
        end
    end
end % classdef
