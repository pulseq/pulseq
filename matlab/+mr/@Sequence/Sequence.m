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
    
    properties(Constant)
        RfRasterTime = 1e-6;
        GradRasterTime = 10e-6;
    end
    
    % Private properties
    %
    properties(GetAccess = public, SetAccess = private)
        definitions     % Optional sequence definitions
        
        blockEvents;    % Event table (references to events)
        rfLibrary;      % Library of RF events
        gradLibrary;    % Library of gradient events
        adcLibrary;     % Library of ADC readouts
        delayLibrary;   % Library of delay events
        shapeLibrary;   % Library of compressed shapes
    end
    
    methods
        
        function obj = Sequence()
            obj.definitions=containers.Map();
            obj.gradLibrary=containers.Map('KeyType','double','ValueType','any');
            obj.shapeLibrary=containers.Map('KeyType','double','ValueType','any');
            obj.rfLibrary=containers.Map('KeyType','double','ValueType','any');
            obj.adcLibrary=containers.Map('KeyType','double','ValueType','any');
            obj.delayLibrary=containers.Map('KeyType','double','ValueType','any');
        end
        
        % See read.m
        read(obj,filename)
        
        % See write.m
        write(obj,filename)
        
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
            setBlock(obj,size(obj.blockEvents,1)+1,varargin{:});
            
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
            
            obj.blockEvents(index,:) = zeros(1,6);
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
                        
                        magShape = mr.compressShape(mag);
                        magId = mr.Sequence.searchLibrary(obj.shapeLibrary,magShape);
                        obj.shapeLibrary(magId) = magShape;
                        
                        phaseShape = mr.compressShape(phase);
                        phaseId = mr.Sequence.searchLibrary(obj.shapeLibrary,phaseShape);
                        obj.shapeLibrary(phaseId) = phaseShape;
                        
                        rf.data = [amplitude magId phaseId event.freqOffset event.phaseOffset];
                        id = mr.Sequence.searchLibrary(obj.rfLibrary,rf);
                        obj.rfLibrary(id) = rf;
                        obj.blockEvents(index,2)=id;
                    case 'grad'
                        amplitude = max(abs(event.waveform));
                        g = event.waveform./amplitude;
                        shape = mr.compressShape(g);
                        shapeId = mr.Sequence.searchLibrary(obj.shapeLibrary,shape);
                        obj.shapeLibrary(shapeId) = shape;
                        
                        grad.data = [amplitude shapeId];
                        grad.type = 'grad';
                        id = mr.Sequence.searchLibrary(obj.gradLibrary,grad);
                        obj.gradLibrary(id) = grad;
                        obj.blockEvents(index,2+event.channel)=id;
                    case 'trap'
                        grad.data = [event.amplitude event.riseTime event.flatTime event.fallTime];
                        grad.type = 'trap';
                        id = mr.Sequence.searchLibrary(obj.gradLibrary,grad);
                        obj.gradLibrary(id) = grad;
                        obj.blockEvents(index,event.channel+2)=id;
                    case 'adc'
                        adc.data = [event.numSamples event.dwell event.delay event.freqOffset event.phaseOffset];
                        id = mr.Sequence.searchLibrary(obj.adcLibrary,adc);
                        obj.adcLibrary(id) = adc;
                        obj.blockEvents(index,6)=id;
                    case 'delay'
                        delay.data = [event.delay];
                        id = mr.Sequence.searchLibrary(obj.delayLibrary,delay);
                        obj.delayLibrary(id) = delay;
                        obj.blockEvents(index,1)=id;
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
            eventInd = obj.blockEvents(index,:);
            
            if eventInd(1)>0
                delay.type = 'delay';
                delay.delay = obj.delayLibrary(eventInd(1)).data;
                block.delay = delay;
            end
            if eventInd(2)>0
                rf.type='rf';
                libData = obj.rfLibrary(eventInd(2)).data;
                
                amplitude = libData(1);
                magShape = libData(2);
                phaseShape = libData(3);
                mag = mr.decompressShape(obj.shapeLibrary(magShape));
                phase = mr.decompressShape(obj.shapeLibrary(phaseShape));
                rf.signal = amplitude*mag.*exp(1j*2*pi*phase);
                rf.t = (1:length(mag))'*mr.Sequence.RfRasterTime;
                
                rf.freqOffset = libData(4);
                rf.phaseOffset = libData(5);
                
                block.rf = rf;
            end
            gradChannels = {'gx','gy','gz'};
            for i=1:length(gradChannels)
                if eventInd(2+i)>0
                    type = obj.gradLibrary(eventInd(2+i)).type;
                    libData = obj.gradLibrary(eventInd(2+i)).data;
                    grad.type = type;
                    grad.channel = i;
                    if strcmp(type,'grad')
                        amplitude = libData(1);
                        shapeId = libData(2);
                        g = mr.decompressShape(obj.shapeLibrary(shapeId));
                        grad.waveform = amplitude*g;
                        grad.t = (1:length(g))'*mr.Sequence.GradRasterTime;
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
                libData = obj.adcLibrary(eventInd(6)).data;
                adc = cell2struct(num2cell(libData),...
                    {'numSamples','dwell','delay','freqOffset','phaseOffset'},2);
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
            for i=1:6
                ax(i)=subplot(3,2,i);
            end
            ax=ax([1 3 5 2 4 6]);
            arrayfun(@(x)hold(x,'on'),ax);
            arrayfun(@(x)grid(x,'on'),ax);
            labels={'ADC','RF mag (Hz)','RF ph (rad)','Gx (Hz/m)','Gy (Hz/m)','Gz (Hz/m)'};
            arrayfun(@(x)ylabel(ax(x),labels{x}),1:6);
            
            tFactorList = [1 1e3 1e6];
            tFactor = tFactorList(strcmp(opt.timeDisp,validTimeUnits));
            xlabel(ax(3),['t (' opt.timeDisp ')']);
            xlabel(ax(6),['t (' opt.timeDisp ')']);
            
            t=0;
            for i=1:size(obj.blockEvents,1)
                block = obj.getBlock(i);
                isValid = t>=opt.timeRange(1) && t<=opt.timeRange(2);
                if ~isempty(block.adc)
                    if isValid
                        adc=block.adc;
                        tt=adc.delay + (0:adc.numSamples-1)*adc.dwell;
                        plot(tFactor*(t+tt),zeros(size(tt)),'rx','Parent',ax(1));
                    end
                end
                if ~isempty(block.rf)
                    if isValid
                        rf=block.rf;
                        tt=rf.t;
                        plot(tFactor*(t+tt),abs(rf.signal),'Parent',ax(2));
                        plot(tFactor*(t+tt),angle(rf.signal),'Parent',ax(3));
                    end
                end
                gradChannels={'gx','gy','gz'};
                for j=1:length(gradChannels)
                    grad=block.(gradChannels{j});
                    if ~isempty(block.(gradChannels{j}))
                        if strcmp(grad.type,'grad')
                            tt=grad.t;
                            waveform=grad.waveform;
                        else
                            tt=cumsum([0 grad.riseTime grad.flatTime grad.fallTime]);
                            waveform=grad.amplitude*[0 1 1 0];
                        end
                        
                        if isValid
                            plot(tFactor*(t+tt),waveform,'Parent',ax(3+j));
                        end
                    end
                end
                t=t+mr.calcDuration(block);
            end
            
            linkaxes(ax(:),'x')
            dispRange = tFactor*[opt.timeRange(1) min(opt.timeRange(2),t)];
            xlim(ax(1),dispRange);
            
            h = zoom;
            setAxesZoomMotion(h,ax(1),'horizontal');
        end
    end
    methods (Static)
        
        function id = searchLibrary(library,dataStruct)
            %searchLibrary Lookup a data structure in the given library.
            %   idx=searchLibrary(lib,data) Return the index of the data in
            %   the library. If the data doesn't exist in the library then
            %   the index for the next new entry is returned.
            %
            %   The dataStruct must have a field .data with numeric valued
            %   array
            %
            %   See also  addBlock
            
            found=0;
            keys=cell2mat(library.keys);
            values=library.values;
            for iL=1:length(keys)
                %                 data=library(keys(iL)).data;
                data=values{iL}.data;
                if length(data)==length(dataStruct.data) && ...
                        norm(data-dataStruct.data)<1e-6
                    id=keys(iL);
                    found=1;
                    break;
                end
            end
            if isempty(keys)
                id=1;
            elseif ~found
                id=max(keys)+1;
            end
        end
        
        
        
    end % Static methods
    
    
    
end % classdef
