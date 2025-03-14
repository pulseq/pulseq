classdef SeqPlot < handle
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

    properties (Access = public)
        f   % figure handle
    end

    properties (Access = private)
        ax % array of plot axes handles
        vLines  % array of vline handles

        hSeq

    end

    properties (Constant = true, Hidden = true)

        % vertical margin (px)
        margin = 6;
        % lower vertical margin (px)
        my1 = 45;
        % left horizontal margin
        mx1 = 70;
        % right horizontal margin
        mx2 = 5;

    end


    methods

        function obj = SeqPlot(seq, varargin)

            validTimeUnits = {'s','ms','us'};
            validLabel = mr.getSupportedLabels();
            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'plot';
                parser.addParamValue('showBlocks',false,@(x)(isnumeric(x) || islogical(x)));
                parser.addParamValue('timeRange',[0 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
                parser.addParamValue('timeDisp',validTimeUnits{1},...
                    @(x) any(validatestring(x,validTimeUnits)));
                parser.addParamValue('label',[]);%,@(x)(isstr(x)));%@(x) any(validatestring(x,validLabel))
                parser.addParamValue('hide',false);%,@(x)(isstr(x)));%@(x) any(validatestring(x,validLabel))
                parser.addParamValue('stacked',false);%,@(x)(isstr(x)));%@(x) any(validatestring(x,validLabel))
                parser.addParamValue('showGuides',true);%,@(x)(isstr(x)));%@(x) any(validatestring(x,validLabel))
            end
            parse(parser,varargin{:});
            opt = parser.Results;
            
            if mr.aux.isOctave()
              if opt.stacked
                warning('Option stacked is not (yet) supported by Octave');
                opt.stacked=false;
              end
              if opt.showBlocks
                warning('Option stacked is not (yet) supported by Octave');
                opt.showBlocks=false;
              end
            end

            obj.f=figure;
            obj.hSeq=seq; % Sequence is a handle-class so copying is cheap...

            set(obj.f, 'Visible', 'off')

            if ~mr.aux.isOctave()
              obj.ax = gobjects(1,6);              
            end
            for i=1:6
                obj.ax(i)=subplot(3,2,i);
            end
            obj.ax=obj.ax([1 3 5 2 4 6]);   % Re-order axes
            arrayfun(@(x)hold(x,'on'),obj.ax);
            arrayfun(@(x)grid(x,'on'),obj.ax);
            labels={'ADC/labels','RF mag (Hz)','RF/ADC ph (rad)','Gx (kHz/m)','Gy (kHz/m)','Gz (kHz/m)'};
            arrayfun(@(x)ylabel(obj.ax(x),labels{x}),1:6);
            
            tFactorList = [1 1e3 1e6];
            tFactor = tFactorList(strcmp(opt.timeDisp,validTimeUnits));
            xlabel(obj.ax(3),['t (' opt.timeDisp ')']);
            xlabel(obj.ax(6),['t (' opt.timeDisp ')']);

            t0=0;
            label_defined=false;
            label_indexes_2plot=[];
            label_legend_2plot=[];
            for i=1:length(validLabel)
                label_store.(validLabel{i})=0;
                if ~isempty(opt.label) && ~isempty(strfind(upper(opt.label),validLabel{i}))
                    label_indexes_2plot=[label_indexes_2plot i];
                    label_legend_2plot=[label_legend_2plot; validLabel{i}];
                end
            end
            if ~isempty(label_indexes_2plot)
                label_colors_2plot=parula(length(label_indexes_2plot)+1); % need +1 because the ADC plot by itself also "eats up" one color
                label_colors_2plot=[label_colors_2plot(end,:); label_colors_2plot(1:end-1,:)]; % we like these colors better ?
            end

            % time format
            switch opt.timeDisp
                case 'us' 
                    timeFormat='%.1f';
                case 'ms'
                    timeFormat='%.4f';
                otherwise 
                    timeFormat='%.7f';
            end
            
            % data cursor callback
            if ~mr.aux.isOctave()
              hDCM = datacursormode(obj.f);
              hDCM.UpdateFcn = @(src, event)DataTipHandler(obj,tFactor,[timeFormat ' ' opt.timeDisp],src,event);
            end

            % time/block range
            timeRange=opt.timeRange;
            blockEdges=[0 cumsum(seq.blockDurations)];
            if opt.blockRange(1)>1 && blockEdges(opt.blockRange(1))>timeRange(1)
                timeRange(1)=blockEdges(opt.blockRange(1));
            end
            if isfinite(opt.blockRange(2)) && opt.blockRange(2)<length(seq.blockDurations) && blockEdges(opt.blockRange(2)+1)<timeRange(2)
                timeRange(2)=blockEdges(opt.blockRange(2)+1);
            end
            % block timings
            blockEdgesInRange=blockEdges(logical((blockEdges>=timeRange(1)).*(blockEdges<=timeRange(2))));
            if strcmp(opt.timeDisp,'us') && ~mr.aux.isOctave()
                for i=1:6
                    xax=get(obj.ax(i),'XAxis');
                    xax.ExponentMode='manual';
                    xax.Exponent=0;
                end
            end
            if opt.showBlocks
                % show block edges in plots
                for i=1:6
                    xax=get(obj.ax(i),'XAxis');
                    xax.TickValues=unique(tFactor.*blockEdgesInRange);
                    set(obj.ax(i),'XTickLabelRotation',90);
                    %xax.MinorTickValues=tFactor.*blockEdgesInRange;
                    %set(obj.ax(i),'XMinorTick', 'on');
                    %set(obj.ax(i),'XMinorGrid', 'on');
                    %set(obj.ax(i),'GridColor',0.8*[1 1 1]);
                    %set(obj.ax(i),'MinorGridColor',0.6*[1 1 1]);
                    %set(obj.ax(i),'MinorGridLineStyle','-');
                end
            end
            %
            gradChannels={'gx','gy','gz'};                    

            % loop through blocks
            for iB=1:length(seq.blockEvents)
                block = seq.getBlock(iB);
                if t0<=timeRange(2)
                    % update the labels / counters even if we are below the display range
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
                end
                isValid = t0+seq.blockDurations(iB)>timeRange(1) && t0<=timeRange(2);
                if isValid
                    if ~isempty(block.adc)
                        adc=block.adc;
                        t=adc.delay + ((0:adc.numSamples-1)'+0.5)*adc.dwell; % according to the information from Klaus Scheffler and indirectly from Siemens this is the present convention (the samples are shifted by 0.5 dwell)
                        p1=plot(tFactor*(t0+t),zeros(size(t)),'rx','Parent',obj.ax(1));
                        if isempty(adc.phaseModulation)
                            adc.phaseModulation=0;
                        end
                        full_freqOffset=adc.freqOffset+adc.freqPPM*1e-6*seq.sys.gamma*seq.sys.B0;
                        full_phaseOffset=adc.phaseOffset+adc.phasePPM*1e-6*seq.sys.gamma*seq.sys.B0;
                        p2=plot(tFactor*(t0+t), angle(exp(1i*(full_phaseOffset+adc.phaseModulation)).*exp(1i*2*pi*t*full_freqOffset)),'b.','MarkerSize',1,'Parent',obj.ax(3)); % plot ADC phase 
                        % labels/counters/flags
                        if label_defined && ~isempty(label_indexes_2plot)
                            set(obj.ax(1),'ColorOrder',label_colors_2plot);
                            label_store_cell=struct2cell(label_store);
                            lbl_vals=[label_store_cell{label_indexes_2plot}];
                            t=t0+adc.delay + (adc.numSamples-1)/2*adc.dwell;
                            p=plot(tFactor*t,lbl_vals,'.','markersize',5,'Parent',obj.ax(1));
                            if ~isempty(label_legend_2plot)
                                legend(obj.ax(1),p,label_legend_2plot,'location','Northwest','AutoUpdate','off');
                                label_legend_2plot=[];
                            end
                        end
                    end
                    if ~isempty(block.rf)
                        rf=block.rf;
                        [tc,ic]=mr.calcRfCenter(rf);
                        sc=rf.signal(ic);
                        if max(abs(diff(rf.t)-rf.t(2)+rf.t(1)))<1e-9 && length(rf.t)>100
                            % homogeneous sampling and long pulses -- use lower time resolution for better display and performance
                            dt=rf.t(2)-rf.t(1);
                            st=round(seq.sys.gradRasterTime/dt);
                            t=rf.t(1:st:end);
                            s=rf.signal(1:st:end);
                            % always include the last point for the accurate display
                            if (t(end)~=rf.t(end))
                                t(end+1)=rf.t(end);
                                s(end+1)=rf.signal(end);
                            end
                        else
                            t=rf.t;
                            s=rf.signal;
                        end
                        sreal=max(abs(imag(s)))/max(abs(real(s)))<1e-6; %all(isreal(s));
                        if abs(s(1))~=0 % fix strangely looking phase / amplitude in the beginning
                            s=[0; s];
                            t=[t(1); t];
                            %ic=ic+1;
                        end
                        if abs(s(end))~=0 % fix strangely looking phase / amplitude in the beginning
                            s=[s; 0];
                            t=[t; t(end)];
                        end
                        full_freqOffset=rf.freqOffset+rf.freqPPM*1e-6*seq.sys.gamma*seq.sys.B0;
                        full_phaseOffset=rf.phaseOffset+rf.phasePPM*1e-6*seq.sys.gamma*seq.sys.B0;
                        if (sreal)
                            p1=plot(tFactor*(t0+t+rf.delay),  real(s),'Parent',obj.ax(2));
                            p2=plot(tFactor*(t0+t+rf.delay),  angle(s.*sign(real(s))*exp(1i*full_phaseOffset).*exp(1i*2*pi*t    *full_freqOffset)), tFactor*(t0+tc+rf.delay), angle(sc*exp(1i*full_phaseOffset).*exp(1i*2*pi*tc*full_freqOffset)),'xb', 'Parent',obj.ax(3));
                        else
                            p1=plot(tFactor*(t0+t+rf.delay),  abs(s),'Parent',obj.ax(2));
                            p2=plot(tFactor*(t0+t+rf.delay),  angle(s*exp(1i*full_phaseOffset).*exp(1i*2*pi*t    *full_freqOffset)), tFactor*(t0+tc+rf.delay), angle(sc*exp(1i*full_phaseOffset).*exp(1i*2*pi*tc*full_freqOffset)),'xb', 'Parent',obj.ax(3));
                        end                        
                    end
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
                            p=plot(tFactor*(t0+t),waveform,'Parent',obj.ax(3+j));
                        end
                    end
                end
                t0=t0+seq.blockDurations(iB);%mr.calcDuration(block);
            end

            % Set axis limits and zoom properties
            dispRange = tFactor*[timeRange(1) min(timeRange(2),t0)];
            arrayfun(@(x)xlim(x,dispRange),obj.ax);
            linkaxes(obj.ax(:),'x')
            if ~mr.aux.isOctave()
              h = zoom(obj.f);
              setAxesZoomMotion(h,obj.ax(1),'horizontal');
            end
            % manually fix the phase vertical scale to +- pi
            ylim(obj.ax(3),[-pi pi]);
            % make Y-axes little bit less tight
            arrayfun(@(x) ylim(x, ylim(x) + 0.03*[-1 1]*sum(ylim(x).*[-1 1])), obj.ax(2:end));
            

            if opt.showGuides
              if mr.aux.isOctave()
                warning('Option showGuides is not implemented in Octave');
              else
                % add vertical lines and make them follow the cursor
                % x-position
                for ii = 1:numel(obj.ax)
                    obj.vLines(ii) = xline(obj.ax(ii), 0, 'r--');
                end
              end
            end

            if opt.stacked
                % vertical stacking is defined in guiResize
                set(obj.f, 'ResizeFcn', @obj.guiResize)
                obj.guiResize()
            end

            if ~opt.hide
                set(obj.f, 'Visible', 'on')
            end

            % do not assign to 'ans' when called without assigned variable
            if nargout == 0
                clear obj
            end
        end

        function guiResize(obj, ~, ~)
            % guiResize()
            %   Is called whenever the figure-shape is changed and makes
            %   sure all UI elements are correctly psotitioned. This
            %   function implements a vertical stacking of the individual
            %   axes.

            nAxes = numel(obj.ax);
            width  = obj.f.Position(3);
            height = obj.f.Position(4);

            axHeight = (height - (nAxes-1)*obj.margin - obj.my1) / nAxes;
            axWidth = width - obj.mx1 - obj.mx2;

            for ii = 1:nAxes
                set(obj.ax(ii), 'units', 'pixels', 'Position', [obj.mx1, height-ii*axHeight-(ii-1)*obj.margin, axWidth, axHeight])
                if ii ~= nAxes
                    set(obj.ax(ii), 'Xlabel', [])
                    set(obj.ax(ii), 'XTickLabel', {})
                end
            end
        end

        function out=DataTipHandler(obj, tfactor, timeFormat, src, event) 
            if ~isa(event,'matlab.graphics.internal.DataTipEvent') || ...
               ~isprop(event, 'Position') || length(event.Position)<2 || ...
               ~isprop(event, 'Target')
                out=[];
                return;
            end
            ax=src.Host.Parent;
            % get the relevant target from the y-axes title
            at=lower(ax.YLabel.String);
            if strcmp(at(1:3),'adc') || ... 
               (strcmp(at(1:6),'rf/adc') && strcmp(event.Target.LineStyle,'none') && strcmp(event.Target.Marker,'.')) % we need to check whether we are dealing with the ADC phase, which is also shown in the same panel as the RF                
                field='adc';
            else
                field=at(1:2);
            end
            % create the custom data tip as tex-formatted cell array of lines
            t=event.Position(1);
            t0=t;
            if isa(event.Target,'matlab.graphics.chart.primitive.Line')
                % for trapezoid gradients the last point may belong to the next block
                t0=event.Target.XData(1); 
            end            
            iB=obj.hSeq.findBlockByTime(t0/tfactor);
            rb=obj.hSeq.getRawBlockContentIDs(iB);
            out={['\bf\color{blue}t:\rm\color{black}' sprintf(timeFormat,t)],...
                 ['\bf\color{blue}Y:\rm\color{black}' num2str(event.Position(2))],...
                 ''};
            if isempty(rb.(field))
                out{3}=['\bf\color{blue}blk:\rm\color{black}' num2str(iB)];
            else
                try
                    switch field(1)
                        case 'a'
                            name = obj.hSeq.adcID2NameMap(rb.(field));
                        case 'r'
                            name = obj.hSeq.rfID2NameMap(rb.(field));
                        otherwise
                            name = obj.hSeq.gradID2NameMap(rb.(field));
                    end
                    out{3}=['\bf\color{blue}blk:\rm\color{black}' num2str(iB) ' \bf\color{blue}' field '\_id:\rm\color{black}' num2str(rb.(field)) ' ''\bf\color{darkGreen}' name '\rm\color{black}'''];
                catch 
                    out{3}=['\bf\color{blue}blk:\rm\color{black}' num2str(iB) ' \bf\color{blue}' field '\_id:\rm\color{black}' num2str(rb.(field))];
                end
            end
            
            % we need to delay the call of the update, otherwise the plot
            % object generates an exception
            t=timer('StartDelay',0e-3,'Period',1e-3,'TimerFcn',@(~,~)updateGuides(obj,t));
            t.start();
        end

        function updateGuides(obj, tPos)
            % updateGuides(tPos)
            %   updates the time-position for all vertical line objects in
            %   all axes 

            for ii = 1:numel(obj.vLines)
                set(obj.vLines(ii), 'Value', tPos);
            end
        end
    end
end

