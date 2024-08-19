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
    %   f=plot(...) Return the new figure handle.
    %

    properties (Access = public)
        f   % figure handle
    end

    properties (Access = private)
        ax % array of plot axes handles
        vLines  % array of vline handles

        t_allValues     % array with all unique timepoints in the sequence

        infoPanel
        hTextInfo

    end

    properties (Constant = true, Hidden = true)

        % vertical margin (px)
        margin = 4;
        % lower vertical margin (px)
        my1 = 45;
        % left horizontal margin
        mx1 = 70;
        % right horizontal margin
        mx2 = 5;

        % height of info panel
        panelHeight = 20;
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
                parser.addOptional('label',[]);%,@(x)(isstr(x)));%@(x) any(validatestring(x,validLabel))
                parser.addOptional('hide',false);%,@(x)(isstr(x)));%@(x) any(validatestring(x,validLabel))
                parser.addOptional('stacked',false);%,@(x)(isstr(x)));%@(x) any(validatestring(x,validLabel))
            end
            parse(parser,varargin{:});
            opt = parser.Results;

            obj.f=figure;

            set(obj.f, 'Visible', 'off')

            obj.ax = gobjects(1,6);
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

            t_allValues = [];

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
                label_colors_2plot=parula(length(label_indexes_2plot)+1); % need +1 because the ADC plot by itself also "eats up" one color
                label_colors_2plot=[label_colors_2plot(end,:); label_colors_2plot(1:end-1,:)]; % we like these colors better ?
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
            if (opt.timeDisp=='us')
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
            %for iB=1:size(seq.blockEvents,1)
            for iB=1:length(seq.blockEvents)
                block = seq.getBlock(iB);
                isValid = t0+seq.blockDurations(iB)>timeRange(1) && t0<=timeRange(2);
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
                        plot(tFactor*(t0+t),zeros(size(t)),'rx','Parent',obj.ax(1));
                        plot(tFactor*(t0+t), angle(exp(1i*adc.phaseOffset).*exp(1i*2*pi*t*adc.freqOffset)),'b.','MarkerSize',1,'Parent',obj.ax(3)); % plot ADC phase
                        obj.t_allValues = [obj.t_allValues; tFactor*(t0+t(:))];
                        if label_defined && ~isempty(label_indexes_2plot)
                            set(obj.ax(1),'ColorOrder',label_colors_2plot);
                            label_store_cell=struct2cell(label_store);
                            lbl_vals=[label_store_cell{label_indexes_2plot}];
                            t=t0+adc.delay + (adc.numSamples-1)/2*adc.dwell;
                            p=plot(tFactor*t,lbl_vals,'.','markersize',5,'Parent',obj.ax(1));
                            obj.t_allValues = [obj.t_allValues; tFactor*t(:)];
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
                        if (sreal)
                            plot(tFactor*(t0+t+rf.delay),  real(s),'Parent',obj.ax(2));
                            plot(tFactor*(t0+t+rf.delay),  angle(s.*sign(real(s))*exp(1i*rf.phaseOffset).*exp(1i*2*pi*t    *rf.freqOffset)), tFactor*(t0+tc+rf.delay), angle(sc*exp(1i*rf.phaseOffset).*exp(1i*2*pi*tc*rf.freqOffset)),'xb', 'Parent',obj.ax(3));
                        else
                            plot(tFactor*(t0+t+rf.delay),  abs(s),'Parent',obj.ax(2));
                            plot(tFactor*(t0+t+rf.delay),  angle(s*exp(1i*rf.phaseOffset).*exp(1i*2*pi*t    *rf.freqOffset)), tFactor*(t0+tc+rf.delay), angle(sc*exp(1i*rf.phaseOffset).*exp(1i*2*pi*tc*rf.freqOffset)),'xb', 'Parent',obj.ax(3));
                        end                        
                        %plot(tFactor*(t0+t+rf.delay),  angle(  exp(1i*rf.phaseOffset).*exp(1i*2*pi*t    *rf.freqOffset)), tFactor*(t0+tc+rf.delay), angle(      exp(1i*rf.phaseOffset).*exp(1i*2*pi*t(ic)*rf.freqOffset)),'xb', 'Parent',obj.ax(3));
                        obj.t_allValues = [obj.t_allValues; tFactor*(t0+t(:)+rf.delay)];
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
                            plot(tFactor*(t0+t),waveform,'Parent',obj.ax(3+j));
                            obj.t_allValues = [obj.t_allValues; tFactor*(t0+t(:))];
                        end
                    end
                end
                t0=t0+seq.blockDurations(iB);%mr.calcDuration(block);
            end

            obj.t_allValues = sort(unique(obj.t_allValues));

            % Set axis limits and zoom properties
            dispRange = tFactor*[timeRange(1) min(timeRange(2),t0)];
            arrayfun(@(x)xlim(x,dispRange),obj.ax);
            linkaxes(obj.ax(:),'x')
            h = zoom(obj.f);
            setAxesZoomMotion(h,obj.ax(1),'horizontal');
            % make Y-axes little bit less tight
            arrayfun(@(x) ylim(x, ylim(x) + 0.03*[-1 1]*sum(ylim(x).*[-1 1])), obj.ax(2:end));

            if opt.stacked

                obj.infoPanel = uipanel( ...
                    'Parent',           obj.f, ...
                    'Units',            'pixels', ...
                    'Position',         [0 0 obj.f.Position(3) obj.panelHeight]);

                obj.hTextInfo = uicontrol( ...
                    'Style',            'text', ...
                    'Parent',           obj.infoPanel, ...
                    'units',            'pixels', ...
                    'Position',         [10 0 300 obj.panelHeight], ...
                    'HorizontalAlignment', 'left', ...
                    'FontUnits',        'normalized', ...
                    'FontSize',         0.8, ...
                    'String',           sprintf('t = %.6f ms',0));

                % add vertical lines and make them follow the cursor
                % x-position
                for ii = 1:numel(obj.ax)
                    obj.vLines(ii) = xline(obj.ax(ii), 0, 'r--');
                end
                set(obj.f, 'WindowButtonMotionFcn', @obj.mouseMovement)

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
            %   Is challed whenever the figure-shape is changed and makes
            %   sure all UI elements are correctly psotitioned. This
            %   function implements a vertical stacking of the individual
            %   axes.

            nAxes = numel(obj.ax);
            width  = obj.f.Position(3);
            height = obj.f.Position(4);

            axHeight = (height - (nAxes-1)*obj.margin - obj.my1 - obj.panelHeight) / nAxes;
            axWidth = width - obj.mx1 - obj.mx2;

            for ii = 1:nAxes
                set(obj.ax(ii), 'units', 'pixels', 'Position', [obj.mx1, height-ii*axHeight-(ii-1)*obj.margin, axWidth, axHeight])
                if ii ~= nAxes
                    set(obj.ax(ii), 'Xlabel', [])
                    set(obj.ax(ii), 'XTickLabel', {})
                end
            end

            set(obj.infoPanel, 'Position', [0 0 width obj.panelHeight])

        end

        function mouseMovement(obj, ~, ~)
            % mouseMovement()
            %   triggered when mouse moves in figure. Identifies axis under
            %   the cursor and calculates cursor-position on x-axis (time).
            
            for ida = 1:numel(obj.ax)
                iteratingAx = obj.ax(ida);
                pAx = get(iteratingAx, 'CurrentPoint');
                if obj.inAxis(iteratingAx, pAx(1, 1), pAx(1, 2))
                    tPos = pAx(1, 1);
                    obj.updateGuides(tPos)
                end
            end
        end

        function b = inAxis(obj, ax, x, y)
            % inAxis(ax, x, y)
            %   checks, whether the point with coordinates (x, y) lies
            %   within the limits of 'ax' and returns a bool

            xlims = get(ax, 'XLim');
            ylims = get(ax, 'YLim');
            if x >= xlims(1) && x <= xlims(2) && y >= ylims(1) && y <= ylims(2)
                b = true;
            else
                b = false;
            end
        end

        function updateGuides(obj, tPos)
            % updateGuides(tPos)
            %   updates the time-position for all vertical line objects in
            %   all axes and updates the info string in the info panel.

            [~, idx] = min(abs(obj.t_allValues - tPos));

            for ii = 1:numel(obj.vLines)
                set(obj.vLines(ii), 'Value', obj.t_allValues(idx))
            end

            obj.hTextInfo.String = sprintf('t = %.6f ms', obj.t_allValues(idx)*1000);

        end

    end
end