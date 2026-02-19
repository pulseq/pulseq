classdef TransformFOV < handle

    properties (Access = public)
        rotation=[]; % 3x3 rotation matrix or empty matrix if none
        translation=[]; % 1x3 vector or empty matrix if none        
    end

    properties (Access = private)
        use_rotation_extension=false; % whether rotations should be explicitly applied to the gradients or implicitly by means of RotationExtension
        prior_phase_cycle=0;
        rotation_quaternion=[]; % rotation quaternion version of the provided rotation matrix or empty matrix if none
        system;
        % high_accuracy;
        labels=struct('NOPOS',0,'NOROT',0,'NOSCL',0);
    end

    methods

        function obj = TransformFOV(varargin)

            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'TransformFOV_f';
                
                addParameter(parser, 'rotation',  [], @(m) (isnumeric(m) && all(size(m)==[3 3])));
                addParameter(parser, 'translation',  [], @(m) (isnumeric(m) && size(m,2)==3 && length(m)==3));
                %addParameter(parser, 'scale',  [], @(m) (isnumeric(m) && size(m,2)==3 && length(m)==3)); % MZ: TODO
                addParamValue(parser, 'transform', [], @(m) (isnumeric(m) && all(size(m)==[4 4])));
                addParameter(parser, 'use_rotation_extension',  false, @(m) (islogical(m) && numel(m)==1));
                addParameter(parser, 'prior_phase_cycle',  0, @(m) (isnumeric(m) && numel(m)==1));
                % addParameter(parser, 'high_accuracy',  false, @(m) (islogical(m) && length(m)==1));
                addParameter(parser, 'system', [], @isstruct); 
            end

            parse(parser, varargin{:});
            opt = parser.Results;

            if ~isempty(opt.transform)
                if ~isempty(opt.rotation) || ~isempty(opt.translation)
                    error('Neither ''translation'' nor ''rotation'' can be provided in combination with the ''transfrom'' option');
                end
                opt.rotation=opt.transform(1:3,1:3);
                off=opt.transform(4,1:3); % TODO: check whether this is indeed the column
                opt.translation=opt.rotation*off; % TODO: check the direction of the roration (or inverse)
            elseif isempty(opt.rotation) && isempty(opt.translation)
                error('At least one transforming parameter needs to be provided');
            end

            obj.rotation=opt.rotation;
            obj.translation=opt.translation;
            obj.prior_phase_cycle=opt.prior_phase_cycle;
            obj.use_rotation_extension=opt.use_rotation_extension;
            % obj.high_accuracy=opt.high_accuracy;

            if obj.use_rotation_extension && ~isempty(obj.rotation) 
                obj.rotation_quaternion=mr.aux.quat.fromRotMat(obj.rotation);
            end

            if isempty(opt.system)
                obj.system=mr.opts();
            else
                obj.system=opt.system;
            end
        end

        function out=applyToBlock(obj, varargin)

            % convert the input into a plain cell array of events
            if ~any(iscell(varargin))
                block_events=varargin;
            else
                block_events={};
                for i=varargin
                    if iscell(i)
                        block_events=[block_events i];
                    else
                        block_events{end+1}=i;
                    end
                end
            end

            % see if we get a block as a single struct and convert it to block_events
            if (isstruct(block_events) && isfield(block_events, 'blockDuration')) || ...
               (iscell(block_events) && ~isempty(block_events) && isstruct(block_events{1}) && isfield(block_events{1}, 'blockDuration'))
                block_events=mr.block2events(block_events);
            end

            % extract various mr events including {rf,adc,gx,gy,gz} from "block_events" input          
            rf = [];
            adc = [];
            grads = cell(1,3);
            other = {};
            rotExtQuaternion = [];
            for i = 1:length(block_events)
                e=block_events{i};
                if length(e)==1 && isstruct(e) && isfield(e, 'type')
                    switch e.type
                        case 'rf'
                            rf = e; % save rf event
                        case 'adc'
                            adc = e; % save adc event
                        case {'trap', 'grad'} % if gradient event, check 'channel'
                            if ~isfield(e, 'channel')
                                error('unspecified gradient channel for the gradient object');
                            end
                            switch e.channel
                                case 'x'
                                    grads{1} = e;
                                case 'y'
                                    grads{2} = e;
                                case 'z'
                                    grads{3} = e;
                                otherwise
                                    error('unsupported gradient channel %s for the gradient object', e.channel);
                            end
                        case 'labelset' %{'labelset', 'labelinc'} % we dont really need 'labelinc', as all labels that are important for us are flags and have no 'inc'
                            for j=1:length(e)
                                switch e(j).label % this switch has only one case on purpose, it is just a lazy way of checking that we deal with a relevant label setting
                                    case {'NOPOS','NOROT','NOSCL'}
                                        obj.labels(e(j).label)=e(j).value;
                                end
                            end
                        case 'rot3D'
                            rotExtQuaternion=e.rotQuaternion;
                        otherwise
                            other{end+1}=e;
                    end
                else
                    other{end+1}=e;
                    %other=[other num2cell(e)]; % e can be an array of structs, and funny enough, num2cell can convert it to a cell array...
                end
            end

            gradRasterTime = obj.system.gradRasterTime;

            %% translation
            if ~isempty(obj.translation)
            % big picture of the algorithm

            % if ~isempty(obj.translation)
            %     phase_cycle_this_block = 0
            %     if NOPOS==0
            %         apply prior_phase_cycle to rf-adc
            %     end
            %     for i=1:3
            %         g = grad{i};
            %         if ~isempty(g)
            %             if translation(i)~=0
            %                 generate piecewise polynomial of the current gradient
            %                 if NOPOS==0
            %                     apply phase and freq offsets to rf-adc by the current gradient
            %                 end
            %                 update phase_cycle_this_block by the current gradient
            %             end
            %         end
            %     end
            %     update prior_phase_cycle =+ phase_cycle_this_block
            % end

                % extract the first and the last time points of possible rf or adc in the block
                if isempty(rf)
                    if isempty(adc)
                        % both ADC and RF are not defined
                        t_end = [];
                        t_start = [];
                    else
                        % only ADC is defined
                        [t_start,t_end] = extract_time(adc);
                    end
                else % RF is defined
                    if isempty(adc)
                        % only RF is defined
                        [t_start,t_end] = extract_time(rf);
                    else
                        % both ADC and RF are defined
                        [t_start_adc,t_end_adc] = extract_time(adc);
                        [t_start_rf,t_end_rf] = extract_time(rf);
                        t_start = min(t_start_adc,t_start_rf);
                        t_end   = max(t_end_adc,t_end_rf);
                    end
                end

                % remove IDs because we will change the objects below
                if ~isempty(rf) && isfield(rf,'id')
                    rf=rmfield(rf,'id');
                end
                if ~isempty(adc) && isfield(adc,'id')
                    adc=rmfield(adc,'id');
                end
                
                % if the current block uses rotation extension and the
                % TransformFOV object is configured not to use the rotation
                % extension then we apply the rotation to the gradients now
                if ~obj.use_rotation_extension && ~isempty(rotExtQuaternion)
                    grads=mr.rotate3D(rotExtQuaternion,grads,'system',obj.sys);
                    rotExtQuaternion=[]; % now that we have applied the current rotation, we can discard it
                end

                % MZ: I think we have to rotate the gradient "backwards" if
                % this block has 'NOROT'. We restore the grads object below
                if obj.labels.NOROT 
                    grads_backup=grads;
                    % MZ: HA! we could rotate obj.translation (or it's copy) in the opposite direction instead
                    % MZ: and, we could use the same mechanism to handle the rotation extention
                    grads=mr.rotate3D(obj.rotation',grads,'system',obj.sys); % MZ: I guess we have to rotate the gradients "back" because we are normally in local logical coordinates, which would be "rotated" if there were NOROT flag
                    % MZ: please check if the above point is correct
                end
                
                % define a temporary parameter for the phase cycle of the current block
                phase_cycle_this_block=0;

                % apply prior_phase_cycle to rf-adc only if NOPOS==0
                if ~obj.labels.NOPOS % do fov positioning for the current block
                    % if there is rf in the current block, adjust its
                    % phase-offset
                    if ~isempty(rf)
                        rf.phaseOffset = rf.phaseOffset + 2*pi * obj.prior_phase_cycle;
                    end
                    % if there is adc in the current block, adjust its
                    % phase-offset
                    if ~isempty(adc)
                        adc.phaseOffset = adc.phaseOffset + 2*pi * obj.prior_phase_cycle;
                    end
                end
                 
                % 1- take all gradintes in the current block and make
                % piecewise polynimials for them
                % 2- only if NOPOS==0 , calculate rf and adc phase and frequency offsets based on
                % the single gradient
                % 3- update the value of 'phase_cycle_this_block'
                % 4- update the value of 'prior_phase_cycle'
                for i=1:3
                    if abs(obj.translation(i))>eps % check whether there is shift in the specific direction (x or y or z)
                        g = grads{i};
                        if ~isempty(g) % check whether there is corresponding gradient in the direction of the shift (gx or gy or gz)
                            % 1- make pp
                            if strcmp(g.type , 'trap') % if g is a simple trapezoid
                                if (abs(g.flatTime)>eps) % interp1 gets confused by triangular gradients (repeating sample)
                                    tt = g.delay+cumsum([0 g.riseTime g.flatTime g.fallTime]);
                                    waveform = g.amplitude*[0 1 1 0];
                                else
                                    if (abs(g.riseTime)>eps && abs(g.fallTime)>eps) % we skip 'empty' gradients
                                        tt = g.delay+cumsum([0 g.riseTime g.fallTime]);
                                        waveform = g.amplitude*[0 1 0];  
                                    else
                                        if abs(g.amplitude)>eps
                                            warning('''empty'' gradient with non-zero magnitude detected');
                                        end
                                    end
                                end
                            else % if g is a extended trapezoid or arbitrary gradient
                                tt = g.delay + g.tt;
                                waveform = g.waveform;                
                            end
            
                            % generate breaks and coefs of the gradient required for
                            % making piecewise polynomial
                            [breaks, coefs, tt_extended, waveform_extended] = generate_breaks_coefs(g, tt, waveform, gradRasterTime, t_start, t_end);
            
                            % make piecewise polynomial for the gradient
                            f_pp = mkpp(breaks, coefs);
                            % integrate it analytically
                            if mr.aux.isOctave()
                              fi_pp = ppint(f_pp);
                            else
                              fi_pp = fnint(f_pp); % MZ: TODO: check whether we need more accurate functions here
                            end
                            
            
                            if ~obj.labels.NOPOS
                                % apply adc or rf phase and freq offset contribution by only one (current) gradient 
                                % (we separate all gradients in the block and apply effect of each one to the offsets)
                                event = {rf,adc};
                                for j=1:length(event)
                                    e = event{j};
                                    if ~isempty(e)
                                        [t_s, t_e] = extract_time(e); % find the first and last time point of the rf or adc event
                                        is_const = is_grad_const(tt_extended, waveform_extended, t_s, t_e); % check whether the gradient is constant during the rf or adc event: 1 means constant
                                        if is_const % in case of constant gradient, we can easily adjust the frequency and phase offset of the rf or ADC event
                                            freq = obj.translation(i) * ppval(f_pp, t_s);
                                            if isfield(e,'t') %e=="rf"
                                                rf.freqOffset  = rf.freqOffset + freq;
                                                phase_cycle = local_frac( accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, t_s, obj.translation(i)) - freq * (t_s-rf.delay) );
                                                rf.phaseOffset = rf.phaseOffset + 2*pi*phase_cycle; 
                                            else %if e=="adc"
                                                adc.freqOffset = adc.freqOffset + freq;
                                                phase_cycle = local_frac( accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, t_s, obj.translation(i)) - freq * (t_s-adc.delay) );
                                                adc.phaseOffset = adc.phaseOffset + 2*pi* phase_cycle;      
                                            end
                                        else 
                                            % in case of non-constant gradient, we should calculate a phase vector for rf or ADC. for rf event we can easily add 
                                            % the phase vector to rf.signal and for ADC event we store the phase vector and use it in image reconstruction
                                            if isfield(e,'t') %e=="rf"
                                                % calculate the frequency at the center
                                                ppval_f_center=ppval(f_pp, rf.delay+rf.center);
                                                freq = obj.translation(i) * ppval_f_center;
                                                rf.freqOffset  = rf.freqOffset + freq;
                                                % ppval_fi_center = ppval(fi_pp, rf.delay+rf.center);
                                                ppval_fi_center_shift = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, rf.delay+rf.center, obj.translation(i));
                                                phase_cycle        = local_frac(ppval_fi_center_shift - freq*rf.center);
                                                rf.phaseOffset = rf.phaseOffset + 2*pi*phase_cycle; 
                                                phase_cycle_vector_tmp = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, rf.t+rf.delay, obj.translation(i));
                                                phase_cycle_vector = local_frac( phase_cycle_vector_tmp -ppval_f_center*(rf.t-rf.center) * obj.translation(i) - ppval_fi_center_shift );
                                                rf.signal = rf.signal .* exp(1i*2*pi*phase_cycle_vector);                                                
                                            else % if e=="adc"
                                                % calculate the frequency at the center
                                                adc_center=0.5*adc.dwell*adc.numSamples;
                                                ppval_f_center=ppval(f_pp, adc.delay+adc_center);
                                                freq = obj.translation(i) * ppval_f_center;
                                                adc.freqOffset  = adc.freqOffset + freq;
                                                % ppval_fi_center=ppval(fi_pp, adc.delay+adc_center);
                                                ppval_fi_center_shift = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, adc.delay+adc_center, obj.translation(i));
                                                % these -0.5 and +0.5 are needed to avoid unnecessay jumps for values that are very close to 0s
                                                phase_cycle        = local_frac( -0.5 + local_frac(0.5 + ppval_fi_center_shift - freq*adc_center) );
                                                adc.phaseOffset = adc.phaseOffset + 2*pi*phase_cycle; 
                                                adc_t=adc.dwell*(0.5:adc.numSamples-0.5);
                                                % these -0.5 and +0.5 are needed to avoid unnecessay jumps for values that are very close to 0s
                                                phase_cycle_vector_tmp = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, adc_t+adc.delay, obj.translation(i));
                                                phase_cycle_vector = local_frac( ( -0.5 + local_frac(0.5 - ppval_f_center*(adc_t-adc_center) ) ) * obj.translation(i) + ppval_fi_center_shift + phase_cycle_vector_tmp );
                                                % phase_cycle_vector = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, adc.dwell*(0:adc.numSamples-1)+adc.delay, obj.translation(i));
                                                if isempty(adc.phaseModulation)
                                                    adc.phaseModulation=2*pi* phase_cycle_vector(:); % store residual adc phase for image reconstruction 
                                                else
                                                    adc.phaseModulation=adc.phaseModulation+2*pi* phase_cycle_vector(:); % store residual adc phase for image reconstruction 
                                                end
                                            end
                                        end
                                    end
                                end
                            end
            
                            % update the phase_cycle_this_block by only the current gradient
                            % phase_cycle_this_block = local_frac( phase_cycle_this_block + local_frac( g.area * obj.translation(i))); 
                            phase_cycle_this_block = local_frac( phase_cycle_this_block + accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, tt_extended(end), obj.translation(i)) ); 
                        end
                    end
                end

                % MZ: now restore the grads object
                if obj.labels.NOROT 
                    grads=grads_backup;
                end

                % now update the phase stored for the next block
                obj.prior_phase_cycle = local_frac( obj.prior_phase_cycle + phase_cycle_this_block ); 
            end

            %% rotation
            if ~isempty(obj.rotation) && ~obj.labels.NOROT
                if obj.use_rotation_extension
                    if isempty(rotExtQuaternion)
                        rotExtQuaternion=obj.rotation_quaternion;
                    else
                        rotExtQuaternion=mr.aux.quat.multiply(rotExtQuaternion, obj.rotation_quaternion); % TODO: check left or right rotation
                    end
                else
                    grads = mr.rotate3D(obj.rotation,grads); 
                end
            end

            %% rotation extension support
            if ~isempty(rotExtQuaternion)
                other{end+1}=mr.makeRotation(rotExtQuaternion);
            end

            %% output
            out=[{rf} {adc} grads(:)' other(:)'];
            out=out(~cellfun(@isempty,out)); % clean up empty
        end

        function seq2 = applyToSeq(obj, seq, varargin)

            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'applyToSeq';
                
                %parser.addRequired('seq');
                parser.addParamValue('sameSeq', false, @islogical); % TODO: add another option for an in-place transform
                parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));
            end
            
            parse(parser, varargin{:});
            opt = parser.Results;

            if ~isfinite(opt.blockRange(2))
                opt.blockRange(2)=length(seq.blockDurations);
            end

            if opt.sameSeq 
                seq2=seq;
            else
                seq2 = mr.Sequence(seq.sys);
                seq2.definitions = seq.definitions; % copy definitions
            end

            obj.labels=struct('NOPOS',0,'NOROT',0,'NOSCL',0);

            for iB=opt.blockRange(1):opt.blockRange(2)
                B=seq.getBlock(iB,opt.sameSeq); % second parameter means 'addIDs'
                B2=obj.applyToBlock(B);
                seq2.addBlock(B2);
            end
        end
    end
end

%% Local Functions

% find the first and last time point of the rf or adc input in the block
function [t_s, t_e] = extract_time(event)
    if strcmp(event.type,'adc')
        t_s = event.delay + event.dwell * 0.5;
        t_e = event.delay + event.dwell * (event.numSamples-0.5);
    elseif strcmp(event.type,'rf')
        t_s = event.delay + event.t(1);
        t_e = event.delay + event.t(end);
    end
end

% check whether gradient is constant during rf or adc time points

% inputs are: t=tt and amp=waveform of the gradient
% and t_start and t_end are the first and last time point of the rf or adc event 
% flag=1 means the gradient is constant
function is_const = is_grad_const(t, amp, t_start, t_end)  
    index_s = find(t <= t_start, 1, 'last');
    index_e = find(t >= t_end, 1, 'first');
    if isempty(index_s)
        index_s=1;
    end
    if isempty(index_e)
        index_e=length(t);
    end
    is_const = all(abs(amp(index_s:index_e) - amp(index_s)) <= 1e-10); % MZ: why this threshold? CHECKME
end


% generate breaks and coefficients for a gradient that is required for
% making piecwise polynomial
function [b, c, tt_extended, waveform_extended] = generate_breaks_coefs(g, tt, waveform, gradRasterTime, t_start, t_end)
if g.type=='grad'
    if abs(g.tt(1)-gradRasterTime/2)<eps % check whether gradient is arbitrary gradient or extended trapezoid
        tt = [tt(1)-gradRasterTime/2; tt(:); tt(end)+gradRasterTime/2]; % if it's arbitrary gradient, we add first and last gradient amplitudes to the waveform
        waveform = [g.first; waveform(:); g.last];
    end
end
tt_extended = tt(:)';
waveform_extended = waveform(:)';

% generate initial breaks and coefs
b = tt_extended;
c = zeros(length(tt)-1,2);
c(:,1) = diff(waveform) ./ diff(tt);
c(find(isnan(c))) = 0;
c(find(isinf(c))) = 0;
c(:,2) = (waveform(1:end-1))';

% MZ: question: why not adding (t_start,0) or (t_end,0) to tt/waveform? %
% MS: if we add 0 at the begining/end of the waveform and then calculate the slope of the first/last time slot, we will get a non-zero slope (if g.first or g.last are not zero). 
% This is because the 0 value is connected to non-zero g.first or g.last and it means there is real gradient between them but should not exist. 
% so it's better to add [0,0] at the begining or end of the coef matrix later on. 
% modify breaks and coefs if the rf or adc event(s) is/are started (or ended)
% before (or after) the first (or last) time point of the gradient
if t_start<tt(1)
    b = [t_start, b];
    c = [0,0;c];
    tt_extended = [t_start, tt_extended];
    waveform_extended = [0, waveform_extended];
end
if t_end>tt(end)
    b = [b, t_end];
    c = [c;0,0];
    tt_extended = [tt_extended, t_end];
    waveform_extended = [waveform_extended, 0];
end
end


% % make piecewise polynomial based on time points and amplitudes of gradient event
% function f_pp = mkpp_linear(t,f)
% breaks = t;
% coefs = zeros(length(t)-1,2);
% coefs(:,1) = diff(f) ./ diff(t);
% coefs(find(isnan(coefs))) = 0;
% coefs(find(isinf(coefs))) = 0;
% coefs(:,2) = (f(1:end-1))';
% f_pp = mkpp(breaks,coefs);
% end


% local function returning the fractional part
function out=local_frac(in)
    out = in-floor(in);
end

% % local mod function
% function out=local_mod(in,m)
% out = zeros(size(in));
% for i=1:length(in)
%     if in(i)>0.5*m
%         z=floor(in(i)/m+0.5);
%         out(i)=in(i)-z*m;
%     elseif in(i)<=-0.5*m
%         z=floor(-in(i)/m+0.5);
%         out(i)=in(i)+z*m;
%     else
%         out(i)=in(i);
%     end
% end
% end


%% Accurate mathematical operations using high resolution times (sum , multiply , mod)
% in this function, we take the time breaks and coefficients of the
% piecewise polynomial of the gradient and desired time point of the rf or ADC event (t) and
% amount of the fov shift. and accurately calculate the area under each time
% slot of the gradient and multiply it by shift and finally calculate the
% mod of that by 1. (phase cycle)
function Mod = accurate_mod_pp(breaks, coefs, t, shift) 
Mod = zeros(size(t));
A = [];
i_breaks = 0;
for c=1:length(t)
    % c
    index = find(breaks <= t(c), 1, 'last');
    t0 = t(c);
    area = 0;

    while index > (i_breaks + 1)
        i_breaks = i_breaks +1;
        AB = coefs(i_breaks,1) * shift * ( (breaks(i_breaks+1)-breaks(i_breaks))^2-(breaks(i_breaks)-breaks(i_breaks))^2 );
        CD = coefs(i_breaks,2) * shift * ( (breaks(i_breaks+1)-breaks(i_breaks))-(breaks(i_breaks)-breaks(i_breaks)) );
        A = [A,local_frac( local_frac(AB) + local_frac(CD) )];
    end
    
    if t0==breaks(i_breaks+1)
        area =  local_frac(sum(A));
    else
        AB_n = coefs(i_breaks+1,1) * shift * ( (t0-breaks(i_breaks+1))^2-(breaks(i_breaks+1)-breaks(i_breaks+1))^2 );
        CD_n = coefs(i_breaks+1,2) * shift * ( (t0-breaks(i_breaks+1))-(breaks(i_breaks+1)-breaks(i_breaks+1)) );
        area =  local_frac( local_frac(sum(A)) + local_frac(AB_n) + local_frac(CD_n) );
    end
    Mod(c) = area;
end
end

