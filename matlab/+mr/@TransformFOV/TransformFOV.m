classdef TransformFOV

    properties (Access = public)
        rotation=[]; % 3x3 rotation matrix or empty matrix if none
        translation=[]; % 1x3 vector or empty matrix if none
    end

    properties (Access = private)
        prior_phase_cycle=0;
        system;
        high_accuracy;
        labels=struct('NOPOS',0,'NOROT',0,'NOSCL',0);
    end

    methods

        function obj = TransformFOV(varargin)

            persistent parser
            if isempty(parser)
                parser = inputParser;
                parser.FunctionName = 'TransformFOV';
                
                addParameter(parser, 'rotation',  [], @(m) (isnumeric(m) && all(size(m)==[3 3])));
                addParameter(parser, 'translation',  [], @(m) (isnumeric(m) && size(m,2)==3 && length(m)==3));
                %addParameter(parser, 'scale',  [], @(m) (isnumeric(m) && size(m,2)==3 && length(m)==3)); % MZ: TODO
                addParameter(parser, 'prior_phase_cycle',  0, @(m) (isnumeric(m) && length(m)==1));
                addParameter(parser, 'high_accuracy',  false, @(m) (islogical(m) && length(m)==1));
                addParameter(parser, 'system', [], @isstruct); 
            end

            parse(parser, varargin{:});
            opt = parser.Results;

            if isempty(opt.rotation) && isempty(opt.translation)
                error('At least one transforming parameter needs to be provided');
            end

            obj.rotation=opt.rotation;
            obj.translation=opt.translation;
            obj.prior_phase_cycle=opt.prior_phase_cycle;
            obj.high_accuracy=opt.high_accuracy;

            if isempty(opt.system)
                obj.system=mr.opts();
            else
                obj.system=opt.system;
            end
        end

        function out=applyToBlock(obj, varargin)

            % convert the input into a plain cell array of events
            if ~any(iscell(varargin))
                events=varargin;
            else
                events={};
                for i=varargin
                    if iscell(i)
                        events=[events i];
                    else
                        events{end+1}=i;
                    end
                end
            end

            % see if we get a block as a single struct and convert it to events
            if (isstruct(events) && isfield(events, 'blockDuration')) || ...
               (iscell(events) && ~isempty(events) && isstruct(events{1}) && isfield(events{1}, 'blockDuration'))
                events=mr.block2events(events);
            end

            % extract various mr events including {rf,adc,gx,gy,gz} from "events" input          
            rf = [];
            adc = [];
            grads = cell(1,3);
            other = {};
            for i = 1:length(events)
                e=events{i};
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
                % extract last time point of possible rf or adc in the block
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

                % MZ: I think we have to rotate the gradient "backwards" if
                % this block has 'NOROT'. We restore the grads object below
                if obj.labels.NOROT 
                    grads_backup=grads;
                    grads=mr.rotate3D(obj.rotation',grads); % MZ: I guess we have to rotate the gradients "back" because we are normally in local logical coordinates, which would be "rotated" if there were NOROT flag
                    % MZ: please check if the above point is correct
                end
            
                phase_cycle_this_block=0;
                        
                % make piecewise polynimials for all gradients in the block and
                % calculate rf and adc phase and frequency offsets based on the gradients 
                for i=1:3
                    if abs(obj.translation(i))>eps % check whether there is shift in the specific direction (x or y or z)
                        g = grads{i};
                        if ~isempty(g) % check whether there is corresponding gradient in the direction of the shift (gx or gy or gz)
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
                            fi_pp = fnint(f_pp); % MZ: TODO: check whether we need more accurate functions here
            
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
                                                if ~obj.high_accuracy % use regular matlab functions for phase calculations
                                                    phase_cycle    = local_frac(local_frac(obj.translation(i)*ppval(fi_pp, t_s)) - freq*(t_s-rf.delay) + obj.prior_phase_cycle);
                                                    rf.phaseOffset = rf.phaseOffset + 2*pi*phase_cycle; 
                                                else % use more accurate implemented functions for phase calculations
                                                    phase_cycle = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, t_s, obj.translation(i)) - multiply_accurate(freq, t_s-rf.delay);
                                                    rf.phaseOffset = rf.phaseOffset + frac_exp2number(multiply_accurate(2*pi, phase_cycle)) + frac_exp2number(multiply_accurate(2*pi, obj.prior_phase_cycle));
                                                end
                                            else %if e=="adc"
                                                adc.freqOffset = adc.freqOffset + freq;
                                                if ~obj.high_accuracy
                                                    phase_cycle = local_frac(local_frac(obj.translation(i)*ppval(fi_pp, t_s)) - freq*(t_s-adc.delay) + obj.prior_phase_cycle);
                                                    adc.phaseOffset = adc.phaseOffset + 2*pi* phase_cycle; 
                                                else
                                                    phase_cycle = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, t_s, obj.translation(i)) - multiply_accurate(freq, t_s-adc.delay);
                                                    adc.phaseOffset = adc.phaseOffset + frac_exp2number(multiply_accurate(2*pi, phase_cycle)) + frac_exp2number(multiply_accurate(2*pi, obj.prior_phase_cycle));
                                                end
                                            end
                                        else 
                                            % in case of non-constant gradient, we should calculate a phase vector for rf or ADC. for rf event we can easily add 
                                            % the phase vector to rf.signal and for ADC event we store the phase vector and use it in image reconstruction
                                            if isfield(e,'t') %e=="rf"
                                                % calculate the frequency at the center
                                                ppval_f_center=ppval(f_pp, rf.delay+rf.center);
                                                freq = obj.translation(i) * ppval_f_center;
                                                rf.freqOffset  = rf.freqOffset + freq;
                                                if ~obj.high_accuracy
                                                    ppval_fi_center=ppval(fi_pp, rf.delay+rf.center);
                                                    phase_cycle        = local_frac(local_frac(obj.translation(i)*ppval_fi_center) - freq*rf.center + obj.prior_phase_cycle);
                                                    rf.phaseOffset = rf.phaseOffset + 2*pi*phase_cycle; 
                                                    phase_cycle_vector = local_frac( (ppval(fi_pp, rf.t+rf.delay)-ppval_fi_center-ppval_f_center*(rf.t-rf.center)) * obj.translation(i));
                                                    rf.signal = rf.signal .* exp(1i*2*pi*phase_cycle_vector);                                                
                                                else
                                                    phase_cycle_vector = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, opt.rf.t-opt.rf.t(1)+opt.rf.delay, opt.translation(i));
                                                    phase_vector = phase_add( zeros(size(phase_cycle_vector)) , phase_cycle_vector );
                                                    opt.rf.signal = opt.rf.signal .* exp(1i*phase_vector); % add phase vector to the signal of rf event
                                                end
                                            else % if e=="adc"
                                                % calculate the frequency at the center
                                                adc_center=0.5*adc.dwell*adc.numSamples;
                                                ppval_f_center=ppval(f_pp, adc.delay+adc_center);
                                                freq = obj.translation(i) * ppval_f_center;
                                                adc.freqOffset  = adc.freqOffset + freq;
                                                if ~obj.high_accuracy
                                                    ppval_fi_center=ppval(fi_pp, adc.delay+adc_center);
                                                    % these -0.5 and +0.5 are needed to avoid unnecessay jumps for values that are very close to 0s
                                                    phase_cycle        = -0.5+local_frac(0.5+local_frac(obj.translation(i)*ppval_fi_center) - freq*adc_center + obj.prior_phase_cycle);
                                                    adc.phaseOffset = adc.phaseOffset + 2*pi*phase_cycle; 
                                                    adc_t=adc.dwell*(0.5:adc.numSamples-0.5);
                                                    % these -0.5 and +0.5 are needed to avoid unnecessay jumps for values that are very close to 0s
                                                    phase_cycle_vector = -0.5+local_frac(0.5+(ppval(fi_pp, adc_t+adc.delay)-ppval_fi_center-ppval_f_center*(adc_t-adc_center) )* obj.translation(i));
                                                    if isempty(adc.phaseModulation)
                                                        adc.phaseModulation=2*pi* phase_cycle_vector(:); % store residual adc phase for image reconstruction 
                                                    else
                                                        adc.phaseModulation=adc.phaseModulation+2*pi* phase_cycle_vector(:); % store residual adc phase for image reconstruction 
                                                    end
                                                else
                                                    phase_cycle_vector = accurate_mod_pp(fi_pp.breaks, fi_pp.coefs, adc.dwell*(0:adc.numSamples-1)+adc.delay, obj.translation(i));
                                                    if isempty(adc.phaseModulation)
                                                        adc.phaseModulation= 2*pi* phase_cycle_vector; % store residual adc phase for image reconstruction 
                                                    else
                                                        adc.phaseModulation=phase_add(adc.phaseModulation,phase_cycle_vector); % store residual adc phase for image reconstruction 
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
            
                            % update the prior phase cycle of the block by only one gradient
                            if ~obj.high_accuracy
                                phase_cycle_this_block = local_frac( phase_cycle_this_block + local_frac( g.area * obj.translation(i))); 
                            else
                                phase_cycle_this_block = mod_accurate( phase_cycle_this_block + mod_accurate(frac_exp2number(multiply_accurate(g.area, obj.translation(i)))) );
                            end
                        end
                    end
                end

                % MZ: now restore the grads object
                if obj.labels.NOROT 
                    grads=grads_backup;
                end

                % now update the phase stored for the next block
                obj.prior_phase_cycle = obj.prior_phase_cycle + phase_cycle_this_block;
            end


            %% rotation
            if ~isempty(obj.rotation) && ~obj.labels.NOROT
                grads = mr.rotate3D(obj.rotation,grads); 
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
        index_e=length(t)
    end
    is_const = all(abs(amp(index_s:index_e) - amp(index_s)) <= 1e-10); % MZ: why this threshold? CHECKME
end


% generate breaks and coefficients for a gradient that is required for
% making piecwise polynomial
function [b, c, tt_extended, waveform_extended] = generate_breaks_coefs(g, tt, waveform, gradRasterTime, t_start, t_end)
if g.type=='grad'
    if g.tt(1)==gradRasterTime/2 % check whether gradient is arbitrary gradient or extended trapezoid
        tt = [tt(1)-gradRasterTime/2, tt, tt(end)+gradRasterTime/2]; % if it's arbitrary gradient, we add first and last gradient amplitudes to the waveform
        waveform = [opt.(g).first, waveform, opt.(g).last];
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

% MZ: question: why not adding (t_start,0) or (t_end,0) to tt/waveform?
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

% local mod function
function out=local_mod(in,m)
out = zeros(size(in));
for i=1:length(in)
    if in(i)>0.5*m
        z=floor(in(i)/m+0.5);
        out(i)=in(i)-z*m;
    elseif in(i)<=-0.5*m
        z=floor(-in(i)/m+0.5);
        out(i)=in(i)+z*m;
    else
        out(i)=in(i);
    end
end
end


%% Accurate mathematical operations using high resolution times (sum , multiply , mod)

% generate fraction and exponent elements of the input
function out = number2frac_exp(in)
[f, e] = log2(in);
out.frac = abs(f);
out.exp = e;
if f<=0
    out.sign = -1;
else
    out.sign = +1;
end
end


% convert fraction and exponent elements of the input to a number
function out = frac_exp2number(in_sfe)
out = in_sfe.sign * in_sfe.frac *2^in_sfe.exp;
end


% accurate mod function
function out = mod_accurate(in) % in should be double
in_sfe = number2frac_exp(in);
frac_digits = 16 - ( floor(log10(abs(in))) + 1 );
if in_sfe.exp == 0
    out = in_sfe.sign * in_sfe.frac;
elseif in_sfe.exp >= 0
    if frac_digits==16
        out = in_sfe.sign * (round(( abs(in) - floor(abs(in)) )*10^(frac_digits-1)))/10^(frac_digits-1);
    else
        out = in_sfe.sign * (round(( abs(in) - floor(abs(in)) )*10^frac_digits))/10^frac_digits;
    end
else
    out = in;
end
end



% accurate multiplication
function out = multiply_accurate(in1, in2)
if isa(in1, 'double')
    in1_sfe = number2frac_exp(in1);
else 
    in1_sfe = in1;
end
if isa(in2, 'double')
    in2_sfe = number2frac_exp(in2);
else
    in2_sfe = in2;
end
out.sign = in1_sfe.sign * in2_sfe.sign;
out.frac = in1_sfe.frac * in2_sfe.frac;
out.exp = in1_sfe.exp + in2_sfe.exp;
end



% accurate division
function out = division_accurate(in1, in2) % in1 / in2
if isa(in1, 'double')
    in1_sfe = number2frac_exp(in1);
else 
    in1_sfe = in1;
end
if isa(in2, 'double')
    in2_sfe = number2frac_exp(in2);
else
    in2_sfe = in2;
end
out.sign = in1_sfe.sign * in2_sfe.sign;
[out.frac, e] = log2(in1_sfe.frac / in2_sfe.frac);
out.exp = e + in1_sfe.exp - in2_sfe.exp;
end



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
    c
    index = find(breaks <= t(c), 1, 'last');
    t0 = t(c);
    area = 0;

    % fast approach
    if index > (i_breaks + 1)
        i_breaks = i_breaks +1;
        AB = multiply_accurate( coefs(i_breaks,1), multiply_accurate( shift, (breaks(i_breaks+1)-breaks(i_breaks))^2-(breaks(i_breaks)-breaks(i_breaks))^2 ) );
        CD = multiply_accurate( coefs(i_breaks,2), multiply_accurate( shift, (breaks(i_breaks+1)-breaks(i_breaks))-(breaks(i_breaks)-breaks(i_breaks)) ) );
        A = [A,mod_accurate( mod_accurate(frac_exp2number(AB)) + mod_accurate(frac_exp2number(CD)) )];
        %A = [A , frac_exp2number(AB) + frac_exp2number(CD)];
    end
    
    AB_n = multiply_accurate( coefs(i_breaks+1,1), multiply_accurate( shift, (t0-breaks(i_breaks+1))^2-(breaks(i_breaks+1)-breaks(i_breaks+1))^2 ) );
    CD_n = multiply_accurate( coefs(i_breaks+1,2), multiply_accurate( shift, (t0-breaks(i_breaks+1))-(breaks(i_breaks+1)-breaks(i_breaks+1)) ) );
    area = mod_accurate ( mod_accurate(sum(A)) + mod_accurate(frac_exp2number(AB_n)) + mod_accurate(frac_exp2number(CD_n)) );
    % area = sum(A) + frac_exp2number(AB_n) + frac_exp2number(CD_n);

    
    % % slow approach
    % for i=1:index
    %     if i<index
    %         % A = multiply__LargSmall_Int_LargSmall_Frac__by__Large_Frac( coefs(i,1) / v^2 , shift ); 
    %         % A = multiply_accurate(coefs(i,1), shift);
    %         % B = number2Int_Frac( u(i+1)^2-u(i)^2 );
    %         % B = number2frac_exp( u(i+1)^2-u(i)^2 );
    %         AB = multiply_accurate( coefs(i,1), multiply_accurate( shift, breaks(i+1)^2-breaks(i)^2 ) );
    %         % C = multiply__LargSmall_Int_LargSmall_Frac__by__Large_Frac( coefs(i,2) / v , shift ); 
    %         % C = multiply_accurate(coefs(i,2), shift);
    %         % D = number2Int_Frac( u(i+1)-u(i) );
    %         % D = number2frac_exp( u(i+1)-u(i) );
    %         CD = multiply_accurate( coefs(i,2), multiply_accurate( shift, breaks(i+1)-breaks(i) ) );
    %         A = [A, mod_accurate(frac_exp2number(AB)) + mod_accurate(frac_exp2number(CD)) ];
    %     elseif i==index
    %         % A = multiply__LargSmall_Int_LargSmall_Frac__by__Large_Frac( coefs(i,1) / v^2 , shift );
    %         % A = multiply_accurate(coefs(i,1), shift);
    %         % B = number2Int_Frac( t0^2-u(i)^2 );
    %         % B = number2frac_exp( t0^2-u(i)^2 );
    %         AB = multiply_accurate( coefs(i,1), multiply_accurate( shift, t0^2-breaks(i)^2 ) );
    %         % C = multiply__LargSmall_Int_LargSmall_Frac__by__Large_Frac( coefs(i,2) / v , shift ); 
    %         % C = multiply_accurate(coefs(i,2), shift);
    %         % D = number2Int_Frac( t0-u(i) );
    %         % D = number2frac_exp( t0-u(i) );
    %         CD = multiply_accurate( coefs(i,2), multiply_accurate( shift, t0-breaks(i) ) );
    %     end
    %     area = mod_accurate ( mod_accurate(area) + mod_accurate(frac_exp2number(AB)) + mod_accurate(frac_exp2number(CD)) );
    %     % area = sum_3_accurate( number2Int_Frac(mod_accurate(area)) , number2Int_Frac(mod_accurate(multiply__LargSmall_Int__by__LargSmall_Int_LargSmall_Frac(B,A))) , number2Int_Frac(mod_accurate(multiply__LargSmall_Int__by__LargSmall_Int_LargSmall_Frac(D,C))) );
    % end

    Mod(c) = area;
end
end



% this function adds phase and phase cycle inputs: out(i) = a(i) + 2*pi*b(i)
function out = phase_add(a , b)  % a is a vector of phase and b is a vector of phase cycles and out is a vector of phase (rf or adc phase vector)
out = zeros(size(b));
for i=1:length(b)
    % out(i) = frac_exp2number(multiply_accurate(2*pi, mod_accurate(a+b(i))));
    out(i) = a(i) + frac_exp2number(multiply_accurate(2*pi, mod_accurate(b(i))));
end
end
