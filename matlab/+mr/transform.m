function [seq2, gw_pp]= transform(seq, varargin)
%transform : create a copy of the provided sequence that is tranformed by 
%    applying an arbitrary rotation or translation or a combination of both. 
%    This function accepts the following optional parameters:
%
%    rotation  : 3x3 rotation matrix
%    offset    : a translation vector in logical Pulseq coordinates [x,y,z]
%    transform : 4x4 homogeneous transfrom matrix containing both the
%                rotation and translation. The translation part is
%                specified in the laboratory (device) coordinates
%    system    : optional MR system description. If not provided system
%                properties of the input sequence object are inherited
%
%   See also  mr.rotate mr.rotate3D

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'transform';
    
    parser.addParamValue('rotation',  [], @(m) (isnumeric(m) && all(size(m)==[3 3])));
    parser.addParamValue('translation',  [], @(m) (isnumeric(m) && all(size(m)==[1 3])));
    parser.addParamValue('offset',    [], @(v) (isnumeric(v) && length(v)==3)); 
    parser.addParamValue('transform', [], @(m) (isnumeric(m) && all(size(m)==[4 4]))); 
    parser.addParamValue('system', [], @isstruct); 
    parser.addParamValue('sameSeq', false, @islogical); 
    parser.addParamValue('blockRange',[1 inf],@(x)(isnumeric(x) && length(x)==2));

    %parser.addParamValue('gw_pp', {}, @(c) (iscell(c) && isequal(size(c), [1, 3]) && all(cellfun(@isstruct, c))));
    parser.addParamValue('gw_pp', {}, @iscell); % the above check is too computationally expensive

end

parse(parser, varargin{:});
opt = parser.Results;

if opt.sameSeq 
    seq2=seq;
else
    if isempty(opt.system)
        seq2 = mr.Sequence(seq.sys);
    else
        seq2 = mr.Sequence(opt.system);
    end
end

if ~isempty(opt.transform)
    if ~isempty(opt.rotation) || ~isempty(opt.translation)
        error('Neither ''translation'' nor ''rotation'' can be provided in combination with the ''transfrom'' option');
    end
    opt.rotation=opt.transform(1:3,1:3);
    off=opt.transform(4,1:3); % TODO: check whether this is indeed the column
    opt.offset=opt.rotation*off; % TODO: check the direction of the roration (or inverse)
elseif isempty(opt.rotation) && isempty(opt.translation)
    error('At least one transforming parameter needs to be provided');
end

if ~isfinite(opt.blockRange(2))
    opt.blockRange(2)=length(seq.blockDurations);
end

if isempty(opt.gw_pp)
    [~,~,~,~,~,~,~,~, gw_pp] = seq.calculateKspacePP('blockRange',opt.blockRange);
else
    gw_pp=opt.gw_pp;
end

% calculate the gradient moment vector from gradient piecewise polynomials 
gm_pp=cell(1,3);
for i=1:3
    gm_pp{i}=fnint(gw_pp{i});
end

% calculate rf and adc phase offsets based on the slope subtraction idea
% 1: seq duration
total_dur=sum(seq.blockDurations);

% 2: calculate the gm_pp value at the final time point to calculate the slope for each direction
gm_total = zeros(1,3);
for i=1:3
    gm_total(i) = ppval(gm_pp{i},total_dur);
end
gw_diff = gm_total/total_dur;

% 3: subtract the calculated slope from the constant part of the gradients
gw_pp_mod = gw_pp;
for i=1:3
    gw_pp_mod{1,i}.coefs(:,2) = gw_pp{1,i}.coefs(:,2) - gw_diff(i);
end

% 4: calculate the modified gradient moments from the modified gradients
gm_pp_mod = cell(1,3);
for i=1:3
    gm_pp_mod{i}=fnint(gw_pp_mod{i});
end

% 5: calculate the rf and adc phase and frequency offsets from the modified gradient moments 
t_blocks = 0;
phase_acc_rf = zeros(1,3);
phase_acc_adc = zeros(1,3);

% loop through all blocks in the original sequence
% * we first apply the translation then rotation

for iB=opt.blockRange(1):opt.blockRange(2)
    B=seq.getBlock(iB,true); % second parameter means 'addIDs'
    if ~isempty(opt.translation)
        if ~isempty(B.rf) && ~isempty(B.gz)
            t_rf = t_blocks + B.rf.delay;
            phase_acc_rf(1) = ppval(gm_pp_mod{1,1},t_rf);  
            phase_acc_rf(2) = ppval(gm_pp_mod{1,2},t_rf);  
            phase_acc_rf(3) = ppval(gm_pp_mod{1,3},t_rf);  
            B.rf.phaseOffset = B.rf.phaseOffset + ...
                local_mod(phase_acc_rf*(opt.translation)',1)*2*pi + angle(exp(1i*2*pi*gw_diff*(opt.translation)'*1e-9)^round(t_rf*1e9));                  
            B.rf.freqOffset = B.rf.freqOffset + ppval(gw_pp{1,3},t_rf+B.rf.t(1)) * opt.translation(3);
            % remove ID because we've changed the RF object'
            B.rf=rmfield(B.rf,'id');
            % preserve the shapeIDs TODO: check if gradient is constant
            % B.rf=rmfield(B.rf,'shapeIDs');
        end
        if ~isempty(B.adc) && ~isempty(B.gx)
            t__adc = t_blocks + B.adc.delay;
            phase_acc_adc(1) = ppval(gm_pp_mod{1,1},t__adc);  
            phase_acc_adc(2) = ppval(gm_pp_mod{1,2},t__adc);  
            phase_acc_adc(3) = ppval(gm_pp_mod{1,3},t__adc); 
            B.adc.phaseOffset = B.adc.phaseOffset + ...
                local_mod(phase_acc_adc*(opt.translation)',1)*2*pi + angle(exp(1i*2*pi*gw_diff*(opt.translation)'*1e-9)^round(t__adc*1e9));
            B.adc.freqOffset = B.adc.freqOffset + ppval(gw_pp{1,1},t__adc+B.adc.dwell/2) * opt.translation(1);
            % remove ID because we've changed the RF object'
            B.adc=rmfield(B.adc,'id');
        end    
        t_blocks = t_blocks + B.blockDuration;
    end

    if ~isempty(opt.rotation)
        B=mr.rotate3D(opt.rotation,B);
        %TODO: check if we have to removeIDs
    end

    % if isempty(opt.translation) 
    %     % without translation RF and ADC objects remain unchanged so they can use the same IDs
    %     for c={'rf', 'adc'}
    %         if ~isempty(rawIDs.(c{1}))
    %             B.(c{1}).id=rawIDs.(c{1});
    %         end
    %     end
    % end
    % if isempty(opt.rotation) 
    %     % without rotation gradient objects remain unchanged so they can use the same IDs
    %     for c={'gx', 'gy', 'gz'}
    %         if ~isempty(rawIDs.(c{1}))
    %             B.(c{1}).id=rawIDs.(c{1});
    %         end
    %     end
    % end
    
    seq2.addBlock(B);
end
end

function out=local_mod(in,m)
    if in>0.5*m
        z=floor(in/m+0.5);
        out=in-z*m;
    elseif in<=-0.5*m
        z=floor(-in/m+0.5);
        out=in+z*m;
    else
        out=in;
    end
end