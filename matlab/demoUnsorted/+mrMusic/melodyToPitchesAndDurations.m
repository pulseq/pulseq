function [pitches, durations] = melodyToPitchesAndDurations(melody, varargin)
% mrMusic.melodyToPitchesAdnDurations : convert *melody* to the 
%     channel-pitch-duration tables

persistent parser
validArticulations={'legato', 'non-legato', 'staccato'};

if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'mrMusic.melodyToPitchesAndDurations';
    parser.addParamValue('timeSignature',4/4,@isnumeric);
    parser.addParamValue('defaultArticulation',validArticulations{1},...
        @(x) any(validatestring(x,validArticulations)));
    %parser.addParamValue('defaultLegato',true,@islogical);
end
parse(parser,varargin{:});
opt = parser.Results;

defaultLegato=strcmp(opt.defaultArticulation,'legato');
timeSignature=opt.timeSignature;

if ~defaultLegato
    % change all 'non-negative' notes to non-legato by making them 0.95 times shorter and introducing a delay of 0.05
    % and fix all 'negative' notes
    if strcmp(opt.defaultArticulation,'non-legato')
        nonLegatoDuration=0.97;
    else
        nonLegatoDuration=0.6;
    end
    for b=1:size(melody,1)
        for c=1:3
            v=melody{b,c};
            vn=[]; 
            for n=1:length(v)
                if real(v(n))>0
                    vn(1,end+1)=real(v(n))+1j*imag(v(n))*nonLegatoDuration;
                    vn(1,end+1)=1j*imag(v(n))*(1-nonLegatoDuration);
                else
                    vn(1,end+1)=abs(real(v(n)))+1j*abs(imag(v(n)));
                end
            end
            melody{b,c}=vn;
        end
    end
end

dur_tol=1e-6; % duration rounding error tolerance

% I know, the code below is very ugly, not vectorized and suboptimal in many ways 
pitches=[];
durations=[];
for b=1:size(melody,1)
    bd=0;
    for c=1:3
        durs{c}=imag(melody{b,c});
        tones{c}=real(melody{b,c})./durs{c}; % restore tone's frequencies
    end
    % until all notes are gone 
    % 1: spit the first note if it needs to be split
    % 2: consume the first note
    % repeat
    c1d=sum(durs{1});
    while ~isempty(durs{1})
        for c=1:3
            if isempty(durs{c})
                if c1d~=timeSignature, warning('duration of chanel 1 in bar %d is %g and not %g (timeSignature)', b,c1d,timeSignature); end
                error('channel %d in bar %d ends too early',c,b);
            end
            durs1(c)=durs{c}(1);
        end
        mind=min(durs1);
        for c=1:3
            if abs(durs1(c)-mind)>=dur_tol
                % need a split
                durs{c}=[mind durs{c}];
                durs{c}(2)=durs{c}(2)-mind;
                tones{c}=[tones{c}(1) tones{c}];
            end
            ft(c)=tones{c}(1);
        end
        assert(abs(durs{1}(1)-durs{2}(1))<dur_tol);
        assert(abs(durs{1}(1)-durs{3}(1))<dur_tol);
        % consume the first duration / tone
        durations(end+1)=mind;
        pitches(:,end+1)=ft;
        for c=1:3
            tones{c}(1)=[];
            durs{c}(1)=[];
        end
        bd=bd+mind;
    end
    for c=2:3
        if ~isempty(durs{c})
            if c1d~=timeSignature, warning('duration of chanel 1 in bar %d is %g and not %g (timeSignature)', b,c1d,timeSignature); end
            error('channel %d in bar %d is longer than chanel 1',c,b);
        end
    end
    if bd~=timeSignature, warning('duration of bar %d is %g and not %g (timeSignature)',b,bd,timeSignature); end
end

