function seq = musicToSequence(seq, pitches, durations, varargin)
% populate the Pulseq sequence baased on the provided pitches and durations

persistent parser
validAxes={'xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx'};

if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'mrMusic.musicToSequence';
    parser.addParamValue('barDurationSeconds',4,@isnumeric);
    parser.addParamValue('pulseqUseWave',false,@islogical);
    parser.addParamValue('addDummyRfAdc',true,@islogical);
    parser.addParamValue('amplitudeScale',1,@isnumeric);
    parser.addParamValue('axesOrder',validAxes{1},...
        @(x) any(validatestring(x,validAxes)));
    %parser.addParamValue('defaultLegato',true,@islogical);
end
parse(parser,varargin{:});
opt = parser.Results;

barDurationSeconds=opt.barDurationSeconds;
pulseqUseWave=opt.pulseqUseWave;
addDummyRfAdc=opt.addDummyRfAdc;

gcs={opt.axesOrder(1), opt.axesOrder(2), opt.axesOrder(3)};
if length(opt.amplitudeScale)==1
    assert(abs(opt.amplitudeScale)<=1);
    ga=opt.amplitudeScale*[1,1,1]*seq.sys.maxGrad*0.95; 
else
    if length(opt.amplitudeScale)~=3
        error('amplitudeScale muste be a scalar or a vector of 3 elements');
    end
    assert(all(abs(opt.amplitudeScale)<=1));
    ga=opt.amplitudeScale*seq.sys.maxGrad*0.95;
end

% dummy MR sequence (FID)
if addDummyRfAdc
    % non-selective pulse 
    rf = mr.makeBlockPulse(pi/2,'Duration',0.1e-3, 'system', seq.sys);
    % ADC event
    adc = mr.makeAdc(256,'Duration',3.2e-3, 'system', seq.sys,'delay',seq.sys.adcDeadTime);
    seq.addBlock(rf,mr.makeDelay(10e-3));
    seq.addBlock(adc,mr.makeDelay(20e-3));
    warning('OFF', 'mr:restoreShape');
end

%% generate the wave 
%sample_rate=44100; %Hz
sample_rate=1/seq.gradRasterTime;
dwell_time=1/sample_rate;
notes=max(size(pitches));
%sound_length=floor(noteLen*notes)+1;
%sound_data=zeros(2,sound_length); %preallocate
pitches=[pitches [0;0;0]]; % add trailing 0
mif=min(pitches(pitches>0)); % minimal frequency for gradient scaling

add=[0 0 0];
for n=1:notes
    noteLen=round(barDurationSeconds/dwell_time*durations(n));
    del=mr.makeDelay(noteLen*dwell_time);
    block={del};
    for c=1:3
        cont=(pitches(c,n)==pitches(c,n+1));
        fs=pitches(c,n)*dwell_time; % frequency in the units of inverse samples
        if fs>0
            if cont 
                thisNoteLen=noteLen;
            else
                thisNoteLen=round(floor((noteLen+add(c))*fs*2)/fs/2)-add(c); % round down to a full number of half-cycles
                if (thisNoteLen<0)
                    thisNoteLen=0; 
                    if ~cont 
                        add(c)=0;
                        continue;
                    end
                end
            end
            a=(mif/pitches(c,n));
            if pulseqUseWave
                % arbitrary gradient - smooth harmonic wave
                noteWave=sin_mix(fs,(1:thisNoteLen)+add(c)-0.5)*a*ga(c);
                first=sin_mix(fs,add(c))*a*ga(c);
                if cont
                    last=sin_mix(fs,thisNoteLen+add(c))*a*ga(c);
                else
                    last=0;
                end
                %cs=(n-1)*noteLen;
                %sound_data(1,(cs+1):(cs+thisNoteLen))=sound_data(1,(cs+1):(cs+thisNoteLen))+0.33*noteWave;
                block{end+1}=mr.makeArbitraryGrad(gcs{c},noteWave,seq.sys,'first',first,'last',last);
            else
                % extended trapezoid - saw tooth wave
                note_data=saw(fs,add(c),(thisNoteLen+add(c)));
                if add(c)==0, note_data(2,1)=0; end;
                if ~cont, note_data(2,end)=0; end;
                block{end+1}=mr.makeExtendedTrapezoid(gcs{c},'times',note_data(1,:)*dwell_time,'amplitudes',ga(c)*a*note_data(2,:),'system',seq.sys);
            end
            if cont
                add(c)=add(c)+noteLen;
            else
                add(c)=0;
            end
        else
            add(c)=0;
        end
    end
    seq.addBlock(block);         
end
%sound_data(2,:)=sound_data(1,:);

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;
if ~ok
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

end % end of the function

%%
function out=sin_mix(f,x) 
    w=[0.6 0.3 0.1];
    m=max(w(1)*sin(2*pi*[0:0.01:1])+w(2)*sin(4*pi*[0:0.01:1])+w(3)*sin(8*pi*[0:0.01:1]));
    w=w/m;
    out=w(1)*sin(2*pi*f*x)+w(2)*sin(4*pi*f*x)+w(3)*sin(8*pi*f*x);
end

function out=saw(f,first,last) 
    points=((floor(first*2*f-1/2):ceil(last*2*f-1/2))+1/2)/2/f;
    assert(points(1)<=first);
    assert(points(end)>=last);
    % round to raster
    pointsr=round(points-first);
    if pointsr(1)==-1, pointsr(1)=-2; end
    if pointsr(2)==1, pointsr(2)=2; end
    if pointsr(end)==ceil(last-first)+1
        pointsr(end)=ceil(last-first)+2;
    end
    lenr=ceil(last-first);
    if pointsr(end-1)==lenr-1
        pointsr(end-1)=lenr-2;
    end
    %
    % figure out whether point1 is +1 or -1
    if points(1)*f-floor(points(1)*f) >= 0.5, o=0; else o=1; end
    values=(-1).^(o+(1:length(points)));
    % 
    if pointsr(1)<0 && pointsr(2)>0
        values(1)= (pointsr(1)*values(1) - pointsr(2)*values(2)) / (pointsr(2)-pointsr(1));
        pointsr(1)=0;
    elseif pointsr(2)==0
        pointsr(1)=[];
        values(1)=[];
    end
    assert(pointsr(1)==0);
    % 
    if pointsr(end-1)<lenr && pointsr(end)>lenr
        values(end)= ((lenr-pointsr(end-1))*values(end) - (lenr-pointsr(end))*values(end-1))/(pointsr(end)-pointsr(end-1));
        pointsr(end)=lenr;
    elseif pointsr(end-1)==lenr
        pointsr(end)=[];
        values(end)=[];
    end
    assert(pointsr(end)==lenr);
    %
    out=[pointsr;values];
end
