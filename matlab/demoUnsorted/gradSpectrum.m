% example script to plot gradient frequency spectrum
% expects: 
%  "seq" object to be already populated 
%  "sys" object to contain system specs
% furthermore, on Siemens you need the *.asc file for your gradient system


dt=sys.gradRasterTime; % time raster
fmax=10000; %10kHz
nwin=5000; % 0.05s
os=3; % frequency oversampling for prettier peaks
ascName=[]; % this disables the display of the system's resonance frequences
%ascName='idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'; % 3T prisma
%ascName='idea/asc/MP_GradSys_P034_X60.asc'; % 3T cima.X

if ischar(ascName)
    ascData=mr.Siemens.readasc(ascName);
end

faxis=(0:(nwin/2-1))/nwin/dt/os;
nfmax=sum(faxis<=fmax);

wave_data=seq.waveforms_and_times();
tmax=max([wave_data{1}(1,end) wave_data{2}(1,end) wave_data{3}(1,end)]);
nt=ceil(tmax/dt);
tmax=nt*dt;

gw=zeros(3,nt);
for i=1:3
    gw(i,:)=interp1(wave_data{i}(1,:),wave_data{i}(2,:),((1:nt)-0.5)*dt,'linear',0);
    % alternative (to be checked in the future)
    % it is actually much more appropriate to calculate the spectrium of
    % the derivative(!) of the gradient wave form and not the waveform
    % itself, at least for the cound/noise of the gradients...
    %gw(i,1:end-1)=diff(interp1(wave_data{i}(1,:),wave_data{i}(2,:),((1:nt)-0.5)*dt,'linear',0));
end

gs=[];

ng=size(gw,1);
for g=1:ng
    x=gw(g,:);
    nx = length(x);

    nx=ceil(nx/nwin)*nwin;
    if nx>length(x) 
        x=[x, zeros(1,nx-length(x))]; % zerofill
    end

    nseg1=nx/nwin; 
    xseg=zeros(nseg1*2-1,nwin*os); 

    xseg(1:2:end,1:nwin)=reshape(x,[nwin,nseg1])';
    if nseg1>1
        xseg(2:2:end,1:nwin)=reshape(x(1+nwin/2:end-nwin/2),[nwin,nseg1-1])';
    end

    xseg_dc=mean(xseg,2);
    xseg=xseg-xseg_dc(:,ones(1,nwin*os));

    if nseg1>1 % WARNING: thisintroduces inconsistency between short and long sequences in term os the peak amplitudes
        cwin=0.5*(1-cos(2*pi*(1:nwin)/nwin));
        xseg(:,1:nwin)=xseg(:,1:nwin).*cwin(ones(size(xseg,1),1),:);
    end

    fseg=abs(fft(xseg,[],2));
    fseg=fseg(:,1:end/2); 
        
    if nseg1>1 
        gs = [gs; mean(fseg.^2).^0.5]; % sos
        %figure; plot(faxis(1:nfmax),sum(fseg(:,1:nfmax).^2).^0.5);
    else
        gs = [gs; abs(fseg)]; % add abs
    end    
end

figure; plot(faxis(1:nfmax),gs(:,1:nfmax));
%%figure; plot(faxis(1:nfmax),sum(gs(:,1:nfmax))); % abs-sum
hold on; plot(faxis(1:nfmax),sum(gs(:,1:nfmax).^2).^0.5); % sos
xlabel('frequency / Hz');

if ischar(ascName)
    hold on; 
    for i=1:length(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency)
        if ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)>0
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i),'-');
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)-ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i)/2,'--');
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)+ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i)/2,'--');
        end
    end
end

legend({'Gx','Gy','Gz','Gtot'});

%nov = floor(nsc/2);
%nff = max(256,2^nextpow2(nsc));

%t = spectrogram(x,hamming(nsc),nov);%,nff);
%t = spectrogram(x,rectwin(nsc),nov);
%maxerr = max(abs(abs(t(:))-abs(s(:))))
