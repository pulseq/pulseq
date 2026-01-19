function [R, Rax, F] = gradSpectrum(obj, FB, fmax, doPlot)
% function [R, Rax, F] = seq.gradSpectrum(FB, [fmax = 3000], [doPlot = true])
%
% Get (and optionally plot) frequency response of a Pulseq sequence.
%
% Input 'FB' is either a struct array with the forbidden bands, 
% or the name of a Siemens ASC file. 
%
% Inputs 
%   seq                   Pulseq sequence object 
%   FB     [num_bands]     Forbidden frequency bands (struct array)
%                         FB(1).freq   center frequency of first forbidden band
%                         FB(1).bw     bandwidth of first forbidden band
%   FB     string          Siemens ASC file name
%   fmax   [1]             Max frequency range
%   doPlot TRUE/false      Plot or just return values
%
% Outputs
%   R      [n]     frequency response (root-sum-of-squares of all axes)
%   Rax    [3 n]   frequency response for individual axes
%   F      [n]     frequency locations (Hz)
%
% Function version of demoUnsorted/gradSpectrum.m

% defaults
if nargin < 4
    doPlot = true;
end
if nargin < 3
    fmax = 3000;
end

% Read ASC file if provided
if ischar(FB)
    ascData=mr.Siemens.readasc(FB);
    clear FB
    for i=1:length(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency)
        if (ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)>0)
            FB(i).freq = ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i);
            FB(i).bw = ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i);
        end
    end
end

% Calculate spectrum/spectrogramm
dt=obj.sys.gradRasterTime; % time raster
nwin=5000; % 0.05s
os=3; % frequency oversampling for prettier peaks

faxis=(0:(nwin/2-1))/nwin/dt/os;
nfmax=sum(faxis<=fmax);

wave_data=obj.waveforms_and_times();
ng=length(wave_data);
tmax=0;
for i=1:ng
    if ~isempty(wave_data{i})
        tmax=max(tmax, wave_data{i}(1,end));
    end
end
if tmax==0
    error('Empty sequence passed to gradSpectrum()');
end
nt=ceil(tmax/dt);
tmax=nt*dt;

gw=zeros(ng,nt);
for i=1:ng
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

    if nseg1>1 % WARNING: this introduces an inconsistency between short and long sequences in term os the peak amplitudes
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

% Define return values
F = faxis(1:nfmax);
Rax = gs(:,1:nfmax);
R = sum(gs(:,1:nfmax).^2).^0.5;

% plot
if ~doPlot
    return;
end

figure; plot(faxis(1:nfmax),gs(:,1:nfmax));
hold on; plot(faxis(1:nfmax),sum(gs(:,1:nfmax).^2).^0.5); % sos
xlabel('frequency / Hz');

% alternative "stained glass" plots
% clr = repmat('rgb', [1 5]);
% if ~isempty(FB)
%     for i=1:length(FB)
%         if FB(i).freq > 0
%             l = FB(i).freq-FB(i).bw/2;
%             r = FB(i).freq+FB(i).bw/2;
%             t = max(R);
%             b = 0;
%             h = fill([l r r l], [b b t t], clr(i));
%             h.FaceAlpha = 0.15;
%         end
%     end
% end

if ~isempty(FB)
    for i=1:length(FB)
        if FB(i).freq > 0
            xline(FB(i).freq,'-');
            xline(FB(i).freq-FB(i).bw/2,'--');
            xline(FB(i).freq+FB(i).bw/2,'--');
        end
    end
end

legend({'Gx','Gy','Gz','Gtot'});

%nov = floor(nsc/2);
%nff = max(256,2^nextpow2(nsc));

%t = spectrogram(x,hamming(nsc),nov);%,nff);
%t = spectrogram(x,rectwin(nsc),nov);
%maxerr = max(abs(abs(t(:))-abs(s(:))))
