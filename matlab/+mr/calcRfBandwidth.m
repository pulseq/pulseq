function [bw,fc,spectrum,f,rfs,t]=calcRfBandwidth(rf, cutoff, df, dt)
%calcRfbandwidth Calculate the spectrum of the RF pulse
%   Returns the bandwidth of the pulse (calculated by a simple FFT, e.g. 
%   pesuming a low-angle approximation) and optionally the spectrum and the
%   frequency axis. The default for the optional parameter 'cutoff' is 0.5 
%   The default for the optional parameter 'dw' is 10 Hz. The default for 
%   the optional parameter 'dt' is 1 us. 
%

if nargin<2
    cutoff=0.5;
end
    
if nargin<3
    df=10; % spectral resolution in Hz
end

if nargin<4
    dt=1e-6; % for now default sampling rate is 1Mhz, it's probably too high 
end

tc=mr.calcRfCenter(rf);

% resample the pulse to a resonable time array
nn=round(1/df/dt);
t=(-floor(nn/2):ceil(nn/2)-1)*dt;

rfs=interp1(rf.t-tc,rf.signal.*exp(1i*(rf.phaseOffset+2*pi*rf.freqOffset*rf.t)),t,'linear',0);
spectrum=fftshift(fft(fftshift(rfs)));
f=(-floor(nn/2):ceil(nn/2)-1)*df;

w1=findFlank(f,spectrum,cutoff);
w2=findFlank(f(end:-1:1),spectrum(end:-1:1),cutoff);

bw=w2-w1;
fc=(w2+w1)/2;

% coarse STE scaling -- we normalize to the max of the spectrum, this works
% better with frequency-shifted pulses than the abs(sum(shape)) -- the
% 0-frequency response; yes, we could take the spectrum at fc but this
% would have a problem for non-symmetric pulses...
%s_ref=max(abs(spectrum));
s_ref=interp1(f,abs(spectrum),fc);
spectrum=sin(2*pi*dt*s_ref)*spectrum/s_ref;

end


function xf=findFlank(x,f,c)

m=max(abs(f));
f=abs(f)/m-c;
[~,i]=find(f>0,1);
if i>1
    f0=f(i-1);
    f1=f(i);
    xf=(f1*x(i-1)-f0*x(i))/(f1-f0); 
else
    xf=x(1);
end

end