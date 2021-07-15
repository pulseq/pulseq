function [bw,spectrum,w]=calcRfBandwidth(rf, cutoff)
%calcRfbandwidth Calculate the spectrum of the RF pulse
%   Returns the bandwidth of the pulse (calculated by a simple FFT, e.g. 
%   pesuming a low-angle approximation) and optionally the spectrum and the
%   frequency axis. The default for the optional parameter 'cutoff' is 0.5 
%

if nargin==1
    cutoff=0.5;
end
    
tc=mr.calcRfCenter(rf);

% resample the pulse to a resonable time array
dw=10; %Hz
dt=1e-6; % for now 1Mhz, it's probably too high 
nn=round(1/dw/dt);
tt=(-floor(nn/2):ceil(nn/2)-1)*dt;

rfs=interp1(rf.t-tc,rf.signal,tt,'linear',0);
spectrum=fftshift(fft(fftshift(rfs)));
w=(-floor(nn/2):ceil(nn/2)-1)*dw;

w1=findFlank(w,spectrum,cutoff);
w2=findFlank(w(end:-1:1),spectrum(end:-1:1),cutoff);

bw=w2-w1;

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