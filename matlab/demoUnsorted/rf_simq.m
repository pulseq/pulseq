function [M_z,M_x_y,F]=rf_simq(pulse,prephase_factor,rephase_factor)
% inspired by the example by Dr. Tony Stoecker
% (https://github.com/stoeckert/mr-simu-example-ismrm19)
% rewritten to quaternions and vectorized for performance by MZ 
% Q in the name stands for both "quick" and "quaternion" 

T     = pulse.t*1e3;                                                       % timesteps axis [ms]
dt    = T(2)-T(1);                                                         % sampling interval
F     = linspace(-2*pi*4,2*pi*4,length(T)*2)';                             % offset frequencies [rad/ms] 

shapea = pulse.signal/max(abs(pulse.signal))*exp(1i*pulse.phaseOffset);
w      = trapz(pulse.t,pulse.signal)*pi*2;
v      = w/abs(trapz(T,shapea));

%a) STA SOLUTION using fft (ignoring pulse phase) -- just debugging
N     = round(numel(F)*.5*pi/dt/max(F))*2;                                 % numuber of points on fft
M_STA = sin(w)*abs(fftshift(fft(shapea,N))) * dt*v/w; 
M_STA = M_STA(N/2+(-numel(F)/2+1:numel(F)/2));                             % extract frequencies of interest

% intialize result vectors
M_ROT=zeros(size(F)); 
Z_ROT=zeros(size(F));
sf=size(F);
q=zeros(sf(1),4);
q(:,1)=1; % init rotation quaternions

% prephaser / left spoiler
W = -F*dt*length(T)*prephase_factor;     % effective field rotation angle
Q = [cos(W/2) zeros(sf) zeros(sf) sin(W/2)];
q=quatmultiply(q,Q);

% RF pulse simulation
for j=1:length(T)
    W = -dt*sqrt(abs(v*shapea(j))^2+F.^2); % effective field rotation angles
    n = dt * [real(v*shapea(j))*ones(sf) imag(v*shapea(j))*ones(sf) F]./abs(W); % effective field rotation axes
    Q = [cos(W/2) sin(W/2).*n];
    q=quatmultiply(q,Q);
end

% rephaser / right spoiler / refocusing pulse
W = -F*dt*length(T)*rephase_factor;     % effective field rotation angle
Q = [cos(W/2) zeros(sf) zeros(sf) sin(W/2)];
q=quatmultiply(q,Q);

%   % summarize results
%   M_ROT(i) = abs(2*Q(1,1)'*Q(2,1));                                        % transv. magn.
%   Z_ROT(i) = Q(1,1)*Q(2,2)+Q(2,1)*Q(1,2);                                  % long. magn.
%   q(1,i)=0.5*real(Q(1,1)+Q(2,2)); % rotation quaternion (for futher anlysis)
%   q(2,i)=0.5*imag(Q(2,1)+Q(1,2)); 
%   q(3,i)=0.5*real(Q(1,2)-Q(2,1));
%   q(4,i)=0.5*imag(Q(1,1)-Q(2,2));

% %plot result
F=F/(2*pi);
% % ,M_ODE(1:2:end),F(1:2:end),'ko'
% figure; plot(F,real(M_STA),F,M_ROT,F,real(Z_ROT),'m','linewidth',1),axis([-4 4 -1.2 1.2])
% grid,legend({'M_tSTA','M_tROT','M_lROT'},'fontsize',12)
% xlabel('frequency offset / kHz','fontsize',12)
% ylabel('magnetisation','fontsize',12)
% title(sprintf('flip angle %2.0f degree',w*180/pi),'fontsize',12)

%figure; plot(F,180/pi*2*acos(q(:,1))); title('flip angle');
%figure; plot(F,180/pi*atan2(q(:,4),(q(:,2).^2+q(:,3).^2).^0.5)); title('theta (90° -> on resonance)');
%figure; plot(F,180/pi*atan2(q(:,3),q(:,4)));title('in-plane angle (RF phase)');

m=zeros(sf(1),4);
m(:,4)=1;
m0rf=quatmultiply(quatconj(q),quatmultiply(m,q));
M_z=m0rf(:,4);
M_x_y=m0rf(:,2)+1i*m0rf(:,3);
figure;plot(F,M_z);axis([-4 4 -1.2 1.2]);title('exitation M_z');
figure;plot(F,abs(M_x_y));axis([-4 4 -.1 1.2]);title('excitation M_x_y');
figure;plot(F,real(M_x_y),F,imag(M_x_y));axis([-4 4 -1.2 1.2]);title('excitation M_x and M_y');
figure;plot(F,angle(M_x_y));axis([-4 4 -3.2 3.2]);title('excitation phase of M_x_y');

figure;plot(F,real(M_STA),F,imag(M_STA));axis([-4 4 -1.2 1.2]);title('M_x and M_y for STA');

cl=13; % convolution length
m=zeros(sf(1),4);
m(:,2)=1;
mxrf=quatmultiply(quatconj(q),quatmultiply(m,q));
mxrf=mxrf(:,2)+1i*mxrf(:,3);
mxrfc=conv(mxrf,ones(cl,1)/cl,'same');
m=zeros(sf(1),4);
m(:,3)=1;
myrf=quatmultiply(quatconj(q),quatmultiply(m,q));
myrf=myrf(:,2)+1i*myrf(:,3);
myrfc=conv(myrf,ones(cl,1)/cl,'same');
%figure;plot(F,sum(mrfc1.^2).^0.5,F,sum(mrfc2.^2).^0.5,F,0.5-0.5*real(Z_ROT));axis([-4 4 -0.1 1.1]); title('refocusing profiles vs ''inversion''');
%figure;plot(F,(0.5*sum(mrfc1.^2,2)+0.5*sum(mrfc2.^2,2)).^0.5,F,0.5-0.5*m0rf(:,4));axis([-4 4 -0.1 1.1]); title('refocusing profiles vs (inverted and scalled) M_z');
figure;plot(F,abs(abs(mxrfc)+1i*abs(myrfc))/2^0.5,F,0.5-0.5*m0rf(:,4));axis([-4 4 -0.1 1.1]); title('refocusing profiles vs (inverted and scalled) M_z');
% refocusing phase
%figure;plot(F,angle(mxrf)/2-(angle(myrf)+pi/2)/2);axis([-4 4 -1.6 1.6]); title('refocusing phase inconsistency');
%figure;plot(F,angle(mxrf)/2+(angle(myrf)+pi/2)/2);axis([-4 4 -1.6 1.6]); title('mean refocusing phase');
figure;plot(F,angle(mxrf),F,(angle(myrf*1i)));axis([-4 4 -3.3 3.3]); title('refocusing phase indicator');
[~,iF_min]=min(abs(F));
%rph=angle(mxrf(iF_min))/2+(angle(myrf(iF_min)+pi/2)/2);% not robust to phase wraps
rph=angle(mxrf(iF_min)+myrf(iF_min)*1i)/2;
mxref=conj(1*exp(-1i*rph))*exp(1i*rph);
myref=conj(1i*exp(-1i*rph))*exp(1i*rph);
figure;plot(F,1-(0.25*abs(mxrf-mxref).^2+0.25*abs(myrf-myref).^2).^0.5);axis([-4 4 -0.1 1.1]); title('refocusing efficiency assuming fixed phase');
figure; plot(F,abs(mxrf+myrf*1i)/2);axis([-4 4 -0.1 1.1]); title('generic refocusing efficiency');
figure; plot(F,angle(mxrf+myrf*1i)/2);axis([-4 4 -0.1 1.1]); title('generic refocusing efficiency phase');
