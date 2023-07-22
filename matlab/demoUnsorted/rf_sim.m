function rf_sim(pulse,prephase_factor,rephase_factor)
% inspired by the example of Dr. Tony Stoecker 

T     = pulse.t*1e3;                                                       % timesteps axis [ms]
dt    = T(2)-T(1);                                                         % sampling interval
F     = linspace(-2*pi*6,2*pi*6,length(T)*2);                                % offset frequencies [rad/ms] 
%w     = pi/2;                                                              % tip-angle [rad]
%shape = @(t) sinc(2*(t-2.5));                                              % sinc RF pulse 
%v     = w / abs(trapz(T,shape(T)));                                        % normalized tip-angle
shapea = pulse.signal/max(abs(pulse.signal))*exp(1i*pulse.phaseOffset);
w      = trapz(pulse.t,pulse.signal)*pi*2;
v      = w/abs(trapz(T,shapea));

%a) STA SOLUTION using fft
N     = round(numel(F)*.5*pi/dt/max(F))*2;                                 % nmuber of points on fft
M_STA = sin(w)*abs(fftshift(fft(shapea,N))) * dt*v/w; 
M_STA = M_STA(N/2+(-numel(F)/2+1:numel(F)/2));                             % extract frequencies of interest

% initialize result vectors
%M_ODE=zeros(size(F)); 
M_ROT=zeros(size(F)); 
Z_ROT=zeros(size(F));
q=zeros(4,length(F));
i=0;                           

for freq = F; i=i+1;                                                      
%%b) ODE SOLUTION                                      
%  bloch = @(t,M) [ freq*M(2);-freq*M(1)+v*shape(t)*M(3);-v*shape(t)*M(2)];
%  [~,M] = ode45(bloch,T,[0;0;1]);                                          % integrate ODE  
%  M_ODE(i)=norm(M(end,1:2));                                               % transv. magn.
%c) ROTATION MATRIX SOLUTION 
  Q=eye(2);                                                                % init Cayley-Klein matrix
  % prephaser / left spoiler
  W = -freq*dt*length(T)*prephase_factor;     % effective field rotation angle
  n = [0 0 1];                               % effective field rotation axis is Z in our case
  a = cos(W/2)-1i*n(3)*sin(W/2);             % Cayley-Klein Parameter (alpha)
  b = (-1i*n(1)+n(2))*sin(W/2);              % Cayley-Klein Parameter (beta)
  Q=[a -b';b a']*Q;                          % "chained" rotation 
  % RF pulse simulation
  for j=1:length(T)%t=T
        t=T(j);
        W = -dt*sqrt(abs(v*shapea(j))^2+freq^2);                            % effective field rotation angle
        n = dt * [v*shapea(j) 0 freq]/abs(W);                               % effective field rotation axis
        a = cos(W/2)-1i*n(3)*sin(W/2);                                     % Cayley-Klein Parameter (alpha)
        b = (-1i*n(1)+n(2))*sin(W/2);                                      % Cayley-Klein Parameter (beta)
        Q=[a -b';b a']*Q;                                                  % chained rotation 
  end
  % rephaser / right spoiler
  W = -freq*dt*length(T)*rephase_factor;     % effective field rotation angle
  n = [0 0 1];                               % effective field rotation axis is Z in our case
  a = cos(W/2)-1i*n(3)*sin(W/2);             % Cayley-Klein Parameter (alpha)
  b = (-1i*n(1)+n(2))*sin(W/2);              % Cayley-Klein Parameter (beta)
  Q=[a -b';b a']*Q;                          % "chained" rotation 

  % summarize results
  M_ROT(i) = abs(2*Q(1,1)'*Q(2,1));                                        % transv. magn.
  Z_ROT(i) = Q(1,1)*Q(2,2)+Q(2,1)*Q(1,2);                                  % long. magn.
  q(1,i)=0.5*real(Q(1,1)+Q(2,2)); % rotation quaternion (for further analysis)
  q(2,i)=0.5*imag(Q(2,1)+Q(1,2)); 
  q(3,i)=0.5*real(Q(1,2)-Q(2,1));
  q(4,i)=0.5*imag(Q(1,1)-Q(2,2));
end

%plot result
F=F/(2*pi);
% ,M_ODE(1:2:end),F(1:2:end),'ko'
figure; plot(F,real(M_STA),F,M_ROT,F,real(Z_ROT),'m','linewidth',1),axis([-4 4 -1.2 1.2])
grid,legend({'M_tSTA','M_tROT','M_lROT'},'fontsize',12)
xlabel('frequency offset / kHz','fontsize',12)
ylabel('magnetisation','fontsize',12)
title(sprintf('flip angle %2.0f degree',w*180/pi),'fontsize',12)

figure; plot(F,180/pi*2*acos(q(1,:))); title('flip angle');
figure; plot(F,180/pi*atan2(q(4,:),(q(2,:).^2+q(3,:).^2).^0.5)); title('theta (90° -> on resonance)');
figure; plot(F,180/pi*atan2(q(3,:),q(4,:)));title('in-plane angle (RF phase)');

m=zeros(4,length(F));
m(4,:)=1;
mrf=quatmultiply(quatconj(q'),quatmultiply(m',q'))';
figure;plot(F,mrf(4,:));axis([-4 4 -1.2 1.2])
figure;plot(F,sum(mrf(2:3,:).^2,1).^0.5);axis([-4 4 -.1 1.2])
figure;plot(F,mrf(2,:),F,mrf(3,:));axis([-4 4 -1.2 1.2])
figure;plot(F,atan2(mrf(3,:),mrf(2,:)));axis([-4 4 -3.2 3.2])

% m=zeros(4,length(F));
% m(2,:)=1;
% mrf=quatmultiply(quatconj(q'),quatmultiply(m',q'))';
% mrfc1=[conv(mrf(2,:),ones(1,3)/3,'same');conv(mrf(3,:),ones(1,3)/3,'same')];
% m=zeros(4,length(F));
% m(3,:)=1;
% mrf=quatmultiply(quatconj(q'),quatmultiply(m',q'))';
% mrfc2=[conv(mrf(2,:),ones(1,3)/3,'same');conv(mrf(3,:),ones(1,3)/3,'same')];
% %figure;plot(F,sum(mrfc1.^2).^0.5,F,sum(mrfc2.^2).^0.5,F,0.5-0.5*real(Z_ROT));axis([-4 4 -0.1 1.1]); title('refocusing profiles vs ''inversion''');
% figure;plot(F,(0.5*sum(mrfc1.^2)+0.5*sum(mrfc2.^2)).^0.5,F,0.5-0.5*real(Z_ROT));axis([-4 4 -0.1 1.1]); title('refocusing profiles vs ''inversion''');
% F;
