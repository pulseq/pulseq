% a basic 3D ZTE/PETRA sequence
% achieves "TE" about 70 us and possibly below

%% high-level sequence parameters

fov=256e-3; dx=2e-3;            % Define FOV and resolution
alpha=4;                        % flip angle
Nr=300;                         % number of readout points (with some oversampling)
R= 8;                           % acceleration/undersampling (for the outer shell)
R_inner= 1;                     % acceleration/undersampling (angular direction of the inner area)
xSpoil=3;%0.6;  %0.6 for 2mm       % the amount of spoiling after the end of the readout (used to ramp to the next point)

Kmax=1/2/dx;
dK=1/fov;

%% more detailed and derived params

rf_duration=10e-6;              % duration of the excitation pulse
ro_duration=300e-6;             % read-out time: controls RO bandwidth and T2-blurring
minRF_to_ADC_time=50e-6;        % the parameter wich defines TE together with ro_discard
%ro_discard=4;                   % how many ADC samples are contaminated by RF switching artifacts and alike
rfSpoilingInc=117;              % RF spoiling increment

% system limits
sys = mr.opts('MaxGrad', 36, 'GradUnit', 'mT/m', ...
    'MaxSlew', 180, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'gamma',42.576e6); % 1H: 42.576e6  23Na: 11.262e6

seq=mr.Sequence(sys);           % Create a new sequence object
seq_sar=mr.Sequence(sys);       % Create an auxillary sequence object for SAR testing

%% create main sequence elements

% %create alpha-degree block pulse 
% rf = mr.makeBlockPulse(alpha*pi/180,'Duration',rf_duration,'system',sys);
%or create alpha-degree gaussian pulse 
rf = mr.makeGaussPulse(alpha*pi/180,'Duration',rf_duration,'timeBwProduct',3,'system',sys);
Tenc=rf_duration/2+minRF_to_ADC_time+ro_duration; %encoding time
Tg=sys.rfDeadTime+rf_duration/2+Tenc+sys.adcDeadTime; % constant gradient time
Tt=ceil((Tenc*(1+xSpoil)-Tg)/sys.gradRasterTime)*sys.gradRasterTime; % transition time
Ag=Kmax/Tenc; % read gradient grdient amplitude
% % Gr=mr.makeTrapezoid('z','Amplitude',Ag,'flatTime',Tg,'riseTime',0); % this graient has no ramps, I wonder if the Matlab mr library would like it...
% Gr=mr.makeExtendedTrapezoid('z','times',[0 Tg],'amplitudes',[Ag Ag]); % this is a constant graient with no ramps
TR=Tg+Tt;

adc=mr.makeAdc(Nr,'Duration', ro_duration, 'Delay', sys.rfDeadTime+rf_duration+minRF_to_ADC_time);

SamplesBookkeeping=[];

%rfbw=1/rf_duration;
rfbw=mr.calcRfBandwidth(rf);
fprintf('Read gradient amplitude %g mT/m, effective "slice thinckess" %f mm\n', 1e3*Ag/sys.gamma, 1/Ag*rfbw*1e3);
FO=50e-3*Ag;
NF=0; % TR is broken now for NF~=0...

%% generate the sampling set on a surface of a sphere
[phi, theta, im]=spherical_samples(Kmax,dK,R);
Ns=length(phi);
SamplesBookkeeping=[SamplesBookkeeping Ns];

%% main ZTE loop
fprintf('Populating ZTE loop (%d TRs)\n', Ns);
tic;
populate_subsequence(sys,seq,rf,adc,phi,theta,im,Ns,Ag,TR,Tt,NF,FO);
toc

%% SPI loop
KstartZTE=Ag*(rf_duration/2+minRF_to_ADC_time+adc.dwell);
nKspi=floor(KstartZTE/dK);
dKspi=KstartZTE/(nKspi+1);
Tenc=rf_duration/2+minRF_to_ADC_time+adc.dwell/2;
fprintf('Populating SPI loop (%d spheres)\n', nKspi);
tic;
for s=nKspi:-1:1
    % generate the sampling set on a surface of a sphere
    [phi, theta, im]=spherical_samples(dKspi*s,dKspi,R_inner); % normally no acceleration
    Ns=length(phi);
    SamplesBookkeeping=[SamplesBookkeeping Ns];
    % update gradients
    Ag=dKspi*s/Tenc; % read gradient grdient amplitude
    % actually create the next sampling sphere
    fprintf('Populating sphere %d (%d TRs)\n',s,Ns);
    fprintf('Effective "slice thinckess" %f mm\n',  1/Ag*rfbw*1e3);
    populate_subsequence(sys,seq,rf,adc,phi,theta,im,Ns,Ag,TR,Tt,NF,FO);
end
% sample the centere of k-space 
seq.addBlock(mr.makeDelay(Tt));
seq.addBlock(rf, adc, mr.makeDelay(TR-Tt));
SamplesBookkeeping=[SamplesBookkeeping 1];
toc
fprintf('Total number of SPI samples: %d; a Cartesian patch would require %d\n', sum(SamplesBookkeeping(2:end)), ceil(nKspi^3*pi*4/3));

% %% OLD!!! SPI loop
% Kspi=ceil(Ag*(rf_duration/2+minRF_to_ADC_time+adc.dwell)/dK);
% rf.delay=rf.delay+Tt;
% Tspi1=mr.calcDuration(rf)+sys.rfRingdownTime;
% adc.delay=adc.delay+Tt-Tspi1;
% Tspi2=TR-Tspi1;
% 
% g=mr.makeTrapezoid('x','system',sys,'area',dK*Kspi,'duration',minRF_to_ADC_time);
% gx=g;gx.channel='x';
% gy=g;gy.channel='y';
% gz=g;gz.channel='z';
% 
% gxr=g;gxr.channel='x';
% gyr=g;gyr.channel='y';
% gzSpoil=mr.makeTrapezoid('x','system',sys,'area',Kmax*(1+xSpoil),'duration',minRF_to_ADC_time);
% 
% fprintf('Populating SPI loop (~%d TRs)\n', round(pi/6*((2*Kspi+1)^3)));
% tic;
% 
% % SPI loop itself
% gxr.amplitude=0;
% gyr.amplitude=0;
% for k=-Kspi:Kspi
%     for j=-Kspi:Kspi
%         for i=-Kspi:Kspi
%             if i^2+j^2+k^2>Kspi^2
%                 continue;
%             end
%             % refocusing and spoiling for the previous TR
%             seq.addBlock( rf, gxr, gyr, gzSpoil, mr.makeDelay(Tspi1));
%             % new TR
%             gx.amplitude=g.amplitude/Kspi*i;
%             gy.amplitude=g.amplitude/Kspi*j;
%             gy.amplitude=g.amplitude/Kspi*k;
%             seq.addBlock( rf, adc, gx, gy, gz, mr.makeDelay(Tspi2));
%             gxr.amplitude=-gx.amplitude;
%             gyr.amplitude=-gy.amplitude;
%         end
%     end
% end
% toc

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%%
%seq.plot('TimeRange',[0 1]);

seq.setDefinition('FOV', [fov fov fov]);
seq.setDefinition('Name', 'petra');
seq.setDefinition('SamplesPerShell', SamplesBookkeeping);

seq.write('zte_petra.seq');
return

%% create an RF-only version of the sequence (e.g. for the SAR or signal evolution testing)

tic;
[total_duration, total_numBlocks]=seq.duration();
for iB=1:total_numBlocks
    b=seq.getBlock(iB);
    bd=seq.blockDurations(iB);
    bs={mr.makeDelay(bd)};
    if ~isempty(b.rf)
        bs{end+1}=b.rf;
    end
    if ~isempty(b.adc)
        bs{end+1}=b.adc;
    end
    seq_sar.addBlock(bs);
end
toc
seq_sar.write('zte_petra_sar.seq');
return

% %% test binary storing
% 
% seq.writeBinary('zte_petra.bin');
% seq_bin=mr.Sequence();          
% seq_bin.readBinary('zte_petra.bin');
% seq_bin.write('zte_petra_bin.seq');
% return

%% visualize the 3D k-space 
tic;
[kfa,~,kf]=seq.calculateKspacePP();
toc

%
figure;plot3(kf(1,:),kf(2,:),kf(3,:));
hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');

%% local functions
function [phi, theta, im]=spherical_samples(Kr,dK,R)
% the number of samples equals the ceil of the area of the sphere divided
% by the area around every sample
Ns=ceil(4*pi*((Kr/dK)^2)/R); 
np=0:(Ns-1);
alpha_gold=pi*(3-sqrt(5));
phi=np*alpha_gold;
%theta=0.5*pi*sqrt(np/Ns); % from Davide Piccini  https://doi.org/10.1002/mrm.22898
theta=acos(1-2*np/(Ns-1));  % from  Anton Semechko (2020). Suite of functions to perform uniform sampling of a sphere (https://github.com/AntonSemechko/S2-Sampling-Toolbox), GitHub. Retrieved October 3, 2020. 

xp=sin(theta).*cos(phi);
yp=sin(theta).*sin(phi);
zp=cos(theta);
%figure; sphere; colormap gray;
%hold on; plot3(xp,yp,zp,'.');

% looking for the optimal interleaving factor
nm=round(Ns/2); % middle of the trajetory
sr=round(sqrt(Ns));
v0=[xp(nm); yp(nm); zp(nm)];
v=[xp(nm+(1:sr)); yp(nm+(1:sr)); zp(nm+(1:sr))];
%figure;plot(vecnorm((v-v0(:,ones(size(v,2),1)))));
[dKm,im]=min(vecnorm((v-v0(:,ones(size(v,2),1)))));

% fprintf('requested dK=%g achieved minimal dK=%g, min acceleration: %g, Ns=%d\n', dK/sqrt(R), dKm*Kr, (dKm*Kr/dK)^2,Ns);

% v=[xp; yp; zp]*Kr/dK/sqrt(R)*100;
% md=zeros(1,Ns);
% d=zeros(1,Ns);
% for i=1:Ns
%     d=vecnorm(v(:,i*ones(1,Ns))-v);
%     d(i)=NaN;
%     md(i)=min(d);
% end
% fprintf('dKmin=%g%%, dKmed=%g%% dKmax=%g%%\n', min(md),median(md),max(md));

end

function populate_subsequence(sys,seq,rf,adc,phi,theta,im,Ns,Ag,TR,Tt,FN, FO)
Azc=Ag*(TR-Tt)/(TR+Tt);  %*0.35;
%Azc=0; % for MoCo we need "no-gradient" event blocks to be able to apply updates

Gr=mr.makeExtendedTrapezoid('z','times',[0 TR-Tt],'amplitudes',[Ag Ag]); % this is a constant graient with no ramps

% pre-ramp the gradient to Azc
if abs(Azc)>eps
    Tpr=max(2,ceil(Azc/sys.maxSlew/sys.gradRasterTime))*sys.gradRasterTime;
    assert(Tpr<=TR);
    % this "dummy TR" doe not have an RF pulse, which is not good when it is called for the inner shells... but otherwise there would be no enough spoiling... 
    seq.addBlock(mr.align('right',mr.makeDelay(TR),mr.makeExtendedTrapezoid('z','system',sys,'times',[0,Tpr],'amplitudes',[0,Azc])));
end

% the loop itself
for j=1:im
    Glast=struct('x',0,'y',0,'z',Azc);
    for i=j:im:Ns
        Gcr=mr.rotate('z',phi(i),mr.rotate('y',theta(i),Gr));
        Gcurr=struct('x',0,'y',0,'z',0);
        for g=1:length(Gcr)
            Gcurr.(Gcr{g}.channel)=Gcr{g}.waveform(1);
        end
        seq.addBlock( ...
            mr.makeExtendedTrapezoid('x','system',sys,'times',[0,Tt],'amplitudes',[Glast.x,Gcurr.x]), ...
            mr.makeExtendedTrapezoid('y','system',sys,'times',[0,Tt],'amplitudes',[Glast.y,Gcurr.y]), ...
            mr.makeExtendedTrapezoid('z','system',sys,'times',[0,Tt],'amplitudes',[Glast.z,Gcurr.z]));
        for f=-FN:FN % 'slices' loop (frequency offsets)
            rf.freqOffset=f*FO;
            rf.phaseOffset=-2*pi*f*FO*rf.t(end)/2;
            seq.addBlock( [{rf, adc}, Gcr] );
        end
        Glast=Gcurr;
    end
    rf_aux=rf;
    rf_aux.delay=rf.delay+Tt;
    if (j==im)
        Azc=0; % on the last interleave we ramp down to 0
    end
    seq.addBlock( rf_aux, ...
        mr.makeExtendedTrapezoid('x','system',sys,'times',[0,TR],'amplitudes',[Glast.x,0]), ...
        mr.makeExtendedTrapezoid('y','system',sys,'times',[0,TR],'amplitudes',[Glast.y,0]), ...
        mr.makeExtendedTrapezoid('z','system',sys,'times',[0,TR],'amplitudes',[Glast.z,Azc])); 
end
end

