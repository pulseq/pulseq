% this is an experimentaal spiral sequence

seq=mr.Sequence();          % Create a new sequence object
fov=256e-3; Nx=128; Ny=128;  % Define FOV and resolution
thickness=3e-3;             % slice thinckness
Nslices=1;
Oversampling=2; % by looking at the periphery of the spiral I would say it needs to be at least 2

% Set system limits
lims = mr.opts('MaxGrad',30,'GradUnit','mT/m',...
    'MaxSlew',140,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6, 'adcDeadTime', 10e-6);  

% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',lims,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

% define k-space parameters
deltak=1/fov;
kRadius = round(Nx/2);
kSamples=round(2*pi*kRadius)*Oversampling;
readoutTime = 4.2e-4;

% calculate a raw Archimedian spiral trajectory
clear ka;
ka(kRadius*kSamples+1)=1i; % init as complex
for c=0:kRadius*kSamples
    r=deltak*c/kSamples;
    a=mod(c,kSamples)*2*pi/kSamples;
    ka(c+1)=r*exp(1i*a);
end
ka=[real(ka); imag(ka)];
% calculate gradients and slew rates
[ga, sa]=mr.traj2grad(ka);

% limit analysis
safety_magrin=0.97; % we need that  otherwise we just about violate the slew rate due to the rounding errors
dt_gcomp=abs(ga)/(lims.maxGrad*safety_magrin)*lims.gradRasterTime;
dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(lims.maxGrad*safety_magrin)*lims.gradRasterTime;
dt_scomp=sqrt(abs(sa)/(lims.maxSlew*safety_magrin))*lims.gradRasterTime;
dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(lims.maxSlew*safety_magrin))*lims.gradRasterTime;

figure;plot([dt_gabs; max(dt_gcomp); dt_sabs; max(dt_scomp)]');

dt_smooth=max([dt_gabs;dt_sabs]);
dt_rough=max([dt_gcomp;dt_scomp]);

% apply the lower limit not to lose the trajectory detail
dt_min=4*lims.gradRasterTime/kSamples; % we want at least 4 points per revolution
dt_smooth0=dt_smooth;
dt_rough0=dt_rough;
dt_smooth(dt_smooth<dt_min)=dt_min;
dt_rough(dt_rough<dt_min)=dt_min;

figure;plot([dt_smooth0; dt_smooth; dt_rough0; dt_rough]');

t_smooth=[0 cumsum(dt_smooth,2)];
t_rough=[0 cumsum(dt_rough,2)];

kopt_smooth=interp1(t_smooth, ka', (0:floor(t_smooth(end)/lims.gradRasterTime))*lims.gradRasterTime)';
kopt_rough=interp1(t_rough, ka', (0:floor(t_rough(end)/lims.gradRasterTime))*lims.gradRasterTime)';

% analyze what we've got
fprintf('duration orig %d us\n', round(1e6*lims.gradRasterTime*length(ka)));
fprintf('duration smooth %d us\n', round(1e6*lims.gradRasterTime*length(kopt_smooth)));
fprintf('duration rough %d us\n', round(1e6*lims.gradRasterTime*length(kopt_rough)));

[gos, sos]=mr.traj2grad(kopt_smooth);
[gor, sor]=mr.traj2grad(kopt_rough);

figure;plot([gos;abs(gos(1,:)+1i*gos(2,:))]');title('gradient with smooth (abs) constraint')
figure;plot([gor;abs(gor(1,:)+1i*gor(2,:))]');title('gradient with rough (component) constraint')

figure;plot([sos;abs(sos(1,:)+1i*sos(2,:))]');title('slew rate with smooth (abs) constraint')
figure;plot([sor;abs(sor(1,:)+1i*sor(2,:))]');title('slew rate with rough (component) constraint')

% Define gradients and ADC events
spiral_grad_shape=gos;
% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',lims,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2);

% calculate ADC
% round-down dwell time to 10 ns
adcTime = lims.gradRasterTime*size(spiral_grad_shape,2);
% actually it is trickier than that: the (Siemens) interpreter sequence 
% per default will try to split the trajectory into segments <=1000 samples
% and every of these segments will have to have duration aligned to the
% gradient raster time
adcSamplesPerSegment=1000; % you may need to play with this number to fill the entire trajectory
adcSamplesDesired=kRadius*kSamples;
adcSegments=round(adcSamplesDesired/adcSamplesPerSegment);
adcSamples=adcSegments*adcSamplesPerSegment;
adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 100 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
if mod(adcSegmentDuration, lims.gradRasterTime)>eps 
    error('ADC segmentation model results in incorrect segment duration');
end
% update segment count
adcSegments=floor(adcTime/adcSegmentDuration);
adcSamples=adcSegments*adcSamplesPerSegment;
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',mr.calcDuration(gzReph));%lims.adcDeadTime);

% extend spiral_grad_shape by repeating the last sample
% this is needed to accomodate for the ADC tuning delay
spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end)];

% readout grad 
gx = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:),'Delay',mr.calcDuration(gzReph));
gy = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:),'Delay',mr.calcDuration(gzReph));

% spoilers
gz_spoil=mr.makeTrapezoid('z',lims,'Area',deltak*Nx*4);
gx_spoil=mr.makeExtendedTrapezoid('x','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(1,end),0]); %todo: make a really good spoiler
gy_spoil=mr.makeExtendedTrapezoid('y','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(2,end),0]); %todo: make a really good spoiler

% because of the ADC alignment requirements the sampling window possibly
% extends past the end of the trajectory (these points will have to be
% discarded in the reconstruction, which is no problem). However, the
% ramp-down parts and the Z-spoiler now have to be added to the readout
% block otherwise there will be a gap inbetween
% gz_spoil.delay=mr.calcDuration(gx);
% gx_spoil.delay=gz_spoil.delay;
% gy_spoil.delay=gz_spoil.delay;
% gx_combined=mr.addGradients([gx,gx_spoil], lims);
% gy_combined=mr.addGradients([gy,gy_spoil], lims);
% gz_combined=mr.addGradients([gzReph,gz_spoil], lims);
 
% Define sequence blocks
for s=1:Nslices
    rf.freqOffset=gz.amplitude*thickness*(s-1-(Nslices-1)/2);
    seq.addBlock(rf,gz);
    seq.addBlock(gzReph,gx,gy,adc);
    seq.addBlock(gx_spoil,gy_spoil,gz_spoil);
    %seq.addBlock(gx_combined,gy_combined,gz_combined,adc);
end

seq.setDefinition('FOV', [fov fov sliceThickness]*1e3);
seq.setDefinition('Name', 'spiral');

seq.write('spiral.seq');   % Output sequence for scanner

% % write the k-space trajectory (now 3D and wothout reconSize)
% save('epi_rs_traj.mat','traj_mat','ktime');%,'reconSize');
% 
seq.plot();             % Plot sequence waveforms
% figure; plot(traj_mat'); % plot k-space trajectory
% figure; plot(traj_mat(1,:),traj_mat(2,:),'b'); % a better plot
% hold;plot(traj_mat(1,:),traj_mat(2,:),'r.');
% 

% new single-function call for trajectory calculation
[ktraj_adc, ktraj, t_excitation, t_refocusing] = seq.calculateKspace();

% plot k-spaces

figure; plot(ktraj'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.');


% % seq.install('siemens');
