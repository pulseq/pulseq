% this is an experimental spiral sequence

fov=256e-3; Nx=64; Ny=Nx;  % Define FOV and resolution
sliceThickness=3e-3;             % slice thinckness
Nslices=4;
adcOversampling=2; % by looking at the periphery of the spiral I would say it needs to be at least 2
phi=pi/4; % orientation of the readout e.g. for interleaving

% Set system limits
sys = mr.opts('MaxGrad',20,'GradUnit','mT/m',...
    'MaxSlew',120,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6, 'adcDeadTime', 10e-6, 'adcSamplesLimit', 8192);  
seq=mr.Sequence(sys);          % Create a new sequence object
%warning('OFF', 'mr:restoreShape'); % restore shape is not compatible with spirals and will throw a warning from each plot() or calcKspace() call

% Create fat-sat pulse 
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
B0=2.89; % 1.5 2.89 3.0
sat_ppm=-3.35;
sat_freq=sat_ppm*1e-6*B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,'dwell',10e-6,...
    'bandwidth',abs(sat_freq),'freqPPM',sat_ppm,'use','saturation');
rf_fs.phasePPM=-2*pi*rf_fs.freqPPM*rf_fs.center; % compensate for the frequency-offset induced phase    

gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create 90 degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(pi/2,'system',sys,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'use','excitation');

% calculate a raw single-shot Archimedian spiral trajectory
% QC: the single-shot Archimedian spiral is a clasical spiral trajectory, whose radius
% (K) is propotional to its Azimuth angle (phi): K(t) = phi(t)/(2*pi*FOV).
% 2025.01.03
% define k-space parameters
deltak = 1/fov ;
kRadius = round(Nx/2) ; % QC: the radius of the spiral, determined by the user-defiled resolutions. it is also the number of circles of the spiral.
kSamples = round(2*pi*kRadius)*adcOversampling ; % QC: the number of samples of the outest circle of the spiral. 2025.01.03
tos_calculation = 10 ; % time oversampling during the trajectory optimization, can be anything, with higher factors giving smoother trajectory but causing slower calculation
gradOversampling = true ; % oversampling of the gradient shape, can be either true or false. QC: should be set to "false" for seq.write_v141. 20250103

clear ka;
ka(kRadius*kSamples+1) = 1i ; % initialize as complex
% QC: the single-shot Archimedian spiral is not really efficient because it has
% a lot of redundant k-space samples. 2025.01.03
for c = 0:kRadius*kSamples*tos_calculation % QC: total number of k-space points, account for oversampling for calculation
    r = deltak*c/kSamples/tos_calculation ;
    a = mod(c,kSamples*tos_calculation)*2*pi/kSamples/tos_calculation ;
    ka(c+1) = r*exp(1i*a) ; % in the polar coordinate, k(t) = r(t) * exp(1i*phi(t))
end
ka = [real(ka); imag(ka)] ; % QC: kx and ky

% calculate gradients and slew rates
dt = sys.gradRasterTime/tos_calculation ; % QC: sampling dwell time with calculation oversampling
[ga, sa] = mr.traj2grad(ka,'RasterTime',dt,'firstGradStepHalfRaster',tos_calculation==1,'conservativeSlewEstimate',true);

% limit analysis
safety_margin = 0.99 ; % we need that, otherwise we just about violate the slew rate due to the rounding errors
dt_gabs = abs(ga(1,:) + 1i*ga(2,:))/(sys.maxGrad*safety_margin)*dt ; % dt in case of decreased g
dt_sabs = sqrt(abs(sa(1,:)+1i*sa(2,:))/(sys.maxSlew*safety_margin))*dt ; % dt in case of decreased slew

%figure;plot([dt_gabs; max(dt_gcomp); dt_sabs; max(dt_scomp)]');title('time stepping defined by gradient and slew-rate');
%%
dt_opt=max([dt_gabs;dt_sabs]) ; % select the larger dt in case of decreased g and s

% apply the lower limit not to lose the trajectory detail
dt_min = 4*sys.gradRasterTime/kSamples/tos_calculation ; % we want at least 4 points per revolution. QC:  at least 4 points per circle. 2025.01.03
dt_opt0 = dt_opt ;
dt_opt(dt_opt<dt_min) = dt_min ;

figure;plot([dt_opt0; dt_opt]');title('combined time stepping');

t_smooth = [0 cumsum(dt_opt,2)] ;

dt_grad = sys.gradRasterTime/(1+gradOversampling) ;
if gradOversampling
    safety_factor_1st_timestep = 0.7 ; % we have to reduce the first time step because otherwise we are likely to exceed the slew rate due to the extreme non-smoothness of the trajectory start
    t_end = t_smooth(end)-(safety_factor_1st_timestep)*dt_grad;
    t_grad = [0 (safety_factor_1st_timestep+(0:floor(t_end/dt_grad)))*dt_grad];
else
    t_end = t_smooth(end)-0.5*dt_grad ;
    t_grad = [0 (0.5+(0:floor(t_end/dt_grad)))*dt_grad] ;
end
kopt=interp1(t_smooth, ka', t_grad)';

% analyze what we've got
fprintf('duration orig %d us\n', round(1e6*dt*length(ka)));
fprintf('duration smooth %d us\n', round(1e6*dt_grad*length(kopt)));

[gos, sos]=mr.traj2grad(kopt,'RasterTime',dt_grad,'firstGradStepHalfRaster',~gradOversampling);

figure;plot([gos;abs(gos(1,:)+1i*gos(2,:))]');
hold on; yline(sys.maxGrad,'--'); title('gradient with the abs constraint');

figure;plot([sos;abs(sos(1,:)+1i*sos(2,:))]');
hold on; yline(sys.maxSlew,'--'); title('slew rate with the abs constraint')

% Define gradients and ADC events
spiral_grad_shape=gos;

% calculate ADC
% round-down dwell time to 10 ns
adcTime = dt_grad*size(spiral_grad_shape,2);
% actually it is trickier than that: the (Siemens) interpreter sequence 
% per default will try to split the trajectory into segments with the number of samples <8192
% and every of these segments will have to have duration aligned to the
% gradient raster time

% adcSamplesPerSegment=1000; % you may need to play with this number to fill the entire trajectory
% adcSamplesDesired=kRadius*kSamples;
% adcSegments=round(adcSamplesDesired/adcSamplesPerSegment);
% adcSamples=adcSegments*adcSamplesPerSegment;
% adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
% adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 100 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
% if mod(adcSegmentDuration, sys.gradRasterTime)>eps 
%     error('ADC segmentation model results in incorrect segment duration');
% end
% % update segment count
% adcSegments=floor(adcTime/adcSegmentDuration);

adcSamplesDesired=kRadius*kSamples; 
adcDwell=round(adcTime/adcSamplesDesired/sys.adcRasterTime)*sys.adcRasterTime; 
adcSamplesDesired=ceil(adcTime/adcDwell);
[adcSegments,adcSamplesPerSegment]=mr.calcAdcSeg(adcSamplesDesired,adcDwell,sys); 

adcSamples=adcSegments*adcSamplesPerSegment;
% we would like to sample the point k=1 with the furst ADC sample (i.e. at 
% t=adcDwell/2), so we advance the ADC and round the delay to the RF raster time
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',round((mr.calcDuration(gzReph)-adcDwell/2)/sys.rfRasterTime)*sys.rfRasterTime);%lims.adcDeadTime);

% extend spiral_grad_shape by repeating the last sample
% this is needed to accomodate for the ADC tuning delay
if ~gradOversampling
    spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end)];
else
    spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end) spiral_grad_shape(:,end)];
    if mod(length(spiral_grad_shape),2)==0
        spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end)]; % the vector length with oversampling must be odd
    end
end

% readout grad 
gx = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:),'Delay',mr.calcDuration(gzReph),'first',0,'last', spiral_grad_shape(1,end),'system',sys,'oversampling',gradOversampling);
gy = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:),'Delay',mr.calcDuration(gzReph),'first',0,'last', spiral_grad_shape(2,end),'system',sys,'oversampling',gradOversampling);

% spoilers
gz_spoil=mr.makeTrapezoid('z',sys,'Area',deltak*Nx*4);
gx_spoil=mr.makeExtendedTrapezoid('x','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(1,end),0],'system',sys); %todo: make a really good spoiler
gy_spoil=mr.makeExtendedTrapezoid('y','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(2,end),0],'system',sys); %todo: make a really good spoiler

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
    seq.addBlock(rf_fs,gz_fs); % fat-sat    
    rf.freqOffset=gz.amplitude*sliceThickness*(s-1-(Nslices-1)/2);
    seq.addBlock(rf,gz);
    seq.addBlock(mr.rotate('z',phi,gzReph,gx,gy,adc,'system',sys));
    seq.addBlock(mr.rotate('z',phi,gx_spoil,gy_spoil,gz_spoil,'system',sys));
end

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'spiral');
seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); % this is important for making the sequence run automatically on siemens scanners without further parameter tweaking

seq.write('spiral.seq');   % Output sequence for scanner

% the sequence is ready, so let's see what we got 
seq.plot();             % Plot sequence waveforms

%% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); title('2D k-space');

% seq.install('siemens');
