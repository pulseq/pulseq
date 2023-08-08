% set system limits
system = mr.opts('MaxGrad', 10, 'GradUnit', 'mT/m', ...
    'MaxSlew', 50, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% high-level sequence parameters 
alpha=90;
sliceThickness=50e-3;
Nx = 2048;
Nrep = 1;
TE=30e-3;
TR = 6000e-3;

rf_90_duration=2.6e-3; % excitation pulse duration
rf_180_duration= 4.6e-3; % refocusing pulse duration

% create and populate the sequence
seq = mr.Sequence(system);              % Create a new sequence object

% Create excitation pulse and gradient
rf_90= mr.makeSLRpulse(pi/2,'duration',rf_90_duration,'timeBwProduct',7.88,'dwell',rf_90_duration/500,'passbandRipple',1,'stopbandRipple',1e-2,'filterType','ms','system',system); 

%Create refocusing pulse and gradient
rf_180_1 = mr.makeAdiabaticPulse('wurst','duration',rf_180_duration,'bandwidth',6000,'dwell',rf_180_duration/500,'n_fac',20,'use','refocusing','system',system); 

timebwproduct_90=8;   % "experimental" value for the used for SLR excitation
timebwproduct_180=22; % "experimental" value for the used refocusing pulse

grad_amplitude_90 = (timebwproduct_90/2.6e-3)/sliceThickness; %BW/thickness
grad_amplitude_180 = (timebwproduct_180/4.6e-3)/sliceThickness; %BW/thickness

% we cannot use mr.calcDuration(rf_90) and mr.calcDuration(rf_180_1) here because this already includes delays and RF dead times
gx = mr.makeTrapezoid('x','flatTime',rf_90_duration,'amplitude',grad_amplitude_90,'system',system);
gy = mr.makeTrapezoid('y','flatTime',rf_180_duration,'amplitude',grad_amplitude_180,'system',system);
gz = mr.makeTrapezoid('z','flatTime',rf_180_duration,'amplitude',grad_amplitude_180,'system',system);
gz_2 = mr.makeTrapezoid('z','flatTime',rf_180_duration,'amplitude',grad_amplitude_180,'system',system);

% refocusing area for the 90-degree pulse (needs to be remembered in the "refocusing budget")
area_90_toRefocus=gx.amplitude*(0.5*gx.fallTime+gx.flatTime-mr.calcRfCenter(rf_90));

% gradient spoiling
spoilMoment=10/sliceThickness;
gzSpoil=mr.makeTrapezoid('z','area',spoilMoment,'system',system);
gxSpoil=mr.makeTrapezoid('x','duration',2*mr.calcDuration(gzSpoil),'area',spoilMoment,'system',system);
gzSpoil_2=mr.makeTrapezoid('z','duration',2*mr.calcDuration(gzSpoil),'area',spoilMoment,'system',system);

% split x gradient at the end of the plato and calculate the rf_90 delay
gx_parts=mr.splitGradientAt(gx,mr.calcDuration(gx)-gx.fallTime+system.rfRingdownTime);
[gx_p1,rf_90,~]=mr.align('right',gx_parts(1),rf_90,mr.makeDelay(system.rfDeadTime+rf_90_duration+system.rfRingdownTime));
gx_parts(2).delay=0; % set the delay to 0 because it will be used in a separate block

% split y gradient at the end of the plato and calculate the rf_180_1 delay
gy_parts=mr.splitGradientAt(gy,mr.calcDuration(gy)-gy.fallTime+system.rfRingdownTime);
rf_180_1.delay=max(mr.calcDuration(gy_parts(1))-rf_180_duration-system.rfRingdownTime, mr.calcDuration(gzSpoil)); 
assert(rf_180_1.delay>=system.rfDeadTime);
gy_p1=gy_parts(1);
gy_p1.delay=rf_180_1.delay+rf_180_duration+system.rfRingdownTime-mr.calcDuration(gy_p1);
gy_parts(2).delay=0; % set the delay to 0 because it will be used in a separate block

% combine y gradient with delays to a new gradient the rf_180_2 delay
rf_180_2=rf_180_1;
rf_180_2.delay=max(mr.calcDuration(gxSpoil),mr.calcDuration(gzSpoil_2));
gy_tmp=mr.splitGradientAt(gy,mr.calcDuration(gy)-gy.fallTime);
gy_tmp(1).delay=rf_180_2.delay+rf_180_duration-mr.calcDuration(gy_tmp(1));
gySpoil=mr.makeExtendedTrapezoidArea('y',gy.amplitude,0,spoilMoment+0.5*gy.amplitude*gy.fallTime,system);
gySpoil.delay=mr.calcDuration(gy_tmp(1));
gy_comb=mr.addGradients({gy_parts(2),gy_tmp(1),gySpoil},'system',system);
gy_comb_parts=mr.splitGradientAt(gy_comb,rf_180_2.delay+rf_180_duration+system.rfRingdownTime);%
gy_comb_parts(2).delay=0;

%
gxSpoil_2=mr.makeTrapezoid('x','area',spoilMoment,'system',system);
rf_180_3=rf_180_1;
rf_180_3.delay=mr.calcDuration(gxSpoil_2);
gz.delay=rf_180_3.delay-gz.riseTime;

%Additional Spoiler gradients 
gzSpoil_semiFinal=mr.makeExtendedTrapezoidArea('z',0,gz.amplitude,spoilMoment,system);
gxSpoil_semiFinal=mr.makeTrapezoid('x','area',spoilMoment+area_90_toRefocus,'system',system); 
gySpoil_semiFinal=mr.makeTrapezoid('y','area',2*spoilMoment,'system',system); 
[gxSpoil_semiFinal,gySpoil_semiFinal,gzSpoil_semiFinal]=mr.align('right',gxSpoil_semiFinal,gySpoil_semiFinal,gzSpoil_semiFinal);

rf_180_4=rf_180_1;
rf_180_4.delay=0;
gz_temp=mr.splitGradientAt(gz_2,gz_2.riseTime);
gz_temp(2).delay=0;
gz_parts=mr.splitGradientAt(gz_temp(2),rf_180_4.delay+rf_180_duration);
gz_parts(1).delay=mr.calcDuration(gzSpoil_semiFinal);
rf_180_4.delay=mr.calcDuration(gzSpoil_semiFinal);

gzSpoil_Final=mr.makeExtendedTrapezoidArea('z',gz.amplitude,0,spoilMoment,system);
gxSpoil_Final=mr.makeTrapezoid('x','area',spoilMoment,'system',system);
gySpoil_Final=mr.makeTrapezoid('y','area',spoilMoment,'system',system);

gzSpoil_Final.delay=rf_180_4.delay+rf_180_duration;
gz_comb=mr.addGradients({gzSpoil_semiFinal,gz_parts(1),gzSpoil_Final},'system',system);

gxSpoil_Final.delay=rf_180_4.delay+rf_180_duration;
gySpoil_Final.delay=rf_180_4.delay+rf_180_duration;

gxSpoil_combi=mr.addGradients({gxSpoil_semiFinal,gxSpoil_Final},'system',system);
gySpoil_combi=mr.addGradients({gySpoil_semiFinal,gySpoil_Final},'system',system);

%timing calculation
lTime1=(rf_90_duration+rf_180_duration)/2+system.rfRingdownTime+mr.calcDuration(gzSpoil);
lTime2=(rf_180_duration+rf_180_duration)/2+system.rfRingdownTime+mr.calcDuration(gxSpoil);
lTime3=(rf_180_duration+rf_180_duration)/2+system.rfRingdownTime+mr.calcDuration(gxSpoil_2);

lTime4=TE/2-lTime2;
lTime5=TE/2-lTime1-lTime3;

%Define ADC events
adc = mr.makeAdc(Nx, 'Dwell', 2e-4, 'system', system);
delayTE1=lTime4-(rf_180_duration+gz.fallTime+mr.calcDuration(gySpoil_semiFinal));
delayTE2=lTime5-(rf_180_duration/2+mr.calcDuration(gySpoil_Final)-gySpoil_Final.delay)-adc.dwell/2;
adc.delay=delayTE2;

% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf_90, gx_p1); 
    seq.addBlock(gzSpoil,gx_parts(2),rf_180_1,gy_p1);
    seq.addBlock(gxSpoil,gzSpoil_2,rf_180_2, gy_comb_parts(1));
    seq.addBlock(gxSpoil_2,gy_comb_parts(2),rf_180_3, gz);
    seq.addBlock(mr.makeDelay(delayTE1));
    seq.addBlock(gxSpoil_combi,gySpoil_combi,gz_comb,rf_180_4);
    seq.addBlock(adc,mr.makeDelay(mr.calcDuration(adc)+system.adcDeadTime));
    
    % this is realy a lazy way of defining the TR delay
    if i==1
        delayTR = TR- seq.duration();
        assert(delayTR>0);
    end
    
    seq.addBlock(mr.makeDelay(delayTR))
end

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
seq.setDefinition('FOV', [sliceThickness sliceThickness sliceThickness]);
seq.setDefinition('Name', 'semiLaser');

seq.write('semiLASER.seq');       % Write to pulseq file

%% do some visualizations

seq.plot('timeDisp','us','showBlocks',1,'timeRange',[0 TE*1.2]);
%seq.plot();             % Plot sequence waveforms

%% trajectory calculation
% [ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();
% t_ktraj=(1:(size(ktraj,2)))*system.gradRasterTime;
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
[ktraj_adc_ofs, t_adc_ofs, ktraj_ofs, t_ktraj_ofs] = seq.calculateKspacePP('gradient_offset',[1.0 .5 -1.0]*1e3); % this will help us to verify the correct echo positions

% plot "k-spaces" (gradient moments)
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
plot(t_ktraj_ofs, ktraj_ofs');
axis([0 TE*2 -250 250]); 
title('gradient moments without and with background gradients');

%% additional timing checks 
% Karl Landheer et al define tau(1:5) and require tau_1+tau_2+tau_3 = tau_2+tau4
tau=diff([t_excitation t_refocusing t_adc(1)]);

if abs(sum(tau)-TE)>5e-5 % we tolerate an error of 1/2 grad rasters
    warning('TE calculation seems to be wrong, check timing!');
end

if abs(sum(tau([1 3 5]))-sum(tau([2 4])))>5e-5 % we tolerate an error of 1/2 grad rasters
    warning('spin echo condition is not fulfilled, check timing!');
end

