% we define here a really crude True-FISP a.k.a. bSSFP sequence
% there is no user control for TR/TE, you just specify the ADC time and the
% RF parameters and the rest is calculated to find the fastest posible
% timing. The sequence intensively uses the extended trapezoid
% functionality to achieve near-optimal timing. Due to the requirement for
% splitting the sequence into blocks the TR is increased by approximately
% 40-60 us (rfDeadTime+adcDeadTime) in comparison to the true minimum TR

seq=mr.Sequence();              % Create a new sequence object
fov=220e-3; Nx=256; Ny=256;     % Define FOV and resolution

% set system limits
% had to slow down ramps and increase adc_duration to avoid stimulation
sys = mr.opts('MaxGrad',36,'GradUnit','mT/m',...
    'MaxSlew',140,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 20e-6);

% ADC duration (controls TR/TE)
adc_dur=2560; %us

% RF parameters 
alpha=40; % deg
thick=4; %mm
rf_dur=500; % us
rf_apo=0.5;
rf_bwt=1;

% Create 'alpha' degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
    'SliceThickness',thick*1e-3,'apodization',rf_apo,'timeBwProduct',rf_bwt,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% now we have to reschuffle gradients to achieve a half-way optimal timing
% new gz will consist of two parts: 
% 1: slice refocusing from the previous TR followed by slice selection 
%    including the plato an a small bit of the ramp-down
% 2: the remainer of the ramp-down and the slice refocusing for the next TR
tp=[0];
am=[0];
tp=[tp gzReph.riseTime+tp(end)];
am=[am gzReph.amplitude];
if gzReph.flatTime>0
    tp=[tp gzReph.flatTime+tp(end)];
    am=[am gzReph.amplitude];
end
tp=[tp gzReph.fallTime+tp(end)];
am=[am 0];
tp=[tp gz.riseTime+tp(end)];
am=[am gz.amplitude];
tp=[tp gz.flatTime+tp(end)];
am=[am gz.amplitude];
t_ramp_part=rf.t(end)-(gz.delay+gz.riseTime+gz.flatTime);
a_ramp_part=gz.amplitude*(gz.fallTime-t_ramp_part)/gz.fallTime;
tp=[tp t_ramp_part+tp(end)];
am=[am a_ramp_part];
gz_1=mr.makeExtendedTrapezoid('z',sys,'times',tp,'amplitudes',am);
rf.delay=tp(end)-rf.t(end);

tp2=[0];
am2=[a_ramp_part];
tp2=[tp2 gz.fallTime-t_ramp_part];
am2=[am2 0];
tp2=[tp2 gzReph.riseTime+tp2(end)];
am2=[am2 gzReph.amplitude];
if gzReph.flatTime>0
    tp2=[tp2 gzReph.flatTime+tp2(end)];
    am2=[am2 gzReph.amplitude];
end
tp2=[tp2 gzReph.fallTime+tp2(end)];
am2=[am2 0];
gz_2=mr.makeExtendedTrapezoid('z',sys,'times',tp2,'amplitudes',am2);

% new gr will consist of two parts: 
% 1: prephaser followed by a part of the read gradient including the 
%    beginning of the ramp-down
% 2: the remainer of the ramp-down and the second "prephaser"
tpr1=[0];
amr1=[0];
tpr1=[tpr1 gxPre.riseTime+tpr1(end)];
amr1=[amr1 gxPre.amplitude];
if gxPre.flatTime>0
    tpr1=[tpr1 gxPre.flatTime+tpr1(end)];
    amr1=[amr1 gxPre.amplitude];
end
tpr1=[tpr1 gxPre.fallTime+tpr1(end)];
amr1=[amr1 0];
tpr1=[tpr1 gx.riseTime+tpr1(end)];
amr1=[amr1 gx.amplitude];
tpr1=[tpr1 gx.flatTime+tpr1(end)];
amr1=[amr1 gx.amplitude];
t_ramp_part=ceil((adc.dwell*adc.numSamples+adc.delay+adc.deadTime)/sys.gradRasterTime)*sys.gradRasterTime - (gx.delay+gx.riseTime+gx.flatTime);
a_ramp_part=gx.amplitude*(gx.fallTime-t_ramp_part)/gx.fallTime;
tpr1=[tpr1 t_ramp_part+tpr1(end)];
amr1=[amr1 a_ramp_part];
gx_1=mr.makeExtendedTrapezoid('x',sys,'times',tpr1,'amplitudes',amr1);
adc.delay=adc.delay+tpr1(end-3);

tpr2=[0];
amr2=[a_ramp_part];
tpr2=[tpr2 gx.fallTime-t_ramp_part+tpr2(end)];
amr2=[amr2 0];
tpr2=[tpr2 gxPre.riseTime+tpr2(end)];
amr2=[amr2 gxPre.amplitude];
if gxPre.flatTime>0
    tpr2=[tpr2 gxPre.flatTime+tpr2(end)];
    amr2=[amr2 gxPre.amplitude];
end
tpr2=[tpr2 gxPre.fallTime+tpr2(end)];
amr2=[amr2 0];
gx_2=mr.makeExtendedTrapezoid('x',sys,'times',tpr2,'amplitudes',amr2);

% adjust delays to align objects
if tpr2(end)>tp(end-2)
    delay_add=tpr2(end)-tp(end-2);
    gz_1.delay=delay_add;
    rf.delay=rf.delay+delay_add;
end

% Calculate timing
pe_dur=min(max(tp2(end),tpr1(end-2)),max(tp(end-3),tpr2(end))); 
TR=mr.calcDuration(gz_1)+mr.calcDuration(gx_1);
TE=TR/2;

% alpha / 2 preparation: shorter TR and half-angle, no PE, no RO
% create 0.5*alpha prep pulse
rf05=rf;
rf05.signal=0.5*rf.signal;
seq.addBlock(rf05,gz_1);
seq.addBlock(gz_2);
% the following delay calculation fails for agressive sequence timing
seq.addBlock(mr.makeDelay(TR/2-mr.calcDuration(gz_1)-mr.calcDuration(gz_2)));

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    gyPre_2 = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',pe_dur); % current PE step
    if i>1
        gyPre_1 = mr.makeTrapezoid('y','Area',-phaseAreas(mod(i+Ny-2,Ny)+1),'Duration',pe_dur); % previous PE step
        seq.addBlock(rf,gz_1, gyPre_1, gx_2);
    else
        seq.addBlock(rf,gz_1);
    end
    %seq.addBlock(gxPre,gyPre,gzReph);
    %seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(gx_1,gyPre_2, gz_2,adc);
    %seq.addBlock(mr.makeDelay(delayTR))
end
% finish the x-grad shape 
seq.addBlock(gx_2);

% check that the calculated TR was reached
% alpha / 2 prep takes 3 blocks
assert(TR==(mr.calcDuration(seq.getBlock(4))+mr.calcDuration(seq.getBlock(5))));

fprintf('Sequence ready\n');
fprintf('TR=%f ms  TE=%f ms\n', TR*1e3, TE*1e3);

seq.write('trufi.seq')       % Write to pulseq file
seq.plot();