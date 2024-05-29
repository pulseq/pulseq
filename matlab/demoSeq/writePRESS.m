system = mr.opts('MaxGrad', 15, 'GradUnit', 'mT/m', ...
                 'MaxSlew', 100, 'SlewUnit', 'T/m/s', ...
                 'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
voxel=[20 30 40]*1e-3; % voxel size
Nx=4096;
Nrep=1;
Ndummy=0;
adcDur=256e-3; 
rfDurEx=3000e-6;
rfDurRef=4000e-6;
TR=400e-3;
TE=120e-3;
spA=0.6e3; % spoiler area in 1/m (=Hz/m*s)
spB=2.0e3; % spoiler area in 1/m (=Hz/m*s)

%% Create slice-selective excitation and refocusing pulses
[rf_ex, g_ex, g_exReph] = mr.makeSincPulse(pi/2,'Duration',rfDurEx,...
    'SliceThickness',voxel(1),'apodization',0.5,'timeBwProduct',8,'system',system);
[rf_ref1, g_ref1] = mr.makeSincPulse(pi,'Duration',rfDurEx,'PhaseOffset',pi/2,...
    'SliceThickness',voxel(2),'apodization',0.6,'timeBwProduct',8,'system',system,'use','refocusing');
[rf_ref2, g_ref2] = mr.makeSincPulse(pi,'Duration',rfDurEx,'PhaseOffset',pi/2,...
    'SliceThickness',voxel(3),'apodization',0.6,'timeBwProduct',8,'system',system,'use','refocusing');
% fix channels for the gradients
g_ex.channel='x';
g_ref1.channel='y';

%% join spoilers with the slice selection pulses of the refocusing gradients
% step 1: create pre-gradient to merge into the plato
g_ref1_pre=mr.makeExtendedTrapezoidArea(g_ref1.channel,0,g_ref1.amplitude,spA,system); 
% step 2: create post-gradient to start at the plato
g_ref1_post=mr.makeExtendedTrapezoidArea(g_ref1.channel,g_ref1.amplitude,0,spA,system);
% step 3: create a composite gradient 
g_refC1=mr.makeExtendedTrapezoid(g_ref1_pre.channel,...
    'times', [g_ref1_pre.tt g_ref1_post.tt+g_ref1_pre.shape_dur+g_ref1.flatTime],...
    'amplitudes',[g_ref1_pre.waveform g_ref1_post.waveform],'system',system);
% see what we've got
figure; plot(g_refC1.tt*1e3,g_refC1.waveform/system.gamma*1e3); title('combined gradient for the first refocunsing pulse');ylabel('mT/m');xlabel('ms');
% same procedure for the second refocusing pulse slice selection
g_ref2_pre =mr.makeExtendedTrapezoidArea(g_ref2.channel,0,g_ref2.amplitude,spB,system);
g_ref2_post=mr.makeExtendedTrapezoidArea(g_ref2.channel,g_ref2.amplitude,0,spB,system);
g_refC2=mr.makeExtendedTrapezoid(g_ref2_pre.channel,...
    'times', [g_ref2_pre.tt g_ref2_post.tt+g_ref2_pre.shape_dur+g_ref2.flatTime],...
    'amplitudes',[g_ref2_pre.waveform g_ref2_post.waveform],'system',system);

%% update RF pulses delays to center them on the central flat parts of the combined gradients
rf_ref1.delay=g_ref1_pre.shape_dur;
rf_ref2.delay=g_ref2_pre.shape_dur;

% now calculate other spoiler gradients
g_spAz1=mr.makeTrapezoid('z','Area',spA,'system',system);
g_spAz2=mr.makeTrapezoid('z','Area',spA,'system',system,'delay',mr.calcDuration(g_spAz1)+g_ref1.flatTime);
g_spAx1=mr.makeTrapezoid('x','Area',spA+g_exReph.area,'system',system);% notice we reduce the area to account for slice refocusing
g_spAx2=mr.makeTrapezoid('x','Area',spA,'system',system,'delay',mr.calcDuration(g_spAz1)+g_ref1.flatTime);
g_spBy1=mr.makeTrapezoid('y','Area',spB,'system',system);
g_spBy2=mr.makeTrapezoid('y','Area',spB,'system',system,'delay',mr.calcDuration(g_spBy1)+g_ref2.flatTime);
g_spBx1=mr.makeTrapezoid('x','Area',spB,'system',system);
g_spBx2=mr.makeTrapezoid('x','Area',spB,'system',system,'delay',mr.calcDuration(g_spBy1)+g_ref2.flatTime);
% combine spoilers to composite gradients
g_spAz=mr.addGradients({g_spAz1,g_spAz2},'system', system);
g_spAx=mr.addGradients({g_spAx1,g_spAx2},'system', system);
g_spBy=mr.addGradients({g_spBy1,g_spBy2},'system', system);
g_spBx=mr.addGradients({g_spBx1,g_spBx2},'system', system);
% update delays in g_refC1, g_refC2, rf_ref1 and rf_ref2 in case g_spAz1 is longer than g_ref1_pre
g_refC1.delay=g_refC1.delay-mr.calcDuration(g_ref1_pre)+mr.calcDuration(g_spAz1);
g_refC2.delay=g_refC2.delay-mr.calcDuration(g_ref2_pre)+mr.calcDuration(g_spBy1);
rf_ref1.delay=rf_ref1.delay-mr.calcDuration(g_ref1_pre)+mr.calcDuration(g_spAz1);
rf_ref2.delay=rf_ref2.delay-mr.calcDuration(g_ref2_pre)+mr.calcDuration(g_spBy1);

%% Define delays and ADC events
delayTE1=1e-3; % this delay allows to shift the spin echo within the ADC window
% we define TE as 2* delay between the centers of the refocusing pulses
% delayTE2=TE/2-(mr.calcDuration(rf_ref1)-mr.calcRfCenter(rf_ref1)-rf_ref1.delay)-g_ref1_post.shape_dur-rf_ref2.delay-mr.calcRfCenter(rf_ref2);
% QC:
delayTE2=TE/2 - rf_ref1.shape_dur/2 - g_ref1_post.shape_dur - rf_ref2.delay - rf_ref2.shape_dur/2 ;
assert(delayTE2>=0);
% we start the ADC object right away after the spoiler
adc = mr.makeAdc(Nx,'Duration',adcDur, 'system', system);

% delayTR=TR-mr.calcDuration(g_ex)-mr.calcDuration(g_refC1)-delayTE1-delayTE2-mr.calcDuration(g_refC2)-mr.calcDuration(adc);
delayTR=TR-max(mr.calcDuration(g_ex), mr.calcDuration(rf_ex))-mr.calcDuration(g_refC1)-delayTE1-delayTE2-mr.calcDuration(g_refC2)-mr.calcDuration(adc);
assert(delayTR>=0);

%% Loop over repetitions and define sequence blocks
for i=(1-Ndummy):Nrep
    seq.addBlock(rf_ex,g_ex);
    seq.addBlock(mr.makeDelay(delayTE1));
    seq.addBlock(rf_ref1,g_refC1,g_spAz,g_spAx);
    seq.addBlock(mr.makeDelay(delayTE2)); 
    seq.addBlock(rf_ref2,g_refC2,g_spBy,g_spBx);
    if i>0 
        seq.addBlock(adc);
    else
        seq.addBlock(mr.makeDelay(mr.calcDuration(adc)));
    end
    seq.addBlock(mr.makeDelay(delayTR));
end

seq.plot('showBlocks',true,'timeDisp','us');

% check whether the timing of the sequence is compatible with the scanner
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.setDefinition('FOV', voxel);
seq.setDefinition('Name', 'press');
seq.write('press.seq')       % Write to pulseq file

%% calculate k-space but only use it to check timing
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('gradient_offset',[0 0 1000]);

figure; plot(t_ktraj,ktraj);title('k-space components as functions of time'); grid on;
