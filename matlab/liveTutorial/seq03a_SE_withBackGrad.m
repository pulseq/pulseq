system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=4096;
Nrep=1;
adcDur=51.2e-3; 
rfDur=1000e-6;
TR=160e-3;
TE=80e-3;
bg=5000; % Hz/m
%todo: vary gradient

% Create non-selective excitation and refocusing pulses 
rf_ex = mr.makeBlockPulse(pi/2,'Duration',rfDur, 'system', system);
rf_ref = mr.makeBlockPulse(pi,'Duration',rfDur, 'system', system, 'use', 'refocusing'); % needed for the proper k-space calculation
    
% Define delays and ADC events
delayTE1=TE/2-mr.calcDuration(rf_ex)/2-mr.calcDuration(rf_ref)/2;
delayTE2=TE/2-mr.calcDuration(rf_ref)+rf_ref.delay+mr.calcRfCenter(rf_ref)-adcDur/2; % this is not perfect, but -adcDur/2/Nx  will break the raster alignment

adc = mr.makeAdc(Nx,'Duration',adcDur, 'system', system, 'delay', delayTE2);

delayTR=TR-mr.calcDuration(rf_ex)-delayTE1-mr.calcDuration(rf_ref);

assert(delayTE1>=0);
assert(delayTE2>=0);
assert(delayTR>=0);

% ramp up the background gradient
ramptime=abs(bg)/system.maxSlew;
ramptime=ceil(ramptime/system.gradRasterTime)*system.gradRasterTime; % round-up to gradient raster 
ramptime=max(ramptime, 3*system.gradRasterTime);                     % and limit it to 3x system.gradRasterTime
seq.addBlock(mr.makeExtendedTrapezoid('x','amplitudes',[0 bg],'times',[0 ramptime]));

% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf_ex,mr.makeExtendedTrapezoid('x','amplitudes',[bg bg],'times',[0 mr.calcDuration(rf_ex)]));
    seq.addBlock(mr.makeDelay(delayTE1),mr.makeExtendedTrapezoid('x','amplitudes',[bg bg],'times',[0 delayTE1])); 
    seq.addBlock(rf_ref,mr.makeExtendedTrapezoid('x','amplitudes',[bg bg],'times',[0 mr.calcDuration(rf_ref)]));
    seq.addBlock(adc,mr.makeDelay(delayTR),mr.makeExtendedTrapezoid('x','amplitudes',[bg bg],'times',[0 delayTR]));   
end

% ramp down the background gradient
seq.addBlock(mr.makeExtendedTrapezoid('x','amplitudes',[bg 0],'times',[0 ramptime]));

seq.plot();

% check whether the timing of the sequence is compatible with the scanner
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.write('se.seq')       % Write to pulseq file
%seq.install('siemens');    % copy to scanner

% calculate and plot k-spaces
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis

% calculate real TE and TR, etc
% reported TR will be slighlty higher for Nrep=1 becuse it uses TA instead
rep = seq.testReport; 
fprintf([rep{:}]); 
