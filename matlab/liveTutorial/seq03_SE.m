system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=4096;
Nrep=1;
adcDur=51.2e-3; 
rfDur=1000e-6;
TR=250e-3;
TE=60e-3;
% todo: change ADC duration

% Create non-selective excitation and refocusing pulses
rf_ex = mr.makeBlockPulse(pi/2,'Duration',rfDur, 'system', system); % for this phantom and this coil I had to reduce flip angle to avoid receiver saturation
rf_ref = mr.makeBlockPulse(pi,'Duration',rfDur, 'system', system, 'use', 'refocusing'); % needed for the proper k-space calculation
    
% Define delays and ADC events
delayTE1=TE/2-mr.calcDuration(rf_ex)/2-mr.calcDuration(rf_ref)/2;
delayTE2=TE/2-mr.calcDuration(rf_ref)+rf_ref.delay+mr.calcRfCenter(rf_ref)-adcDur/2; % this is not perfect, but -adcDur/2/Nx  will break the raster alignment

adc = mr.makeAdc(Nx,'Duration',adcDur, 'system', system, 'delay', delayTE2);

delayTR=TR-mr.calcDuration(rf_ex)-delayTE1-mr.calcDuration(rf_ref);

assert(delayTE1>=0);
assert(delayTE2>=0);
assert(delayTR>=0);

% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf_ex);
    seq.addBlock(mr.makeDelay(delayTE1)); 
    seq.addBlock(rf_ref);
    seq.addBlock(adc,mr.makeDelay(delayTR));  
end

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

%% calculate k-space but only use it to check the TE calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

assert(abs(t_refocusing-t_excitation-TE/2)<1e-6); % check that the refocusing happens at the 1/2 of TE
assert(abs(t_adc(Nx/2)-t_excitation-TE)<adc.dwell); % check that the echo happens as close as possible to the middle of the ADC elent
