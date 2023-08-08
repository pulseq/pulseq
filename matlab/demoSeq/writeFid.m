system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=256;
Nrep=16;

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',0.1e-3, 'system', system);

% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',3.2e-3, 'system', system,'delay',system.adcDeadTime);
delayTE=20e-3;
delayTR=5000e-3;
%
assert(delayTE>=mr.calcDuration(rf));
assert(delayTR>=mr.calcDuration(adc));
% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf,mr.makeDelay(delayTE));
    seq.addBlock(adc,mr.makeDelay(delayTR));
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
seq.setDefinition('Name', 'fid');
seq.write('fid.seq')       % Write to pulseq file
