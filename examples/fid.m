system = mr.opts('rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq = mr.Sequence(system);              % Create a new sequence object
Nx = 256;
Nrep = 16;

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2, 'Duration', 0.1e-3, 'system', system);

% Define delays and ADC events
adc = mr.makeAdc(Nx, 'Duration', 3.2e-3, 'system', system,'delay', system.adcDeadTime);
delayTE = 20e-3;
delayTR = 1000e-3;

% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf);
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(adc,mr.makeDelay(mr.calcDuration(adc)));
    seq.addBlock(mr.makeDelay(delayTR))
end

seq.write('fid.seq');       % Write to pulseq file
parsemr('fid.seq');