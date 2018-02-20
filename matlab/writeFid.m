system = mr.opts('rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=256;
Nrep=1;

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',0.1e-3, 'system', system);

% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',3.2e-3, 'system', system);
delayTE=20e-3;
delayTR=1000e-3;

% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf);
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(adc);
    seq.addBlock(mr.makeDelay(delayTR))
end

seq.write('fid.seq')       % Write to pulseq file