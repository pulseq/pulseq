system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq = mr.Sequence(system);              % Create a new sequence object
Nx = 1024;

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2, 'Duration', 0.1e-3, 'system', system);

% Define delays and ADC events
adc = mr.makeAdc(Nx, 'Duration', 0.32, 'system', system,'delay', system.adcDeadTime);
delayTE = 5e-3;

seq.addBlock(rf);
seq.addBlock(mr.makeDelay(delayTE));
seq.addBlock(adc);

seq.setDefinition('Name', 'fid');
seq.write('fiddisp.seq');       % Write to pulseq file
