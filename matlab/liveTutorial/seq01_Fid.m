system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=4096;
Nrep=1;
adcDur=51.2e-3; 
rfDur=500e-6;
TR=5000e-3;
TE=10e-3; % the used coil has a very long switching time!
% todo1: change flip angle to reduce SNR
% todo2: change repetitions to increase SNR

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',rfDur, 'system', system);
    
% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',adcDur, 'system', system, 'delay', TE-rfDur/2-system.rfRingdownTime);

delayTR=TR-mr.calcDuration(rf);
assert(delayTR>=0);

% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf);
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

seq.write('fid.seq')       % Write to pulseq file
%seq.install('siemens');    % copy to scanner

% optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport; 
fprintf([rep{:}]); 
