system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=4096;
Nrep=16;
Ndummy=4

RfDuration=0.3e-3;
AdcDuration=256e-3;

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration', RfDuration, 'system', system,'use','excitation');

% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration', AdcDuration, 'system', system,'delay',system.adcDeadTime);
delayTE=20e-3;
delayTR=10000e-3;
%
assert(delayTE>=mr.calcDuration(rf));
assert(delayTR>=mr.calcDuration(adc));
% Loop over repetitions and define sequence blocks
for i=(1-Ndummy):Nrep
    seq.addBlock(rf,delayTE);
    if (i>0)
        seq.addBlock(adc,delayTR);
    else
        seq.addBlock(delayTR);
    end
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
