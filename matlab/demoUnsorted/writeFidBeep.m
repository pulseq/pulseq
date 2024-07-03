system = mr.opts('MaxGrad',15,'GradUnit','mT/m',...
                 'MaxSlew',100,'SlewUnit','T/m/s',...
                 'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=4096;
Nrep=1;

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',0.3e-3, 'system', system);

% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',512e-3, 'system', system,'delay',system.adcDeadTime);
delayTE=20e-3;
delayTR=1000e-3;
%
assert(delayTE>=mr.calcDuration(rf));
assert(delayTR>=mr.calcDuration(adc));
% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf,mr.makeDelay(delayTE));
    seq.addBlock(adc,mr.makeDelay(delayTR));
end

% we use the mr.Music library for the beep-beep-beep part
mrMusic.init;
melody = {
    % bar 1
    [c1/8 o/8 o/4      o/4      o/4],...
    [o/4      c1/8 o/8 o/4      o/4],...
    [o/4      o/4      c1/8 o/8 o/4] };

% add beeps
[pitches, durations] = mrMusic.melodyToPitchesAndDurations(melody,'timeSignature',4/4);
seq = mrMusic.musicToSequence(seq, pitches, durations, 'barDurationSeconds', 2, 'pulseqUseWave', false, 'addDummyRfAdc', false);

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
seq.setDefinition('Name', 'fid-beep');
seq.write('fid-beep.seq')       % Write to pulseq file
