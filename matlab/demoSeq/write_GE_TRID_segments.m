% -------------------------------------------------------------------------
% NOTE:
% This example sequence is NOT physically meaningful and is NOT intended
% for actual MRI experiments. Its sole purpose is to demonstrate the
% automatic handling of GE-specific TRID segment labels in Pulseq.
%
% The key idea is that the user only needs to call
%       seq.addTRID('segment_name')
% before the first block of each logical sequence segment.
%
% The numeric TRID IDs required by GE scanners are assigned automatically
% by the Sequence object based on the first occurrence of each label name.
% This removes the need for the sequence programmer to manually plan,
% track, or maintain a global TRID ID table, which is especially helpful
% for large, modular, or loop-based sequence constructions.
%
% When flag_ge = false, TRID labels are silently ignored, allowing the same
% source code to be used across different vendor platforms.
% -------------------------------------------------------------------------

% version: Maximilian Gram; University of Wuerzburg; 20.02.2026 

%% init system and seq object
clear

% set system limits
sys = mr.opts('MaxGrad', 40,  'GradUnit', 'mT/m', ...
              'MaxSlew', 180, 'SlewUnit', 'T/m/s', ... 
              'rfDeadTime',          1e-4, ...
              'rfRingdownTime',      1e-4, ...
              'adcDeadTime',         4e-5, ...
              'adcRasterTime',       2e-6, ...
              'rfRasterTime',        2e-6, ...
              'gradRasterTime',      4e-6, ...
              'blockDurationRaster', 4e-6, ...
              'flag_ge', 1 ); % IMPORTANT: 0/1 -> defines if GE specific TRID labels are added

% init sequence object
seq = mr.Sequence(sys);

%% calculate some sequence objects

% ----- Inversion Recovery Segment -----
rfInv   = mr.makeSincPulse(180*pi/180, sys, 'Duration', 10e-3, 'timeBwProduct', 8, 'use', 'inversion'); % Global inversion-like RF pulse
gzCrush = mr.makeTrapezoid('z', sys, 'Area', 2000, 'Duration', 3e-3); % Crusher gradient (arbitrary area)
tRec    = [5 25 100 250] *1e-3; % Recovery delays

% ----- Fat Suppression Segment -----
rfFat      = mr.makeGaussPulse(110*pi/180, sys, 'Duration', 8e-3, 'bandwidth', 3.45*1e-6*sys.B0*sys.gamma, 'freqOffset', -3.45*1e-6*sys.B0*sys.gamma, 'use', 'saturation'); % simple fat selective pulse
gzFatCrush = mr.makeTrapezoid('z', sys, 'Area', 1000, 'Duration', 2e-3); 

% Slice-selective excitation RF (sinc) + its slice gradient is inside rfExc.gz
[rfExc, gz] = mr.makeSincPulse(10*pi/180, sys, 'Duration', 2e-3, 'SliceThickness', 5e-3, 'apodization', 0.5, 'timeBwProduct', 4, 'use', 'excitation');

% Simple readout gradient and ADC (arbitrary timing)
gx  = mr.makeTrapezoid('x', sys, 'FlatTime', 2e-3, 'FlatArea', 10);
adc = mr.makeAdc(128, sys, 'Duration', gx.flatTime, 'Delay', gx.riseTime);

%% add sequence objects and TRID labels to seq object

for j = 1 : numel(tRec)

    seq.addTRID(['inv_prep_' num2str(j)]); % TRID label with auto-incrementing
    seq.addBlock(rfInv);
    seq.addBlock(gzCrush);
    seq.addBlock(mr.makeDelay(tRec(j)));

    seq.addTRID('fat_suppression'); % TRID label with fixed increment
    seq.addBlock(rfFat);
    seq.addBlock(gzFatCrush);

    for k = 1:8
        seq.addTRID('readout'); % TRID label with fixed increment
        seq.addBlock(rfExc, gz);
        seq.addBlock(gx, adc);
    end

end

%% print TRID label list and history, plot sequence, write seq file

disp('TRID history (order of calls):');
disp(seq.tridHistory);

disp('TRID id -> name (first occurrence defines numeric ID):');
disp(seq.tridId2Name);

seq.plot('stacked', true, 'showGuides', false);

seq.write('demo_trid_simple.seq');
