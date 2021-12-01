gamma = 42.57E6;

sequencerRasterTime = 1/(122.88E6); % make sure all times are a multiple of sequencer raster time

grad_interval = ceil(10E-6/sequencerRasterTime)*sequencerRasterTime;

rf_interval = ceil(1E-6/sequencerRasterTime)*sequencerRasterTime;


fov=10e-3; Nx=200; Ny=80;   % Define FOV and resolution

Ndummy = 2;

TR=5; % [s]     

ETL=8;

ESP=10e-3;

oversampling_factor = 4;

sliceThickness = 10;

use_slice = 0;

sp_amplitude = 2000; % spoiler area in 1/m (=Hz/m*s)

sp_duration = ceil(0.5E-3/sequencerRasterTime)*sequencerRasterTime;

requested_gxFlatTime = 3e-3;  % = adc read time [s]

requested_rf90duration = 0.1e-3;

% put the dwell time and rf90duration on a 2*sequencerRasterTime raster, so that

% mr.calcDuration(gx)/2 and rf90duration/2 are still on the sequencer raster 

dwellTime = ceil(requested_gxFlatTime/(Nx*oversampling_factor)/(2*sequencerRasterTime))*(2*sequencerRasterTime);

gxFlatTime = dwellTime * Nx * oversampling_factor;

rf90duration=ceil(requested_rf90duration/(2*sequencerRasterTime))*(2*sequencerRasterTime);

rf180duration=2*rf90duration;

% set system limits

maxGrad = 400; % [mT/m], value for tabletop coils and gpa fhdo

rfDeadTime = ceil(500e-6/sequencerRasterTime)*sequencerRasterTime; % [us], minicircuits PA needs 500 us to turn on

adcDeadTime = 0;

sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', 800, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, 'adcDeadTime', adcDeadTime, ...
    'rfRasterTime', rf_interval, 'gradRasterTime', grad_interval);

seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses, 500 us delay for tx gate

rf90duration=0.1e-3;

if use_slice == 1

    [rf90, gs] = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
        'PhaseOffset', 0, 'sys', sys, 'SliceThickness', sliceThickness);

    gs.channel='z'; % change it to X because we want sagittal orientation

else

    rf90 = mr.makeBlockPulse(pi/2, 'duration', rf180duration,...
        'PhaseOffset', 0, 'sys', sys);    

end

rf180 = mr.makeBlockPulse(pi, 'duration', rf90duration*2,...
    'PhaseOffset', pi/2, 'use','refocusing', 'sys',sys);

% Define other gradients and ADC events

deltak=1/fov;

kWidth=deltak*Nx;

kHeight=deltak*Ny;

gx = mr.makeTrapezoid('x','FlatArea',kWidth,'FlatTime',gxFlatTime,'sys',sys);

fprintf('Sequence bandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov);

fprintf('Pixelbandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov/Nx);

gx.delay = 0; % assumes rfDeadTime > gx.riseTime !!

gxPre = mr.makeTrapezoid('x','Area',gx.area/2,'Duration',gx.flatTime/2,'sys',sys);

g_sp = mr.makeTrapezoid('x','Area',sp_amplitude,'Duration',sp_duration,'system',sys);

gy_area = kHeight/2;

gy = mr.makeTrapezoid('y','Area',gy_area,'Duration',gx.flatTime/2,'sys',sys);

Ny = round(Ny);

adc = mr.makeAdc(round(oversampling_factor*Nx),'Duration',gx.flatTime,'Delay',gx.riseTime,'sys',sys);

% Calculate timing

delayTE3 = 0.5E-3;

delayTE = round((ESP/2 - (mr.calcDuration(rf90) - rf90.delay)/2 ...
    - mr.calcDuration(gxPre) - mr.calcDuration(g_sp) ...
    - rf180.delay - (mr.calcDuration(rf180) - rf180.delay)/2)/sequencerRasterTime)*sequencerRasterTime;

delayTE1 = round((ESP/2 -  mr.calcDuration(gx)/2 - gx.flatTime/2 ...
    - mr.calcDuration(g_sp) - rf180.delay - (mr.calcDuration(rf180) - rf180.delay)/2)/sequencerRasterTime)*sequencerRasterTime;

delayTE2 = round((ESP/2 - (mr.calcDuration(rf180) - rf180.delay)/2 ...
    - mr.calcDuration(gx)/2  -  mr.calcDuration(g_sp) -  mr.calcDuration(gy))/sequencerRasterTime)*sequencerRasterTime;

delayTR = TR - ETL*ESP -rf90.delay -(mr.calcDuration(rf90) - rf90.delay)/2 - mr.calcDuration(gx)/2;

fprintf('delay1: %.3f ms \ndelay2: %.3f ms \n',delayTE1*1E3,delayTE2*1E3)

phase_factor = linspace(-1,1,Ny);

n = 1;

while n <= Ny

    if n>Ny

        break

    end    

    if use_slice == 1

        seq.addBlock(rf90, gs);

    else

        seq.addBlock(rf90);

    end

    seq.addBlock(mr.makeDelay(delayTE));

    seq.addBlock(gxPre);

    for m=1:ETL  

        if n>Ny

            break

        end

        gy = mr.makeTrapezoid('y','Area',gy_area*phase_factor(n),'Duration',gx.flatTime/2,'sys',sys);    

        gy_rev = mr.makeTrapezoid('y','Area',-gy_area*phase_factor(n),'Duration',gx.flatTime/2,'sys',sys);    

        if m ~= 1

            seq.addBlock(mr.makeDelay(delayTE1));

        end

        seq.addBlock(g_sp);    

        seq.addBlock(rf180);

        seq.addBlock(g_sp);    

        seq.addBlock(mr.makeDelay(delayTE2));

        seq.addBlock(gy);        

        seq.addBlock(gx,adc);

        seq.addBlock(gy_rev);

        n = n+1;

    end

    seq.addBlock(mr.makeDelay(delayTR));

end

%% prepare sequence export

seq.setDefinition('Name', 'tse_2d');

seq.setDefinition('FOV', [fov fov]);

seq.setDefinition('ESP [s]', ESP);

seq.setDefinition('ELT', ETL);

seq.setDefinition('TR', TR);

seq.setDefinition('Nx', Nx);

seq.setDefinition('Ny', Ny);

seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);

seq.setDefinition('grad_t', grad_interval*1E6);

seq.setDefinition('tx_t', rf_interval*1E6);

seq.setDefinition('SliceThickness', sliceThickness);

seq.plot();

seq.write('tabletop_tse_pulseq.seq')       % Write to pulseq file

parsemr('tabletop_tse_pulseq.seq');

%% check whether the timing of the sequence is compatible with the scanner
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% do some visualizations

seq.plot();             % Plot sequence waveforms

% trajectory calculation
[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

% plot k-spaces
time_axis=(1:(size(ktraj,2)))*lims.gradRasterTime;
figure; plot(time_axis, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
%axis off;

