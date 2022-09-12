% Create a new sequence object
seq=mr.Sequence();

% Define FOV and resolution
fov = 256e-3;
sliceThickness = 5e-3;
Nx = 128;
Ny = Nx;

% Define sequence parameters
TE = 8e-3;
TR = 16e-3;
alpha=30;

% set system limits
sys = mr.opts('MaxGrad',25,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6);

% Create slice selection alpha-pulse and corresponding gradients
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180, 'Duration', 4e-3,...
    'SliceThickness', sliceThickness, 'apodization', 0.5,'timeBwProduct', 4, ...
    'system' ,sys);

% Define other gradients and ADC events
deltak = 1/fov; % Pulseq toolbox defaults to k-space units of m^-1
gx = mr.makeTrapezoid('x', 'FlatArea', Nx*deltak, 'FlatTime', 6.4e-3);
adc = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
gxPre = mr.makeTrapezoid('x', 'Area', -gx.area/2, 'Duration', 2e-3);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% Calculate timing
delayTE = round((TE - mr.calcDuration(gxPre) - mr.calcDuration(gz)/2 ...
                    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = round((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
                    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
% Loop over phase encodes and define sequence blocks
for i=1:Ny
    seq.addBlock(rf, gz);
    gyPre = mr.makeTrapezoid('y', 'Area', phaseAreas(i), 'Duration', 2e-3);
    seq.addBlock(gxPre, gyPre, gzReph);
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(gx, adc);
    seq.addBlock(mr.makeDelay(delayTR));
end

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (~ok)
    fprintf('Timing check failed! Sequence probably will not run on the scanner.\n'); 
end

% export definitions
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'DEMO_gre0'); % if submitting a sequence please write your name to the Name field of the definition section

seq.write('DEMO_grep0.seq')       % Write to pulseq file

seq.plot('timeRange', [0 2*TR])

return % stop here as the following cells are slower

%% plot gradients to check for gaps and optimality of the timing
gw=seq.gradient_waveforms();
figure; plot(gw'); % plot the entire gradient shape

%% new single-function call for trajectory calculation
[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

% plot k-spaces
time_axis=(1:(size(ktraj,2)))*sys.gradRasterTime;
figure; plot(time_axis, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport;
fprintf([rep{:}]);

%% listen to the sequence
seq.sound();
