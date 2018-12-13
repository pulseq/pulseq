% step
% 0 ... Basic sequence
% 1 ... Add spoiler in read, phase and slice (vary spoiler?)
% 2 ... Refocus in phase
% 3 ... Vary RF phase quasi-randomly
% 4 ... Make receiver phase follow transmitter phase
step = 0;

% Create a new sequence object
seq=mr.Sequence();

% Define FOV and resolution
fov = 256e-3;
sliceThickness = 5e-3;
Nx = 128;
Ny = Nx;

% Define sequence parameters
TE = 10e-3;
TR = 20e-3;
alpha=30;

% set system limits
sys = mr.opts('MaxGrad',25,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6);

% Create 20 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180, 'Duration', 4e-3,...
    'SliceThickness', sliceThickness, 'apodization', 0.5,'timeBwProduct', 4, ...
    'system' ,sys);

% Define other gradients and ADC events
deltak = 1/fov;
gx = mr.makeTrapezoid('x', 'FlatArea', Nx*deltak, 'FlatTime', 6.4e-3);
adc = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
gxPre = mr.makeTrapezoid('x', 'Area', -gx.area/2, 'Duration', 2e-3);
gzReph = mr.makeTrapezoid('z', 'Area', -gz.area/2, 'Duration', 2e-3);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;


% Calculate timing
delayTE = round((TE - mr.calcDuration(gxPre) - mr.calcDuration(gz)/2 ...
                    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = round((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
                    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;

if step > 0
    spoilArea=4*gx.area();
    % Add spoilers in read, refocus in phase and spoiler in slice
    gxPost = mr.makeTrapezoid('x', 'Area', spoilArea, 'system', sys); % we pass 'system' here to calculate shortest time gradient
    gyPost = mr.makeTrapezoid('y', 'Area', spoilArea, 'system', sys);
    gzPost = mr.makeTrapezoid('z', 'Area', spoilArea, 'system', sys);
end

if step > 1
    gyPost = mr.makeTrapezoid('y', 'Area', -max(phaseAreas(:)), 'Duration', 2e-3);
end

if step > 0
    delayTR = delayTR - mr.calcDuration(gxPost, gyPost, gzPost);
end

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    if step > 2
        % Vary RF phase quasi-randomly
        rand_phase = mod(117*(i^2 + i + 2), 360)*pi/180;
        [rf, gz] = mr.makeSincPulse(20*pi/180, 'Duration', 4e-3,...
                                    'SliceThickness', 5e-3, ...
                                    'apodization', 0.5, ...
                                    'timeBwProduct', 4, ...
                                    'system', sys, ...
                                    'phaseOffset', rand_phase);
    end
    seq.addBlock(rf, gz);
    gyPre = mr.makeTrapezoid('y', 'Area', phaseAreas(i), 'Duration', 2e-3);
    seq.addBlock(gxPre, gyPre, gzReph);
    seq.addBlock(mr.makeDelay(delayTE));
    if step > 3
        % Make receiver phase follow transmitter phase
        adc = mr.makeAdc(Nx, 'Duration', gx.flatTime,...
                         'Delay', gx.riseTime,...
                         'phaseOffset', rand_phase);
    end
    seq.addBlock(gx, adc);
    if step > 1
        gyPost = mr.makeTrapezoid('y', 'Area', -gyPre.area, 'Duration', 2e-3);
    end
    if step > 0
        % Add spoilers in read and slice and may be in phase
        seq.addBlock(gxPost, gyPost, gzPost);
    end
    seq.addBlock(mr.makeDelay(delayTR));
end

seq.setDefinition('FOV', [fov fov sliceThickness]*1e3);
seq.setDefinition('Name', 'gre');

seq.write(['DEMO_gre' num2str(step) '.seq'])       % Write to pulseq file

seq.plot('timeRange', [0 2*TR])

% do not run the rest of the script automatically
return

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

%% listen to the sequence
seq.sound();
