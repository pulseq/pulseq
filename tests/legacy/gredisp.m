% Define FOV, resolution and other sequence parameters
fov = 256e-3; slThck = 5e-3;
Nx = 32; Ny = Nx;
TE = 8e-3; TR = 22e-3; alpha=15;

% set system limits
sys = mr.opts('MaxGrad',20,'GradUnit','mT/m',...
    'MaxSlew',100,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6, 'setAsDefault', true);

% Create a new sequence object
seq=mr.Sequence(sys);

% Create slice selective alpha-pulse and corresponding gradients
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180, 'Duration', 1e-3,...
    'SliceThickness', slThck, 'apodization', 0.5,'timeBwProduct', 4,'use','excitation');
rf=mr.makeBlockPulse(alpha*pi/180, 'Duration', 1e-3,'use','excitation'); % replace the RF pulse by a block pulse to make the sequence listing shorter

% Define other gradients and ADC events
gx = mr.makeTrapezoid('x', 'FlatArea', Nx/fov, 'FlatTime', 6.4e-3);
adc = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
gxPre = mr.makeTrapezoid('x', 'Area', -gx.area/2, 'Duration', 2e-3);
gyPre = mr.makeTrapezoid('y', 'Area', Ny/2/fov, 'Duration', 2e-3);

% Calculate encoding tables and timing
phaseAreas = (0:Ny-1)/(Ny/2)-1;
delayTE = round((TE - mr.calcDuration(gxPre) - gz.flatTime/2 - gz.fallTime ...
                    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = round((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
                    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    seq.addBlock(rf, gz);    
    seq.addBlock(gxPre, mr.scaleGrad(gyPre,phaseAreas(i)), gzReph);
    seq.addBlock(delayTE);
    seq.addBlock(gx, adc);
    seq.addBlock(delayTR);
end

% export definitions
seq.setDefinition('FOV', [fov fov slThck]);
seq.setDefinition('Name', 'gre');

seq.write('gredisp.seq')       % Write to pulseq file
