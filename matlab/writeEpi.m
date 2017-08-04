seq=mr.Sequence();              % Create a new sequence object
fov=220e-3; Nx=64; Ny=64;       % Define FOV and resolution

% Set system limits
lims = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s');  

% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',lims,'Duration',3e-3,...
    'SliceThickness',3e-3,'apodization',0.5,'timeBwProduct',4);

% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
readoutTime = 3.2e-4;
gx = mr.makeTrapezoid('x',lims,'FlatArea',kWidth,'FlatTime',readoutTime);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);

% Pre-phasing gradients
preTime=8e-4;
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2-deltak/2,'Duration',preTime);
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak,'Duration',preTime);

% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur);

% Define sequence blocks
seq.addBlock(rf,gz);
seq.addBlock(gxPre,gyPre,gzReph);
for i=1:Ny
    seq.addBlock(gx,adc);           % Read one line of k-space
    seq.addBlock(gy);               % Phase blip
    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
end

seq.write('epi.seq');   % Output sequence for scanner
seq.plot();             % Plot sequence waveforms