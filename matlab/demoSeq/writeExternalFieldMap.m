seq=mr.Sequence();              % Create a new sequence object
fov=220e-3; Nx=256; Ny=256;     % Define FOV and resolution

blipAmp = 0.2;


% Create 20 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(20*pi/180,'Duration',4e-3,...
    'SliceThickness',5e-3,'apodization',0.5,'timeBwProduct',4);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',6.4e-3);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',2e-3);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',2e-3);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% Calculate timing
delayTE=10e-3 - mr.calcDuration(gxPre) - mr.calcDuration(rf)/2 ...
    - mr.calcDuration(gx)/2;
delayTR=40e-3 - mr.calcDuration(gxPre) - mr.calcDuration(rf) ...
    - mr.calcDuration(gx) - delayTE;

% Loop over phase encodes and define sequence blocks
for b=1:1
    % Define gradient blip on external gradient channel
    g1 = mr.makeTrapezoid('1','Duration',2e-3,'Amplitude',blipAmp*b);
    
    for i=1:Ny
        seq.addBlock(rf,gz);
        gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',2e-3);
        seq.addBlock(gxPre,gyPre,gzReph,g1);
        seq.addBlock(mr.makeDelay(delayTE));
        seq.addBlock(gx,adc);
        seq.addBlock(mr.makeDelay(delayTR))
    end
end
seq.write('fieldmap.seq')       % Write to pulseq file