seq=mr.Sequence();              % Create a new sequence object
fov=[190e-3 190e-3 190e-3];     % Define FOV
Nx=64; Ny=Nx; Nz=Nx;            % Define FOV and resolution
Tread=3.2e-3;
Tpre=3e-3;
riseTime=400e-6;
Ndummy=50;
sys=mr.opts('maxGrad',20,'gradUnit','mT/m','riseTime',riseTime,...
    'RFdeadTime',10e-6,'ADCdeadTime',10e-6);


% Create non-selective pulse
rf = mr.makeBlockPulse(8*pi/180,sys,'Duration',0.2e-3);

% Define other gradients and ADC events
deltak=1./fov;
gx = mr.makeTrapezoid('x',sys,'FlatArea',Nx*deltak(1),'FlatTime',Tread);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'Duration',Tpre);
gxSpoil = mr.makeTrapezoid('x',sys,'Area',gx.area,'Duration',Tpre);
areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);

% Calculate timing
TE=10e-3;
TR=40e-3;
delayTE = TE - mr.calcDuration(rf)/2 - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gx)/2;
delayTR = TR - mr.calcDuration(rf) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - mr.calcDuration(gxSpoil) - delayTE;

% Drive magnetization to steady state
for iY=1:Ndummy
    % RF
    rf.phaseOffset = mod(117*(iY^2+iY+2)*pi/180,2*pi);
    seq.addBlock(rf);
    % Gradients
    gyPre = mr.makeTrapezoid('y','Area',areaY(floor(Ny/2)),'Duration',Tpre);
    gyReph = mr.makeTrapezoid('y','Area',-areaY(floor(Ny/2)),'Duration',Tpre);
    seq.addBlock(gxPre,gyPre);
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(gx);
    seq.addBlock(gyReph,gxSpoil);
    seq.addBlock(mr.makeDelay(delayTR))
end

% Make trapezoids for inner loop to save computation
for iY=1:Ny
    gyPre(iY) = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',Tpre);
    gyReph(iY) = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',Tpre);
end

% Loop over phase encodes and define sequence blocks
for iZ=1:Nz
    gzPre = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',Tpre);
    gzReph = mr.makeTrapezoid('z','Area',-areaZ(iZ),'Duration',Tpre);
    for iY=1:Ny
        % RF spoiling
        rf.phaseOffset = mod(117*(iY^2+iY+2)*pi/180,2*pi);
        adc.phaseOffset = rf.phaseOffset;
        
        % Excitation
        seq.addBlock(rf);
        
        % Encoding
        seq.addBlock(gxPre,gyPre(iY),gzPre);
        seq.addBlock(mr.makeDelay(delayTE));
        seq.addBlock(gx,adc);
        seq.addBlock(gyReph(iY),gzReph,gxSpoil);
        seq.addBlock(mr.makeDelay(delayTR))
    end
end

% Visualise sequence and output for execution
seq.plot('TimeRange',[Ndummy+1 Ndummy+3]*TR)
seq.write('external.seq');
