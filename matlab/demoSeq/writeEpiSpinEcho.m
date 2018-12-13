seq=mr.Sequence();              % Create a new sequence object
fov=256e-3; Nx=64; Ny=64;       % Define FOV and resolution

% Set system limits
lims = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m',...
    'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  

% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,lims,'Duration',3e-3,...
    'SliceThickness',3e-3,'apodization',0.5,'timeBwProduct',4);

% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
readoutTime = 3.2e-4;
gx = mr.makeTrapezoid('x',lims,'FlatArea',kWidth,'FlatTime',readoutTime);
adc = mr.makeAdc(Nx,lims,'Duration',gx.flatTime,'Delay',gx.riseTime);

% Pre-phasing gradients
preTime=8e-4;
%gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2-deltak/2,'Duration',preTime);
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);
%gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak,'Duration',preTime);
% we need no minus for in-plane prephasers because of the spin-echo (position reflection in k-space)
gxPre = mr.makeTrapezoid('x',lims,'Area',gx.area/2-deltak/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',lims,'Area',Ny/2*deltak,'Duration',preTime);

% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur);

% Refocusing pulse with spoiling gradients
rf180 = mr.makeBlockPulse(pi,lims,'Duration',500e-6,'use','refocusing');
gzSpoil = mr.makeTrapezoid('z',lims,'Area',gz.area*2,'Duration',3*preTime);

% Calculate delay time
TE=60e-3;
durationToCenter = (Nx/2+0.5)*mr.calcDuration(gx) + Ny/2*mr.calcDuration(gy);
rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);
rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
delayTE1=TE/2 - mr.calcDuration(gz) + rfCenterInclDelay - preTime - mr.calcDuration(gzSpoil) - rf180centerInclDelay;
delayTE2=TE/2 - mr.calcDuration(rf180) + rf180centerInclDelay - mr.calcDuration(gzSpoil) - durationToCenter;

% Define sequence blocks
seq.addBlock(rf,gz);
seq.addBlock(gxPre,gyPre,gzReph); 
seq.addBlock(mr.makeDelay(delayTE1));
seq.addBlock(gzSpoil);
seq.addBlock(rf180);
seq.addBlock(gzSpoil);
seq.addBlock(mr.makeDelay(delayTE2));
for i=1:Ny
    seq.addBlock(gx,adc);           % Read one line of k-space
    seq.addBlock(gy);               % Phase blip
    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
end
seq.addBlock(mr.makeDelay(1e-4));

seq.write('epi_se.seq');   % Output sequence for scanner
seq.plot();             % Plot sequence waveforms

%% calculate trajectory 
[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

%% plot k-spaces
time_axis=(1:(size(ktraj,2)))*lims.gradRasterTime;
figure; plot(time_axis, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis

figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display

%% sanity checks
TE_check=(t_refocusing(1)-t_excitation(1))*2;
fprintf('intended TE=%.03f ms, actual spin echo TE=%.03fms\n', TE*1e3, TE_check*1e3); 
