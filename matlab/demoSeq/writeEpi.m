% this is a demo low-performance EPI sequence;
% it doesn't use ramp-samping and is only good for educational purposes.
%
seq=mr.Sequence();              % Create a new sequence object
fov=220e-3; Nx=64; Ny=64;       % Define FOV and resolution
thickness=3e-3;                 % slice thinckness
Nslices=3;

% Set system limits
lims = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
               'MaxSlew',130,'SlewUnit','T/m/s', ...
               'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6);


% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',lims,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,...
    'use', 'excitation');

% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
dwellTime = 4e-6; % I want it to be divisible by 2
readoutTime = Nx*dwellTime;
flatTime=ceil(readoutTime*1e5)*1e-5; % round-up to the gradient raster
gx = mr.makeTrapezoid('x',lims,'Amplitude',kWidth/readoutTime,'FlatTime',flatTime);
adc = mr.makeAdc(Nx,'Duration',readoutTime,'Delay',gx.riseTime+flatTime/2-(readoutTime-dwellTime)/2);

% Pre-phasing gradients
preTime=8e-4;
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2,'Duration',preTime); % removed -deltak/2 to aligh the echo between the samples
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak,'Duration',preTime);

% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',dur);

% Define sequence blocks
% seq.addBlock(mr.makeDelay(1)); % older scanners like Trio may need this
                                 % dummy delay to keep up with timing
for s=1:Nslices
    rf.freqOffset=gz.amplitude*thickness*(s-1-(Nslices-1)/2);
    seq.addBlock(rf,gz);
    seq.addBlock(gxPre,gyPre,gzReph);
    for i=1:Ny
        seq.addBlock(gx,adc);           % Read one line of k-space
        seq.addBlock(gy);               % Phase blip
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
    end
    if s==1
        TR_1slice=seq.duration; % note the actual TR per slice
    end
end

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Plot sequence waveforms
seq.plot();

seq.plot('stacked',1,'timeRange',[0 TR_1slice], 'timeDisp','ms'); % niceer plot for the 1st sclice

%% trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.');

seq.write('epi.seq');   % Output sequence for scanner
% seq.sound(); % simulate the seq's tone
