%% set system limits and other options

sys = mr.opts('MaxGrad', 24, 'GradUnit', 'mT/m', ...
              'MaxSlew', 140, 'SlewUnit', 'T/m/s', ...
              'blockDurationRaster', 10e-6, 'gradRasterTime', 10e-6, ...
              'rfRingdownTime', 20e-6, 'rfDeadTime', 1100e-6,'rfRasterTime', 10e-6, ...
              'adcRasterTime', 1e-10, 'adcDeadTime', 350e-6);
% 'gradRasterTime', 6.4e-6

seq=mr.Sequence(sys);           % Create a new sequence object
fov=256e-3; Nx=256; Ny=256;     % Define FOV and resolution
alpha=15;                       % flip angle
sliceThickness=3e-3;            % slice
TR=12e-3;                       % repetition time TR
TE=5e-3;                        % echo time TE  
%TE=[7.38 9.84]*1e-3;           % alternatively give a vector here to have multiple TEs (e.g. for field mapping)

% pTx related options
num_tx = 8;                     % Number of Tx channels
slice_delay = 2;                % Delay in seconds between slices
make_power_monitoring_scan = false;     % Make a fancy pattern, as used in the GIF of ISMRM abstract & presentation

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
roDuration=3.2e-3;              % ADC duration

%% Create Pulseq objects and calculate sequence timing

% Create fat-sat pulse 
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
% B0=2.89; % 1.5 2.89 3.0
% sat_ppm=-3.45;
% sat_freq=sat_ppm*1e-6*B0*lims.gamma;
% rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims,'Duration',8e-3,...
%     'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
% gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys);
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys, 'Delay', (sys.adcDeadTime-gx.riseTime));

adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',sys.adcDeadTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',1e-3,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',1e-3,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
gyPre = mr.makeTrapezoid('y','Area',max(abs(phaseAreas)),'Duration',mr.calcDuration(gxPre),'system',sys);
peScales=phaseAreas/gyPre.area;
        

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTE>=0));
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

%%

% Eye matrix that uses (only) one of the Tx channels (in order)
eye_tx = eye(num_tx);

if make_power_monitoring_scan
    % Make a more interesting (power usage) pattern
    tx_pattern = eye_tx;
    tx_pattern((end+1):(end+7), :) = eye_tx(7:-1:1, :); % Go back again
    
    tx_pattern((end):(end+3), :) = eye_tx(1:4, :) + eye_tx(8:-1:5, :); % Collapse in the centre!
else
    % Otherwise only do a slice with single Tx channels
    tx_pattern = eye_tx;
end

% Prepare the RF pulses
ptx = {}; id = {}; shapeids = {};
for i = 1:size(tx_pattern, 1)

    ptx{i} = rf;

    ptx{i}.signal = reshape(rf.signal' .* tx_pattern(i, :), 1, []);
    ptx{i}.t = repmat(rf.t, 1, num_tx);

    %[id{i} shapeids{i} ] = seq.registerRfEvent(ptx{i});
end



%%

% Loop over the different transmit patterns / pulses
for cur_tx = 1:size(tx_pattern, 1)
    rf_phase=0;
    rf_inc=0;
    % Loop over phase encodes and define sequence blocks
    for cur_y=1:Ny
        for c=1:length(TE)
            %seq.addBlock(rf_fs,gz_fs); % fat-sat
            ptx{cur_tx}.phaseOffset=rf_phase/180*pi;
            adc.phaseOffset=rf_phase/180*pi;
            rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase=mod(rf_phase+rf_inc, 360.0);
            %
            seq.addBlock(ptx{cur_tx},gz);
            seq.addBlock(gxPre,mr.scaleGrad(gyPre,peScales(cur_y)),gzReph);
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gx,adc);
            %gyPre.amplitude=-gyPre.amplitude;
            seq.addBlock(mr.makeDelay(delayTR(c)),gxSpoil,mr.scaleGrad(gyPre,-peScales(cur_y)),gzSpoil)
        end
    end
    seq.addBlock(mr.makeDelay(slice_delay))
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

%% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'gre2d_pTxSingleChan_8Tx_forDario');

seq.write('gre2d_pTxSingleChan_8Tx_forDario.seq')       % Write to pulseq file

%seq.install('siemens');

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5]*TR);

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport;
fprintf([rep{:}]);

