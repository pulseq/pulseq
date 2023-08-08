% this is a very naiive and non-optimized cardiac cine GRE sequence 

seq=mr.Sequence();              % Create a new sequence object
fov=256e-3; Nx=128; Ny=Nx;      % Define FOV and resolution
alpha=5;                        % flip angle
sliceThickness=5e-3;            % slice
%TE=[7.38 9.84]*1e-3;            % give a vector here to have multiple TEs (e.g. for field mapping)
TE=4.92e-3;
TR=9e-3;                        % only a single value for now

% cardiac features
phases = 8;
hearbeats = 15; % odd numbers of heartbeats / segments work better

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
rf_duration = 2e-3;
adc_duration = 3.2e-3;
pre_duration = 1e-3;

% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% Create fat-sat pulse 
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
% B0=2.89; % 1.5 2.89 3.0
% sat_ppm=-3.45;
% sat_freq=sat_ppm*1e-6*B0*lims.gamma;
% rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims,'Duration',8e-3,...
%     'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
% gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% define the trigger to play out
trig=mr.makeTrigger('physio1','duration', 2000e-6); % duration after
%trig=mr.makeTriggerPulse('osc0','duration', 4100e-6); % possible channels: 'osc0','osc1','ext1'
trig_out=mr.makeDigitalOutputPulse('ext1','duration', 100e-6,'delay',500e-6); % possible channels: 'osc0','osc1','ext1'


% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',rf_duration,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',adc_duration,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',pre_duration,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',pre_duration,'system',sys);

lines_per_segment = round(Ny/hearbeats);
Ns=ceil(Ny/lines_per_segment);
Ny=Ns*lines_per_segment; % it can be that because of the rounding above we measure few more k-space lines...
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
% now reverse the order in every second segment
phaseAreasSeg=reshape(phaseAreas,lines_per_segment,Ns);
phaseAreasSeg(:,2:2:end)=phaseAreasSeg(end:-1:1,2:2:end);
phaseAreas=phaseAreasSeg(:);

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

fprintf('the sequence will acquire %d lines per segment resulting in a temporal resolution of %g ms per phase\n', lines_per_segment, TR*lines_per_segment*1e3);
fprintf('cardiac acquisition window is: %g ms\n', TR*phases*lines_per_segment*1e3);

rf_phase=0;
rf_inc=0;

% Loop over phase encodes and define sequence blocks
for s=1:Ns
    seq.addBlock(trig); % wait for the cardiac trigger
    for p=1:phases
        for l=1:lines_per_segment
            % restore the line counter
            i=(s-1)*lines_per_segment+l;
            %seq.addBlock(rf_fs,gz_fs); % fat-sat
            rf.phaseOffset=rf_phase/180*pi;
            adc.phaseOffset=rf_phase/180*pi;
            rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase=mod(rf_phase+rf_inc, 360.0);
            %
            seq.addBlock(rf,gz,trig_out);
            gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',pre_duration,'system',sys);
            seq.addBlock(gxPre,gyPre,gzReph);
            if delayTE>0 
                seq.addBlock(mr.makeDelay(delayTE));
            end
            seq.addBlock(gx,adc);
            gyPre.amplitude=-gyPre.amplitude;
            seq.addBlock(mr.makeDelay(delayTR),gxSpoil,gyPre,gzSpoil)
        end
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

%% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'cine-gre');

seq.write('cine_gre.seq')       % Write to pulseq file

%seq.install('siemens');
return

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5*TR]);

% new single-function call for trajectory calculation
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

