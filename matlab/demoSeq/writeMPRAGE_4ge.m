% GE-related variables/calculations are in the 'GE' struct,
% used here as a namespace to avoid confusion.
% All other variables are Pulseq/Siemens specific, or general/common.

ex.mode = 'hard';  % 'we' = 1-1 binomial water excitation; 'hard' = rectangular pulse

% set system limits (slew rate 130 and max_grad 30 work on Prisma)
% Reduce slew to make the inner TR about the same as the GE scan
sys = mr.opts('MaxGrad', 24, 'GradUnit', 'mT/m', ...
    'MaxSlew', 80, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'B0', 123/128*3, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq = mr.Sequence(sys);           % Create a new sequence object
fov = [192 256 256]*1e-3;         % Define FOV and resolution
N = [192 256 256];                % matrix sizes
alpha = 8;                        % flip angle
%ro_dur=5017.6e-6; % BW=200Hz/pix
%load GE   % see bw.m. % TODO: bw.m is poorly named
GE = getGEparams;
ro_dur = GE.ro_dur;
%ro_os=2;                        % readout oversampling
ro_os = GE.decimation;
ro_spoil = 3;                    % additional k-max excursion for RO spoiling
TI = 1.07;                       % to match the ABCD protocol. Was: 1.1
TRout=2.5;

% TE & TR in the inner loop are as short as possible derived from the above parameters and the system specs.
% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
rfLen=100e-6;
ax=struct; % encoding axes
ax.d1='z'; % the fastest dimension (readout)
ax.d2='x'; % the second-fastest dimension (the inner pe loop)
ax.d3=setdiff('xyz',[ax.d1 ax.d2]); % automatically set the slowest dimension
ax.n1=strfind('xyz',ax.d1);
ax.n2=strfind('xyz',ax.d2);
ax.n3=strfind('xyz',ax.d3);

% Create alpha-degree hard or binominal water excitation pulse.
% Make duration multiple of sys.gradRasterTime so it passes seq.checkTiming().
% Pad duration to multiple of both 4us (GE raster) and 10us (sys.gradRasterTime)
% Design on GE raster time, then interpolate to sys.rfRasterTime.
fatChemShift = 3.5; % ppm
gamma = 4.2576e7;   % Hz/Tesla
gamG = 4.2576e3;     % Hz/Gauss
fatOffres = gamma*sys.B0*fatChemShift*1e-6;  % fat resonance frequency offset (Hz)
ex.raster = GE.raster; % s 
ex.nhard = 10; % number of waveform samples in each alpha/2 hard pulse
ex.hard = [(alpha/2/360) / (gamG * ex.nhard * ex.raster) * ones(ex.nhard,1)];
if strcmp(ex.mode, 'we')
    ex.ngap = round(1/fatOffres/2/ex.raster) - ex.nhard; 
else
    ex.ngap = 0;
end
tmpWav = [ex.hard; zeros(ex.ngap,1); ex.hard];  % Gauss; 4us raster
ex.dur = pulsegeq.roundtoraster(ex.raster*length(tmpWav), sys.gradRasterTime*2);
ex.n = round(ex.dur/ex.raster);
tmpWav = [tmpWav; zeros(ex.n-length(tmpWav),1)];
ex.signal = pulsegeq.rf2pulseq(tmpWav, ex.raster, sys);  % Gauss -> Hz; 4us -> 1us.
%rf = mr.makeBlockPulse(alpha*pi/180, sys, 'Duration',rfLen); 
rf = mr.makeArbitraryRf(ex.signal, alpha/180*pi, 'system', sys);

% Create adiabatic inversion pulse
rf180 = mr.makeAdiabaticPulse('hypsec', sys, 'Duration', 10.24e-3, 'dwell',1e-5);

% Define other gradients and ADC events
deltak=1./fov;
gro = mr.makeTrapezoid(ax.d1, ...
    'Amplitude', N(ax.n1)*deltak(ax.n1)/ro_dur, ...
    'FlatTime', ceil(ro_dur/sys.gradRasterTime)*sys.gradRasterTime, ...
    'system',sys);
adc = mr.makeAdc(N(ax.n1)*ro_os, ...
    'Duration', ro_dur,...
    'Delay', gro.riseTime, ...
    'system',sys);
groPre = mr.makeTrapezoid(ax.d1, ...
    'Area', -gro.amplitude*(adc.dwell*(adc.numSamples/2+0.5)+0.5*gro.riseTime),...
    'system',sys); % the first 0.5 is necessary to acount for the Siemens sampling in the center of the dwell periods
gpe1 = mr.makeTrapezoid(ax.d2, ...   
    'Area', deltak(ax.n2)*(N(ax.n2)/2),...  
    'system',sys);  % maximum PE1 gradient
gpe2 = mr.makeTrapezoid(ax.d3, ...
    'Area', deltak(ax.n3)*(N(ax.n3)/2), ...
    'system',sys);  % maximum PE2 gradient
gslSp = mr.makeTrapezoid(ax.d3, ...
    'Area', max(deltak.*N)*4, ...  % spoil with 4x cycles per voxel
    'Duration', 10e-3, ...
    'system',sys);  

% we cut the RO gradient into two parts for the optimal spoiler timing
[gro1,groSp]=mr.splitGradientAt(gro,gro.riseTime+gro.flatTime);
% gradient spoiling
if ro_spoil>0
    groSp = mr.makeExtendedTrapezoidArea(gro.channel, gro.amplitude, 0, deltak(ax.n1)/2*N(ax.n1)*ro_spoil, sys);
end

% calculate timing of the fast loop 
% we will have two blocks in the inner loop:
% 1: RF 
% 2: prewinder,phase enconding + readout + spoilers/rewinders
[groPre,~] = mr.align('right',groPre,mr.makeDelay(mr.calcDuration(gpe1,gpe2)-gro.riseTime));
gro1.delay = mr.calcDuration(groPre);
groSp.delay = mr.calcDuration(gro1);
adc.delay = gro1.delay+gro.riseTime;
gro1 = mr.addGradients({gro1,groPre,groSp},'system',sys);
gpe1c = mr.addGradients({gpe1, ...
    mr.makeTrapezoid(ax.d2, 'Area', -gpe1.area, 'duration', groSp.shape_dur, 'delay', groSp.delay, 'system',sys)});
gpe2c = mr.addGradients({gpe2, ...
    mr.makeTrapezoid(ax.d3, 'Area', -gpe2.area, 'duration', groSp.shape_dur, 'delay', groSp.delay, 'system',sys)});
TRinner = mr.calcDuration(rf)+mr.calcDuration(gro1); % we'll need it for the TI delay

% peSteps -- control reordering
pe1Steps = ((0:N(ax.n2)-1)-N(ax.n2)/2)/N(ax.n2)*2;
pe2Steps = ((0:N(ax.n3)-1)-N(ax.n3)/2)/N(ax.n3)*2;

% TI calc
TIdelay=round((TI-(find(pe1Steps==0)-1)*TRinner-(mr.calcDuration(rf180)-mr.calcRfCenter(rf180)-rf180.delay)-rf.delay-mr.calcRfCenter(rf))/sys.blockDurationRaster)*sys.blockDurationRaster;
TRoutDelay=TRout-TRinner*N(ax.n2)-TIdelay-mr.calcDuration(rf180);

% pre-register objects that do not change while looping
gslSp.id=seq.registerGradEvent(gslSp);
gro1.id=seq.registerGradEvent(gro1);
[~, gpe1c.shapeIDs]=seq.registerGradEvent(gpe1c);
[~, gpe2c.shapeIDs]=seq.registerGradEvent(gpe2c);
[~, rf.shapeIDs]=seq.registerRfEvent(rf); % the phase of the RF object will change, therefore we only per-register the shapes 
[rf180.id, rf180.shapeIDs]=seq.registerRfEvent(rf180); % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create .mod files and modules.txt needed for execution on GE scanners

% Write inversion pulse to inversion.mod
GE.rf180.t = GE.raster/2:GE.raster:max(rf180.t);
tmp = interp1(rf180.t, rf180.signal, GE.rf180.t, 'linear', 'extrap');   % interpolate to GE raster time (4us)
%tmp(isnan(tmp)) = 0;   % due to interp1
%GE.nChop = [0 0]; % trim this many samples at beginning/end of RF waveform
%tmp = [zeros(GE.nChop(1)+2,1); tmp(:); zeros(GE.nChop(2)+2,1)]; % must start + end with zero 
GE.rf180.signal = [0; tmp(:)/gamG; 0]           % Gauss
GE.rf180.signal = toppe.makeGElength(GE.rf180.signal);  % enforce 4-sample boundary
toppe.writemod(GE.sys, ...
    'ofname', 'inversion.mod', ...
    'rf', GE.rf180.signal);

% Write spoil.mod (only used after inversion pulse)
GE.spoil.mxs = 8;  % Gauss/cm/ms. Lower to reduce PNS.
res = fov(1)/N(1)*100;   % spatial resolution (cm)
GE.spoil.nSpoilCycles = 4;
gCrush = toppe.utils.makecrusher(GE.spoil.nSpoilCycles, res, GE.sys, 0, ...
    GE.spoil.mxs, GE.sys.maxGrad);
toppe.writemod(GE.sys, 'gy', gCrush, 'ofname', 'spoil.mod');

% Write binomial water excitation pulse to tipdown.mod
% Note that waveform must start and end with 0.
GE.fatOffres = gamma*128/128*3*fatChemShift*1e-6;
if strcmp(ex.mode, 'we')
    GE.ex.ngap = round(1/GE.fatOffres/2/ex.raster) - ex.nhard;
else
    GE.ex.ngap = 0;
end
%tmpWav = [zeros(GE.nChop(1)+2,1); ex.hard; zeros(GE.ex.ngap,1); ex.hard; zeros(GE.nChop(2)+2,1)]; % Gauss; 4us raster 
GE.ex.signal = toppe.makeGElength([0; ex.hard; ex.hard; 0]);  % enforce 4-sample boundary
toppe.writemod(GE.sys, ...
    'ofname', 'tipdown.mod', ...
    'rf', GE.ex.signal);

% Create readout.mod 
% This file contains the x/y/z readout gradient waveforms at full amplitude;
% the gradient amplitude during scanning are then set dynamically
% as specified in scanloop.txt.
% Here we use the helper function 'makegre' to design the waveforms,
% but that's not a requirement.
% Note: In TOPPE, all waveforms are 'arbitrary' (using a Pulseq term), 
% i.e., there is no special handling of 'trapezoids' as in Pulseq.
% Reduce slew during design to reduce PNS.
GE.zres = fov(ax.n3)/N(ax.n3) * 1e2;   % cm
tmpSys = toppe.systemspecs('maxSlew', min(15, GE.sys.maxSlew), 'slewUnit', 'Gauss/cm/ms', ... 
    'maxGrad', GE.sys.maxGrad, 'gradUnit', 'Gauss/cm');
[GE.ro.wav, GE.pey.wav, GE.pez.wav] = toppe.utils.makegre(fov(ax.n1)*1e2, ... % logical y, z
    N(ax.n1), GE.zres, tmpSys, ...
    'oprbw', GE.oprbw, ...
    'ncycles', ro_spoil/2 + 0.5, ...  % number of cycles of spoiling along readout
    'ofname', []);  % don't write gradients to file
toppe.writemod(GE.sys, ...
    'ofname', 'readout.mod', ...
    'gx', GE.pez.wav, 'gy', GE.pey.wav, 'gz', GE.ro.wav, ...
    'nChop', [64 64]);  % turn off ADC during at beginning and end, to make room for ADC dead time

% create modules.txt file for TOPPE
% If placing in a folder other than /usr/g/bin/,
% you must provide the full filename path.
GE.modFileText = ['' ...
'Total number of unique cores\n' ...
'4\n' ...
'fname  duration(us)    hasRF?  hasDAQ?\n' ...
'inversion.mod\t0\t1\t0\n' ...   % entries are tab-separated   
'spoil.mod\t0\t0\t0\n' ...
'tipdown.mod\t0\t1\t0\n' ...
'readout.mod\t0\t0\t1' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, GE.modFileText);
fclose(fid);

% Calculate delays for TOPPE to achieve desired TI and volume TR (TRout)
% toppe.getTRtime contructs the sequence from the TOPPE files in the
% current directory (*.mod, modules.txt, scanloop.txt)
toppe.write2loop('setup', GE.sys, 'version', 4);
toppe.write2loop('inversion.mod', GE.sys);
toppe.write2loop('spoil.mod', GE.sys);
toppe.write2loop('finish', GE.sys);
GE.inversion.mindur = toppe.getTRtime(1, 2, GE.sys);  % sec
toppe.write2loop('setup', GE.sys, 'version', 4);
toppe.write2loop('tipdown.mod', GE.sys);
toppe.write2loop('readout.mod', GE.sys);
toppe.write2loop('finish', GE.sys);
GE.inner.tr = toppe.getTRtime(1, 2, GE.sys);   % sec
GE.inner.dur = GE.inner.tr * N(ax.n2);
GE.inversion.delay = TI - GE.inversion.mindur - GE.inner.dur/2; % sec
GE.outer.delay = TRout - GE.inversion.mindur - GE.inversion.delay - GE.inner.dur;

%% done creating .mod files and modules.txt for TOPPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define GRAPPA undersampling pattern along slowest dimension (outer pe loop)
nACS = 32;  % number of auto-calibration lines (dense sampling in center)
R = 2;      % undersampling outside ACS region
S = 0*(1:N(ax.n3));
S(1:R:end) = 1;
S( (end/2-nACS/2):(end/2+nACS/2-1) ) = 1;
J = find(S == 1);

% intialize scanloop.txt file for TOPPE
toppe.write2loop('setup', GE.sys, 'version', 4);

% Scan loop
GE.view = 1;  
for j = J  % J(1:6)  % 1:N(ax.n3) 
    for ib = 1:40
        fprintf('\b');
    end
    fprintf('Writing scanloop: pe %d of %d', j, N(ax.n3));

    % inversion pulse, spoiler, and delay
    seq.addBlock(rf180);
    seq.addBlock(mr.makeDelay(TIdelay),gslSp);

    % for TOPPE
    toppe.write2loop('inversion.mod', GE.sys, ...
        'RFamplitude', 1.0);
    toppe.write2loop('spoil.mod', GE.sys, ...
        'textra', GE.inversion.delay*1e3); % ms

    rf_phase=0;
    rf_inc=0;

    % pre-register the PE gradients that repeat in the inner loop
    gpe2cj=mr.scaleGrad(gpe2c,pe2Steps(j));
    gpe2cj.id=seq.registerGradEvent(gpe2cj);

    for i=1:N(ax.n2)
        % excitation and readout (for .seq file)
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        seq.addBlock(rf);
        seq.addBlock(adc,gro1,mr.scaleGrad(gpe1c,pe1Steps(i)),gpe2cj);

        % excitation and readout (for GE)
        % GE data is stored in 'slice', 'echo', and 'view' indeces
        toppe.write2loop('tipdown.mod', GE.sys, ...
            'RFphase', rf_phase/180*pi); 
        GE.textra = (i == N(ax.n2)) * GE.outer.delay*1e3; % ms
        toppe.write2loop('readout.mod', GE.sys, ...
            'Gamplitude', [-pe1Steps(i) -pe2Steps(j) 1.0]', ... % pe1 = logical z = physical x; pe2 = logical y = physical y
            'DAQphase', rf_phase/180*pi, ...
            'textra', GE.textra, ...
            'slice', i, 'view', GE.view);  

        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
    end
    
    GE.view = GE.view + 1;

    seq.addBlock(mr.makeDelay(TRoutDelay));
end

% Add noise scans.
% Do this last, since receive gain for GE is based on signal from 
% beginning of sequence.
% First insert ~5s pause to allow magnetization/system to settle.
% Set delay of adc block so duration is multiple of sys.gradRasterTime.
adcDur = adc.numSamples * adc.dwell; % s
adcDurNew = pulsegeq.roundtoraster(adcDur, sys.gradRasterTime);
adc.delay = adc.delay + (adcDurNew - adcDur);
seq.addBlock(mr.makeDelay(5)); % sec
toppe.write2loop('spoil.mod', GE.sys, ...
    'Gamplitude', [0 0 0]', ...
    'textra', 5e3-2);  % ms
for i = 1:N(ax.n2)  
    seq.addBlock(adc);
end
for j = 1:2  % so max slice and view is even (for GE)
    for i = 1:N(ax.n2)  
        toppe.write2loop('readout.mod', GE.sys, ...
        'Gamplitude', [0 0 0]', ... % turn off gradients 
        'slice', i, 'view', GE.view);  
    end
    GE.view = GE.view + 1;
end

% finalize scanloop.txt
toppe.write2loop('finish', GE.sys);   

fprintf('. Sequence ready\n');

%save GE GE  % need GE.sys to plot sequence (for exact timing)

%% check whether the timing of the Pulseq sequence is correct
fprintf('Checking Pulseq timing... ');
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Pulseq timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% plot, etc
seq.plot('TimeRange',[0 TRout*2]); 
figure; toppe.plotseq(1,400,GE.sys); pause(1);

%%
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'mprage');

seq.write('mprage.seq')       % Write to pulseq file
%seq.install('siemens');

%% GE: Create seqstamp file and bundle scan files
toppe.preflightcheck('toppe1.entry', 'seqstamp.txt', GE.sys);
system('tar cf mprage.tar toppe1.entry modules.txt scanloop.txt *.mod seqstamp.txt');

return;



%% visualize the 3D k-space (only makes sense for low-res, otherwise one sees nothing)
%if Nx<=32
    tic;
    [kfa,ta,kf]=seq.calculateKspacePP();
    toc
    figure;plot3(kf(1,:),kf(2,:),kf(3,:));
    hold on;plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
%end



%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slew rate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 