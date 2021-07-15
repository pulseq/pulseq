% a very basic UTE-like sequence, without ramp-sampling, ramp-RF and other
% tricks yet. Achieves TE in the range of 300-400 us

% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 100, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq=mr.Sequence(sys);           % Create a new sequence object
fov=250e-3; Nx=256;             % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=3e-3;            % slice
TR=10e-3;                       % TR
Nr=128;                         % number of radial spokes
delta= 2* pi / Nr;              % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it
ro_duration=2.56e-3;            % read-out time: controls RO bandwidth and T2-blurring
ro_os=2;                        % oversampling
ro_asymmetry=1;                 % 0: fully symmetric 1: half-echo
minRF_to_ADC_time=50e-6;        % the parameter wich defines TE (together with the RO asymmetyry)

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment

% Create alpha-degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',1e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',2,...
    'centerpos',1,'system',sys);

% Align RO assymmetry to ADC samples
Nxo=round(ro_os*Nx);
ro_asymmetry = round(ro_asymmetry*Nxo/2)/Nxo*2; % check whether we need to use 2Nx or so...
% Define other gradients and ADC events
deltak=1/fov/(1+ro_asymmetry);
ro_area=Nx*deltak;
gx = mr.makeTrapezoid('x','FlatArea',ro_area,'FlatTime',ro_duration,'system',sys);
adc = mr.makeAdc(Nxo,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-(gx.area-ro_area)/2 -gx.amplitude*adc.dwell/2 - ro_area/2*(1-ro_asymmetry),'system',sys);

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',0.2*Nx*deltak,'system',sys);

% Calculate timing
%ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
%    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
% calculate actual achieved TE
TE = gz.fallTime + mr.calcDuration(gxPre,gzReph)+gx.riseTime + adc.dwell*Nxo/2*(1-ro_asymmetry);
delayTR=ceil((TR - mr.calcDuration(gxPre,gzReph) - mr.calcDuration(gz) ...
    - mr.calcDuration(gx))/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTR>=mr.calcDuration(gxSpoil)));

fprintf('TE= %d us\n', round(TE*1e6));

if mr.calcDuration(gzReph) > mr.calcDuration(gxPre) 
    gxPre.delay=mr.calcDuration(gzReph) - mr.calcDuration(gxPre);
end

rf_phase=0;
rf_inc=0;

for i=1:Nr
    for c=1:2
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        % UTE: alternate GZ
        gz.amplitude=-gz.amplitude;
        gzReph.amplitude=-gzReph.amplitude;
        %
        seq.addBlock(rf,gz);
        phi=delta*(i-1);
        gpc=gxPre;   gps=gxPre;   gpc.amplitude=gxPre.amplitude*cos(phi);   gps.amplitude=gxPre.amplitude*sin(phi);   gps.channel='y';
        grc=gx;      grs=gx;      grc.amplitude=gx.amplitude*cos(phi);      grs.amplitude=gx.amplitude*sin(phi);      grs.channel='y';
        gsc=gxSpoil; gss=gxSpoil; gsc.amplitude=gxSpoil.amplitude*cos(phi); gss.amplitude=gxSpoil.amplitude*sin(phi); gss.channel='y';
        seq.addBlock(gpc,gps,gzReph);
        seq.addBlock(grc,grs,adc);
        seq.addBlock(gsc,gss,mr.makeDelay(delayTR));
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

%%
seq.plot();

%% plot gradients to check for gaps and optimality of the timing
gw=seq.waveforms_and_times();
figure; plot(gw{1}(1,:),gw{1}(2,:),gw{2}(1,:),gw{2}(2,:),gw{3}(1,:),gw{3}(2,:)); % plot the entire gradient shape

%% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'ute');

seq.write('ute.seq');       % Write to pulseq file

%seq.install('siemens');
