% a basic UTE-like sequence
% achieves "TE" below 100 us

seq=mr.Sequence();              % Create a new sequence object
fov=250e-3; Nx=256;             % Define FOV and resolution
alpha=30;                       % flip angle
sliceThickness=3e-3;            % slice
TR=10e-3;                       % TR
Nr=321;                         % number of radial spokes
Ndummy=20;                      % number of dummy scans
delta= 2* pi / Nr;              % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it
rf_duration=1.0e-3;             % duration of the excitation pulse
ro_duration=1.6e-3;             % read-out time: controls RO bandwidth and T2-blurring
ro_os=2;                        % oversampling
minRF_to_ADC_time=70e-6;        % the parameter wich defines TE together with ro_discard
ro_discard=0;                   % dummy ADC samples to discard (due to ADC filter 
ro_spoil=2;                     % extend RO to achieve spoiling

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment

% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 0e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

%% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',rf_duration,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',2,...
    'centerpos',1,'system',sys);

% resample the RF pulse to the ramp
gza=[0 1 1 0];
gzt=cumsum([0 gz.riseTime gz.flatTime gz.fallTime]);
gzas_0=interp1(gzt+gz.delay,gza,rf.t+rf.delay);
rft_1=[sys.rfRasterTime:sys.rfRasterTime:rf_duration+0.5*gz.fallTime];
gzas_1=interp1(gzt+gz.delay,gza,rft_1+rf.delay+gz.fallTime*0.5);
kzs_0=cumsum(gzas_0);
kzs_1=cumsum(gzas_1);
kzs_0=kzs_0-max(kzs_0);
kzs_1=kzs_1-max(kzs_1);
rfs_1=interp1(kzs_0,rf.signal,kzs_1);
%figure; plot(kzs_1, abs(rfs_1),'r.');
%hold; plot(kzs_0, abs(rf.signal),'b-');
rf.t=rft_1;
rf.signal=rfs_1.*gzas_1;
gz.flatTime=gz.flatTime-gz.fallTime*0.5; % oops, we can get off gradient raster here, FIXME

% Align RO assymmetry to ADC samples
Nxo=round(ro_os*Nx);
% Define other gradients and ADC events
deltak=1/fov/2;
ro_area=Nx*deltak;
gx = mr.makeTrapezoid('x','FlatArea',ro_area,'FlatTime',ro_duration,'system',sys);
adc = mr.makeAdc(Nxo,'Duration',gx.flatTime,'system',sys);

% ro-spoiling
gx.flatTime=gx.flatTime*ro_spoil;

% Calculate timing
%ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
%    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
% calculate actual achieved TE
TE = ceil((minRF_to_ADC_time + adc.dwell*ro_discard)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gz) ...
    - mr.calcDuration(gx) - TE)/seq.gradRasterTime)*seq.gradRasterTime;

fprintf('TE= %d us\n', round(TE*1e6));

% set up timing
gx.delay=mr.calcDuration(gz)+TE;
adc.delay=gx.delay-adc.dwell*ro_discard;

rf_phase=0;
rf_inc=0;

for i=(-Ndummy):Nr
    for c=1:2
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        % UTE: alternate GZ
        gz.amplitude=-gz.amplitude;
        %
        phi=delta*(i-1);
        grc=gx;      grs=gx;      grc.amplitude=gx.amplitude*cos(phi);      grs.amplitude=gx.amplitude*sin(phi);      grs.channel='y';
        if (i>0)
            seq.addBlock(rf,gz,grc,grs,adc);
        else
            seq.addBlock(rf,gz,grc,grs);
        end
        seq.addBlock(mr.makeDelay(delayTR));
    end
end

seq.plot();

%% plot gradients to check for gaps and optimality of the timing
gw=seq.gradient_waveforms();
figure; plot(gw'); % plot the entire gradient shape

%% trajectory calculation
[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

% plot k-spaces
time_axis=(1:(size(ktraj,2)))*sys.gradRasterTime;
figure; plot(time_axis, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%
seq.setDefinition('FOV', [fov fov sliceThickness]*1e3);
seq.setDefinition('Name', 'ute_rs');

seq.write('ute_rs.seq');       % Write to pulseq file

%seq.install('siemens');
