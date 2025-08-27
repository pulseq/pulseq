% set system limits (slew rate 130 and max_grad 30 work on Prisma)
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq=mr.Sequence(sys);           % Create a new sequence object
fov=240e-3; Nx=240;             % Define FOV and resolution
alpha=5;                        % flip angle
sliceThickness=6e-3;            % slice
Nr=256;                         % number of radial spokes
Ndummy=10;                      % number of dummy scans
delta= pi / Nr;                 % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it
ro_dur=1200e-6;                 % RO duration
ro_os=2;                        % readout oversampling
ro_spoil=0.5;                   % additional k-max excursion for RO spoiling
sl_spoil=1.5;                   % spoil area compared to the slice thickness

% TE & TR are as short as possible derived from the above parameters and
% the system specs below

% more in-depth parameters
rfSpoilingInc=84;              % RF spoiling increment
rfTimeBwProd = 2;              % time-bandwidth product for the RF pulse

% Create alpha-degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',400e-6,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',rfTimeBwProd,'system',sys,...
    'use','excitation');
% gradient spoiling in slice direction
if sl_spoil>0
    sp_area_needed=sl_spoil/sliceThickness*rfTimeBwProd-gz.flatArea/2;
    %gzSpoil=mr.makeTrapezoid('z','Area',sp_area_needed,'system',sys,'Delay',gx.riseTime+gx.flatTime);
    gzSpoil=mr.makeExtendedTrapezoidArea('z', 0, gz.amplitude, sp_area_needed, sys);
    gz=mr.makeExtendedTrapezoid('z','times',[0, gz.flatTime,gz.flatTime+gz.fallTime],'amplitudes',[gz.amplitude,gz.amplitude,0],'system',sys);
    gz.delay=gzSpoil.shape_dur;
    gz=mr.addGradients({gz, gzSpoil}, 'system', sys);
    rf.delay=max(gzSpoil.shape_dur,sys.rfDeadTime);
    gz.delay=rf.delay-gzSpoil.shape_dur;
end
% join rephaser
gzReph.delay=mr.calcDuration(gz);
gzComb=mr.addGradients({gz, gzReph}, 'system', sys);
[gzSpoil,gz]=mr.splitGradientAt(gzComb,rf.delay-sys.rfDeadTime);
gz.delay=gz.delay-mr.calcDuration(gzSpoil);
rf.delay=sys.rfDeadTime;
clear gzComb;

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','Amplitude',Nx*deltak/ro_dur,'FlatTime',ceil(ro_dur/sys.gradRasterTime)*sys.gradRasterTime,'system',sys);
adc = mr.makeAdc(Nx*ro_os,'Duration',ro_dur,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.amplitude*(ro_dur/Nx/ro_os*(Nx*ro_os/2-0.5)+0.5*gx.riseTime),'system',sys); % 0.5 is necessary to acount for the Siemens sampling in the center of the dwell periods
%
[gxPre,~]=mr.align('right', gxPre, 'right', gz);
addDelay=mr.calcDuration(rf)-gxPre.delay;
if addDelay>0
    gxPre.delay = gxPre.delay+ceil(addDelay/sys.gradRasterTime)*sys.gradRasterTime;
end

% gradient spoiling in slice direction
if ro_spoil>0
    %ro_spoil_area=(gx.area-gx.flatArea)/2;
    ro_add_time=ceil(((gx.area/Nx*(Nx/2+1)*ro_spoil)/gx.amplitude)/sys.gradRasterTime)*sys.gradRasterTime;
    gx.flatTime=gx.flatTime+ro_add_time; % careful, areas stored in the object are now wrong
end

% join slice spoiler with the slice selection
%if (rf.delay>mr.calcDuration()) no, this does not work to be really optimal we need a new function with start, stop and area
% could be done with mr.makeExtendedTrapezoidArea()

% Calculate timing
% TODO: just calculate actual TE and TR here

% start the sequence
rf_phase=0;
rf_inc=0;
TR=0;
seq.addBlock(gzSpoil);
[gzSpoil,~,~]=mr.align('right',gzSpoil,'left',gx,mr.makeDelay(ceil((adc.delay+adc.dwell*adc.numSamples)/sys.gradRasterTime)*sys.gradRasterTime+gzSpoil.shape_dur));
for i=(1-Ndummy):Nr
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);
    %
    phi=delta*(i-1);
    seq.addBlock(mr.rotate('z',phi,rf,gz,gxPre));
    if (i>0)
        seq.addBlock(mr.rotate('z',phi,gx,adc,gzSpoil));
    else
        seq.addBlock(mr.rotate('z',phi,gx,gzSpoil));
    end
    if TR<=0
        TR=seq.duration;
    end
end
% dummy slice select to ramp down to 0 Z spoil gradient
seq.addBlock(gz);

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% plot, etc
seq.plot();
%return;
%% trajectory calculation
%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'gre_rad');

seq.write('fast_gre_rad.seq')       % Write to pulseq file

%seq.install('siemens');

return;
%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slew rate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 



