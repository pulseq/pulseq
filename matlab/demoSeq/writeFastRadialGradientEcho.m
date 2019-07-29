seq=mr.Sequence();              % Create a new sequence object
fov=256e-3; Nx=128;             % Define FOV and resolution
alpha=5;                       % flip angle
sliceThickness=4e-3;            % slice
Nr=128;                         % number of radial spokes
Ndummy=0;                      % number of dummy scans
delta=pi / Nr;                  % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it
ro_dur=600e-6;                  % RO duration
ro_os=2;                        % readout oversampling
ro_spoil=0;%0.5;                   % additional k-max excursion for RO spoiling
sl_spoil=0;%4;                     % spoil area compared to the slice thickness

% more in-depth parameters
rfSpoilingInc=0;              % RF spoiling increment

% set system limits
sys = mr.opts('MaxGrad', 35, 'GradUnit', 'mT/m', ...
    'MaxSlew', 180, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% Create alpha-degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',400e-6,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',2,'system',sys);
gzReph.delay=mr.calcDuration(gz);
gzComb=mr.addGradients({gz, gzReph}, 'system', sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',ro_dur,'system',sys);
adc = mr.makeAdc(Nx*ro_os,'Duration',ro_dur,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.flatArea/Nx*floor(Nx/2)-(gx.area-gx.flatArea)/2,'system',sys);
%
gxPre=mr.align('right', gxPre, 'right', gzComb);
addDelay=mr.calcDuration(rf)-gxPre.delay;
if addDelay>0
    gxPre.delay = gxPre.delay+ceil(addDelay/sys.gradRasterTime)*sys.gradRasterTime;
end

% gradient spoiling
if sl_spoil>0
    sp_area_needed=sl_spoil/sliceThickness-gz.area/2;
    gzSpoil=mr.makeTrapezoid('z','Area',sp_area_needed,'system',sys,'Delay',gx.riseTime+gx.flatTime);
else
    gzSpoil=[];
end

if ro_spoil>0
    ro_spoil_area=(gx.area-gx.flatArea)/2;
    ro_add_time=ceil(((gx.area/Nx*(Nx/2+1)*(1+ro_spoil))/gx.amplitude)/sys.gradRasterTime)*sys.gradRasterTime;
    gx.flatTime=gx.flatTime+ro_add_time; % careful, areas stored in the object are now wrong
end

% join slice spoiler with the slice selection
%if (rf.delay>mr.calcDuration()) no, this does not work to be really optimal we need a new function with start, stop and area

% Calculate timing
% TODO: just calculate actual TE and TR here

% start the sequence
rf_phase=0;
rf_inc=0;
TR=0;
for i=(1-Ndummy):Nr
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);
    %
    phi=delta*(i-1);
    seq.addBlock(mr.rotate('z',phi,rf,gzComb,gxPre));
    if (i>0)
        seq.addBlock(mr.rotate('z',phi,gx,adc,gzSpoil));
    else
        seq.addBlock(mr.rotate('z',phi,gx,gzSpoil));
    end
    if TR<=0
        TR=seq.duration;
    end
end

seq.plot();
%return;
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
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'gre_rad');

seq.write('gre_rad.seq')       % Write to pulseq file

%seq.install('siemens');

return;
%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 

%out = function makeTransitionArea(dur, area, endampl, sys)
% assume we know ampl (plato amplitude)
% tr1=ampl/sys.maxSlew;
% tr2=abs(ampl-endampl)/sys.maxSlew; =adampl/sys.maxSlew
% tp=dur-tr1-tr2;
% area=(tr1/2+tp)*ampl        + (ampl+endampl)*tr2/2
%     =(dur-tr1/2)*ampl + (endampl-ampl)*tr2/2
%     =dur*ampl - ampl^2/2/sys.maxSlew + (endampl-ampl)*adampl/2/sys.maxSlew
% depending on the sign of (ampl-endampl) there are two options for the
% second order equation. There are 0,1 or two solutions if attemting to
% invert it (and two options above)
%end


