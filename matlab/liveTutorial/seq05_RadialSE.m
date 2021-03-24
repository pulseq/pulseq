system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
adcDur=51.2e-3; 
rfDur1=3e-3;
rfDur2=1e-3;
TR=200e-3;
TE=60e-3;
spA=1000; % spoiler area in 1/m (=Hz/m*s)

sliceThickness=3e-3;            % slice
fov=250e-3; Nx=256;             % Define FOV and resolution
Nr=256;                         % number of radial spokes
Ndummy=0;                      % number of dummy scans
delta=2*pi / Nr;                % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it


% Create 90 degree slice selection pulse and gradient and even the
% refocusing gradient; we will not use it however but will bubtract its
% area from the first spoiler
[rf_ex, gs, gsr] = mr.makeSincPulse(pi/2,system,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.4,'timeBwProduct',4);
gs.channel='x'; % change it to X because we want sagittal orientation

% Create non-selective refocusing pulse
rf_ref = mr.makeBlockPulse(pi,'Duration',rfDur2, 'system', system, 'use', 'refocusing'); % needed for the proper k-space calculation
    
% calculate spoiler gradient
g_sp1=mr.makeTrapezoid('x','Area',spA-gsr.area,'system',system);
rf_ref.delay=max(mr.calcDuration(g_sp1),rf_ref.delay);

g_sp2=mr.makeTrapezoid('x','Area',spA,'system',system);

% Define delays and ADC events
delayTE1=TE/2-(mr.calcDuration(gs)-mr.calcRfCenter(rf_ex)-rf_ex.delay)-rf_ref.delay-mr.calcRfCenter(rf_ref);
delayTE2=TE/2-mr.calcDuration(rf_ref)+rf_ref.delay+mr.calcRfCenter(rf_ref)-adcDur/2; % this is not perfect, but -adcDur/2/Nx  will break the raster alignment
assert(delayTE2>mr.calcDuration(g_sp1));

deltak=1/fov;
gr = mr.makeTrapezoid('z',system,'FlatArea',Nx*deltak,'FlatTime',adcDur);
adc = mr.makeAdc(Nx,system,'Duration',adcDur,'delay',delayTE2);
gr.delay=delayTE2-gr.riseTime;

grPre = mr.makeTrapezoid('z',system,'Area',gr.area/2+deltak/2,'Duration',delayTE1);

delayTR=TR-mr.calcDuration(rf_ex)-delayTE1-mr.calcDuration(rf_ref);

assert(delayTE1>=0);
assert(delayTE2>=0);
assert(delayTR>=0);

% Loop over repetitions and define sequence blocks
for i=(1-Ndummy):Nr
    seq.addBlock(rf_ex,gs);
    %seq.addBlock(grPre); 
    seq.addBlock(mr.rotate('x',delta*(i-1),grPre)); 
    seq.addBlock(rf_ref,g_sp1);
    if (i>0)
        %seq.addBlock(adc,gr,g_sp1,mr.makeDelay(delayTR)); 
        seq.addBlock(mr.rotate('x',delta*(i-1),adc,gr,g_sp2,mr.makeDelay(delayTR)));  
    else
        %seq.addBlock(gr,g_sp1,mr.makeDelay(delayTR));  
        seq.addBlock(mr.rotate('x',delta*(i-1),gr,g_sp2,mr.makeDelay(delayTR)));  
    end
end

seq.plot();%('showBlocks',1);

% check whether the timing of the sequence is compatible with the scanner
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.write('se_radial.seq')       % Write to pulseq file
%seq.install('siemens');    % copy to scanner

% calculate k-space but only use it to check timing
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay',[0 0 0]*1e-6); % play with anisotropic trajectory delays -- zoom in to see the trouble ;-)

if Ndummy==0
    assert(abs(t_refocusing(1)-t_excitation(1)-TE/2)<1e-6); % check that the refocusing happens at the 1/2 of TE
    assert(abs(t_adc(Nx/2)-t_excitation(1)-TE)<adc.dwell); % check that the echo happens as close as possible to the middle of the ADC elent
end

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(3,:),'.'); % and sampling points on the kz-axis
figure; plot(ktraj(3,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;plot(ktraj_adc(3,:),ktraj_adc(2,:),'r.'); % plot the sampling points
