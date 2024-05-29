% ACR spin-echo sequence for quality control
% transversal T1 series

system = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6,...
    'adcDeadTime', 20e-6, 'B0', 2.89 ... % this is Siemens' 3T
);

seq = mr.Sequence(system) ;              % Create a new sequence object
adcDur = 2*2.56e-3 ; 
disp(['readout bandwidht = ', num2str(1/adcDur), ' Hz/pixel']) ;
rfDur1 = 4e-3 ;
rfDur2 = 4e-3 ;
TR = 500e-3 ;
TE = 20e-3 ; % 20ms still works with the chosen parameters & sysyem props
spAx = 2000 ;
spAy = 2000 ;
spAz = 2000 ; % spoiler area in 1/m (=Hz/m*s) % MZ: need 5000 for my oil phantom
ro_os = 2 ;
sliceThickness = 5e-3 ;            % slice thickness
sliceGap = 5e-3 ;               % slice gap
Nslices = 11 ;
fov=250e-3 ; Nx = 256 ;             % Define FOV and resolution
Ny = Nx ;                          % number of radial spokes
Ndummy = 16 ;                    % number of dummy scans
Nrep = 2 ;                      % averages

% Create 90 degree slice selection pulse and gradient and even the
% refocusing gradient; we will not use the latter however but will subtract 
% its area from the first (left) spoiler
[rf_ex, gz, gzr] = mr.makeSLRpulse(pi/2,'duration',rfDur1,'SliceThickness',sliceThickness,...
    'timeBwProduct',5,'dwell',rfDur2/500,'passbandRipple',1,'stopbandRipple',1e-2,...
    'filterType','ms','system',system,'use','excitation', 'PhaseOffset' ,0);

% Create non-selective refocusing pulse
[rf_ref, g_ref] =  mr.makeSLRpulse(pi,'duration',rfDur2,'SliceThickness',sliceThickness,...
    'timeBwProduct',5,'dwell',rfDur2/500,'passbandRipple',1,'stopbandRipple',1e-2,...
    'filterType','ms','system',system,'use','refocusing', 'PhaseOffset' ,pi/2); 

gamma_H1 = 42.58 ; % [MHz/T]
rf_ex_peak = max(abs(rf_ex.signal))/gamma_H1 ; % [uT]
rf_ref_peak = max(abs(rf_ref.signal))/gamma_H1 ; % [uT]
disp(['The peak rf_ex amplitude = ', num2str(rf_ex_peak), ' uT']) ;
disp(['The peak rf_ref amplitude = ', num2str(rf_ref_peak), ' uT']) ;

% join spoilers with the slice selection pulses of the refocusing gradients
g_ref_pre = mr.makeExtendedTrapezoidArea(g_ref.channel, 0, g_ref.amplitude, spAz+gzr.area, system) ;
g_ref_post = mr.makeExtendedTrapezoidArea(g_ref.channel, g_ref.amplitude, 0, spAz, system) ;
g_refC = mr.makeExtendedTrapezoid(g_ref_pre.channel, ...
    'times', [g_ref_pre.tt g_ref_post.tt+g_ref_pre.shape_dur+g_ref.flatTime],...
    'amplitudes', [g_ref_pre.waveform g_ref_post.waveform], 'system', system) ;
rf_ref.delay = g_ref_pre.shape_dur ;
% calculate spoiler gradients in x- and y-directions
% Gx spoilers before and after refocusing gradient
g_SPx_pre = mr.makeTrapezoid('x', 'Area', spAx, 'system', system) ;
g_SPx_post = mr.makeTrapezoid('x', 'Area', spAx, 'system', system, 'delay', mr.calcDuration(g_SPx_pre)+g_ref.flatTime) ;
g_SPx = mr.addGradients({g_SPx_pre, g_SPx_post},'system', system) ;
% Gy spoilers before and after refocusing gradient
g_SPy_pre = mr.makeTrapezoid('y', 'Area', spAy, 'system', system) ;
g_SPy_post = mr.makeTrapezoid('y', 'Area', spAy, 'system', system, 'delay', mr.calcDuration(g_SPy_pre)+g_ref.flatTime) ;
g_SPy = mr.addGradients({g_SPy_pre, g_SPy_post},'system', system) ;
% update delays
if (mr.calcDuration(g_ref_pre) > mr.calcDuration(g_SPx_pre, g_SPy_pre))
    g_SPx.delay = g_SPx.delay + mr.calcDuration(g_ref_pre) - mr.calcDuration(g_SPx_pre) ;
    g_SPy.delay = g_SPy.delay + mr.calcDuration(g_ref_pre) - mr.calcDuration(g_SPy_pre) ;
else
    g_refC.delay = g_refC.delay - mr.calcDuration(g_ref_pre) + mr.calcDuration(g_SPx_pre, g_SPy_pre) ;
    rf_ref.delay = rf_ref.delay - mr.calcDuration(g_ref_pre) + mr.calcDuration(g_SPx_pre, g_SPy_pre) ;
end

% Define delays and ADC events
deltak = 1/fov ;
gr = mr.makeTrapezoid('x', system, 'FlatArea', Nx*deltak, 'FlatTime', ceil(adcDur/system.gradRasterTime)*system.gradRasterTime) ;
adc = mr.makeAdc(Nx*ro_os, system, 'Duration', adcDur, 'delay', gr.riseTime) ;
disp(['ADC dwell time = ', num2str(adc.dwell*1e6), ' us']) ;
grPredur = 2e-3 ; % use a fixed time to make this gradient visible on the plot
grPre = mr.makeTrapezoid('x', 'system', system, 'Area', gr.area/2+deltak/2, 'Duration', grPredur) ;  % we need this "deltak/2" because of the ADC sampling taking place in the middle of the dwell time
phaseAreas = ((0:Ny-1)-Ny/2)*deltak ;
PEscale = phaseAreas / max(abs(phaseAreas)) ;
gy = mr.makeTrapezoid('y','Area', max(abs(phaseAreas)), 'Duration', mr.calcDuration(grPre),'system', system) ;
gx_spoil = mr.makeTrapezoid('x','Area', spAx, 'Duration', mr.calcDuration(gy),'system', system) ;
gz_spoil = mr.makeTrapezoid('z','Area', spAz, 'Duration', mr.calcDuration(gy),'system', system) ;

% slice positions
slicePositions = (sliceThickness + sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2) ;
slicePositions=slicePositions([1:2:Nslices 2:2:Nslices]) ; % reorder slices for an interleaved acquisition (optional)
%slicePositions=slicePositions([1:3:Nslices 2:3:Nslices 3:3:Nslices]); % reorder slices for an interleaved acquisition (optional)

delayTE1 = TE/2 - ( mr.calcDuration(gz) - rf_ex.shape_dur/2 - rf_ex.delay) ...
    - mr.calcDuration(grPre) ...
    - mr.calcDuration(g_refC) + g_ref.flatTime/2 + mr.calcDuration(g_ref_post) ;
delayTE2 = TE/2 - g_ref.flatTime/2 - mr.calcDuration(g_ref_post) - mr.calcDuration(gr)/2 ;
delayTR = TR - Nslices * ( rf_ex.delay + rf_ex.shape_dur/2 + TE + mr.calcDuration(gr)/2 + mr.calcDuration(gy)) ;
delayTR_1slice = ceil(delayTR/Nslices/system.blockDurationRaster) * system.blockDurationRaster ;

assert(delayTE1 >= 0) ;
assert(delayTE2 >= 0) ;
assert(delayTR >= 0) ;

% change orientation to match the siemens product sequence
% reverse the polarity of all gradients in readout direction (Gz)
grPre.amplitude = -grPre.amplitude ;
g_SPx = mr.scaleGrad(g_SPx, -1) ;
gr.amplitude = -gr.amplitude ;
gx_spoil.amplitude = -gx_spoil.amplitude ;

seq.addBlock(mr.makeLabel('SET','REP', 0)) ;
for r=1:Nrep
    seq.addBlock(mr.makeLabel('SET','LIN', 0) ) ;
    for i=(1-Ndummy):Ny
        % loop over slices
        seq.addBlock(mr.makeLabel('SET','SLC', 0)) ;
        for s=1:Nslices
            rf_ex.freqOffset = gz.amplitude * slicePositions(s) ;
            rf_ex.phaseOffset = -2*pi*rf_ex.freqOffset*mr.calcRfCenter(rf_ex) ; % compensate for the slice-offset induced phase, unit radian (2*pi radians in a circle)
            rf_ref.freqOffset = g_ref.amplitude * slicePositions(s) ;
            rf_ref.phaseOffset = pi/2 - 2*pi*rf_ref.freqOffset*mr.calcRfCenter(rf_ref) ;
            seq.addBlock(rf_ex, gz) ;
            if (i>0) % semi-negative index -- dummy scans
                seq.addBlock(grPre, mr.scaleGrad(gy, PEscale(i))) ;
            else
                seq.addBlock(grPre) ;
            end
            seq.addBlock(mr.makeDelay(delayTE1)) ;
            seq.addBlock(rf_ref, g_refC, g_SPx, g_SPy) ;
            seq.addBlock(mr.makeDelay(delayTE2));
            if (i>0) % semi-negative index -- dummy scans
                seq.addBlock(adc, gr) ;
                seq.addBlock(mr.makeLabel('INC', 'SLC', 1)) ;
                seq.addBlock(gx_spoil, mr.scaleGrad(gy, PEscale(i)), gz_spoil) ;
            else
                seq.addBlock(gr) ;
                seq.addBlock(gx_spoil, gz_spoil) ;
            end
            seq.addBlock(mr.makeDelay(delayTR_1slice)) ;
        end
        % seq.addBlock(mr.makeDelay(delayTR)) ;
        if (i>0)
            seq.addBlock(mr.makeLabel('INC', 'LIN', 1)) ;
        end
    end
    seq.addBlock(mr.makeLabel('INC','REP', 1)) ;
end

% show the first non-dummy TR with the block structure
seq.plot('showBlocks', 1, 'timeRange', TR*(Ndummy+[0 2]), 'timeDisp', 'us') ;

% check whether the timing of the sequence is compatible with the scanner
[ok, error_report] = seq.checkTiming ;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% evaluate label settings more specifically
%seq.plot('timeRange', [0 32]*TRout, 'TimeDisp', 'ms', 'Label', 'LIN');
adc_lbl=seq.evalLabels('evolution','adc');
figure; plot(adc_lbl.REP);
hold on; plot(adc_lbl.SLC);
plot(adc_lbl.LIN) ;
legend('REF','SLC', 'LIN');
title('evolution of labels/counters');
%%
seq.setDefinition('Name', 'QA_T1');
seq.setDefinition('FOV', [fov fov max(slicePositions)-min(slicePositions)+sliceThickness]);
% the following definitions have effect in conjunction with LABELs 
seq.setDefinition('SlicePositions', slicePositions);
seq.setDefinition('SliceThickness', sliceThickness);
seq.setDefinition('SliceGap', sliceGap);
seq.setDefinition('ReadoutOversamplingFactor', ro_os) ;

seq.write('QA_T1.seq')       % Write to pulseq file
%seq.install('siemens');    % copy to scanner

% calculate k-space but only use it to check timing
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay',[0 0 0]*1e-6); % play with anisotropic trajectory delays -- zoom in to see the trouble ;-)

if Ndummy==0 % if nDummy equals 0 we can do additional checks
    assert(abs(t_refocusing(1)-t_excitation(1)-TE/2)<1e-6); % check that the refocusing happens at the 1/2 of TE
    assert(abs(t_adc(Nx/2)-t_excitation(1)-TE)<adc.dwell); % check that the echo happens as close as possible to the middle of the ADC elent
end

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kz-axis
title('k-vector components as functions of time'); xlabel('time /s'); ylabel('k-component /m^-^1');
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('2D k-space trajectory'); xlabel('k_x /m^-^1'); ylabel('k_y /m^-^1');
