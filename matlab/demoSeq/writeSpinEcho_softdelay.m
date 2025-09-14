% spin-echo sequence with adjustable TE and TR using the soft delay
% extension

system = mr.opts('MaxGrad',24,'GradUnit','mT/m',...
    'MaxSlew',50,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6,...
    'adcDeadTime', 20e-6, 'B0', 2.89 ... % this is Siemens' 3T
);
vendor = 'siemens' ;
seq = mr.Sequence(system) ;              % Create a new sequence object
adcDur = 2*2.56e-3 ; 
disp(['readout bandwidth = ', num2str(1/adcDur), ' Hz/pixel']) ;
rfDur1 = 3e-3 ;
rfDur2 = 8.8e-3 ;
TR = 1400e-3 ;
TE = 19e-3 ; % 19ms still works with the chosen parameters & system props
spAx = 0 ;
spAy = 0 ;
spAz = 1500 ; % spoiler area in 1/m (=Hz/m*s) % MZ: need 5000 for my oil phantom
ro_os = 2 ;
sliceThickness = 5e-3 ;            % slice thickness
sliceGap = 5e-3 ;               % slice gap
Nslices = 2 ;
fov=250e-3 ; Nx = 256 ;             % Define FOV and resolution
Ny = Nx ;                          % number of radial spokes
Ndummy = 1 ;                    % number of dummy scans
Nrep = 1 ;                      % averages
% correction factors RF
sth_ex = 1 ;
sth_ref = 1.25 ;
max_TE = 120e-3 ;

% Create 90 degree slice selection pulse and gradient and even the
% refocusing gradient; we will not use the latter however but will subtract 
% its area from the first (left) spoiler
[rf_ex, gz, gzr] = mr.makeSLRpulse(pi/2,'duration',rfDur1,'SliceThickness',sliceThickness*sth_ex,...
    'timeBwProduct',5,'dwell',rfDur1/500,'passbandRipple',1,'stopbandRipple',1e-2,...
    'filterType','ms','system',system,'use','excitation', 'PhaseOffset' ,pi/2); % MZ: other RF cycle

% Create selective refocusing pulse
[rf_ref, g_ref] =  mr.makeSLRpulse(pi,'duration',rfDur2,'SliceThickness',sliceThickness*sth_ref,...
    'timeBwProduct',6,'dwell',rfDur2/500,'passbandRipple',1,'stopbandRipple',1e-2,...
    'filterType','ms','system',system,'use','refocusing', 'PhaseOffset' ,0); % MZ: other RF cycle

% check RF profile alighnment
[M_z90,M_xy90,F2_90]=mr.simRf(rf_ex);
sl_th_90=mr.aux.findFlank(F2_90(end:-1:1)/gz.amplitude,abs(M_xy90(end:-1:1)),0.5)-mr.aux.findFlank(F2_90/gz.amplitude,abs(M_xy90),0.5);
[M_z180,M_xy180,F2_180,ref_eff]=mr.simRf(rf_ref);
sl_th_180=mr.aux.findFlank(F2_180(end:-1:1)/g_ref.amplitude,abs(ref_eff(end:-1:1)),0.5)-mr.aux.findFlank(F2_180/g_ref.amplitude,abs(ref_eff),0.5);
figure; plot(F2_90/gz.amplitude*1000,abs(M_xy90)); title('RF profiles'); xlabel('through-slice pos, mm');
hold on; plot(F2_180/g_ref.amplitude*1000,abs(ref_eff)); legend('ex (Mxy)','ref (ref-eff)');
fprintf('slice thicknes 90-degree excitation pulse: %.3f mm\n',sl_th_90*1e3);
fprintf('slice thicknes 180-degree refocusing pulse: %.3f mm\n',sl_th_180*1e3);

rf_ex_peak = max(abs(rf_ex.signal))/system.gamma ; % [uT]
rf_ref_peak = max(abs(rf_ref.signal))/system.gamma ; % [uT]
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

grPredur = mr.calcDuration(g_ref_post) ; % use a fixed time to make this gradient visible on the plot
grPre = mr.makeTrapezoid('x', 'system', system, 'Area', (gr.area/2+deltak/2), 'Duration', grPredur) ;  % we need this "deltak/2" because of the ADC sampling taking place in the middle of the dwell time
phaseAreas = ((0:Ny-1)-Ny/2)*deltak ;
PEscale = phaseAreas / max(abs(phaseAreas)) ;
gyPre = mr.makeTrapezoid('y','Area', -max(abs(phaseAreas)), 'Duration', grPredur,'system', system ) ;

gyPost = mr.makeTrapezoid('y','Area', max(abs(phaseAreas)), 'Duration', grPredur,'system', system) ;
gx_spoil = mr.makeTrapezoid('x','Area', spAx,'system', system ) ; %, 'Duration', mr.calcDuration(gy)
gz_spoil = mr.makeTrapezoid('z','Area', spAz,'system', system ) ; %, 'Duration', mr.calcDuration(gy)

% slice positions
slicePositions = (sliceThickness + sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2) ;
slicePositions = slicePositions([1:2:Nslices 2:2:Nslices]) ; % reorder slices for an interleaved acquisition (optional)
%slicePositions=slicePositions([1:3:Nslices 2:3:Nslices 3:3:Nslices]); % reorder slices for an interleaved acquisition (optional)

delayTE1 = TE/2 - ( mr.calcDuration(gz) - rf_ex.center - rf_ex.delay) ...
    - rf_ref.center - rf_ref.delay - mr.calcDuration(grPre,gyPre);
delayTE2 = TE/2 - ( mr.calcDuration(g_refC, g_SPx, g_SPy) - rf_ref.center - rf_ref.delay ) - mr.calcDuration(gr)/2 ;
delayTR = TR - Nslices * ( rf_ex.delay + rf_ex.shape_dur/2 + TE + mr.calcDuration(gr)/2 + mr.calcDuration(gyPost, gx_spoil, gz_spoil)+max_TE-TE) ;
delayTR_1slice = ceil(delayTR/Nslices/system.blockDurationRaster) * system.blockDurationRaster ;

delayTE1 = round(delayTE1 / system.gradRasterTime) * system.gradRasterTime ;
delayTE2 = round(delayTE2 / system.gradRasterTime) * system.gradRasterTime ;
delayTR_1slice = round(delayTR_1slice / system.gradRasterTime) * system.gradRasterTime ;

assert(delayTE1 >= 10e-6); % for softDelays it needs to be positive
assert(delayTE2 >= 10e-6); % for softDelays it needs to be positive
assert(delayTR_1slice > 0) ;

% change orientation to match the siemens product sequence
% reverse the polarity of all gradients in readout direction (Gz)
if vendor(1) == 's' || vendor(1) == 'S' % if vendor is Siemens
    grPre = mr.scaleGrad(grPre, -1) ;
    g_SPx = mr.scaleGrad(g_SPx, -1) ;
    gr = mr.scaleGrad(gr, -1) ;
    gx_spoil = mr.scaleGrad(gx_spoil, -1) ;
elseif vendor(1) == 'g' || vendor(1) == 'G' % if vendor is GE
    gyPre = mr.scaleGrad(gyPre, -1) ;
    gyPost = mr.scaleGrad(gyPost, -1) ;
    g_SPy = mr.scaleGrad(g_SPy, -1) ;
else
    disp('Please define your vendor informaiton (Siemens or GE)') ;
    return ;
end

% seq.addBlock(mr.makeLabel('SET','REP', 0)) ;
for r=1:Nrep
    seq.addBlock(mr.makeLabel('SET','LIN', 0) ) ;
    for i=(1-Ndummy):Ny
        % loop over slices
        seq.addBlock(mr.makeLabel('SET','SLC', 0)) ;
        for s=1:Nslices
            % MZ: we alternate RF phase for the RF excitation and ADC with the i-counter
            adc.phaseOffset = mod(i,2)*pi;
            rf_ex.freqOffset = gz.amplitude * slicePositions(s) ;
            rf_ref.freqOffset = g_ref.amplitude * slicePositions(s) ;
            if vendor(1)=='g' || vendor(1) == 'G'
                rf_ex.freqOffset = -rf_ex.freqOffset ;
                rf_ref.freqOffset = -rf_ref.freqOffset ;
            end
            rf_ex.phaseOffset = mod(i,2)*pi+pi/2 -2*pi*rf_ex.freqOffset*mr.calcRfCenter(rf_ex) ; % compensate for the slice-offset induced phase, unit radian (2*pi radians in a circle)
            rf_ref.phaseOffset = -2*pi*rf_ref.freqOffset*mr.calcRfCenter(rf_ref) ;
            
            seq.addBlock(rf_ex, gz) ;
            % softdelay TE1
            % seq.addBlock(mr.makeDelay(delayTE1)) ;
            seq.addBlock(delayTE1, mr.makeSoftDelay(0,'TE','offset',delayTE1-TE/2, 'factor', 2) ) ;
            if (i>0) % semi-negative index -- dummy scans
                seq.addBlock(grPre, mr.scaleGrad(gyPre, PEscale(i)));
                seq.addBlock(rf_ref, g_refC, g_SPx, g_SPy) ;
            else
                seq.addBlock(grPre);
                seq.addBlock(rf_ref, g_refC, g_SPx, g_SPy) ;
            end
            % softdelay TE2
            % seq.addBlock(mr.makeDelay(delayTE2));
            seq.addBlock(delayTE2, mr.makeSoftDelay(0,'TE','offset',delayTE2-TE/2, 'factor', 2) ) ;
            if (i>0) % semi-negative index -- dummy scans
                seq.addBlock(adc, gr) ;
                seq.addBlock(mr.makeLabel('INC', 'SLC', 1)) ;
                seq.addBlock(gx_spoil, mr.scaleGrad(gyPost, PEscale(i)), gz_spoil) ;
            else
                seq.addBlock(gr) ;
                seq.addBlock(gx_spoil, gz_spoil) ;
            end
            % softdelay for resolving TE/TR dependency -- quite a tricky thing to calculate... but checkTiming() helps to find out if it was right
            seq.addBlock(max_TE-TE, mr.makeSoftDelay(0,'TE','factor',-1,'offset',max_TE));
            if r==1 && i==(1-Ndummy) && s == 1
                durPerSlc = sum(seq.blockDurations) ;
            end
            % seq.addBlock(mr.makeDelay(delayTR_1slice)) ;
            seq.addBlock(delayTR_1slice,mr.makeSoftDelay(1,'TR','offset',-durPerSlc,'factor',Nslices));
        end
        % seq.addBlock(mr.makeDelay(delayTR)) ;
        if (i>0)
            seq.addBlock(mr.makeLabel('INC', 'LIN', 1)) ;
        end
    end
    seq.addBlock(mr.makeLabel('INC','REP', 1)) ;
end

% show the first non-dummy TR with the block structure
%seq.plot('showBlocks', 1, 'timeRange', TR*(Ndummy+[0 2]), 'timeDisp', 'us') ;
seq.plot('showBlocks', 1, 'timeRange', TR*(Ndummy+[0 1]), 'timeDisp', 'us', 'stacked', 1) ;
%% check whether the timing of the sequence is compatible with the scanner
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
seq.setDefinition('ReceiverGainHigh',1);

seq.write('se_softdelay.seq')       % Write to pulseq file
%seq.install('siemens');    % copy to scanner

%% calculate k-space but only use it to check timing
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP;%('blockRange',[1,150]);
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay',[0 0 0]*1e-6); % play with anisotropic trajectory delays -- zoom in to see the trouble ;-)

%if Ndummy==0 % if nDummy equals 0 we can do additional checks
%    assert(abs(t_refocusing(1)-t_excitation(1)-TE/2)<1e-6); % check that the refocusing happens at the 1/2 of TE
%    assert(abs(t_adc(Nx/2)-t_excitation(1)-TE)<adc.dwell); % check that the echo happens as close as possible to the middle of the ADC elent
%end

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kz-axis
title('k-vector components as functions of time'); xlabel('time /s'); ylabel('k-component /m^-^1');
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('2D k-space trajectory'); xlabel('k_x /m^-^1'); ylabel('k_y /m^-^1');

% check TE and TR mannually

t_echo=t_adc(round(Nx*ro_os/2)+1);
i_ex=find(t_excitation<t_echo,1,'last');
i_ref=find(t_refocusing<t_echo,1,'last');

fprintf('Spin echo is produced at %g us after the excitation pulse\n', 1e6*(t_refocusing(i_ref)-t_excitation(i_ex))*2);
fprintf('Echo column is sampled at %g us after the excitation pulse\n', 1e6*(t_echo-t_excitation(i_ex)));

fprintf('TR per slice %g ms\n',(t_excitation(2)-t_excitation(1))*1e3);

return

%% apply soft delay values to see how the sequence timing is changing

seq.applySoftDelay('TE',30e-3);
seq.plot('showBlocks', 1, 'timeRange', TR*(Ndummy+[0 1]), 'timeDisp', 'us', 'stacked', 1) ;

seq.applySoftDelay('TR',2000e-3);
seq.plot('showBlocks', 1, 'timeRange', TR*(Ndummy+[0 1]), 'timeDisp', 'us', 'stacked', 1) ;

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slew rate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 
