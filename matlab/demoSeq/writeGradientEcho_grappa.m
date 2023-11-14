%% ISMRM virtual meeting. Nov. 16, 2023. Qingping Chen
% Implement GRAPPA acceleration to the 2D GRE sequence with optimized 
% spoiler and accelerated computation. Use LABEL extension to produce raw
% data reconstructable by the integrated ICE/Gadgetron reconstructio on the
% scanner
clear all; close all; clc ;

%% set system limits and parameters
sys = mr.opts('MaxGrad', 22, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6) ;

seq = mr.Sequence(sys) ;           % Create a new sequence object
fov = 256e-3 ; Nx = 256 ; Ny = 256 ;     % Define FOV and resolution
alpha = 10 ;                       % flip angle
sliceThickness = 3e-3 ;            % slice
sliceGap = 1e-3 ;
TR = 20e-3 ;                       % repetition time TR
TE = 6e-3 ;                        % echo time TE
roDuration = 5.12e-3 ;              % ADC duration

%% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

%% Define other gradients and ADC events
deltak = 1/fov ;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys) ;
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys) ;
gxPre = mr.makeTrapezoid('x','Area',-gx.amplitude*(adc.dwell*(adc.numSamples/2+0.5)+0.5*gx.riseTime),...
    'Duration',1e-3,'system',sys) ;
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',1e-3,'system',sys) ;
phaseAreas = ((0:Ny-1) - Ny/2) * deltak ;
gy = mr.makeTrapezoid('y','Area',max(abs(phaseAreas)),'Duration',mr.calcDuration(gxPre),'system',sys) ;
peScales = phaseAreas/gy.area ;
clear gyPre gyReph ;
for iY = 1:Ny
    gyPre(iY) = mr.scaleGrad(gy, peScales(iY)) ;
    gyReph(iY) = mr.scaleGrad(gy, -peScales(iY)) ;
end
%% gradient spoiling
% we cut the RO gradient into two parts for the optimal spoiler timing
[gx1, ~] = mr.splitGradientAt(gx, gx.riseTime + gx.flatTime) ;
% gradient spoiling
gxSpoil = mr.makeExtendedTrapezoidArea(gx.channel, gx.amplitude, 0, 2*Nx*deltak, sys) ;
gxSpoil.delay = mr.calcDuration(gx1) ;
gx_add = mr.addGradients({gx1, gxSpoil},'system',sys) ;

gzSpoil = mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys) ;
gzSpoil.delay = max(mr.calcDuration(gx1), gzSpoil.delay) ;

for iY = 1:Ny
    gyReph(iY).delay = max(mr.calcDuration(gx1), gyReph(iY).delay) ;
end

%% Calculate timing
delayTE = ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime) * seq.gradRasterTime ;
delayTR = ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx_add, gzSpoil, gyReph) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime ;
assert(delayTE >= 0 ) ;
assert(delayTR >= 0 ) ;

%% implement GRAPPA pattern
% set ACS lines for GRAPPA simulation (fully sampled central k-space
% region)
accelFactorPE = 2 ;
ACShw = 12 ; % GRAPPA ACS half width i.e. here 28 lines are ACS
PEsamp_u = 1:accelFactorPE:Ny ; % undersampling by every alternate line. Phase encoding in Y direction
PEsamp_ACS = Ny/2-ACShw+1 : Ny/2+ACShw ; % GRAPPA autocalibration lines
PEsamp = union(PEsamp_u, PEsamp_ACS) ; % actually sampled lines
nPEsamp = length(PEsamp) ; % number of actually sampled
PEsamp_INC = diff([PEsamp, PEsamp(end)]) ;

%% label setting
% Set PAT scan flag
% Mdh.setPATRefScan; Mdh.setPATRefAndImaScan
lblSetRefScan = mr.makeLabel('SET','REF', true) ;
lblSetRefAndImaScan = mr.makeLabel('SET','IMA', true) ;
lblResetRefScan = mr.makeLabel('SET','REF', false) ;
lblResetRefAndImaScan = mr.makeLabel('SET','IMA', false) ;

%% accelerate computations
% preregister constant objects to accelerate computations
% this is not necessary, but accelerates the sequence creation
gxPre.id = seq.registerGradEvent(gxPre) ;
gx1.id = seq.registerGradEvent(gx1) ;
gzSpoil.id = seq.registerGradEvent(gzSpoil) ;
% the phase of the RF object will change, therefore we only per-register the shapes
% [~, rf.shapeIDs] = seq.registerRfEvent(rf) ; 
rf.id = seq.registerRfEvent(rf) ; % NO GO EXAMPLE!!!
for iY = 1:Ny
    gyPre(iY).id = seq.registerGradEvent(gyPre(iY)) ;
    gyReph(iY).id = seq.registerGradEvent(gyReph(iY)) ;
end

lblSetRefScan.id = seq.registerLabelEvent(lblSetRefScan) ;
lblSetRefAndImaScan.id = seq.registerLabelEvent(lblSetRefAndImaScan) ;
lblResetRefScan.id = seq.registerLabelEvent(lblResetRefScan) ;
lblResetRefAndImaScan.id = seq.registerLabelEvent(lblResetRefAndImaScan) ;

%% Loop over phase encodes and define sequence blocks
tic ;
for count=1:nPEsamp
    % set GRAPPA labels
    if ismember(PEsamp(count),PEsamp_ACS)
        if ismember(PEsamp(count),PEsamp_u)
            seq.addBlock(lblSetRefAndImaScan, lblSetRefScan) ;
        else
            seq.addBlock(lblResetRefAndImaScan, lblSetRefScan) ;
        end
    else
        seq.addBlock(lblResetRefAndImaScan, lblResetRefScan) ;
    end
    % RF spoiling (vary RF phase pseudo-randomly)
    rand_phase = mod(117*(count^2 + count + 2), 360) * pi/180 ;
    rf.phaseOffset = rand_phase ;
    adc.phaseOffset = rand_phase ;
    %
    seq.addBlock(rf, gz) ;
    seq.addBlock(gxPre, gyPre(PEsamp(count)), gzReph) ;
    seq.addBlock(mr.makeDelay(delayTE)) ;
    seq.addBlock(gx_add, adc, gzSpoil, gyReph(PEsamp(count)) ) ;
    seq.addBlock(mr.makeDelay(delayTR) ) ;
    seq.addBlock(mr.makeLabel('INC','LIN', PEsamp_INC(count))) ;
end
toc ;

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
seq.setDefinition('FOV', [fov fov sliceThickness]) ;
seq.setDefinition('Name', 'gre_p2') ;
% the following definitions have effect in conjunction with LABELs 
seq.setDefinition('SlicePositions', 0) ;
seq.setDefinition('SliceThickness', sliceThickness) ;
seq.setDefinition('SliceGap', sliceGap) ;
seq.write('gre_p2.seq')       % Write to pulseq file

%% evaluate label settings more specifically
adc_lbl = seq.evalLabels('evolution','adc') ;
figure ; 
hold on ;
plot(adc_lbl.LIN) ;
plot(adc_lbl.REF) ;
plot(adc_lbl.IMA) ;

legend('LIN','REF', 'IMA') ;
title('evolution of labels/counters') ;

%% plot sequence and k-space diagrams
seq.plot('timeRange', [0 5]*TR) ;

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b') ; % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points



