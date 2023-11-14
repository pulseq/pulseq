% this is a demo GRE sequence, which uses LABEL extension to produce raw
% data reconstuctable by the integrated image reconstruction on the scanner

% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq=mr.Sequence(sys);         % Create a new sequence object

fov=256e-3; Nx=64; Ny=Nx; % Define FOV and resolution
alpha=10;                  % flip angle
thickness=3e-3;            % slice
Nslices=1;
sliceGap = 1e-3 ;
TR=30e-3; 
TE=4.3e-3;
%TE=[7.38 9.84 15]*1e-3;            % alternatively give a vector here to have multiple TEs (e.g. for field mapping)

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
roDuration=3.2e-3;              % ADC duration

% Create fat-sat pulse 
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
% B0=2.89; % 1.5 2.89 3.0
% sat_ppm=-3.45;
% sat_freq=sat_ppm*1e-6*B0*lims.gamma;
% rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims,'Duration',8e-3,...
%     'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
% gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',1e-3,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',1e-3,'system',sys);
phaseAreas = -((0:Ny-1)-Ny/2)*deltak; % phase area should be Kmax for clin=0 and -Kmax for clin=Ny... strange

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/thickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTE>=0));
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

rf_phase=0;
rf_inc=0;

% set ACS lines for GRAPPA simulation (fully sampled central k-space
% region)
accelFactorPE = 2 ;
ACShw = 12 ; % GRAPPA ACS half width i.e. here 28 lines are ACS
PEsamp_u = 1:accelFactorPE:Ny ; % undersampling by every alternate line. Phase encoding in Y direction
PEsamp_ACS = Ny/2-ACShw+1 : Ny/2+ACShw ; % GRAPPA autocalibration lines
PEsamp = union(PEsamp_u, PEsamp_ACS) ; % actually sampled lines
nPEsamp = length(PEsamp) ; % number of actually sampled
PEsamp_INC = diff([PEsamp, PEsamp(end)]) ;
% PEsamp indexes the actually sampled lines to the encoded k-space line number. 
% For example, if there were just regular factor 2 undersampling 
% (with no ACS lines), PEsamp would have length 128 and be [1 3 5 ... 255].
% With ACS lines, the elements of PEsamp are separated by 2 near the k-space
% edges, and by 1 in the central ACS region.

% all LABELS / counters an flags are automatically initialized to 0 in the beginning, no need to define initial 0's  
% so we will just increment LIN after the ADC event (e.g. during the spoiler)

% Set PAT scan flag
% Mdh.setPATRefScan; Mdh.setPATRefAndImaScan
lblSetRefScan = mr.makeLabel('SET','REF', true) ;
lblSetRefAndImaScan = mr.makeLabel('SET','IMA', true) ;
lblResetRefScan = mr.makeLabel('SET','REF', false) ;
lblResetRefAndImaScan = mr.makeLabel('SET','IMA', false) ;

%lblSetRefScan.id=seq.registerLabelEvent(lblSetRefScan);
%lblSetRefAndImaScan.id=seq.registerLabelEvent(lblSetRefAndImaScan);
%lblResetRefScan.id=seq.registerLabelEvent(lblResetRefScan);
%lblResetRefAndImaScan.id=seq.registerLabelEvent(lblResetRefAndImaScan);

% Add noise scans.
seq.addBlock(mr.makeLabel('SET', 'LIN', 0)) ;
seq.addBlock(adc, mr.makeLabel('SET', 'NOISE', true),lblResetRefScan,lblResetRefAndImaScan) ;
seq.addBlock(mr.makeLabel('SET', 'NOISE', false)) ;

if ismember(PEsamp(1),PEsamp_ACS)
    if ismember(PEsamp(1),PEsamp_u)
        seq.addBlock(mr.makeLabel('SET', 'LIN', PEsamp(1)-1), lblSetRefAndImaScan, lblSetRefScan) ;
    else
        seq.addBlock(mr.makeLabel('SET', 'LIN', PEsamp(1)-1), lblResetRefAndImaScan, lblSetRefScan) ;
    end
else
    seq.addBlock(mr.makeLabel('SET', 'LIN', PEsamp(1)-1), lblResetRefAndImaScan, lblResetRefScan) ; % Set first PE PAT scan line
end

% slice positions
slicePositions=(thickness+sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2);
slicePositions=slicePositions([1:2:Nslices 2:2:Nslices]); % reorder slices for an interleaved acquisition (optional)

rf.freqOffset=gz.amplitude*thickness*(1-1-(Nslices-1)/2);
% loop over phase encodes and define sequence blocks
for count = 1:nPEsamp
    if ismember(PEsamp(count),PEsamp_ACS)
        if ismember(PEsamp(count),PEsamp_u)
            seq.addBlock(lblSetRefAndImaScan, lblSetRefScan) ;
        else
            seq.addBlock(lblResetRefAndImaScan, lblSetRefScan) ;
        end
    else
        seq.addBlock(lblResetRefAndImaScan, lblResetRefScan) ;
    end
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);
    %
    seq.addBlock(rf,gz);
    gyPre = mr.makeTrapezoid('y','Area',phaseAreas(PEsamp(count)),'Duration',mr.calcDuration(gxPre),'system',sys);
    seq.addBlock(gxPre,gyPre,gzReph);
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(gx,adc);
    gyPre.amplitude=-gyPre.amplitude;

    seq.addBlock(mr.makeDelay(delayTR),gxSpoil,gyPre,gzSpoil);
    seq.addBlock(mr.makeLabel('INC','LIN', PEsamp_INC(count)));
end
seq.addBlock(mr.makeLabel('SET', 'LIN', 0)) ;

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
seq.setDefinition('FOV', [fov fov max(slicePositions)-min(slicePositions)+thickness]);
seq.setDefinition('Name', 'gre_gt');
% the following definitions have effect in conjunction with LABELs 
seq.setDefinition('SlicePositions', slicePositions);
seq.setDefinition('SliceThickness', thickness);
seq.setDefinition('SliceGap', sliceGap);

seq.write('gre_gt.seq')       % Write to pulseq file

%seq.install('siemens');
%return

%% evaluate label settings more specifically
%seq.plot('timeRange', [0 32]*TRout, 'TimeDisp', 'ms', 'Label', 'LIN');
adc_lbl=seq.evalLabels('evolution','adc');
figure; plot(adc_lbl.REF);
hold on; plot(adc_lbl.LIN);plot(adc_lbl.IMA) ; plot(adc_lbl.NOISE);
plot(adc_lbl.IMA) ;
legend('REF','LIN', 'IMA', 'NOISE');
title('evolution of labels/counters');

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 256]*TR, 'TimeDisp', 'ms', 'label', 'lin');

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
title('k-space components as functions of time');
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('2D k-space');

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
return ;
rep = seq.testReport;
fprintf([rep{:}]);

