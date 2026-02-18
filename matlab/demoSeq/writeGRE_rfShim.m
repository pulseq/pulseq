% this is a demo GRE sequence, which uses LABEL extension to produce raw
% data reconstuctable by the integrated image reconstruction on the scanner

% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq=mr.Sequence(sys);         % Create a new sequence object

fov=224e-3; Nx=256; Ny=Nx; % Define FOV and resolution
alpha=15;                  % flip angle
thickness=5e-3;            % slice
Nslices=1;
Nrep=4;
TR=20e-3; 
TE=4.3e-3;
%TE=[4.38 6.84 12]*1e-3;            % alternatively give a vector here to have multiple TEs (e.g. for field mapping)

% more in-depth parameters
% RF spoiling increment = 84° for smoother transient decay, https://doi.org/10.1002/mrm.1910350216, 169° for diffusion independent rf spoiling in steady-state https://doi.org/10.1371/journal.pone.0324455
rfSpoilingInc=84;              % RF spoiling increment
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
    'SliceThickness',thickness,'apodization',0.42,'timeBwProduct',4,'system',sys,...
    'use','excitation');

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

% accelerate sequence calculation for objects that do not change un the loops
gz.id=seq.registerGradEvent(gz);
gxPre.id=seq.registerGradEvent(gxPre);
gzReph.id=seq.registerGradEvent(gzReph);
gx.id=seq.registerGradEvent(gx);
gxSpoil.id=seq.registerGradEvent(gxSpoil);
gzSpoil.id=seq.registerGradEvent(gzSpoil);

% RF spoiling
rf_phase=0;
rf_inc=0;

% all LABELS / counters an flags are automatically initialized to 0 in the beginning, no need to define initial 0's  
% so we will just increment LIN after the ADC event (e.g. during the spoiler)

seq.addBlock(mr.makeLabel('SET','REV', 1)); % left-right swap fix (needed for 1.4.0 and later)

% loop over repetitions
for r=1:Nrep
    seq.addBlock(mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','SLC', 0)); % needed to make it compatible to multiple REPs
    % loop over slices
    for s=1:Nslices
        rf.freqOffset=gz.amplitude*thickness*(s-1-(Nslices-1)/2);
        % loop over phase encodes and define sequence blocks
        for i=1:Ny
          for c=1:length(TE)
            rf.phaseOffset=rf_phase/180*pi;
            adc.phaseOffset=rf_phase/180*pi;
            rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase=mod(rf_phase+rf_inc, 360.0);
            %
            if r==1 
                seq.addBlock(rf,gz);
            else
                seq.addBlock(rf,gz,mr.makeRfShim([1, exp(1i*2*pi/(Nrep)*(r-1))]));
            end
            gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(gxPre),'system',sys);
            seq.addBlock(gxPre,gyPre,gzReph);
            seq.addBlock(delayTE(c));
            seq.addBlock(gx,adc);
            gyPre.amplitude=-gyPre.amplitude;
            spoilBlockContents={mr.makeDelay(delayTR(c)),gxSpoil,gyPre,gzSpoil}; % here we demonstrate the technique to combine variable counter-dependent content into the same block
            if c~=length(TE)
                spoilBlockContents=[spoilBlockContents {mr.makeLabel('INC','ECO', 1)}];
            else
                if length(TE)>1
                    spoilBlockContents=[spoilBlockContents {mr.makeLabel('SET','ECO', 0)}];
                end
                if i~=Ny
                    spoilBlockContents=[spoilBlockContents {mr.makeLabel('INC','LIN', 1)}];
                else
                    spoilBlockContents=[spoilBlockContents {mr.makeLabel('SET','LIN', 0), mr.makeLabel('INC','SLC', 1)}];
                end
            end
            seq.addBlock(spoilBlockContents{:});
          end
        end
    end
    seq.addBlock(mr.makeLabel('INC','REP', 1));
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

%% prepare sequence export
seq.setDefinition('FOV', [fov fov thickness*Nslices]);
seq.setDefinition('Name', 'gre_rfshim');
seq.setDefinition('TE', TE) ; % 3 TE values
seq.setDefinition('TR', TR) ;
seq.write('gre_rfshim.seq')       % Write to pulseq file

%seq.install('siemens');
%return

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 32]*TR, 'TimeDisp', 'ms', 'Label', 'LIN,SLC'); % looks like there is abug in the coloring of the multiple labels in the plot

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

%% evaluate label settings more specifically

lbls=seq.evalLabels('evolution','adc');
lbl_names=fieldnames(lbls);
figure; hold on;
for n=1:length(lbl_names)
    plot(lbls.(lbl_names{n}));
end
legend(lbl_names(:));
title('evolution of labels/counters/flags');
xlabel('adc number');

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport;
fprintf([rep{:}]);

