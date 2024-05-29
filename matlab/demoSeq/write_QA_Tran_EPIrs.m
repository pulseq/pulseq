% Quality control for fMRI stability measures
% scan plane: straight transversal, TR = 2000 ms, TE = 30 ms
% this is an experimentaal high-performance EPI sequence
% which uses split gradients to overlap blips with the readout
% gradients combined with ramp-samping

% Set system limits
sys = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6,...
    'adcDeadTime', 10e-6, 'B0', 2.89 ... % this is Siemens' 3T
) ;  

seq = mr.Sequence(sys) ;      % Create a new sequence object
fov = 220e-3 ; Nx = 64 ; Ny = Nx ;  % Define FOV and resolution
thickness = 4e-3 ;            % slice thinckness in mm
sliceGap = 1e-3 ;             % slice gap im mm
Nslices = 27 ;
Nrep = 200 ;
TR = 2 ;
pe_enable = 1 ;               % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
ro_os = 2 ;                   % oversampling factor (in contrast to the product sequence we don't really need it)
% scanner = 'prisma' ;
%if (scanner(1) == 'c' || scanner(1) == 'C') % cimax
 %   readoutTime=580e-6 ;%seq4.2e-4;        % this controls the readout bandwidth
%elseif (scanner(1) == 'p' || scanner(1) == 'P') % prisma
%readoutTime=580e-6 ;
%elseif (scanner(1) == 't' || scanner(1) == 'T') % trio
%    readoutTime=580e-6 ;
%else
readoutTime = 520e-6;%770e-6 ; % default value
%    disp('no scanner selected, use default value of readout time!') ;
%end

readoutBW = 1/readoutTime ; % readout bandwidth
disp(['Readout bandwidth = ', num2str(readoutBW), ' Hz/Px']) ;
partFourierFactor = 1 ;       % partial Fourier factor: 1: full sampling 0: start with ky=0
Navigator = 3 ;		   % navigator echoes for ghost supprerssion

% Create fat-sat pulse 
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*sys.B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,'dwell',10e-6,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase    
gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',0.1/1e-4); % spoil up to 0.1mm
% Create 90 degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(pi/2,'system',sys,'Duration',2e-3,...
    'SliceThickness',thickness,'apodization',0.42,'timeBwProduct',4,'use','excitation');
gamma_H1 = 42.58 ; % [MHz/T]
rf_fs_peak = max(abs(rf_fs.signal))/gamma_H1 ; % [uT]
rf_peak = max(abs(rf.signal))/gamma_H1 ; % [uT]
disp(['The peak rf_fs amplitude = ', num2str(rf_fs_peak), ' uT']) ;
disp(['The peak rf amplitude = ', num2str(rf_peak), ' uT']) ;
% define the output trigger to play out with every slice excitatuion
trig=mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'

% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;

% Phase blip in shortest possible time
blip_dur = ceil(2*sqrt(deltak/sys.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
% the split code below fails if this really makes a trpezoid instead of a triangle...
gy = mr.makeTrapezoid('y',sys,'Area',-deltak,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center
%gy = mr.makeTrapezoid('y',lims,'amplitude',deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area 
extra_area=blip_dur/2*blip_dur/2*sys.maxSlew; % check unit!;
gx = mr.makeTrapezoid('x',sys,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
% actual sampled area = whole area - rampup deadarea - rampdown deadarea
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx.amplitude=gx.amplitude/actual_area*kWidth; % rescale amplitude to make sampled area = kWidth
gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2); % udpate parameters relative to amplitude
gx.flatArea = gx.amplitude*gx.flatTime;
ESP = 1e3 * mr.calcDuration(gx) ; % echo spacing, ms
disp(['echo spacing = ', num2str(ESP), ' ms']) ;
% calculate ADC
% % we use ramp sampling, so we have to calculate the dwell time and the
% % number of samples, which are will be qite different from Nx and
% % readoutTime/Nx, respectively. 
% adcDwellNyquist=deltak/gx.amplitude/ro_os;
% % round-down dwell time to 100 ns
% adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
% adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
% simple Siemens-like calculations with the fixed oversampling
assert(ro_os>=2) ;
adcSamples=Nx*ro_os ;
adcDwell=floor(readoutTime/adcSamples*1e7)*1e-7;
disp(['ADC bandwidth = ', num2str(1/adcDwell/1000), ' kHz']) ;
fprintf('Actual RO oversampling factor is %g, Siemens recommends it to be above 1.3\n', deltak/gx.amplitude/adcDwell)
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
disp(['ADC dwell time = ', num2str(adc.dwell*1e6), ' us']) ;
% realign the ADC with respect to the gradient
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us 
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

% FOV positioning requires alignment to grad. raster... -> TODO

% split the blip into two halves and produce a combined synthetic gradient
gy_parts = mr.splitGradientAt(gy, blip_dur/2, sys);
[gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, sys);

% pe_enable support
gy_blipup.waveform=gy_blipup.waveform*pe_enable;
gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;

% phase encoding and partial Fourier
Ny_pre=round(partFourierFactor*Ny/2-1); % PE steps prior to ky=0, excluding the central line
Ny_post=round(Ny/2+1); % PE lines after the k-space center including the central line
Ny_meas=Ny_pre+Ny_post;

% Pre-phasing gradients
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2);
gyPre = mr.makeTrapezoid('y',sys,'Area',Ny_pre*deltak);
[gxPre,gyPre,gzReph]=mr.align('right',gxPre,'left',gyPre,gzReph);
% relax the PE prepahser to reduce stimulation
gyPre = mr.makeTrapezoid('y',sys,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre,gzReph));
gyPre.amplitude=gyPre.amplitude*pe_enable;

% slice positions
slicePositions=(thickness+sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2);
slicePositions=slicePositions([1:2:Nslices 2:2:Nslices]); % reorder slices for an interleaved acquisition (optional)
%slicePositions=slicePositions([1:3:Nslices 2:3:Nslices 3:3:Nslices]); % reorder slices for an interleaved acquisition (optional)
% Define sequence blocks
TR_1slice = mr.calcDuration(gz_fs) + mr.calcDuration(gz) + mr.calcDuration(gzReph)+...
    Navigator*mr.calcDuration(gx) + mr.calcDuration(gyPre) + ...
    Ny_meas*mr.calcDuration(gx) ;
TRdelay = TR - TR_1slice * Nslices ;
TRdelay_perSlice = ceil(TRdelay / Nslices / sys.blockDurationRaster) * sys.blockDurationRaster ;
assert(TRdelay_perSlice>=0 ) ;

TE = rf.shape_dur/2 + rf.ringdownTime + mr.calcDuration(gzReph)+...
    Navigator*mr.calcDuration(gx) + mr.calcDuration(gyPre) + ...
    Ny_meas/2*mr.calcDuration(gx) - mr.calcDuration(gx)/2;
disp(['TR = ', num2str(TR), ' s', ', TE = ', num2str(1000*TE), ' ms']) ;
% if (scanner(1) == 't' || scanner(1) == 'T') % trio
%     seq.addBlock(mr.makeDelay(1)) ; % older scanners like Trio may need this % dummy delay to keep up with timing
% end

% change orientation to match the siemens product sequence
% reverse the polarity of all gradients in readout direction (Gx)
%gxPre.amplitude = -gxPre.amplitude ;
%gx.amplitude = -gx.amplitude ;
% reverse the polarity of all gradients in slice encoding direction (Gz)?
% gz_fs.amplitude = -gz_fs.amplitude ;
% gz.amplitude = -gz.amplitude ;
tic ;
seq.addBlock(mr.makeLabel('SET','REP', 0)) ;
for r=1:Nrep
    disp(['current repetition = ', num2str(r), '/', num2str(Nrep)]) ;
    seq.addBlock(mr.makeLabel('SET', 'SLC', 0) ) ;
    for s=1:Nslices
        seq.addBlock(rf_fs, gz_fs) ;
        rf.freqOffset=gz.amplitude*slicePositions(s);
        rf.phaseOffset=-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
        seq.addBlock(rf,gz,trig);
        if Navigator>0
            seq.addBlock(gxPre,gzReph, ...
                mr.makeLabel('SET','NAV',1),...
                mr.makeLabel('SET','LIN', floor(Ny/2))) ; % k-space center line
            for n=1:Navigator
                seq.addBlock(gx,adc, ...
                    mr.makeLabel('SET','REV', gx.amplitude>0), ...
                    mr.makeLabel('SET','SEG', gx.amplitude<0), ...
                    mr.makeLabel('SET','AVG', floor(n/2)));
                gx.amplitude = -gx.amplitude ;   % Reverse polarity of read gradient
            end
            seq.addBlock(gyPre,  ...
                mr.makeLabel('SET','LIN', -1), ...
                mr.makeLabel('SET','NAV', 0), ...
                mr.makeLabel('SET','AVG', 0) );% lin/nav/avg reset
        else
            seq.addBlock(gxPre,gyPre,gzReph,...
                mr.makeLabel('SET','LIN', -1), ...
                mr.makeLabel('SET','NAV', 0), ...
                mr.makeLabel('SET','AVG', 0) );% lin/nav/avg reset
        end

        for i=1:Ny_meas
            lrev=mr.makeLabel('SET','REV', gx.amplitude>0 ) ;
            lseg=mr.makeLabel('SET','SEG', gx.amplitude<0 ) ;
            llin=mr.makeLabel('INC','LIN', 1) ;
            if i==1 % first phase encoding step
                seq.addBlock(gx, gy_blipup, adc, lrev, lseg, llin) ; % Read the first line of k-space with a single half-blip at the end
            elseif i==Ny_meas % last phase encoding step
                seq.addBlock(gx, gy_blipdown, adc, lrev, lseg, llin) ; % Read the last line of k-space with a single half-blip at the beginning
            else % phase encoding steps in between
                seq.addBlock(gx, gy_blipdownup, adc, lrev, lseg, llin) ; % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
            end
    	    gx.amplitude = -gx.amplitude ;   % Reverse polarity of read gradient
        end
        seq.addBlock(mr.makeLabel('INC','SLC', 1)) ;
        if rem(Navigator+Ny_meas,2)~=0 % if the polarity of gx is not the same as original one
            gx.amplitude = -gx.amplitude;
        end
       seq.addBlock(mr.makeDelay(TRdelay_perSlice)) ;
    end
    seq.addBlock(mr.makeLabel('INC','REP', 1)) ;
end
toc ;
%% check whether the timing of the sequence is correct
if (Nrep<20)
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
end

% %% do some visualizations
% 
% seq.plot('Label', 'SEG,LIN,SLC');             % Plot all sequence waveforms
% 
% seq.plot('timeDisp','us','showBlocks',1,'timeRange',[0 25e-3]); %detailed view
% 
% rf.freqOffset=0;
% rf.phaseOffset=0;
% [rf_bw,rf_f0,rf_spectrum,rf_w]=mr.calcRfBandwidth(rf);
% figure;plot(rf_w,abs(rf_spectrum));
% title('Excitation pulse profile (low-angle approximation)');
% xlabel('Frequency, Hz');
% xlim(3*[-rf_bw rf_bw]);
% 
% %% trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
% 
% % plot k-spaces
% figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
% hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
% title('k-space vector components as functions of time');
% 
% figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
% hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
% axis('equal'); % enforce aspect ratio for the correct trajectory display
% title('k-space trajectory (k_x/k_y)');
% 
% figure; plot(t_slicepos, slicepos, '*');
% title('slice position (vector components) as a function or time');
% %axis off;

%% prepare the sequence output for the scanner
seq.setDefinition('Name', 'QA_epi'); 
seq.setDefinition('FOV', [fov fov max(slicePositions)-min(slicePositions)+thickness]);
% the following definitions have effect in conjunction with LABELs 
seq.setDefinition('SlicePositions', slicePositions);
seq.setDefinition('SliceThickness', thickness);
seq.setDefinition('SliceGap', sliceGap);
% readout without oversampling
%seq.setDefinition('IntendedReadoutSamples',Nx); % number of samples without oversamping
% readout gridding parameters
seq.setDefinition('ReceiverGainHigh',1);
seq.setDefinition('ReadoutOversamplingFactor',ro_os);
seq.setDefinition('TargetGriddedSamples',adc.numSamples); % number of samples after gridding (with oversamping)
seq.setDefinition('TrapezoidGriddingParameters', [gx.riseTime gx.flatTime gx.fallTime adc.delay-gx.delay adc.duration]); % rise,flat,fall,adc_delay,adc_dur

seq.write('QA_epi.seq'); 

% seq.install('siemens');

% seq.sound(); % simulate the seq's tone

%% evaluate label settings
adc_lbl=seq.evalLabels('evolution','adc');
figure; plot(adc_lbl.SLC);
hold on; plot(adc_lbl.LIN);
plot(adc_lbl.SEG); plot(adc_lbl.AVG);plot(adc_lbl.REP);
legend('slc','lin','seg','avg','rep');
title('evolution of labels/counters');
return;

%% another manual pretty plot option for gradients

lw=1;
%gw=seq.gradient_waveforms();
wave_data=seq.waveforms_and_times(true); % also export RF
gwm=max(abs([wave_data{1:3}]'));
rfm=max(abs([wave_data{4}]'));
ofs=2.05*gwm(2);

% plot "axes"
figure; 
axis_clr=[0.5,0.5,0.5];
plot([-0.01*gwm(1),1.01*gwm(1)],[0 0]*ofs,'Color',axis_clr,'LineWidth',lw/5); hold on; 
plot([-0.01*gwm(1),1.01*gwm(1)],[1 1]*ofs,'Color',axis_clr,'LineWidth',lw/5);
plot([-0.01*gwm(1),1.01*gwm(1)],[2 2]*ofs,'Color',axis_clr,'LineWidth',lw/5);
plot([-0.01*gwm(1),1.01*gwm(1)],[3 3]*ofs,'Color',axis_clr,'LineWidth',lw/5);

% plot the RF waveform
plot(wave_data{4}(1,:), abs(wave_data{4}(2,:))/rfm(2)*gwm(2)*0.75+3*ofs,'k','LineWidth',lw); 

% plot the entire gradient waveforms
plot(wave_data{3}(1,:), wave_data{3}(2,:)+2*ofs,'Color',[0,0.5,0.3],'LineWidth',lw); 
plot(wave_data{2}(1,:), wave_data{2}(2,:)+1*ofs,'r','LineWidth',lw);
plot(wave_data{1}(1,:), wave_data{1}(2,:),'b','LineWidth',lw);
t_adc_gr=t_adc+0.5*seq.gradRasterTime; % we have to shift the time axis because it is otherwise adpted to the k-space, which is a one-sided integration of the trajectory
gwr_adc=interp1(wave_data{1}(1,:), wave_data{1}(2,:),t_adc_gr);
plot(t_adc_gr,gwr_adc,'b.','MarkerSize',5*lw); % and sampling points on the kx-axis

xlim([-0.03*gwm(1),1.03*gwm(1)]);

set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slew rate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 

return
%% create a smoothly rotating 3D k-space plot
[kfa,ta,kf]=seq.calculateKspacePP();

figure;plot3(kf(1,:),-kf(3,:),kf(2,:));
hold on;plot3(kfa(1,:),-kfa(3,:),kfa(2,:),'r.');
set(gca,'visible','off'); % hide axes
set(gca, 'CameraViewAngle',get(gca, 'CameraViewAngle')); % freeze the view
kabsmax=max(abs(kf)')';
kxyabsmax=max(kabsmax([1 3]));
kxyzabsmax=max(kabsmax);
%axis([-kxyabsmax kxyabsmax -kxyabsmax kxyabsmax -kabsmax(2) kabsmax(2)])
axis(0.01*[-kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax -kxyzabsmax kxyzabsmax])
%s1=1.2;    
%axis([ -kabsmax(2)*s1 kabsmax(2)*s1  -kabsmax(2)*s1 kabsmax(2)*s1 min(kf(1,:)) kabsmax(1)]);
[caz,cel] = view;
folder='kspace3d';
mkdir(folder);
for caz_add=0:5:359 
    view(caz+caz_add,cel);
    drawnow;
    %print( '-r100', '-dpng', [folder '/frame_' num2str(caz_add,'%03d') '.png']);
    % use convert frame_???.png -gravity center -crop 300x300+0+0 +repage -delay 0.1 -loop 0 kspace_gre3d.gif
    % to create a GIF movie
end

