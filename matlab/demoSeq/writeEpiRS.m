% this is an experimentaal high-performance EPI sequence
% which uses split gradients to overlap blips with the readout
% gradients combined with ramp-samping

% Set system limits
sys = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6,...
    'B0', 2.89 ... % this is Siemens' 3T
);  

seq=mr.Sequence(sys);      % Create a new sequence object
fov=256e-3; Nx=64; Ny=Nx;  % Define FOV and resolution
thickness=4e-3;            % slice thinckness
Nslices=5;

pe_enable=1;               % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
ro_os=1;                   % oversampling factor (in contrast to the product sequence we don't really need it)
readoutTime=4.2e-4;        % this controls the readout bandwidth
partFourierFactor=1;       % partial Fourier factor: 1: full sampling 0: start with ky=0

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
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx.amplitude=gx.amplitude/actual_area*kWidth;
gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
gx.flatArea = gx.amplitude*gx.flatTime;

% calculate ADC
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively. 
adcDwellNyquist=deltak/gx.amplitude/ro_os;
% round-down dwell time to 100 ns
adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
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

% Define sequence blocks
%seq.addBlock(mr.makeDelay(1)); % older scanners like Trio may need this
                                % dummy delay to keep up with timing
for s=1:Nslices
    seq.addBlock(rf_fs,gz_fs);
    rf.freqOffset=gz.amplitude*thickness*(s-1-(Nslices-1)/2);
    rf.phaseOffset=-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
    seq.addBlock(rf,gz,trig);
    seq.addBlock(gxPre,gyPre,gzReph);
    for i=1:Ny_meas
        if i==1
            seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
        elseif i==Ny_meas
            seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
        else
            seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
        end 
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
    end
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

%% do some visualizations

seq.plot();             % Plot all sequence waveforms

seq.plot('timeDisp','us','showBlocks',1,'timeRange',[0 25e-3]); %detailed view

[rf_bw,rf_spectrum,rf_w]=mr.calcRfBandwidth(rf);
figure;plot(rf_w,abs(rf_spectrum));
title('Excitation pulse profile (low-angle approximation)');
xlabel('Frequency, Hz');
xlim(3*[-rf_bw rf_bw]);

% trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
title('k-space vector components as functions of time');

figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
axis('equal'); % enforce aspect ratio for the correct trajectory display
title('k-space trajectory (k_x/k_y)');

figure; plot(t_slicepos, slicepos, '*');
title('slice position (vector components) as a function or time');
%axis off;

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

%% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov thickness]);
seq.setDefinition('Name', 'epi');

seq.write('epi_rs.seq'); 

% seq.install('siemens');

% seq.sound(); % simulate the seq's tone

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


