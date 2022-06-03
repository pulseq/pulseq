% this is a demo low-performance EPI sequence;
% it doesn't use ramp-samping and is only good for educational purposes.
% in addition, it demonstrated how LABEL extension can be used to set data
% header values, which can be used either in combination with integrated
% image reconstruction or to guide the off-line reconstruction tools
% 
seq=mr.Sequence();              % Create a new sequence object
fov=220e-3; Nx=96; Ny=Nx;       % Define FOV and resolution
thickness=3e-3;                 % slice thinckness
Nslices=7;
Nreps = 4;
Navigator = 3;

% Set system limits
lims = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
               'MaxSlew',130,'SlewUnit','T/m/s', ...
               'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6);


% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',lims,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

% define the trigger to play out
trig=mr.makeTrigger('physio1','duration', 2000e-6); % duration after

% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
dwellTime = 4e-6; % I want it to be divisible by 2
readoutTime = Nx*dwellTime;
flatTime=ceil(readoutTime*1e5)*1e-5; % round-up to the gradient raster
gx = mr.makeTrapezoid('x',lims,'Amplitude',kWidth/readoutTime,'FlatTime',flatTime);
adc = mr.makeAdc(Nx,'Duration',readoutTime,'Delay',gx.riseTime+flatTime/2-(readoutTime-dwellTime)/2);

% Pre-phasing gradients
preTime=8e-4;
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2,'Duration',preTime); % removed -deltak/2 to aligh the echo between the samples
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',lims,'Area',+Ny/2*deltak,'Duration',preTime); % phase area should be Kmax for clin=0 and -Kmax for clin=Ny... strange

% Phase blip in shortest possible time
dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6)*10e-6;
gy = mr.makeTrapezoid('y',lims,'Area',-deltak,'Duration',dur); % phase area should be Kmax for clin=0 and -Kmax for clin=Ny... strange

gz_spoil=mr.makeTrapezoid('z',lims,'Area',deltak*Nx*4);

% Define sequence blocks
for r=1:Nreps
    seq.addBlock(trig, mr.makeLabel('SET','SLC', 0)); 
    for s=1:Nslices
        rf.freqOffset=gz.amplitude*thickness*(s-1-(Nslices-1)/2);
        rf.phaseOffset=-rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
        seq.addBlock(rf,gz);
        seq.addBlock(gxPre,gzReph, ...
                     mr.makeLabel('SET','NAV',1),...
                     mr.makeLabel('SET','LIN', round(Ny/2)));
        for n=1:Navigator
            seq.addBlock(gx,adc, ...
                         mr.makeLabel('SET','REV', gx.amplitude<0), ...
                         mr.makeLabel('SET','SEG', gx.amplitude<0), ...
                         mr.makeLabel('SET','AVG',n==3));
            if (n~=Navigator)
                seq.addBlock(mr.makeDelay(mr.calcDuration(gy))); % we need this dummy blip pulse to maintain ientical RO gradient timing anf the correspnding eddy currents  
            end
            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
        end
        %seq.addBlock(gxPre,gyPre,gzReph);
        seq.addBlock(gyPre, ...
                     mr.makeLabel('SET','LIN', 0), ...
                     mr.makeLabel('SET','NAV', 0), ...
                     mr.makeLabel('SET','AVG', 0) );% lin/nav/avg reset
        for i=1:Ny
            seq.addBlock(mr.makeLabel('SET','REV', gx.amplitude<0), ...
                         mr.makeLabel('SET','SEG', gx.amplitude<0));
            seq.addBlock(gx,adc);   % Read one line of k-space
            seq.addBlock(gy,mr.makeLabel('INC','LIN', 1));  % Phase blip
            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
        end
        seq.addBlock(gz_spoil,mr.makeDelay(0.1),mr.makeLabel('INC','SLC', 1)); 
        if rem(Navigator+Ny,2)~=0
            gx.amplitude = -gx.amplitude;
        end
    end
    seq.addBlock(mr.makeLabel('INC','REP', 1)); % zero-duration block, but it's OK
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
seq.setDefinition('Name', 'epi_lbl');

seq.write('epi_label.seq');   % Output sequence for scanner

return

%% plots, etc.
seq.plot('TimeRange',[0 0.1], 'TimeDisp', 'ms', 'Label', 'SEG,LIN,SLC');
% trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.');

% seq.sound(); % simulate the seq's tone
% test read-write
%seq.read('epi_label.seq');
%seq.write('epi_label2.seq');
