% this is an experimentaal high-performance EPI sequence
% which uses split gradients to overlap blips with the readout
% gradients combined with ramp-samping

seq=mr.Sequence();         % Create a new sequence object
fov=220e-3; Nx=64; Ny=64;  % Define FOV and resolution
thickness=3e-3;            % slice thinckness
Nslices=16;
traj_recon_delay=-7.5e-6; % adjust this parameter to supress ghosting (negative allowed)
                           % -1.75e-6 makes the plot look pretty, on our
                           % TRIO -7.5e-6 results in a reasonable ghost supression
pe_enable=1;               % a flag to quickly disable phase encoding (1/0)

% Set system limits
lims = mr.opts('MaxGrad',32,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6);  

% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',lims,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
readoutTime = 4.2e-4;

% Phase blip in shortest possible time
blip_dur = ceil(2*sqrt(deltak/lims.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
gy = mr.makeTrapezoid('y',lims,'Area',deltak,'Duration',blip_dur);

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area 
extra_area=blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
gx = mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx.amplitude=gx.amplitude/actual_area*kWidth;
gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
gx.flatArea = gx.amplitude*gx.flatTime;

% calculate ADC
% round-down dwell time to 10 ns
adcTime = floor((gx.riseTime+gx.flatTime+gx.fallTime-blip_dur)/Nx*1e8)*1e-8*Nx;
adc = mr.makeAdc(Nx,'Duration',adcTime,'Delay',blip_dur/2);

% split the blip into two halves and produnce a combined synthetic gradient
gy_parts = mr.splitGradient(gy, lims);
gy_blipup=gy_parts(1);
gy_blipdown=gy_parts(3);
gy_blipup.delay=gx.riseTime+gx.flatTime+gx.fallTime-blip_dur/2;
gy_blipdown.delay=0;
gy_blipdownup=mr.addGradients([gy_blipdown, gy_blipup], lims);

% pe_enable support
gy_blipup.waveform=gy_blipup.waveform*pe_enable;
gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;

% Pre-phasing gradients
preTime=8e-4;
gxPre = mr.makeTrapezoid('x',lims,'Area',-gx.area/2-deltak/2,'Duration',preTime);
gzReph = mr.makeTrapezoid('z',lims,'Area',-gz.area/2,'Duration',preTime);
gyPre = mr.makeTrapezoid('y',lims,'Area',-Ny/2*deltak*pe_enable,'Duration',preTime);

% Define sequence blocks
for s=1:Nslices
    rf.freqOffset=gz.amplitude*thickness*(s-1-(Nslices-1)/2);
    seq.addBlock(rf,gz);
    seq.addBlock(gxPre,gyPre,gzReph);
    for i=1:Ny
        if i==1
            seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
        else if i==Ny
            seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
        else
            seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
        end; end;
        gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
    end
    % generate waveforms and k-space trajectory for the 1st slice
    if s==1
        grad_waveforms=seq.gradient_waveforms();
        % detect excitation
        [a,i]=max(abs(rf.signal));
        t_start=rf.delay+rf.t(i);
        i_start=round(t_start/lims.gradRasterTime);
        % zerofill waveforms prior to the excitation 
        % (we prefer this to trimming because of the common time axis)
        grad_waveforms(:,1:(i_start))=0; 
        % kind of a 1/2 pixel shift ...
        grad_waveforms(:,i_start+1)=0.5*grad_waveforms(:,i_start+1); 
        
        grad_wavelen=size(grad_waveforms,2);
        % calculate gradient integrals
        grad_moments=cumsum(grad_waveforms,2)*lims.gradRasterTime;        
        % compute ADC sampling point times
        ktime=[];
        t0=0;
        for i=1:length(seq.blockEvents)
            block=seq.getBlock(i);
            if ~isempty(block.adc)
                ktime=[ktime, ((0:(block.adc.numSamples-1))*block.adc.dwell + block.adc.delay + t0 + traj_recon_delay)];
            end
            t0=t0+mr.calcDuration(block);
        end
        % sample the moments at the ADC time points
        traj_mat=interp1((0:(grad_wavelen-1))*lims.gradRasterTime, grad_moments', ktime)';
    end
end

seq.write('epi_rs.seq');   % Output sequence for scanner

% write the k-space trajectory (now 3D and wothout reconSize)
save('epi_rs_traj.mat','traj_mat','ktime');%,'reconSize');

seq.plot();             % Plot sequence waveforms
figure; plot(traj_mat'); % plot k-space trajectory
figure; plot(traj_mat(1,:),traj_mat(2,:),'b'); % a better plot
hold;plot(traj_mat(1,:),traj_mat(2,:),'r.');

% seq.install('siemens');
