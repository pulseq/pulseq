% very basic and crude recon for single-shot EPI with ramp-sampling 
%
% needs mapVBVD in the path

%% Load the latest file from a dir
path='../IceNIH_RawSend/'; % directory to be scanned for data files
%path='~/Dropbox/shared/data/siemens/';
%path='~/Dropbox/shared/data/siemens/demo_epi/';

pattern='/*.dat';

D=dir([path pattern]);
[~,I]=sort([D(:).datenum]);
data_file_path=[path D(I(end)).name]; % use end-1 to reconstruct the second-last data set, etc.
%%
twix_obj = mapVBVD(data_file_path);

%% Load sequence from file (optional)

seq_file_path = [data_file_path(1:end-3) 'seq'];

seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%[ktraj_adc, ktraj, t_excitation, t_refocusing] = seq.calculateKspace();
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal');

%% define raw data

if iscell(twix_obj)
    rawdata = double(twix_obj{end}.image.unsorted());
else
    rawdata = double(twix_obj.image.unsorted());
end

%% if necessary re-tune the trajectory delay to supress ghosting
traj_recon_delay=0e-6;%3.23e-6;%-1e-6;%3.90e-6;%-1.03e-6; % adjust this parameter to supress ghosting (negative allowed) (our trio -1.0e-6, prisma +3.9e-6; avanto +3.88)
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace('trajectory_delay', traj_recon_delay);
%ktraj_adc_nodelay=seq.calculateKspace('trajectory_delay', 10e-6);

%% automatic detection of the measurement parameters (FOV, matrix size, etc)
nADC = size(rawdata, 1);
k_last=ktraj_adc(:,end);
k_2last=ktraj_adc(:,end-nADC);
delta_ky=k_last(2)-k_2last(2);
fov=1/abs(delta_ky);
Ny_post=round(abs(k_last(2)/delta_ky));
if k_last(2)>0
    Ny_pre=round(abs(min(ktraj_adc(2,:))/delta_ky));
else
    Ny_pre=round(abs(max(ktraj_adc(2,:))/delta_ky));
end
Nx=2*max([Ny_post,Ny_pre]);
Ny=Nx;
Ny_sampled=Ny_pre+Ny_post+1;
% for Ny_sampled=1:Ny 
%     if abs((ktraj_adc(2,end-nADC*(Ny_sampled-1))-ktraj_adc(2,end-nADC*Ny_sampled))-delta_ky)>1e-6
%         break;
%     end
% end

%% manually (re-)define FOV and resolution
%fov=256e-3; 
%Nx=64; Ny=Nx; Ny_sampled=Ny; 

%DWI example:
%Nx=112; Ny=Nx; Ny_sampled=98; fov=224e-3; 

%SE EPI RS example:
%Nx=64; Ny=Nx; Ny_sampled=56; fov=250e-3; 

%% classical phase correction / trajectory delay calculation 
%  here we assume we are dealing with the calibration data
data_odd=ifftshift(ifft(ifftshift(rawdata(:,:,1:2:end),1)),1);
data_even=ifftshift(ifft(ifftshift(rawdata(end:-1:1,:,2:2:end),1)),1);
cmplx_diff=data_even.*conj(data_odd);
cmplx_slope=cmplx_diff(2:end,:,:).*conj(cmplx_diff(1:end-1,:,:));
mslope_phs=angle(sum(cmplx_slope(:)));
dwell_time=(t_adc(nADC)-t_adc(1))/(nADC-1);
measured_traj_delay=mslope_phs/2/2/pi*nADC*dwell_time;
fprintf('measured trajectory delay (assuming it is a calibration data set) is %g s\n', measured_traj_delay);
fprintf('type this value in the section above and re-run the script\n');
% we do not calculate the constant phase term here because it depends on
% the definitions of the center of k-space and image-space 

%% analyze the trajecotory, resample the data
% here we expect rawdata ktraj_adc loaded (and having the same dimensions)
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);

kxmin=min(ktraj_adc(1,:));
kxmax=max(ktraj_adc(1,:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions
ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);

data_resampled=zeros(length(kxx), nCoils, nAcq);
ktraj_resampled=zeros(nD, length(kxx), nAcq);
t_adc_resampled=zeros(length(kxx), nAcq);
for a=1:nAcq
    for c=1:nCoils
        data_resampled(:,c,a)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a),kxx,'spline',0);
    end
    ktraj_resampled(1,:,a)=kxx;
    for d=2:nD
        ktraj_resampled(d,:,a)=interp1(ktraj_adc2(1,:,a),ktraj_adc2(d,:,a),kxx,'linear',NaN);
    end
    t_adc_resampled(:,a)=interp1(ktraj_adc2(1,:,a),t_adc2(:,a),kxx,'linear',NaN);
end

figure;imagesc(squeeze(abs(data_resampled(:,1,:)))');axis('square');

%% in some cases because of the incorrectly calculated trajectory phase correction may be needed
%  one such case is the use of the frequency shift proportional to gradient
%  in combination with the gradient delay and FOV offset in the RO direction
%  this calculation is best done with the calibration data, but also seems
%  to work with the actual image data

% here we assume we are dealing with the calibration data
data_odd=ifftshift(ifft(ifftshift(data_resampled(:,:,1:2:end),1)),1);
data_even=ifftshift(ifft(ifftshift(data_resampled(:,:,2:2:end),1)),1);
cmplx_diff1=data_even.*conj(data_odd);
cmplx_diff2=data_even(:,:,1:end-1).*conj(data_odd(:,:,2:end));
mphase1=angle(sum(cmplx_diff1(:)));
mphase2=angle(sum(cmplx_diff2(:)));
mphase=angle(sum([cmplx_diff1(:); cmplx_diff2(:)]));

%%
pc_coef=0;
%pc_coef=mphase1/2/pi;

data_pc=data_resampled;
for c=1:nCoils
    for i=1:size(data_resampled,1)
        data_pc(i,c,:)=squeeze(data_resampled(i,1,:)).*exp(1i*2*pi*pc_coef*mod((1:size(data_pc,3))',2));
    end
end
figure;imagesc(squeeze(angle(data_pc(:,1,:)))');axis('square');

%% reshape for multiple slices or repetitions
n4 = nAcq / Ny_sampled;
data_pc = reshape(data_pc, [size(data_pc,1),nCoils,Ny_sampled,n4]);

%% display results

%figure;imagesc(squeeze(abs(data_resampled)));axis('square');

data_xky=ifftshift(ifft(ifftshift(data_pc,1)),1);

if Ny_sampled~=Ny
    data_xky1=zeros(size(data_pc,1),nCoils,Ny,n4);
    data_xky1(:,:,(1:Ny_sampled)+Ny-Ny_sampled,:)=data_xky;
    data_xky=data_xky1;
end

%figure;imagesc(abs(squeeze(data_xky(:,1,:,1)))');axis('square');colormap('gray');

% strictly speaking we have to do an intensity compensation here due to the
% convolution at the interp1() step, but for now we ignore it...

data_xy=ifftshift(ifft(ifftshift(data_xky,3),[],3),3);

%%
if delta_ky>0
    figure;imab(sqrt(squeeze(sum(abs(data_xy(:,:,end:-1:1,:).^2),2))));axis('square');colormap('gray'); 
else
    figure;imab(sqrt(squeeze(sum(abs(data_xy(:,:,:,:).^2),2))));axis('square');colormap('gray');
end
