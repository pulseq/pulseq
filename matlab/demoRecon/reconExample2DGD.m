% very basic and crude non-Cartesian recon using griddata()
%
% needs mapVBVD in the path

%% Load the latest file from a dir
path='../IceNIH_RawSend/'; % directory to be scanned for data files
%path='~/Dropbox/shared/data/siemens/';
pattern='*.dat';

D=dir([path pattern]);
[~,I]=sort([D(:).datenum]);
data_file_path=[path D(I(end)).name]; % use end-1 to reconstruct the second-last data set, etc.

%%
twix_obj = mapVBVD(data_file_path);

%% Load sequence from file (optional)

seq_file_path = [data_file_path(1:end-3) 'seq'];

traj_recon_delay=0e-6;%1.75e-6; % adjust this parameter to potentially improve resolution & geometric accuracy. It can be calibrated by inverting the spiral revolution dimension and making two images match. for our Prisma and a particular trajectory we found 1.75e-6

seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');
[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace('trajectory_delay', traj_recon_delay);
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal');

%% Define FOV and resolution and simple off-resonance frequency correction 

fov=256e-3; Nx=128; Ny=Nx; 
deltak=1/fov;
os=2; % oversampling factor (we oversample both in image and k-space)
offresonance=0; % global off-resonance in Hz

%%
if iscell(twix_obj)
    data_unsorted = double(twix_obj{end}.image.unsorted());
else
    data_unsorted = double(twix_obj.image.unsorted());
end
rawdata = permute(data_unsorted, [1,3,2]);
rawdata = reshape(rawdata, [size(rawdata,1)*size(rawdata,2),size(rawdata,3)]);
channels=size(rawdata,2);

for c=1:channels
    rawdata(:,c) = rawdata(:,c) .* exp(-1i*2*pi*t_adc*offresonance);
end

%% here we expect Nx, Ny, deltak to be set already
% and rawdata ktraj_adc loaded (and having the same dimensions)

kxm=round(os*os*Nx/2);
kym=round(os*os*Ny/2);

[kyy,kxx] = meshgrid(-kxm:(kxm-1), -kym:(kym-1));
kyy=-kyy*deltak/os;
kxx=kxx*deltak/os;

kgd=zeros([size(kxx) channels]);
for c=1:channels
    kgd(:,:,c)=griddata(ktraj_adc(1,:),ktraj_adc(2,:),rawdata(:,c),kxx,kyy,'cubic'); % we swap the order ind invert one sign to account for Matlab's strange column/line convention
end
kgd(isnan(kgd))=0;

figure;imagesc(log(abs(kgd(:,:,1))));axis('square');

igd=ifftshift(ifft2(ifftshift(kgd())));

Nxo=round(Nx*os);
Nyo=round(Ny*os);
Nxs=round((size(igd,1)-Nxo)/2);
Nys=round((size(igd,2)-Nyo)/2);
igdc = igd((Nxs+1):(Nxs+Nxo),(Nys+1):(Nys+Nyo),:);
figure;imab(abs(igdc));colormap('gray');
%axis('equal');

%% Sum of squares combination
sos=abs(sum(igdc.^2,ndims(igdc)).^(1/2));
sos=sos./max(sos(:));
figure;imab(sos);colormap('gray');
%imwrite(sos, ['img_combined.png'])

