% very basic generic non-Cartesian recon using BART
% assumes raw data *.mat (or *.dat) fiels are stored together with the 
% *.seq Pulseq files withe the identical(!) names
%
% optionally needs mapVBVD in the path (if raw data have not yet been converted to *.mat)
% requires BART, modify or remove the line below depending on your path settings
bart_path='~/pulseq_home/bart/bart-0.6.00';
addpath([bart_path '/matlab']);
cur_dir=pwd;
cd(bart_path);
setenv('TOOLBOX_PATH', pwd);
cd(cur_dir);
%% Load the latest file from a dir
path='../IceNIH_RawSend/'; % directory to be scanned for data files
%path='/data/zte_petra/';
%path='~/Dropbox/shared/data/siemens/';
%path='~/20211025-AMR/data';

pattern='*.seq';

if path(end)~=filesep, path=[path filesep]; end

D=dir([path pattern]);
[~,I]=sort([D(:).datenum]);
seq_file_path = [path D(I(end-0)).name]; % use end-1 to reconstruct the second-last data set, etc.

%% alternatively just provide the path to the .seq file

%seq_file_path = '../interpreters/siemens/data_example/gre_example.seq'
%seq_file_path = [path '2020-11-10-073628.seq'];

%% Load sequence from the file with the same name as the raw data
fprintf(['loading `' seq_file_path '´ ...\n']);
seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');

%% keep basic filename without the extension
[p,n,e] = fileparts(seq_file_path);
basic_file_path=fullfile(p,n);

%% load the raw data file
data_file_path= [basic_file_path '.mat']; % try to load a matlab file with raw data...
try
    fprintf(['loading `' data_file_path '´ ...\n']);
    data_unsorted = load(data_file_path);
    if isstruct(data_unsorted)
        fn=fieldnames(data_unsorted);
        assert(length(fn)==1); % we only expect a single variable
        data_unsorted=data_unsorted.(fn{1});
    end
catch
    data_file_path= [basic_file_path '.dat']; % now try to load a raw data file...
    fprintf(['falling back to `' data_file_path '´ ...\n']);
    twix_obj = mapVBVD(data_file_path);
    if iscell(twix_obj)
        data_unsorted = twix_obj{end}.image.unsorted();
    else
        data_unsorted = twix_obj.image.unsorted();
    end
    seqHash_twix=twix_obj.hdr.Dicom.tSequenceVariant;
    if length(seqHash_twix)==32
        fprintf(['raw data contain pulseq-file signature ' seqHash_twix '\n']);
    end
    clear twix_obj
end

%% calculate k-space trajectory
traj_recon_delay=0*1e-6;%1.75e-6; % adjust this parameter to potentially improve resolution & geometric accuracy. It can be calibrated by inverting the spiral revolution dimension and making two images match. for our Prisma and a particular trajectory we found 1.75e-6
samples_to_mask=0; % each ADC may contain damaged samples in the begining
fprintf('calculating k-space trajectory ...');

%[ktraj_adc, t_adc, ktraj, t, t_excitation, t_refocusing]=seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace('trajectory_delay', traj_recon_delay);
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);

% optionally transform or patch the trajectory
ktraj_adc(2,:)=-ktraj_adc(2,:); % seems to be needed for the correct orientaton
%ktraj_adc(3,:)=0;
fprintf(' done\n');

figure; plot3(ktraj_adc(1,:),ktraj_adc(2,:),ktraj_adc(3,:),'r.');
axis('equal');
drawnow;

fov=seq.getDefinition('FOV');
if isempty(fov)
    fov=[256 256 256]; % default FOV 
    fprintf('WARNING: no FOV defined in the pulseq file, assuming FOV of %g x %g x %g m, matrix size will be wrong if these are incorrect\n', fov);
end

%fov=fov*1.5;
%% raw data preparation

adc_len=size(data_unsorted,1);
channels=size(data_unsorted,2);
readouts=size(data_unsorted,3);

%% average repetitions 
na=numel(data_unsorted)/channels/numel(t_adc); 
if (na>1)
    fprintf('averaging over %d acquisitions ...\n',na);
    data_unsorted=reshape(sum(reshape(data_unsorted,[adc_len,channels,readouts/na,na]),4),[adc_len,channels,readouts/na]);
    %data_unsorted=reshape(data_unsorted,[adc_len,channels,readouts/na,na]);
    %data_unsorted=reshape(data_unsorted(:,:,:,1),[adc_len,channels,readouts/na]);
    readouts=readouts/na;
end

%% mask data if requested
if samples_to_mask>0
    ktraj_adc=reshape(ktraj_adc,[3,adc_len,readouts]);
    ktraj_adc(:,1:samples_to_mask,:)=[];
    t_adc=reshape(t_adc,[adc_len,readouts]);
    t_adc(1:samples_to_mask,:)=[];
    data_unsorted(1:samples_to_mask,:,:)=[];
    adc_len=adc_len-samples_to_mask;
    ktraj_adc=reshape(ktraj_adc,[3,adc_len*readouts]);
    t_adc=reshape(t_adc,[1,adc_len*readouts]);
end

%% data preparation and petra/zte specific stuff
rawdata = permute(data_unsorted, [1,3,2]);

if strcmp('petra',seq.getDefinition('Name'))
    fprintf('performing special data handling for petra ...\n');
    SamplesPerShell=seq.getDefinition('SamplesPerShell');
    waspi_samples=1; % 1 16 32 tested; BART (nlinv) does not seem to like waspi (background noise goes up)
    if ~isempty(SamplesPerShell)
        SamplesSPI=sum(SamplesPerShell(2:end));
        data_mask=boolean([ones(SamplesPerShell(1),adc_len); ...
                          ones(SamplesSPI,waspi_samples), zeros(SamplesSPI,adc_len-waspi_samples)])';
    else
        data_mask=boolean(ones(adc_len,readouts));
    end    
    ktraj_adc=ktraj_adc(:,data_mask(:));

    rdtmp=reshape(rawdata, [adc_len*readouts,channels]);
    rawdata=rdtmp(data_mask(:),:);
    clear rdtmp;
    % reshape to 4D (1 samples 1 channels) as BART needs it
    rawdata=reshape(rawdata,[1 size(rawdata,1) 1 channels]);

    kyz_plane=abs(ktraj_adc(1,:))<0.9; % *(1/fov(1)) -- for BART-like trajectory not needed
    figure; plot(ktraj_adc(2,kyz_plane),ktraj_adc(3,kyz_plane),'r.','MarkerSize',0.5); axis equal;title('Kyz plane');
    %figure; plot3(kspace(kxz_plane,1),kspace(kxz_plane,2),kspace(kxz_plane,3),'r.'); axis equal;
    drawnow;
elseif strcmp('ute_rs',seq.getDefinition('Name'))
    fprintf('performing special data handling for ute ...\n');
    % average every 2nd spoke because of the half-rf excitation
    ktraj_adc=reshape(ktraj_adc,[3,adc_len,readouts]);
    ktraj_adc=ktraj_adc(:,:,1:2:end-1);
    %rawdata=rawdata(:,1:2:end-1,:);
    rawdata=rawdata(:,1:2:end-1,:)+rawdata(:,2:2:end,:);
    readouts=readouts/2;
    ktraj_adc=reshape(ktraj_adc,[3,adc_len*readouts]);
    % make it 2D
    ktraj_adc(3,:)=0;
    % reshape raw data to 4D (1 samples 1 channels) as BART needs it
    rawdata=reshape(rawdata,[1 adc_len*readouts 1 channels]);
else
    % just a general reshape to 4D (1 samples 1 channels) as BART needs it
    rawdata = reshape(rawdata, [1 adc_len*readouts 1 channels]); 
end

%% BART expects the trajectory in units of 1/FOV, Pulseq's trajectory is in inverse meters, need to rescale
for i=1:3, ktraj_adc(i,:)=ktraj_adc(i,:)*fov(i); end

% left-right swap (tested on TRIO)
ktraj_adc(1,:)=-ktraj_adc(1,:);

%% detect and clean-up 2D trajectory
kmax=max(abs(ktraj_adc'));
klim_almost_zero=1e-3; % why? 
if any(kmax<klim_almost_zero)
    fprintf('2D trajectory detected\n');
    ktraj_adc(kmax<klim_almost_zero,:)=0;
end

%% calculate signal power in respective channels
datapwr=squeeze(sum(rawdata.*conj(rawdata),2));
[~,pwrsort]=sort(-datapwr);

%% call bart for the actual recon -- first try inverse gridding
igrid = bart('nufft -i -t -lowmem1', ktraj_adc, rawdata); % inverse interative gridding (-i) % -lowmem1 may be of help with large data
save([basic_file_path '_recon_igrid'],'igrid');
% display (needs imab for displaying multiple slices/partitions)
figure; imab(sum(abs(igrid).^2,4).^0.5); colormap gray; 
title('BART recon, inverse gridding'); drawnow;

%% call bart for the actual recon -- now try iterative non-linear inversion which also uses multiple coils in a parallel imaging way
nlinv = bart('nlinv -n -d4 -m1 -i16 -t', ktraj_adc, rawdata); % there is a -R<n> key, but I did not notice any effect
%nlinv = bart('nlinv -n -d4 -m1 -i16 -t', ktraj_adc, rawdata(:,:,:,pwrsort(1:3))); % only include N most important channels (to save time and/or RAM)
save([basic_file_path '_recon_nlinv'],'nlinv');
figure; imab(abs(nlinv)/max(abs(nlinv(:)))); colormap gray;
title('BART recon, nonlinear inversion'); drawnow;

%% call bart for the pics algorithm with compressed sensing (with dummy coil sensitivities)
pics = bart('pics -S -l2 -r5 -t', ktraj_adc, rawdata, ones(size(igrid))); % inverse interative gridding (-i) % -lowmem1 may be of help with large data
%pics_l1 = bart('pics -S -l1 -i100 -s0.0002 -r0.003 -t', ktraj_adc, rawdata, ones(size(igrid))); % inverse interative gridding (-i) % -lowmem1 may be of help with large data
% tested r 0.002 0.004 and 0.01 (l1, l1s or ++, l1ss or +++)
save([basic_file_path '_recon_pics'],'pics');
% display (needs imab for displaying multiple slices/partitions)
figure; imab(sum(abs(pics).^2,4).^0.5); colormap gray; 
title('BART recon, pics-l2'); drawnow;

% %%
% pics_l2 = bart('pics -S -l2 -r5 -t', ktraj_adc, rawdata, ones(size(igrid))); % inverse interative gridding (-i) % -lowmem1 may be of help with large data
% save([basic_file_path '_recon_picsL2'],'pics_l2');
% % display (needs imab for displaying multiple slices/partitions)
% figure; imab(sum(abs(pics_l2).^2,4).^0.5); colormap gray; drawnow;

%% 3D viewer 
% img3d=abs(pics_l1(:,:,end:-1:1));
% %img3d=imag(pics);
% %img3d(img3d<0)=0;
% % img3d=sum(abs(igrid).^2,4).^0.5;
% 
% s=size(img3d);
% [xx,yy,zz]=ndgrid(-s(1)/2:s(1)/2-1,-s(2)/2:s(2)/2-1,-s(3)/2:s(3)/2-1);
% img3d(vecnorm([xx(:),yy(:),zz(:)]')>min(s)/2)=0;
% figure; imab(permute(img3d(116:140,:,:),[2,3,1])); colormap gray;
% volumeViewer(img3d);
