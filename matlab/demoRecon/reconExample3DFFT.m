% Reconstruction of 2D Cartesian Pulseq data
% provides an example on how data reordering can be detected from the MR
% sequence with almost no additional prior knowledge
%
% it loads Matlab .mat files with the rawdata in the format 
%     adclen x channels x readouts
% if Matlab .mat file not available it attempt to load Siemens .dat (which needs mapVBVD in the path)
% but first it seeks the accompanying .seq file with the same name to interpret
%     the data

%% Load the latest file from the specified directory
path='../IceNIH_RawSend/'; % directory to be scanned for data files
%path='/data/Dropbox/ismrm2021pulseq_liveDemo/dataLive/Vienna_7T_Siemens'; % directory to be scanned for data files

pattern='*.seq';
D=dir([path filesep pattern]);
[~,I]=sort([D(:).datenum]);
seq_file_path=[path filesep D(I(end-0)).name]; % use end-1 to reconstruct the second-last data set, etc...
                                                % or replace I(end-0) with I(1) to process the first dataset, I(2) for the second, etc...
%seq_file_path='../interpreters/siemens/data_example/gre_example.seq'

% keep basic filename without the extension
[p,n,e] = fileparts(seq_file_path);
basic_file_path=fullfile(p,n);

% try loading Matlab data
data_file_path=[basic_file_path '.mat'];
if isfile(data_file_path)
    fprintf(['loading `' data_file_path '´ ...\n']);
    data_unsorted = load(data_file_path);
    if isstruct(data_unsorted)
        fn=fieldnames(data_unsorted);
        assert(length(fn)==1); % we only expect a single variable
        data_unsorted=data_unsorted.(fn{1});
    end
else
    % revert to Siemens .dat file
    data_file_path=[basic_file_path '.dat'];
    fprintf(['loading `' data_file_path '´ ...\n']);
    twix_obj = mapVBVD(data_file_path);
    if iscell(twix_obj)
        data_unsorted = twix_obj{end}.image.unsorted();
    else
        data_unsorted = twix_obj.image.unsorted();
    end
end
[adc_len,channels,readouts]=size(data_unsorted);

%% Load sequence from file 
fprintf(['loading `' seq_file_path '´ ...\n']);
seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); title('2D kx/ky k-space trajectory');

figure; plot(ktraj(1,:),ktraj(3,:),'b',...
             ktraj_adc(1,:),ktraj_adc(3,:),'r.'); % another 2D plot
axis('equal'); title('2D kx/kz k-space trajectory');

%% Analyze the trajectory data (ktraj_adc)
fprintf('analyzing the k-space trajectory ...\n');
k_extent=max(abs(ktraj_adc),[],2);
k_scale=max(k_extent);
k_threshold=k_scale/5000;

% detect unused dimensions and delete them
if any(k_extent<k_threshold)
    ktraj_adc(k_extent<k_threshold,:)=[]; % delete rows
    k_extent(k_extent<k_threshold)=[];
end

% detect dK, k-space reordering and repetitions (or slices, etc)
kt_sorted=sort(ktraj_adc,2);
dk_all=kt_sorted(:,2:end)-kt_sorted(:,1:(end-1));
dk_all(dk_all<k_threshold)=NaN;
dk_min=min(dk_all,[],2);
dk_max=max(dk_all,[],2);
dk_all(dk_all-dk_min(:,ones(1,size(dk_all,2)))>k_threshold)=NaN;
dk_all_cnt=sum(isfinite(dk_all),2);
dk_all(~isfinite(dk_all))=0;
dk=sum(dk_all,2)./dk_all_cnt;
dk(~isfinite(dk))=0;
[~,k0_ind]=min(sum(ktraj_adc.^2,1));
kindex=round((ktraj_adc-ktraj_adc(:,k0_ind*ones(1,size(ktraj_adc,2))))./dk(:,ones(1,size(ktraj_adc,2))));
kindex(~isfinite(kindex))=0;
kindex_min=min(kindex,[],2);
kindex_mat=kindex-kindex_min(:,ones(1,size(ktraj_adc,2)))+1;
kindex_end=max(kindex_mat,[],2);
sampler=zeros(kindex_end');
repeat=zeros(1,size(ktraj_adc,2));
for i=1:size(kindex_mat,2)
    if (size(kindex_mat,1)==3)
        ind=sub2ind(kindex_end,kindex_mat(1,i),kindex_mat(2,i),kindex_mat(3,i));
    else
        ind=sub2ind(kindex_end,kindex_mat(1,i),kindex_mat(2,i)); 
    end
    repeat(i)=sampler(ind);
    sampler(ind)=repeat(i)+1;
end
if (max(repeat(:))>0)
    kindex=[kindex;(repeat+1)];
    kindex_mat=[kindex_mat;(repeat+1)];
    kindex_end=max(kindex_mat,[],2);
end
%figure; plot(kindex(1,:),kindex(2,:),'.-');

%% sort the k-space data into the data matrix
% the incoming data order is [kx coils acquisitions]
data_coils_last = permute(data_unsorted, [1, 3, 2]);
data_coils_last = reshape(data_coils_last, [adc_len*readouts, channels]);

data=zeros([kindex_end' channels]);
if (size(kindex,1)==3)
    for i=1:size(kindex,2)
        data(kindex_mat(1,i),kindex_mat(2,i),kindex_mat(3,i),:)=data_coils_last(i,:);
    end
else
    for i=1:size(kindex,2)
        data(kindex_mat(1,i),kindex_mat(2,i),:)=data_coils_last(i,:);
    end
end

if size(kindex,1)==3
    nImages=size(data,3);
else
    nImages=1;
    data=reshape(data, [size(data,1) size(data,2) 1 size(data,3)]); % we need a dummy images/slices dimension
end

%figure; imab(log(abs(data))); title('k-space data');

%% Reconstruct coil images

images = zeros(size(data));
%figure;

for ii = 1:channels
    images(:,:,:,ii) = fftshift(fft2(fftshift(data(:,:,:,ii)))); % 1.4.0. does not need inversion of the read direction
    %for ni = 1:nImages
        %tmp = abs(images(:,:,ni,ii));
        %tmp = tmp./max(tmp(:));
        %imwrite(tmp, ['img_coil_' num2str(ii) '.png'])
    %end
end

% Phase images (possibly channel-by-channel and echo-by-echo)
%figure;imab(angle(images));colormap('jet');
%figure;imab(abs(images));colormap('gray');

%% Image display with optional sum of squares combination
% figure;
% if nCoils>1
%     sos=abs(sum(images.*conj(images),ndims(images))).^(1/2);
%     sos=sos./max(sos(:));    
%     imab(sos);
%     %imwrite(sos, ['img_combined.png']
% else
%     imab(abs(images));
% end
% colormap('gray');

%% 3D image now

images3D=fftshift(fft(fftshift(images,3),[],3),3);
if channels>1
    sos3D=abs(sum(images3D.*conj(images3D),ndims(images3D))).^(1/2);    
    im3D=sos3D./max(sos3D(:));
else
    im3D=abs(images3D);
end

figure;imab(im3D);colormap('gray');title('all reconstructed image(s)');
figure;imab(im3D(:,:,end/2+1));colormap('gray'); title('central partition');
