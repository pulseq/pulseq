% Reconstruction of 2D Cartesian Pulseq data
% provides an example on how data reordering can be detected from the MR
% sequence with almost no additional prior knowledge
%
% needs mapVBVD in the path

%% Load the latest file from a dir
path='../IceNIH_RawSend/'; % directory to be scanned for data files
%path='~/Dropbox/shared/data/siemens/';
%path='~/Dropbox/shared/data/siemens/demo_gre/';
pattern='*.dat';

D=dir([path pattern]);
[~,I]=sort([D(:).datenum]);
%
data_file_path=[path D(I(end)).name];  % use end-1 to reconstruct the second-last data set, etc.
%%
twix_obj = mapVBVD(data_file_path);

%% Load sequence from file (optional, you may still have it in memory)

seq_file_path = [data_file_path(1:end-3) 'seq'];

seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');
[ktraj_adc, ktraj, t_excitation, t_refocusing] = seq.calculateKspace();
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal');

%% Analyze the trajectory data (ktraj_adc)
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
[~,k0_ind]=min(sum(ktraj_adc.^2,1));
kindex=round((ktraj_adc-ktraj_adc(:,k0_ind*ones(1,size(ktraj_adc,2))))./dk(:,ones(1,size(ktraj_adc,2))));
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

%% sort in the k-space data
if iscell(twix_obj)
    data_unsorted = twix_obj{end}.image.unsorted();
else
    data_unsorted = twix_obj.image.unsorted();
end

% the incoming data order is [kx coils acquisitions]
data_coils_last = permute(data_unsorted, [1, 3, 2]);
nCoils = size(data_coils_last, 3);

data_coils_last=reshape(data_coils_last, [size(data_coils_last,1)*size(data_coils_last,2), nCoils]);

data=zeros([kindex_end' nCoils]);
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

%figure; imab(data);

%% Reconstruct coil images

images = zeros(size(data));
%figure;

for ii = 1:nCoils
    %images(:,:,:,ii) = fliplr(rot90(fftshift(fft2(fftshift(data(:,:,:,ii))))));
    images(:,:,:,ii) = fftshift(fft2(fftshift(data(end:-1:1,:,:,ii))));
    %subplot(2,2,ii);
    %imshow(abs(images(:,:,ii)), []);
    %title(['RF Coil ' num2str(ii)]);
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
figure;
if nCoils>1
    sos=abs(sum(images.*conj(images),ndims(images))).^(1/2);
    sos=sos./max(sos(:));    
    imab(sos);
    %imwrite(sos, ['img_combined.png']
else
    imab(abs(images));
end
colormap('gray');

%% reconstruct field map (optional)

if size(images,3)==2
    cmplx_diff=images(:,:,2,:).*conj(images(:,:,1,:));
    phase_diff_image=angle(sum(cmplx_diff,4));
    figure;
    imab(phase_diff_image);colormap('jet');
end

%% gif movie export

scale=1.2;
filename='fftReconMovie';
fps=5;

if size(images,3)>4

    clm=gray(256);
    
    if nCoils>1
        imex=sos;
    else
        imex=abs(images);
        imex=imex/max(imex(:));
    end
    
    imex=permute(imex(:,end:-1:1,:),[2,1,3]);

    imind=uint8(scale*imex*(size(clm,1)-1)+1);
    imwrite(imind(:,:,1),clm, [filename, '.gif'],'DelayTime',1/fps,'Loopcount',inf);
    for i=2:size(imind,3),
        imwrite(imind(:,:,i),clm, [filename, '.gif'],'DelayTime',1/fps, 'WriteMode','append');
    end

end
