% very basic FID data handling example
%
% needs mapVBVD in the path

%% Load the latest file from a dir
path='../IceNIH_RawSend/'; % directory to be scanned for data files
pattern='*.dat';
D=dir([path pattern]);
[~,I]=sort([D(:).datenum]);
data_file_path=[path D(I(end-0)).name]; % use end-1 to reconstruct the second-last data set, etc.

twix_obj = mapVBVD(data_file_path);

%% raw data preparation
if iscell(twix_obj)
    data_unsorted = double(twix_obj{end}.image.unsorted());
else
    data_unsorted = double(twix_obj.image.unsorted());
end

channels=size(data_unsorted,2);
adc_len=size(data_unsorted,1);
readouts=size(data_unsorted,3);

rawdata = permute(data_unsorted, [1,3,2]);

%% Load sequence from file 
seq = mr.Sequence();              % Create a new sequence object
seq_file_path = [data_file_path(1:end-3) 'seq'];
seq.read(seq_file_path,'detectRFuse');

%% we just want t_adc
[~,t_adc,~,~,t_excitation]=seq.calculateKspacePP();

t_excitation_adc=t_excitation((end-readouts+1):end); % remove dummy/prep scans
t_adc=reshape(t_adc,[adc_len,readouts]);
t_e=t_adc-t_excitation_adc(ones(1,adc_len),:);

% plot raw data
figure; plot(t_e, abs(rawdata(:,:))); title('raw signal(s)'); xlabel('time /s');
hold on; plot(t_e(1),0); % trick to force Y axis scaling

if (readouts>1)
    figure; plot(abs(rawdata(1,:))); title('time evolution (first FID point)'); 
    dataavg=squeeze(mean(rawdata,2));
    figure; plot(t_e, abs(dataavg)); title('averaged FIDs');
end

% calculate frequency axis
f=zeros(size(t_e));
adc_dur=(t_adc(end,1)-t_adc(1,1))/(adc_len-1)*adc_len;
f(:)=1/adc_dur;
f=cumsum(f,1)-0.5/adc_dur*adc_len;
% plot spectra
figure; plot(f, abs(ifftshift(ifft(ifftshift(rawdata(:,:),1)),1))); title('spectr(um/a)'); xlabel('frequency /Hz');
