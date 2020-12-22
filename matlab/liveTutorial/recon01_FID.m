% very basic FID data handling example
%
% needs mapVBVD in the path

%% Load the latest file from a dir
path='../IceNIH_RawSend/'; % directory to be scanned for data files
pattern='*.dat';
D=dir([path pattern]);
[~,I]=sort([D(:).datenum]);
data_file_path=[path D(I(end)).name]; % use end-1 to reconstruct the second-last data set, etc.

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

%% we just want t_adc
[ktraj_adc, t_adc, ktraj, t, t_excitation, t_refocusing]=seq.calculateKspacePP();

t_adc=reshape(t_adc,[adc_len,readouts]);
t_e=t_adc-t_excitation(ones(1,adc_len),:);

% plot raw data
figure; plot(t_e, abs(rawdata(:,:))); title('raw FIDs');
hold on; plot(t_e(1),0); % tric to force Y axis scaling

if (readouts>1)
    figure; plot(abs(rawdata(1,:))); title('time evolution (first FID point)'); 
    dataavg=squeeze(mean(rawdata,2));
    figure; plot(t_e, abs(dataavg)); title('averaged FIDs');
end