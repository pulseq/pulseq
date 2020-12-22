% very basic generic non-Cartesian recon using BART
% assumes raw data *.dat fiels are stored together with the *.seq Pulseq
% files withe the identical(!) names
%
% needs mapVBVD in the path
% modify or remove the line below depending on zour path settings
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
pattern='*.dat';

D=dir([path pattern]);
[~,I]=sort([D(:).datenum]);
data_file_path=[path D(I(end)).name]; % use end-1 to reconstruct the second-last data set, etc.

%% alternatively just provide the path to the .dat file

%data_file_path='../interpreters/siemens/data_example/gre_example.dat'
%data_file_path=[path '2020-11-10-073628.dat'];


%% load the raw data file
twix_obj = mapVBVD(data_file_path);

%% Load sequence from the file with the same name as the raw data
seq_file_path = [data_file_path(1:end-3) 'seq'];

seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');

%% calculate k-space trajectory
traj_recon_delay=0e-6;%1.75e-6; % adjust this parameter to potentially improve resolution & geometric accuracy. It can be calibrated by inverting the spiral revolution dimension and making two images match. for our Prisma and a particular trajectory we found 1.75e-6

%[ktraj_adc, t_adc, ktraj, t, t_excitation, t_refocusing]=seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace('trajectory_delay', traj_recon_delay);

% optionally transform or patch the trajectory
ktraj_adc(2,:)=-ktraj_adc(2,:); % seems to be needed for the correct orientaton
ktraj_adc(3,:)=0;

%% raw data preparation
if iscell(twix_obj)
    data_unsorted = double(twix_obj{end}.image.unsorted());
else
    data_unsorted = double(twix_obj.image.unsorted());
end

channels=size(data_unsorted,2);
adc_len=size(data_unsorted,1);
readouts=size(data_unsorted,3);
runs=readouts*adc_len/size(ktraj_adc,2);

rawdata = permute(data_unsorted, [1,3,2]);
rawdata = reshape(rawdata, [1 adc_len*readouts channels]); 

if runs>1 
    % average over runs
    readouts=readouts/runs;
    rawdata = reshape(squeeze(sum(reshape(rawdata, [adc_len*readouts runs channels]),2)),[1 adc_len*readouts channels]);
end

%% call bart for the actual recon
%igrid = bart('nufft -i -t', ktraj_adc/pi/(2^0.5), rawdata(:,:,1));
nlinv = bart('nlinv -n -d4 -m1 -i12 -t', ktraj_adc/pi/(2^0.5), rawdata(:,:,1)); % no idea why we seem to need /(2^0.5)

%% display (needs imab for displaying multiple slices)
figure; 
%imab(abs(igrid));
imab(abs(nlinv));
colormap gray;
