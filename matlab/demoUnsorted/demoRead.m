%% Read a sequence into MATLAB
%
% The |Sequence| class provides an implementation of the _open file
% format_ for MR sequences described here: http://pulseq.github.io/specification.pdf
%
% This example demonstrates parsing an MRI sequence stored in this format,
% accessing sequence parameters and visualising the sequence.

%% Read a sequence file
% A sequence can be loaded from the open MR file format using the |read|
% method.
%seq_name='gre.seq'; 
%seq_name='epi_rs.seq';
seq_name='spiral.seq';
%seq_name='idea/external.seq';
%seq_name='trufi.seq'; 
%seq_name='epi_se.seq';
%seq_name='../tests/gre.seq';
%seq_name='../tests/epi_rs.seq';
%seq_name='/home/zaitsev/pulseq_home/open_sequence/ufr-collab/tmp/gre.seq';
%seq_name='historical_format_tests/gre_example_120.seq';
%seq_name='/home/zaitsev/range_software/pulseq/matlab/user_seq/test592_MSE_vFA_T2-45ms_nEchos-1_nSlices-1_TE-8ms_TR-4133ms_FOV-256mm_Nx-262_Ny-1_gm-32_sm-120_sTexc-3mm_sTrefoc-9mm_acc_noAcc.seq'
%sys = mr.opts('B0', 2.89); % we need system here if we want 'detectRFuse' to detect fat-sat pulses
sys = mr.opts('B0', 2.89, 'MaxGrad', 22, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 10e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);
seq2=mr.Sequence(sys);
%seq.read(seq_name,'detectRFuse');
seq2.read(seq_name);

%%
[ok, error_report]=seq2.checkTiming ;

if (ok)
    fprintf('Timing check passed successfully\n') ;
else
    fprintf('Timing check failed! Error listing follows:\n') ;
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq2.plot('timeDisp','us','showBlocks',1,'timeRange',[0 40e-3]); %detailed view
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq2.calculateKspacePP();
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
axis('equal'); % enforce aspect ratio for the correct trajectory display
title('k-space trajectory (k_x/k_y)');

return
%%
rf_ph=[];
rf_freq=[];
adc_ph=[];
adc_freq=[];
for iB=1:length(seq2.blockDurations)
    b=seq2.getBlock(iB);
    if ~isempty(b.rf)
        rf_ph(end+1)=b.rf.phaseOffset;
        rf_freq(end+1)=b.rf.freqOffset;
    end
    if ~isempty(b.adc)
        adc_ph(end+1)=b.adc.phaseOffset;
        adc_freq(end+1)=b.adc.freqOffset;
    end
end

%%
figure; plot(adc_ph);
figure; plot(rf_ph);

%% sanity check to see if the reading and writing are consistent
% seq.write('read_test.seq');
% system(['diff -s -u ' seq_name ' read_test.seq'],'-echo');

%% Access sequence parameters and blocks
% Parameters defined with in the |[DEFINITIONS]| section of the sequence file
% are accessed with the |getDefinition| method. These are user-specified
% definitions and do not effect the execution of the sequence.
seqName=seq2.getDefinition('Name')

%% calculate and display real TE, TR as well as slew rates and gradient amplitudes

rep = seq2.testReport; 
fprintf([rep{:}]); 

%%
% Sequence blocks are accessed with the |getBlock| method. As shown in the
% output the first block is a selective excitation block and contains an RF
% pulse and gradient and on the z-channel.
b1=seq2.getBlock(1)

%%
% Further information about each event can be obtained by accessing the
% appropriate fields of the block struct. In particular, the complex RF
% signal is stored in the field |signal|.
rf=b1.rf

figure;
subplot(211);
plot(rf.t,abs(rf.signal));
ylabel('RF magnitude');
subplot(212);
plot(1e3*rf.t,angle(rf.signal))
ylabel('RF phase'); xlabel('t (ms)');

%%
% The next three blocks contain: three gradient events; a delay; and
% readout gradient with ADC event, each with corresponding fields defining
% the details of the events.
b2 = seq2.getBlock(2);
b3 = seq2.getBlock(3);
b4 = seq2.getBlock(4);
b2.gx
%b3.delay
b4.adc


%% Plot the sequence 
% Visualise the sequence using the |plot| method of the class. This creates
% a new figure and shows ADC, RF and gradient events. The axes are linked
% so zooming is consistent. In this example, a simple gradient echo sequence
% for MRI is displayed.
seq2.plot()
return
%%
% The details of individual pulses are not well-represented when the entire
% sequence is visualised. Interactive zooming is helpful here.
% Alternatively, a time range can be specified.  An additional parameter
% also allows the display units to be changed for easy reading. 
% Further, the handle of the created figure can be returned if required.
fig=seq2.plot('TimeRange',[0 16e-3],'timeDisp','ms')

%% Modifying sequence blocks
% In addition to loading a sequence and accessing sequence blocks, blocks
% can be modified. In this example, a Hamming window is applied to the
% first RF pulse of the sequence and the flip angle is changed to 45
% degrees. The remaining RF pulses are unchanged. 

rf2=rf;
duration=rf2.t(end);
t=rf2.t-duration/2;                                 % Centre time about 0
alpha=0.5;
BW=4/duration;                                      % time bandwidth product = 4
window = (1.0-alpha+alpha*cos(2*pi*t/duration));    % Hamming window
signal = window.*sinc(BW*t);

% Normalise area to achieve 2*pi rotation
signal=signal./(seq2.rfRasterTime*sum(real(signal)));

% Scale to 45 degree flip angle
rf2.signal=signal.*45/360;

b1.rf=rf2;
seq2.setBlock(1,b1);

%% second check to see what we have changed
% seq.write('read_test2.seq');
% system(['diff -s -u ' seq_name ' read_test2.seq'],'-echo');


%%
% The amplitude of the first rf pulse is reduced due to the reduced
% flip-angle. Notice the reduction is not exactly a factor of two due to
% the windowing function.
amp1_in_Hz = max(abs(seq2.getBlock(1).rf.signal))
amp2_in_Hz = max(abs(seq2.getBlock(6).rf.signal))
