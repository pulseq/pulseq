%% Generate a sequence in MATLAB and write to files
% 
% The |Sequence| class provides functionality to create magnetic
% resonance sequences (MRI or NMR) from basic building blocks.
%
% This provides an implementation of the open file format for MR sequences
% described here: http://pulseq.github.io/specification.pdf
%
% This example performs the following steps:
% 
% # Create slice selective RF pulse for imaging.
% # Create readout gradient and phase encode strategy.
% # Loop through phase encoding and generate sequence blocks.
% # Write the sequence to an open file format suitable for execution on a
% scanner.
% 

%%% Instantiation and gradient limits
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence();

%%%
% The gradient limits can be specified in various units _mT/m_, _Hz/cm_, or
% _Hz/m_. However the limits will be stored internally in units of _Hz/m_
% for amplitude and _Hz/m/s_ for slew.
gradOpts = mr.opts('MaxGradAmp',30,'AmpUnit','mT/m',...
    'MaxGradSlew',170,'SlewUnit','T/m/s')


%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
fov=200e-3;
Nx=64; Ny=64;
sliceThickness=3e-3;

%%% Slice selection
% Key concepts in the sequence description are *blocks* and *events*.
% Blocks describe a group of events that are executed simulaneously. This
% heirachial structure means that one event can be used in multiple blocks,
% a common occurance in MR sequences, particularly in imaging sequences.
%
% First, a slice selective RF pulse (and corresponding slice gradient) can
% be generated using the |makeSincPulse| function.
%
flip=90*pi/180;
[rf, gz] = mr.makeSincPulse(flip,'Duration',3e-3,'SliceThickness',5e-3,...
    'apodization',0.5,'timeBwProduct',4,gradOpts);

figure
plot(rf.t,real(rf.signal))


%%% Readout gradient
% To define the remaining encoding gradients we need to calculate the
% $k$-space sampling. The Fourier relationship
%
% $$\Delta k = \frac{2\pi}{FOV}$$
% 
% Therefore the area of the readout gradient is $n\Delta k$.
deltak=2*pi/fov;
kWidth = Nx*deltak;
readoutTime = 3.2e-3;
gx = mr.makeTrapezoid('x','FlatArea',kWidth,'FlatTime',readoutTime,gradOpts);

adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime);

%%% Phase encoding
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used. Furthermore rephasing of the slice
% select gradient is required.
gxPre = mr.makeTrapezoid('x','Area',gx.area/2,'Duration',2e-3,gradOpts);
gzReph = mr.makeTrapezoid('x','Area',gz.area/2,'Duration',2e-3,gradOpts);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

%%% Calculate timing
delayTE=20e-3 - mr.calcDuration(gxPre) - mr.calcDuration(rf)/2 ...
    - mr.calcDuration(gx)/2;
delayTR=50e-3 - mr.calcDuration(gxPre) - mr.calcDuration(rf) ...
    - mr.calcDuration(gx) - delayTE;
delay1 = mr.makeDelay(delayTE);
delay2 = mr.makeDelay(delayTR);

%% Define sequence blocks
% Next, the blocks are put together to form the sequence

for i=1:Ny
    seq.addBlock(rf,gz);
    gyPre = mr.makeTrapezoid('x','Area',phaseAreas(i),'Duration',2e-3,gradOpts);
    seq.addBlock(gxPre,gyPre,gzReph);
    seq.addBlock(delay1);
    seq.addBlock(gx,adc);
    seq.addBlock(delay2)
end

%% Write to file
% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
seq.write('external.seq')

%%
% Display the first few lines of the output file
s=fileread('external.seq');
disp(s(1:309))

