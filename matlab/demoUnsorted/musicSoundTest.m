%% init the music variables

mrMusic.init;

%% test with a harmonic scale generated explicitly
%  we don't use the 'melody' mechanism here and create pitches directly
scale = real([ c1 cis1 d1 dis1 e1 f1 fis1 g1 gis1 a1 ais1 h1 ]);
pitches= [scale/4 scale/2 scale scale*2 scale*4 real(c3)*2; ...
          scale/4 scale/2 scale scale*2 scale*4 real(c3)*2; ...
          scale/4 scale/2 scale scale*2 scale*4 real(c3)*2];
% pitches= real([0  0  0  0  0  0  0  0  0 0  0  0  0  0  0  0  0  ; ...
%                0  0  0  0  0  0  0  0  0 0  0  0  0  0  0  0  0  ; ...
%                c1 d1 e1 f1 g1 a1 h1 c2 0 c2 h1 a1 g1 f1 e1 d1 c1 ]);

durations=1/4 * ones(1,size(pitches,2));

%% Pulseq sequence 

% never use full gradient performance because it is almost impossible to
% avoid the mechanical resonances of the gradient system completely
sys = mr.opts('MaxGrad',18,'GradUnit','mT/m',...
    'MaxSlew',160,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6 ...
);  
seq=mr.Sequence(sys);      % Create a new sequence object

pulseqUseWave=false; % use the "UseWave" option with case as it is really demanding both on the Pulseq environment and the scanner (shape memory)

seq = mrMusic.musicToSequence(seq, pitches, durations, 'barDurationSeconds', 1, 'pulseqUseWave', pulseqUseWave);

%% output
if pulseqUseWave
    seq.setDefinition('Name', 'scale');
    seq.write('scale.seq');
else
    seq.setDefinition('Name', 'scale1');
    seq.write('scale1.seq');
end

return
%% play

seq.sound();

return
%% optional slow step for checking whether we are staying within slew rate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 

