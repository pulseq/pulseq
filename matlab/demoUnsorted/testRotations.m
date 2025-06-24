% create two sequence objects
seq1=mr.Sequence; % object 1 (traditional, non-extension)
seq2=mr.Sequence; % object 2 (new, rotation extension) 

gx = mr.makeTrapezoid('x', 'area', 500);
gy = mr.makeExtendedTrapezoid('y', 'times', [0 100 300 400]*1e-6, 'amplitudes', [0 1 -1 0]*8e5);
gz = mr.makeArbitraryGrad('z',(1-cos((0:40)*pi/20))*4e5, 'first', 0, 'last', 0); 

% need to scale down all gradients to enable free rotations without exceeding limits
gx=mr.scaleGrad(gx,0.55);
gy=mr.scaleGrad(gy,0.55);
gz=mr.scaleGrad(gz,0.55);

% sanity check: no rotations
seq1.addBlock(gx,gy,gz);
seq2.addBlock(gx,gy,gz);

% simple polar angle
seq1.addBlock(mr.rotate('z', 20*pi/180, gx,gy,gz));
seq2.addBlock(gx,gy,gz,mr.makeRotation(20*pi/180));

% same test with rotation matrix on seq1 and a different angle
seq1.addBlock(mr.rotate3D(rotz(-35), gx,gy,gz));
seq2.addBlock(gx,gy,gz,mr.makeRotation(-35*pi/180));

% same with rotation matrix on both seq1 and seq2
seq1.addBlock(mr.rotate3D(rotz(35), gx,gy,gz));
seq2.addBlock(gx,gy,gz,mr.makeRotation(rotz(35)));

% rotation by the azimuthal angle only 
seq1.addBlock(mr.rotate('x', -20*pi/180, gx,gy,gz));
seq2.addBlock(gx,gy,gz,mr.makeRotation(0,-20*pi/180));

% same with rotation matrix on seq1
seq1.addBlock(mr.rotate3D(rotx(40), gx,gy,gz));
seq2.addBlock(gx,gy,gz,mr.makeRotation(0,40*pi/180));

% both azimuthal and polar angles
seq1.addBlock(mr.rotate('z',55*pi/180,mr.rotate('x', 70*pi/180, gx,gy,gz)));
seq2.addBlock(gx,gy,gz,mr.makeRotation(55*pi/180,70*pi/180));

% same with the rotation matrix on the extension side
seq1.addBlock(mr.rotate('z',75*pi/180,mr.rotate('x', -65*pi/180, gx,gy,gz)));
seq2.addBlock(gx,gy,gz,mr.makeRotation(rotz(75)*rotx(-65)));

% same with the rotation matrix on the conventional side
seq1.addBlock(mr.rotate3D(rotz(-80)*rotx(-35), gx,gy,gz));
seq2.addBlock(gx,gy,gz,mr.makeRotation(-80*pi/180,-35*pi/180));

% rotation about Y axis
seq1.addBlock(mr.rotate3D(roty(-50), gx,gy,gz));
seq2.addBlock(gx,gy,gz,mr.makeRotation([0 1 0],-50*pi/180));

% test waveforms
w1=seq1.waveforms_and_times();
w2=seq2.waveforms_and_times();

figure; plot(w1{1}(1,:),w1{1}(2,:),w1{2}(1,:),w1{2}(2,:),w1{3}(1,:),w1{3}(2,:));
hold on; plot(w2{1}(1,:),w2{1}(2,:),'--',w2{2}(1,:),w2{2}(2,:),'--',w2{3}(1,:),w2{3}(2,:),'--');
title('waveform comparisons');

%figure; plot(w1{1}(1,:),w1{1}(2,:)-w2{1}(2,:),w1{2}(1,:),w1{2}(2,:)-w2{2}(2,:),w1{3}(1,:),w1{3}(2,:)-w2{3}(2,:));
%title('waveform differences');

% k-space tests
[~,~,k1,t1]=seq1.calculateKspacePP();
[~,~,k2,t2]=seq2.calculateKspacePP();

figure; plot(t1,k1(1,:),t1,k1(2,:),t1,k1(3,:));
hold on; plot(t2,k2(1,:),'--',t2,k2(2,:),'--',t2,k2(3,:),'--');
title('k-space comparisons');
