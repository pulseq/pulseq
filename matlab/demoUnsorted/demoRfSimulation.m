%% Simulate some RF pulses and analyze their properties
%
% This example demonstrates the use of mr.simRf() function.

%% Create a system object
% we need it here to create fat-sat pulses
sys = mr.opts('B0', 2.89); 
%seq=mr.Sequence(sys);

%% 30 degree slice selective SINC pulse 
rf30_sinc = mr.makeSincPulse(pi/6,'system',sys,'Duration',3e-3,'use','excitation',...
    'PhaseOffset',pi/2,'apodization',0.3,'timeBwProduct',4);

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf30_sinc);
[M_z,M_xy,F2]=mr.simRf(rf30_sinc);

figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
axis([f0-2*bw, f0+2*bw, -0.1, 1.2]);
legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('STA vs. simulation, flip angle 30°');

figure; plot(F2,real(M_xy),F2,imag(M_xy));
axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
legend({'M_xSIM','M_ySIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Real and imag. parts of transverse magnetisation, 30° flip');

%% 90 degree slice selective SINC pulse 
rf90_sinc = mr.makeSincPulse(pi/2,'system',sys,'Duration',3e-3,'use','excitation',...
    'PhaseOffset',pi/2,'apodization',0.6,'timeBwProduct',8);

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf90_sinc);
[M_z,M_xy,F2]=mr.simRf(rf90_sinc);

%%
figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
axis([f0-2*bw, f0+2*bw, -0.1, 1.2]);
legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('STA vs. simulation, flip angle 90°');

figure; plot(F2,atan2(abs(M_xy),M_z)/pi*180);
axis([f0-2*bw, f0+2*bw, -5, 100]);
xlabel('frequency offset / Hz');
ylabel('flip ange [°]');
legend({'SINC'});
grid on;
title('Achieved flip angle for the nominal 90° flip');

figure; plot(F2,real(M_xy),F2,imag(M_xy));
axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
legend({'M_xSIM','M_ySIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Real and imag. parts of transverse magnetisation, 90° flip');

%%
figure; plot(F2,angle(M_xy));
axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
xlabel('frequency offset / Hz');
ylabel('phase');
title('Phase of transverse magnetisation, 90° flip'); 

M_xy_masked=M_xy;
M_xy_masked(F2<f0-bw/2)=0;
M_xy_masked(F2>f0+bw/2)=0;

i_phase_slope=sum(M_xy_masked(2:end).*conj(M_xy_masked(1:end-1)));
i_phase_slope=1i*angle(i_phase_slope)/(F2(2)-F2(1));
i_phase_offset=sum(M_xy_masked.*exp(-i_phase_slope*F2));
i_phase_offset=i_phase_offset/abs(i_phase_offset);

hold on; plot(F2,angle(exp(i_phase_slope*F2)*i_phase_offset),'--');

xline(f0,'-');
xline(f0-bw/2,'--');
xline(f0+bw/2,'--');

legend({'SINC-phase','linear fit'});
fprintf('SINC90 rf center error: %g (%g %%)\n', abs(i_phase_slope)/2*pi/rf90_sinc.shape_dur/10, 100/0.5*abs(i_phase_slope)/2*pi/rf90_sinc.shape_dur/10); % no idea where this 10 comes from

%% 90 degree slice selective SLR pulse 
rf_90slr= mr.makeSLRpulse(pi/2,'duration',3e-3,'timeBwProduct',4,'PhaseOffset',pi/2,'use','excitation',...
    'passbandRipple',1,'stopbandRipple',1e-2,'filterType','ms','system',sys); 

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf_90slr);
[M_z,M_xy,F2]=mr.simRf(rf_90slr);

figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
axis([f0-2*bw, f0+2*bw, -0.1, 1.2]);
legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('SLR: STA vs. simulation, flip angle 90°');

figure; plot(F2,real(M_xy),F2,imag(M_xy));
axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
legend({'M_xSIM','M_ySIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('SLR: Real and imag. parts of transverse magnetisation, 90° flip');

figure; plot(F2,angle(M_xy));
axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
xlabel('frequency offset / Hz');
ylabel('phase');
legend({'SLR'});
title('SLR: Phase of transverse magnetisation, 90° flip'); 

figure; plot(F2,atan2(abs(M_xy),M_z)/pi*180);
axis([f0-2*bw, f0+2*bw, -5, 100]);
xlabel('frequency offset / Hz');
ylabel('flip ange [°]');
legend({'SLR'});
grid on;
title('SLR: Achieved flip angle for the nominal 90° flip');

%% 60 degree slice selective SINC pulse 
rf60_block = mr.makeBlockPulse(pi/3,'system',sys,'Duration',0.5e-3,'use','excitation','PhaseOffset',pi/2);

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf60_block);
[M_z,M_xy,F2]=mr.simRf(rf60_block);

figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
axis([f0-2*bw, f0+2*bw, -0.1, 1.2]);
legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('STA vs. simulation, 60° hard pulse');

figure; plot(F2,real(M_xy),F2,imag(M_xy));
axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
legend({'M_xSIM','M_ySIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Real and imag. parts of transverse magnetisation, 60° hard pulse');

%% fat-sat pulse 
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*sys.B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase    

[M_z,M_xy,F2]=mr.simRf(rf_fs);

figure; plot(F2,real(M_xy),F2,imag(M_xy),F2,M_z);
axis([sat_freq-900, sat_freq+900, -1.2, 1.2]);
legend({'M_x','M_y','M_z'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Simulation, Gaussian fat-sat pulse');

%% 180 degree slice selective SINC pulse 
rf180_sinc = mr.makeSincPulse(pi,'system',sys,'Duration',4e-3,'use','refocusing',...
    'apodization',0.3,'timeBwProduct',6);

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf180_sinc);
[M_z,M_xy,F2,ref_eff]=mr.simRf(rf180_sinc);

% figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
% axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
% legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
% xlabel('frequency offset / Hz');
% ylabel('magnetisation');
% title('STA vs. simulation, flip angle 180°');
% 
% figure; plot(F2,real(M_xy),F2,imag(M_xy));
% axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
% legend({'M_xSIM','M_ySIM'});
% xlabel('frequency offset / Hz');
% ylabel('magnetisation');
% title('Real and imag. parts of transverse magnetisation, 180° flip');
% 
% figure; plot(F2,angle(M_xy));
% axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
% xlabel('frequency offset / Hz');
% ylabel('phase');
% title('Phase of transverse magnetisation, 180° flip'); 

figure; plot(F2,atan2(abs(M_xy),M_z)/pi*180);
axis([f0-2*bw, f0+2*bw, -5, 190]);
xlabel('frequency offset / Hz');
ylabel('flip ange [°]');
legend({'SINC'});
grid on;
title('Achieved flip angle for the nominal 180° flip');

figure; plot(F2,abs(ref_eff)); 
axis([f0-2*bw, f0+2*bw, -0.1, 1.1]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'SINC'});
title('refocusing efficiency'); 

figure; plot(F2,angle(ref_eff)); 
axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'SINC'});
title('refocusing efficiency phase (~2x RF phase)'); 

%% 180 degree slice selective SLR pulse 
rf180_slr= mr.makeSLRpulse(pi,'duration',4e-3,'timeBwProduct',6,'use','refocusing','filterType','ms','system',sys); 

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf180_slr);
[M_z,M_xy,F2_slr,ref_eff_slr]=mr.simRf(rf180_slr);

figure; plot(F2_slr,atan2(abs(M_xy),M_z)/pi*180);
axis([f0-2*bw, f0+2*bw, -5, 190]);
xlabel('frequency offset / Hz');
ylabel('flip ange [°]');
legend({'SLR'});
grid on;
title('Achieved flip angle for the nominal 180° flip');

figure; plot(F2_slr,abs(ref_eff_slr)); 
axis([f0-2*bw, f0+2*bw, -0.1, 1.1]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'SLR'});
title('refocusing efficiency'); 

figure; plot(F2_slr,angle(ref_eff_slr)); 
axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
title('SLR: refocusing efficiency phase (~2x RF phase)'); 

figure; plot(F2,abs(ref_eff),F2_slr,abs(ref_eff_slr)); 
axis([f0-2*bw, f0+2*bw, -0.1, 1.1]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'SINC','SLR'});
title('refocusing efficiency: SINC vs SLR'); 

%% adiabatic pulse 
rf180_ad = mr.makeAdiabaticPulse('wurst','duration',4.6e-3,'bandwidth',6000,'n_fac',20,'use','inversion','system',sys); 

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf180_ad);
[M_z,M_xy,F2,ref_eff,ref_mx,ref_my]=mr.simRf(rf180_ad,-0.5);

figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('STA vs. simulation, adiabatic pulse');

figure; plot(F2,real(M_xy),F2,imag(M_xy));
axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
legend({'M_xSIM','M_ySIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Real and imag. parts of transverse magnetisation, adiabatic pulse');

figure; plot(F2,angle(M_xy));
axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
xlabel('frequency offset / Hz');
ylabel('phase');
legend({'WURST'});
title('Phase of transverse magnetisation, adiabatic pulse'); 

figure; plot(F2,atan2(abs(M_xy),M_z)/pi*180);
axis([f0-2*bw, f0+2*bw, -5, 190]);
xlabel('frequency offset / Hz');
ylabel('flip ange [°]');
legend({'WURST'});
grid on;
title('Achieved flip angle adiabatic pulse');

figure; plot(F2,abs(ref_eff)); 
axis([f0-2*bw, f0+2*bw, -0.1, 1.1]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'WURST'});
title('refocusing efficiency'); 

figure; plot(F2,angle(ref_eff)); 
axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'WURST'});
title('refocusing efficiency phase (~2x RF phase)'); 

%% spoiling simulation for the same pulse used and refocusing pulse

spoiling_factor=5; % area of the left/righ spoiler; reasonable range 1..10
cl=13; % convolution length to simulate intravoxel dephasing

[M_z,M_xy,F2,ref_eff,mxrf,myrf]=mr.simRf(rf180_ad,spoiling_factor,spoiling_factor); 

mxrfc=conv(mxrf,ones(cl,1)/cl,'same');
myrfc=conv(myrf,ones(cl,1)/cl,'same');

figure;plot(F2,abs(abs(mxrfc)+1i*abs(myrfc))/2^0.5,F2,abs(ref_eff),F2,0.5-0.5*M_z,'--'); 
xlabel('frequency offset / Hz');
ylabel('signal');
legend({'spoiling','ref.eff.','.5-.5*M_z'});
title('signal with spoiling vs refocusing efficiency'); 

%% investigare RF center shift as a function of clip angle (is it an artifact?)

alphas=[5:5:150];
rfce=[];
for a=alphas
    rfAlpha = mr.makeSincPulse(a/180*pi,'system',sys,'Duration',3e-3,'use','excitation','apodization',0.3,'timeBwProduct',4);
    %rfAlpha = mr.makeGaussPulse(a*pi/180,'system',sys,'Duration',3e-3,'timeBwProduct',8,'use','excitation');

    [bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rfAlpha);
    [M_z,M_xy,F2]=mr.simRf(rfAlpha);
    
    M_xy_masked=M_xy;
    M_xy_masked(F2<f0-bw/2)=0;
    M_xy_masked(F2>f0+bw/2)=0;
    
    i_phase_slope=sum(M_xy_masked(2:end).*conj(M_xy_masked(1:end-1)));
    i_phase_slope=1i*angle(i_phase_slope)/(F2(2)-F2(1));
    
    rfce=[rfce,abs(i_phase_slope)/2*pi/rf90_sinc.shape_dur/10]; % no idea where this 10 comes from
end

figure;plot(alphas,rfce/0.5*100); title('excess gradient refocusing needed in %');
xticks([0:30:alphas(end)]);
yticks([0:2:10]);
grid on;
xlabel('flip angle / °');
ylabel('refocusing / %');

