%% Simulate some RF pulses and analyze their properties
%
% This example demonstrates the use of mr.simRf() function.

%% Create a system object
% we need it here to create fat-sat pulses
sys = mr.opts('B0', 2.89); 
%seq=mr.Sequence(sys);
thickness_mm=5;

%% 30 degree slice selective SINC pulse 
[rf30_sinc, gz] = mr.makeSincPulse(pi/6,'system',sys,'Duration',3e-3,'use','excitation',...
    'PhaseOffset',pi/2,'apodization',0.4,'timeBwProduct',4,'SliceThickness',thickness_mm*1e-3);

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

sl_th=mr.aux.findFlank(F2(end:-1:1)/gz.amplitude,M_xy(end:-1:1),0.5)-mr.aux.findFlank(F2/gz.amplitude,M_xy,0.5);
figure; plot(F2/gz.amplitude*1000,abs(M_xy),'LineWidth',1.5); title('simulated slice profile, 30° flip, SINC'); xlabel('through-slice pos, mm');
hold on; yline(0.25,'-.'); xline([-0.5 0.5]*thickness_mm,'--'); legend({'slice profile','half-amplitude line','desired thickness'});
fprintf('actual slice thickness : %.3f mm\n',sl_th*1e3);

%% 90 degree slice selective SINC pulse 
[rf90_sinc, gz] = mr.makeSincPulse(pi/2,'system',sys,'Duration',3e-3,'use','excitation',...
    'PhaseOffset',pi/2,'apodization',0.4,'timeBwProduct',4,'SliceThickness',thickness_mm*1e-3);

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf90_sinc);
[M_z,M_xy,F2]=mr.simRf(rf90_sinc);

%
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

%
sl_th=mr.aux.findFlank(F2(end:-1:1)/gz.amplitude,M_xy(end:-1:1),0.5)-mr.aux.findFlank(F2/gz.amplitude,M_xy,0.5);
figure; plot(F2/gz.amplitude*1000,abs(M_xy),'LineWidth',1.5); title('simulated slice profile, 90° flip, SINC'); xlabel('through-slice pos, mm');
hold on; yline(0.5,'-.'); xline([-0.5 0.5]*thickness_mm,'--'); legend({'slice profile','half-amplitude line','desired thickness'});
fprintf('actual slice thickness : %.3f mm\n',sl_th*1e3);

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
[rf_90slr, gz] = mr.makeSLRpulse(pi/2,'duration',3e-3,'timeBwProduct',4,'PhaseOffset',pi/2,'use','excitation',...
    'passbandRipple',1,'stopbandRipple',1e-2,'filterType','ms','system',sys,'SliceThickness',thickness_mm*1e-3); 

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

%
sl_th=mr.aux.findFlank(F2(end:-1:1)/gz.amplitude,M_xy(end:-1:1),0.5)-mr.aux.findFlank(F2/gz.amplitude,M_xy,0.5);
figure; plot(F2/gz.amplitude*1000,abs(M_xy),'LineWidth',1.5); title('simulated slice profile, 90° flip, SLR'); xlabel('through-slice pos, mm');
hold on; yline(0.5,'-.'); xline([-0.5 0.5]*thickness_mm,'--'); legend({'slice profile','half-amplitude line','desired thickness'});
fprintf('actual slice thickness : %.3f mm\n',sl_th*1e3);

%% 60 degree BLOCK pulse 
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
fs_dur= 8e-3;
fs_bw_mul=1.4;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',fs_dur,...
    'bandwidth',fs_bw_mul*abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase    

[M_z,M_xy,F2]=mr.simRf(rf_fs);

figure; plot(F2,real(M_xy),F2,imag(M_xy),F2,M_z);
axis([sat_freq-900, sat_freq+900, -1.2, 1.2]);
legend({'M_x','M_y','M_z'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Simulation, Gaussian fat-sat pulse');

%% SLR fat-sat pulse 
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*sys.B0*sys.gamma;
fs_dur= 12e-3; % duration of 8 ms is sufficient for Gauss but is insufficinet for SLR
fs_bw_mul=1.2;
rf_fs = mr.makeSLRpulse(110*pi/180,'duration',fs_dur,'timeBwProduct',fs_bw_mul*abs(sat_freq)*fs_dur,'freqOffset',sat_freq,'use','saturation',...
    'passbandRipple',1,'stopbandRipple',1e-2,'filterType','ms','system',sys); 
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase    

[M_z,M_xy,F2]=mr.simRf(rf_fs);

figure; plot(F2,real(M_xy),F2,imag(M_xy),F2,M_z);
axis([sat_freq-900, sat_freq+900, -1.2, 1.2]);
legend({'M_x','M_y','M_z'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Simulation, SLR fat-sat pulse');


%% 180 degree slice selective SINC pulse 
[rf180_sinc, gz] = mr.makeSincPulse(pi,'system',sys,'Duration',4e-3,'use','refocusing',...
    'apodization',0.3,'timeBwProduct',6,'SliceThickness',thickness_mm*1e-3);

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

sl_th=mr.aux.findFlank(F2(end:-1:1)/gz.amplitude,ref_eff(end:-1:1),0.5)-mr.aux.findFlank(F2/gz.amplitude,ref_eff,0.5);
figure; plot(F2/gz.amplitude*1000,abs(ref_eff),'LineWidth',1.5); title('simulated slice profile, 180° flip, SINC'); xlabel('through-slice pos, mm');
hold on; yline(0.5,'-.'); xline([-0.5 0.5]*thickness_mm,'--'); legend({'slice profile','half-amplitude line','desired thickness'});
xlim([-10 10]); ylim([0 1.05]);
fprintf('actual slice thickness : %.3f mm\n',sl_th*1e3);

figure; plot(F2,angle(ref_eff)); 
axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'SINC'});
title('refocusing efficiency phase (~2x RF phase)'); 

%% 180 degree slice selective SLR pulse 
[rf180_slr, gz]= mr.makeSLRpulse(pi,'duration',4e-3,'timeBwProduct',6,'use','refocusing','filterType','ms','system',sys,'SliceThickness',thickness_mm*1e-3); 

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
%%
sl_th=mr.aux.findFlank(F2_slr(end:-1:1)/gz.amplitude,ref_eff_slr(end:-1:1),0.5)-mr.aux.findFlank(F2_slr/gz.amplitude,ref_eff_slr,0.5);
figure; plot(F2_slr/gz.amplitude*1000,abs(ref_eff_slr),'LineWidth',1.5); title('simulated slice profile, 180° flip, SLR'); xlabel('through-slice pos, mm');
hold on; yline(0.5,'-.'); xline([-0.5 0.5]*thickness_mm,'--'); legend({'slice profile','half-amplitude line','desired thickness'});
xlim([-10 10]); ylim([0 1.05]);
fprintf('actual slice thickness : %.3f mm\n',sl_th*1e3);
%%
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

%% investigare the RF center shift and the effective slice thickness as a function of flip angle

alphas=[5:5:180];
rfce=[];
sl_th_Mxy=[];
sl_th_Mz=[];
th_nom=5e-3;
for a=alphas
    [rfAlpha, gz] = mr.makeSincPulse(a/180*pi,'system',sys,'Duration',3e-3,'use','excitation','apodization',0.4,'timeBwProduct',4,'SliceThickness',th_nom);
    %[rfAlpha, gz] = mr.makeSLRpulse(a/180*pi,'system',sys,'Duration',3e-3,'SliceThickness',th_nom,'timeBwProduct',4,'dwell',5e-6,'passbandRipple',1,'stopbandRipple',1e-2,'filterType','ms','use','excitation');
    %[rfAlpha, gz] = mr.makeGaussPulse(a*pi/180,'system',sys,'Duration',3e-3,'timeBwProduct',8,'use','excitation','SliceThickness',th_nom);

    [bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rfAlpha);
    [M_z,M_xy,F2]=mr.simRf(rfAlpha);
    
    M_xy_masked=M_xy;
    M_xy_masked(F2<f0-bw/2)=0;
    M_xy_masked(F2>f0+bw/2)=0;
    
    i_phase_slope=sum(M_xy_masked(2:end).*conj(M_xy_masked(1:end-1)));
    i_phase_slope=1i*angle(i_phase_slope)/(F2(2)-F2(1));
    
    rfce(end+1)=abs(i_phase_slope)/2*pi/rf90_sinc.shape_dur/10; % no idea where this 10 comes from
    
    sl_th_Mxy(end+1)=mr.aux.findFlank(F2(end:-1:1)/gz.amplitude,M_xy(end:-1:1),0.5)-mr.aux.findFlank(F2/gz.amplitude,M_xy,0.5); % thinkness of the excite slice
    sl_th_Mz(end+1)=mr.aux.findFlank(F2(end:-1:1)/gz.amplitude,1-M_z(end:-1:1),0.5)-mr.aux.findFlank(F2/gz.amplitude,1-M_z,0.5); % related to the thickness of the slice if used as a refocusing pulse
end
% fix phase fit results for too large flip angles 
rfce(alphas>140)=NaN;

figure;plot(alphas,rfce/0.5*100); title('excess gradient refocusing needed in %');
xticks([0:30:alphas(end)]);
yticks([0:2:10]);
grid on;
xlabel('flip angle / °');
ylabel('grad. refocusing / %');
%
figure; plot(alphas,sl_th_Mxy/th_nom*100); title('slice thickness in %');
hold on; plot(alphas,sl_th_Mz/th_nom*100);
legend('Mxy', 'Mz','Location','northwest');
xticks([0:30:alphas(end)]);
%yticks([0:2:10]);
grid on;
xlabel('flip angle / °');
ylabel('slice thickness / %');

