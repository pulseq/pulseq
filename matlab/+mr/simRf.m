function [M_z,M_xy,F,ref_eff,ref_mx,ref_my]=simRf(rf,rephase_factor,prephase_factor) 
%simRf Simulate an RF pulse with the given pulse shape.
%   [M_z,M_xy,ref_eff,F]=simRf(pulse,prephase_factor,rephase_factor) 
%   Performs a rapid RF pulse simulation based on the rotation formalism.
%   The algorithm is additionally optimized by using quaternions to
%   represent rotations (may require some Matlab toolboxes). 
%   The compulsory parameter 'rf' is the Pulseq RF pulse. Optional
%   parameter 'rephase_factor' is needed in several cases e.g. to correclty 
%   visualize the phase of the magnetization for slice-selective
%   excitation. Another optional parameter 'prephase_factor' is an 
%   experimental parameter useful for simulating refocusing pulses or
%   spoiling needed. 
%
%   The implementation was inspired by the example by Dr. Tony Stoecker
%   (https://github.com/stoeckert/mr-simu-example-ismrm19)
%   The algorithm was rewritten to quaternions and vectorized for 
%   performance by MZ 
%

bw_mul=4;    % simulation bandwidth (multiplier of the pulse bandwidth)
df=1;        % spectral resolution [Hz]
dt=10e-6;    % (re-)sampling interval

if nargin < 2
    if isfield(rf,'use') && strcmp(rf.use,'refocusing')
        rephase_factor = 0;
    else
        rephase_factor = -(rf.shape_dur-mr.calcRfCenter(rf))/rf.shape_dur;
    end
end

if nargin < 3
    prephase_factor = 0;
end
        
[bw,f0,spectrum,FF,rfs,tt]=mr.calcRfBandwidth(rf,0.5,10,10e-6);

T     = (1:round(rf.shape_dur/dt))*dt-0.5*dt;                           % timesteps axis [s]
F     = 2*pi*linspace(f0-bw_mul*bw/2,f0+bw_mul*bw/2,bw/df)';               % offset frequencies [rad/s] 

shapea = interp1(rf.t, 2*pi*rf.signal.*exp(1i*(rf.phaseOffset+2*pi*rf.freqOffset*rf.t)),T,'linear',0);

% intialize result vectors
M_ROT=zeros(size(F)); 
Z_ROT=zeros(size(F));
sf=size(F);
q=zeros(sf(1),4);
q(:,1)=1; % init rotation quaternions

% prephaser / left spoiler
W = -F*dt*length(T)*prephase_factor;     % effective field rotation angle
Q = [cos(W/2) zeros(sf) zeros(sf) sin(W/2)];
q=quatmultiply(q,Q);

% RF pulse simulation
for j=1:length(T)
    W = -dt*sqrt(abs(shapea(j))^2+F.^2); % effective field rotation angles
    n = dt * [real(shapea(j))*ones(sf) imag(shapea(j))*ones(sf) F]./abs(W); % effective field rotation axes
    Q = [cos(W/2) sin(W/2).*n];
    q=quatmultiply(q,Q);
end

% rephaser / right spoiler / refocusing pulse
W = -F*dt*length(T)*rephase_factor;     % effective field rotation angle
Q = [cos(W/2) zeros(sf) zeros(sf) sin(W/2)];
q=quatmultiply(q,Q);

% export results
F=F/(2*pi);
m=zeros(sf(1),4);

% excitation: start with M0=M_z
m(:,4)=1;
m0rf=quatmultiply(quatconj(q),quatmultiply(m,q));
M_z=m0rf(:,4);
M_xy=m0rf(:,2)+1i*m0rf(:,3);

% refocusing: start both with M0=M_x and them M0=M_y
m=zeros(sf(1),4);
m(:,2)=1;
ref_mx=quatmultiply(quatconj(q),quatmultiply(m,q));
ref_mx=ref_mx(:,2)+1i*ref_mx(:,3);
m=zeros(sf(1),4);
m(:,3)=1;
ref_my=quatmultiply(quatconj(q),quatmultiply(m,q));
ref_my=ref_my(:,2)+1i*ref_my(:,3);
ref_eff=(ref_mx+ref_my*1i)/2;
