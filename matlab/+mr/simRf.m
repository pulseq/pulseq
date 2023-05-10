function [Mz_z,Mz_xy,F,ref_eff,Mx_xy,My_xy]=simRf(rf,rephase_factor,prephase_factor) 
%simRf Simulate an RF pulse with the given pulse shape.
%   [Mz_z,Mz_xy,F,ref_eff,Mx_xy,My_xy]=simRf(pulse,prephase_factor,rephase_factor) 
%   Performs a rapid RF pulse simulation based on the rotation formalism.
%   The algorithm is optimized by using quaternions to represent rotations. 
%   The compulsory parameter 'rf' is the Pulseq RF pulse. Optional
%   parameter 'rephase_factor' is needed in several cases e.g. to correclty 
%   visualize the phase of the magnetization for slice-selective
%   excitation. Another optional parameter 'prephase_factor' is an 
%   experimental parameter useful for simulating refocusing pulses or
%   spoiling needed. 
%   Return values: 
%     Mz_z,Mz_xy:  z and xy comnponents of the magnetisation after the pulse
%                  assuming the unit magnetization was aligned with z before
%                  the pulse. Useful for assessing excitation RF pulses. 
%     F:           frequency axis in Hz 
%     ref_eff:     Refocusing efficiency of the pulse as a complex value.
%                  Magnitude of ref_eff seems to closely follow Mz_z. Phase
%                  of ref_eff is related to the effective phase of the RF
%                  pulse, e.g. the axis of the planar flip.
%     Mx_xy,My_xy: xy magnetizations after the RF pulse assuming the unit
%                  magnetization was aligned with x or y axis prior to the
%                  pulse, respectively. Useful for detailed analyses of
%                  refocusing pulses.
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
        
[bw,f0]=mr.calcRfBandwidth(rf,0.5,df*10,dt);

% adapt time stepping -- just some compromizes -- we stick to dt~1/bw/50
if bw>4e3
    dt=5e-6;
    if bw>1e4
        dt=2e-6;
        if bw>20000
            dt=1e-6;
        end
    end
end

T     = (1:round(rf.shape_dur/dt))*dt-0.5*dt;                           % timesteps axis [s]
F     = 2*pi*linspace(f0-bw_mul*bw/2,f0+bw_mul*bw/2,bw/df)';            % offset frequencies [rad/s] 

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
q=quat_multiply(q,Q);

% RF pulse simulation
for j=1:length(T)
    W = -dt*sqrt(abs(shapea(j))^2+F.^2); % effective field rotation angles
    n = dt * [real(shapea(j))*ones(sf) imag(shapea(j))*ones(sf) F]./abs(W); % effective field rotation axes
    Q = [cos(W/2) sin(W/2).*n];
    q=quat_multiply(q,Q);
end

% rephaser / right spoiler / refocusing pulse
W = -F*dt*length(T)*rephase_factor;     % effective field rotation angle
Q = [cos(W/2) zeros(sf) zeros(sf) sin(W/2)];
q=quat_multiply(q,Q);

% export results
F=F/(2*pi);
m=zeros(sf(1),4);

% excitation: start with M0=M_z
m(:,4)=1;
m0rf=quat_multiply(quat_conj(q),quat_multiply(m,q));
Mz_z=m0rf(:,4);
Mz_xy=m0rf(:,2)+1i*m0rf(:,3);

% refocusing: start both with M0=M_x and them M0=M_y
m=zeros(sf(1),4);
m(:,2)=1;
Mx_xy=quat_multiply(quat_conj(q),quat_multiply(m,q));
Mx_xy=Mx_xy(:,2)+1i*Mx_xy(:,3);
m=zeros(sf(1),4);
m(:,3)=1;
My_xy=quat_multiply(quat_conj(q),quat_multiply(m,q));
My_xy=My_xy(:,2)+1i*My_xy(:,3);
ref_eff=(Mx_xy+My_xy*1i)/2;
end 

function qout = quat_multiply( q, r )
%  quat_multiply: Calculate the product of two quaternions.

% Calculate vector portion of quaternion product
% vec = s1*v2 + s2*v1 + cross(v1,v2)
vec = [q(:,1).*r(:,2) q(:,1).*r(:,3) q(:,1).*r(:,4)] + ...
         [r(:,1).*q(:,2) r(:,1).*q(:,3) r(:,1).*q(:,4)]+...
         [ q(:,3).*r(:,4)-q(:,4).*r(:,3) ...
           q(:,4).*r(:,2)-q(:,2).*r(:,4) ...
           q(:,2).*r(:,3)-q(:,3).*r(:,2)];

% Calculate scalar portion of quaternion product
% scalar = s1*s2 - dot(v1,v2)
scalar = q(:,1).*r(:,1) - q(:,2).*r(:,2) - ...
             q(:,3).*r(:,3) - q(:,4).*r(:,4);

qout = [scalar  vec];
end
       
function q = quat_conj( q ) 
%  quat_conj Calculate the conjugate of a quaternion.
q(:,2:4) = -q(:,2:4);
end
