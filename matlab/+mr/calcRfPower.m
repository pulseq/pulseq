function [total_energy, peak_pwr, rf_rms]=calcRfPower(rf, dt)
%calcRfPower : Calculate the relative** power of the RF pulse
%   Returns the (relative) energy of the pulse expressed in the units of 
%   RF amplitude squared multiplied by time, e.g. in Pulseq these are
%   Hz * Hz * s = Hz. Sounds strange, but is true. 
%   Also returns peak power [Hz^2] and RMS B1 amplitude [Hz].
%   Optional parameter 'dt' can be used to specify sampling of RF pulses 
%   defined on a variable raster (currently only block pulses)
%   ** Note: the power and rf amplitude calculated by this function is 
%   relative as it is calculated in units of Hz^2 or Hz. The rf amplitude 
%   can be converted to T by dividing the resulting value by gamma. 
%   Correspondingly, The power can be converted to mT^2*s by dividing
%   the given value by gamma^2. Nonetheless, the absolute SAR is related to
%   the electric field, so the further scaling coeficient is both tx-coil-
%   dependent (e.g. depends on the coil design) and also subject-dependent 
%   (e.g. depends on the reference voltage). 

if nargin<2
    dt=1e-6; % for now default sampling rate is 1Mhz
end 

% resample the pulse to a resonable time array
nn=round(rf.shape_dur/dt);
t=((0:(nn-1))+0.5)*dt;
rfs=interp1(rf.t,rf.signal,t,'linear',0);
% TODO: avoid resampling of already sampled pulses (above)

rfs_sq=rfs.*conj(rfs);
total_energy=sum(rfs_sq)*dt;
if nargout>1 
    peak_pwr=max(rfs_sq);
    if nargout>2
        rf_rms=sqrt(total_energy/rf.shape_dur);
    end
end

