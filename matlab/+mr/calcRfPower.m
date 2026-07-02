function [total_energy, peak_pwr, rf_rms]=calcRfPower(rf, dt)
%calcRfPower Compute relative energy, peak power, and RMS amplitude of an RF pulse.
%
%   PURPOSE
%     Computes the relative energy, peak power, and RMS B1 amplitude of an
%     RF pulse by resampling its waveform onto a uniform time grid and
%     integrating |rf|^2. Returns 'relative' quantities in units of Hz^2
%     and Hz, not absolute SI power, because Pulseq RF amplitudes are
%     stored in Hz (gamma-scaled). Used for pulse-level sanity checks and
%     as a building block for the sequence-level mr.Sequence/calcRfPower.
%
%   SIGNATURES
%     total_energy                     = mr.calcRfPower(rf)         % default dt=1 us
%     [total_energy, peak_pwr]         = mr.calcRfPower(rf)
%     [total_energy, peak_pwr, rf_rms] = mr.calcRfPower(rf)
%     [...]                            = mr.calcRfPower(rf, dt)     % override sampling step
%
%     The two arguments must be passed positionally and in order; they
%     cannot be given as name/value pairs. Outputs beyond nargout are not
%     computed.
%
%   INPUTS
%     rf  [required]  struct, RF event struct from mr.makeSincPulse,
%                     mr.makeBlockPulse, mr.makeGaussPulse, mr.makeArbitraryRf,
%                     mr.makeSLRpulse, or mr.makeAdiabaticPulse. Must have
%                     fields .t (seconds), .signal (complex Hz), and .shape_dur
%                     (seconds).
%     dt  [optional]  double, resampling step in seconds. Default: 1e-6.
%
%   OUTPUT
%     total_energy  double, Hz (= Hz^2 * s), integral of |rf|^2 over the pulse duration
%     peak_pwr      double, Hz^2, maximum of |rf|^2 over the resampled waveform
%     rf_rms        double, Hz, RMS B1 amplitude, sqrt(total_energy/rf.shape_dur)
%
%   NOTES
%     - Outputs are 'relative': amplitude is in Hz, so total_energy is in
%       Hz^2 * s and peak_pwr in Hz^2. To convert to SI: divide rf_rms by
%       sys.gamma to get Tesla; divide total_energy by sys.gamma^2 to get
%       Tesla^2 * s. Absolute SAR requires further scaling by tx-coil and
%       subject-specific factors (reference voltage, coil design) and is
%       not computed here.
%     - The pulse is resampled onto a uniform grid of step dt with bin
%       midpoints at ((0:nn-1)+0.5)*dt where nn = round(rf.shape_dur/dt).
%       Samples outside rf.t are filled with 0 by linear extrapolation.
%       For pulses already sampled on a uniform raster (sinc, Gauss, SLR,
%       arbitrary, adiabatic) the default dt=1e-6 oversamples mildly; for
%       block pulses (rf.t has only the two endpoints) resampling is
%       required.
%     - peak_pwr depends on dt because the resampled waveform may not hit
%       the original peak exactly; for narrow pulses, reduce dt to tighten
%       the estimate.
%     - rf_rms uses rf.shape_dur, which excludes any rf.delay or
%       post-pulse ringdown. For SAR-relevant duty cycle over a sequence
%       use mr.Sequence/calcRfPower instead.
%     - rf.freqOffset, rf.phaseOffset, and rf.freqPPM are ignored: the
%       calculation works on |rf.signal|^2, which is invariant under
%       frequency or phase modulation.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     rf  = mr.makeSincPulse(pi/2, 'Duration', 3e-3, ...
%                            'SliceThickness', 3e-3, 'apodization', 0.5, ...
%                            'timeBwProduct', 4, 'system', sys);
%     [total_energy, peak_pwr, rf_rms] = mr.calcRfPower(rf);
%     fprintf('energy = %.3g Hz, peak = %.3g Hz^2, rms = %.3g Hz\n', ...
%             total_energy, peak_pwr, rf_rms);
%
%   SEE ALSO
%     mr.calcRfBandwidth, mr.calcRfCenter, mr.simRf, mr.makeSincPulse,
%     mr.makeBlockPulse, mr.makeGaussPulse, mr.makeArbitraryRf,
%     mr.makeSLRpulse, mr.makeAdiabaticPulse

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

