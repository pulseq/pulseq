function [bw,fc,spectrum,f,rfs,t]=calcRfBandwidth(rf, cutoff, df, dt)
%calcRfBandwidth Compute the bandwidth and spectrum of an RF pulse.
%
%   PURPOSE
%     Estimates the bandwidth of an RF pulse using a simple FFT under the
%     low-angle (small-tip) approximation. Also returns the pulse's center
%     frequency, complex spectrum, and the corresponding frequency and
%     resampled time axes. Typically used for plotting excitation profile
%     magnitude vs. frequency and for sanity-checking RF pulses produced by
%     mr.makeSincPulse, mr.makeGaussPulse, mr.makeBlockPulse,
%     mr.makeArbitraryRf, mr.makeSLRpulse, and mr.makeAdiabaticPulse.
%
%   SIGNATURES
%     bw                            = mr.calcRfBandwidth(rf)                    % default cutoff=0.5 (FWHM), df=10 Hz, dt=1 us
%     [bw,fc]                       = mr.calcRfBandwidth(rf)                    % also return center frequency
%     [bw,fc,spectrum,f]            = mr.calcRfBandwidth(rf)                    % also return spectrum and frequency axis
%     [bw,fc,spectrum,f,rfs,t]      = mr.calcRfBandwidth(rf)                    % also return resampled RF waveform and time axis
%     [...] = mr.calcRfBandwidth(rf, cutoff)                                    % override fractional cutoff
%     [...] = mr.calcRfBandwidth(rf, cutoff, df)                                % override spectral resolution
%     [...] = mr.calcRfBandwidth(rf, cutoff, df, dt)                            % override time sampling step
%
%     Bandwidth is the frequency width at which |spectrum| drops to
%     cutoff*max(|spectrum|), with linear interpolation at the first and
%     last crossings. All four positional arguments must be given in order;
%     they cannot be passed as name/value pairs.
%
%   INPUTS
%     rf      [required]  struct, RF event struct from mr.makeSincPulse,
%                         mr.makeBlockPulse, mr.makeGaussPulse, mr.makeArbitraryRf,
%                         mr.makeSLRpulse, or mr.makeAdiabaticPulse. Must have
%                         fields .t (seconds), .signal (complex Hz),
%                         .freqOffset (Hz), .freqPPM, .phaseOffset (radians),
%                         and .center (seconds).
%     cutoff  [optional]  double, fractional threshold for bandwidth measurement,
%                         dimensionless in (0,1]. 0.5 gives FWHM. Default: 0.5.
%     df      [optional]  double, spectral resolution in Hz. Default: 10.
%     dt      [optional]  double, time sampling step in seconds. Default: 1e-6.
%
%   OUTPUT
%     bw        double, Hz, bandwidth at the given cutoff threshold
%     fc        double, Hz, center frequency (midpoint of the two threshold crossings)
%     spectrum  complex vector, small-tip excitation response, scaled by
%                 sin(2*pi*dt*s_ref)/s_ref with s_ref = |spectrum(fc)|
%     f         double vector, Hz, frequency axis, step df, length round(1/df/dt)
%     rfs       complex vector, RF waveform resampled onto the uniform
%                 time grid t, with rf.freqOffset, rf.freqPPM, and
%                 rf.phaseOffset folded into the phase modulation
%     t         double vector, seconds, uniform time axis centered on
%                 rf.center, step dt, span 1/df
%
%   NOTES
%     - Computed from the RF envelope's Fourier transform (small-tip
%       equivalent). The Bloch-simulated bandwidth grows modestly at
%       large flip angles. For exact slice-profile analysis use mr.simRf.
%     - If rf.freqPPM is non-zero, the function calls mr.opts() to fetch
%       sys.gamma and sys.B0 for the PPM-to-Hz conversion, and emits a
%       warning reminding the caller to set system defaults via
%       mr.opts('setAsDefault', true). The system defaults are used
%       silently otherwise.
%     - The FFT length is round(1/df/dt) (default 1e5 points). Reducing df
%       or dt below their defaults increases compute time quadratically in
%       the product, not in either one alone.
%     - The spectrum is normalized so that |spectrum(fc)| approximates the
%       on-resonance small-tip excitation amplitude rather than an
%       unnormalized FFT magnitude.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     rf  = mr.makeSincPulse(pi/2, 'Duration', 3e-3, ...
%                            'SliceThickness', 3e-3, 'apodization', 0.5, ...
%                            'timeBwProduct', 4, 'system', sys);
%     [bw, fc, spectrum, f] = mr.calcRfBandwidth(rf);
%     fprintf('BW = %.1f Hz, center = %.2f Hz\n', bw, fc);
%     figure; plot(f, abs(spectrum)); xlim(3*[-bw bw]);
%     xlabel('Frequency, Hz'); title('Excitation pulse profile');
%
%   SEE ALSO
%     mr.calcRfCenter, mr.calcRfPower, mr.simRf, mr.makeSincPulse,
%     mr.makeGaussPulse, mr.makeBlockPulse, mr.makeArbitraryRf,
%     mr.makeSLRpulse, mr.makeAdiabaticPulse
%

if nargin<2
    cutoff=0.5;
end
    
if nargin<3
    df=10; % spectral resolution in Hz
end

if nargin<4
    dt=1e-6; % for now default sampling rate is 1Mhz, it's probably too high 
end

if abs(rf.freqPPM)>eps
    warning('mr.calcRfBandwidth((): relying on the system properties, like B0 and gamma, stored in the global environment by callimg mr.lims(''setAsDefault'',true)');
    sys=mr.opts();
    full_freqOffset=rf.freqOffset+rf.freqPPM*1e-6*sys.gamma*sys.B0;
else
    full_freqOffset=rf.freqOffset;
end

tc=rf.center;

% resample the pulse to a resonable time array
nn=round(1/df/dt);
t=(-floor(nn/2):ceil(nn/2)-1)*dt;

rfs=interp1(rf.t-tc,rf.signal.*exp(1i*(rf.phaseOffset+2*pi*full_freqOffset*rf.t)),t,'linear',0);
spectrum=fftshift(fft(fftshift(rfs)));
f=(-floor(nn/2):ceil(nn/2)-1)*df;

w1=mr.aux.findFlank(f,spectrum,cutoff);
w2=mr.aux.findFlank(f(end:-1:1),spectrum(end:-1:1),cutoff);

bw=w2-w1;
fc=(w2+w1)/2;

% coarse STE scaling -- we normalize to the max of the spectrum, this works
% better with frequency-shifted pulses than the abs(sum(shape)) -- the
% 0-frequency response; yes, we could take the spectrum at fc but this
% would have a problem for non-symmetric pulses...
%s_ref=max(abs(spectrum));
s_ref=interp1(f,abs(spectrum),fc);
spectrum=sin(2*pi*dt*s_ref)*spectrum/s_ref;

end

