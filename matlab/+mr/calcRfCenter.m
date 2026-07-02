function [tc, ic, fi]=calcRfCenter(rf)
%calcRfCenter Calculate the effective center time of an RF pulse.
%
%   PURPOSE
%     Returns the time point of the effective rotation of an RF pulse,
%     used to align echo timing, set off-resonance phase compensation, and
%     compute TE/TI delays. For shaped pulses this is the peak of the RF
%     amplitude; for block pulses it is the midpoint of the constant
%     plateau. Zero-padding in rf.signal is treated as part of the shape.
%     The rf.delay field is NOT included in the returned time, so callers
%     that need an absolute time within a block typically compute
%     rf.delay + mr.calcRfCenter(rf).
%
%   SIGNATURES
%     tc            = mr.calcRfCenter(rf)    % time of center, seconds
%     [tc, ic]      = mr.calcRfCenter(rf)    % also integer sample index
%     [tc, ic, fi]  = mr.calcRfCenter(rf)    % also fractional offset
%
%     If rf has a .center field (set by the modern RF constructors
%     mr.makeSincPulse, mr.makeBlockPulse, mr.makeGaussPulse,
%     mr.makeArbitraryRf, mr.makeSLRpulse, mr.makeAdiabaticPulse), tc is
%     taken directly from rf.center and ic is the nearest index in rf.t.
%     Otherwise tc is computed from the amplitude peak of rf.signal.
%
%   INPUTS
%     rf  [required]  struct, RF event struct. Must have fields .t (seconds, time
%                     axis on the RF raster) and .signal (complex Hz, waveform).
%                     If field .center (seconds) is present it is used directly;
%                     otherwise the function detects the peak of abs(rf.signal).
%                     Typically produced by mr.makeSincPulse, mr.makeBlockPulse,
%                     mr.makeGaussPulse, mr.makeArbitraryRf, mr.makeSLRpulse,
%                     or mr.makeAdiabaticPulse.
%
%   OUTPUT
%     tc   double, seconds, time of the RF center relative to the start
%            of the RF shape (rf.delay not included)
%     ic   integer, 1-based index into rf.t / rf.signal of the sample
%            nearest to tc
%     fi   double, dimensionless, fractional offset in [-0.5, 0.5] from
%            rf.t(ic) toward the previous (negative) or next (positive)
%            sample, normalized by the local raster step. 0 when tc lies
%            exactly on rf.t(ic) (within 1 ns)
%
%   NOTES
%     - When rf.center is absent, the peak detector treats samples within
%       0.001% of max(abs(rf.signal)) as part of the same plateau and
%       returns the midpoint. This is what makes block pulses (constant
%       amplitude) yield a center at the middle of the pulse rather than
%       at the first sample.
%     - rf.delay is intentionally excluded from tc. If a sequence places
%       an RF event after a delay, the absolute time within the block is
%       rf.delay + tc.
%     - Returned tc is on the RF raster (rf.t edges); fi captures the
%       sub-raster offset when rf.center does not coincide with a sample.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s', ...
%                   'rfDeadTime', 100e-6);          % bumps rf.delay to 100 us
%     rf  = mr.makeSincPulse(pi/2, 'Duration', 3e-3, ...
%                            'SliceThickness', 3e-3, 'system', sys);
%     % off-resonance phase compensation for an off-center slice:
%     rf.freqOffset = sys.gamma * 1e-3 * 5e-3;        % 5 mm offset
%     rf.phaseOffset = -2*pi*rf.freqOffset*mr.calcRfCenter(rf);
%     % absolute time of the RF center within its block:
%     tCenter = rf.delay + mr.calcRfCenter(rf);       % 100 us + 1.5 ms = 1.6 ms
%
%   SEE ALSO
%     mr.calcRfBandwidth, mr.calcRfPower, mr.makeSincPulse,
%     mr.makeBlockPulse, mr.makeGaussPulse, mr.makeArbitraryRf,
%     mr.makeSLRpulse, mr.makeAdiabaticPulse
%

%     % detect zero-padding
%     last=length(rf.signal);
%     for first=1:last
%         if abs(rf.signal(first))>eps
%             break;
%         end
%     end
%     for last=last:-1:first
%         if abs(rf.signal(last))>eps
%             break;
%         end
%     end

%     rfmax=max(abs(rf.signal(first:last)));
%     ipeak=find(abs(rf.signal(first:last))>=rfmax-eps);

    if isfield(rf,'center')
        tc=rf.center;
        [~,ic]=min(abs(rf.t-tc));
    else

        % we detect the excitation peak and if it is a plato we take its center
        rfmax=max(abs(rf.signal));
        ipeak=find(abs(rf.signal)>=rfmax*0.99999);
        tc=(rf.t(ipeak(1))+rf.t(ipeak(end)))/2;
        ic=ipeak(round(end/2));
    end

    ft=tc-rf.t(ic);
    if ic<length(rf.t) && ft>1e-9 % 1 ns
        fi=ft/(rf.t(ic+1)-rf.t(ic));
    elseif ic>1 && ft<-1e-9 % -1 ns
        fi=ft/(rf.t(ic)-rf.t(ic-1));
    else
        fi=0;
    end

%     % detect the excitation peak (this code is far from being ideal...)
%     rfmin=min(abs(rf.signal(first:last))); % pure max check fails for the block pulse!!!
%     [rfmax,ic]=max(abs(rf.signal(first:last))); 
%     if (rfmax-rfmin)<=eps
%         ic=round((last-first+1)/2); % we take the center of the pulse for block pulses
%         tc=(rf.t(first)+rf.t(last))/2;
%     else
%         tc=rf.t(first-1+ic);
%     end
end