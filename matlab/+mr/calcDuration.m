function duration=calcDuration(varargin)
%calcDuration Calculate the duration of an event, set of events, or block.
%
%   PURPOSE
%     Compute the time, in seconds, that one or more sequence events
%     occupy. Used pervasively when laying out blocks: TE/TR delay
%     calculations, gradient alignment, and consistency checks between
%     block timing and event timing all call this function.
%
%   SIGNATURES
%     dur = mr.calcDuration(event)              % single event struct
%     dur = mr.calcDuration(e1, e2, ...)        % max over multiple events
%     dur = mr.calcDuration(block)              % a block struct from seq.getBlock(n)
%     dur = mr.calcDuration({e1, e2, ...})      % cell array of events (auto-unwrapped)
%
%     For multiple events the return value is the maximum of each event's
%     duration (events in a block run concurrently, not sequentially).
%     Unknown event types are silently ignored and contribute 0.
%
%   INPUTS
%     varargin  [required]  One or more arguments, each one of:
%                             - event struct with a .type field set to one of
%                               'delay', 'rf', 'grad', 'trap', 'adc', 'output',
%                               'trigger'
%                             - a block struct (i.e. has a .rf field; typically
%                               obtained from mr.Sequence/getBlock). Only a single
%                               block struct may be passed.
%                             - a numeric scalar interpreted as a blockDuration
%                               field (seconds). Used internally when block2events
%                               expands a block struct; rarely passed directly by
%                               user code.
%
%   OUTPUT
%     duration  double, seconds. The maximum event duration encountered.
%               Returned as 0 if all arguments are unknown event types.
%
%   ERRORS
%     - MATLAB:assertion:failed: a numeric blockDuration argument is
%       smaller than the maximum event duration encountered before it
%       in argument order. Indicates an inconsistent block whose declared
%       blockDuration is shorter than its longest event.
%     - 'Only a single block structure can be added' (from mr.block2events):
%       more than one argument was a block struct.
%     - 'Index exceeds array bounds' (from mr.block2events): called with
%       no arguments. Always pass at least one event.
%
%   NOTES
%     - Per-event duration formulas (all in seconds):
%         delay   : event.delay
%         rf      : event.delay + event.shape_dur + event.ringdownTime
%         grad    : event.delay + event.shape_dur (arbitrary gradient)
%         trap    : event.delay + event.riseTime + event.flatTime + event.fallTime
%         adc     : event.delay + event.numSamples*event.dwell + event.deadTime
%         output  : event.delay + event.duration
%         trigger : event.delay + event.duration
%     - For 'rf' and 'grad' events the function relies on a precomputed
%       event.shape_dur field; it does not re-derive duration from the
%       waveform samples.
%     - For 'adc', the deadTime addend is the post-acquisition dead time
%       captured at construction (adc.deadTime), not the pre-acquisition
%       delay (which is already included via adc.delay).
%     - Unknown .type values are silently skipped. A struct without a
%       .type field whose value matches no case will not raise.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     Nx = 256;  fov = 256e-3;  deltak = 1/fov;
%     gx    = mr.makeTrapezoid('x', sys, 'FlatArea', Nx*deltak, 'FlatTime', 6.4e-3);
%     gxPre = mr.makeTrapezoid('x', sys, 'Area', -gx.area/2, 'Duration', 1e-3);
%     adc   = mr.makeAdc(Nx, sys, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
%     % Longest event in this readout block
%     dur   = mr.calcDuration(gxPre, gx, adc);
%     % TE delay rounded up to gradient raster
%     TE = 10e-3;
%     delayTE = ceil((TE - mr.calcDuration(gx)/2)/sys.gradRasterTime)*sys.gradRasterTime;
%
%   SEE ALSO
%     mr.block2events, mr.makeDelay, mr.Sequence/getBlock,
%     mr.Sequence/addBlock

duration=0;

% Convert block structure to cell array of events
varargin=mr.block2events(varargin);

% Loop over events and calculate maximum duration
for i=1:length(varargin)
    events = varargin{i};
    for j=1:length(events)
        event=events(j);
        if isnumeric(event) % this is the "blockDuration" field
            assert(duration<=event);
            duration=event;
            continue;
        end
        switch event.type
            case 'delay'
                duration=max(duration, event.delay);
            case 'rf'
                duration=max(duration, event.delay + event.shape_dur + event.ringdownTime);  % rf has now a field called 'shape_dur' because it is impossible to calculate duration here in a general case...
            case 'grad'
%                 % duration=max(duration, event.t(end) + event.delay ); 
%                 % MZ: we need to increase the last timestamp by one gradient 
%                 % raster time because otherwise the duration of the gradint 
%                 % containing one sample will be zero
%                 % however, we do not have access to the gradient raster time
%                 % here, so we opt of a hack, but it will actually fail for the
%                 % gradint containg one sample...
%                 % event.t(2) - event.t(1) gives us one gradient raster time
%                 duration=max(duration, event.t(end) + event.t(2) - event.t(1) + event.delay ); 
                duration=max(duration, event.delay+event.shape_dur); % shaped gradient has now a field called 'shape_dur' because it is impossible to calculate duration here in a general case...
            case 'adc'
                duration=max(duration, event.delay + ...
                             event.numSamples*event.dwell + event.deadTime);
            case 'trap'
                duration=max(duration, event.delay + event.riseTime + ...
                             event.flatTime + event.fallTime);
    %             duration=max(duration, event.riseTime+ event.flatTime +event.fallTime);
            case {'output','trigger'} 
                duration=max(duration, event.delay + event.duration);
        end
    end
end

end
