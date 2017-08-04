function duration=calcDuration(varargin)
%calcDuration Calculate the duration of an event or block.
%   dur=calcDuration(e1,e2,...) Determine the maximum duration of the
%   provided events
%
%   dur=calcDuration(b) Calculate the duration of the block structure

% Convert block structure to cell array of events
varargin=mr.block2events(varargin);

% Loop over events and calculate maximum duration
duration=0;
for i=1:length(varargin)
    event = varargin{i};
    switch event.type
        case 'delay'
            duration=max(duration,event.delay);
        case 'rf'
            duration=max(duration,event.t(end)+event.deadTime+event.ringdownTime);
        case 'grad'
            duration=max(duration,event.t(end));
        case 'adc'
            adcTime = event.delay + event.numSamples*event.dwell + event.deadTime;
            duration=max(duration,adcTime);
        case 'trap'
            duration=max(duration,event.riseTime+event.flatTime+event.fallTime);
    end
end

end