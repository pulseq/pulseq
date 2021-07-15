function duration=calcDuration(varargin)
%calcDuration Calculate the duration of an event or block.
%   dur=calcDuration(e1,e2,...) Determine the maximum duration of the
%   provided events
%
%   dur=calcDuration(b) Calculate the duration of the block structure

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
