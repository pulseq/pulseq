function [ is_ok, text_error ] = checkTiming( system, varargin )
%checkTiming(sys, objects, ...) 
%   Function checks whether timing of the specified evets are aligned to
%   the corresponding raster

    total_dur=mr.calcDuration(varargin{:});
    is_ok=div_check(total_dur,system.gradRasterTime);
    if (is_ok)
        text_error=[];
    else
        text_error=['total duration:' num2str(total_dur*1e6) 'us' ];
    end
    for i=1:length(varargin)
        e=varargin{i};
        ok=true;
        if isfield(e, 'type') && (strcmp(e.type,'adc') || strcmp(e.type,'rf'))
            raster=system.rfRasterTime;
        else
            raster=system.gradRasterTime;
        end
        if isfield(e, 'delay')
            if ~div_check(e.delay,raster)
                ok=false;
            end
        end
        if isfield(e, 'type') && strcmp(e.type,'trap')
            if ~div_check(e.riseTime, system.gradRasterTime) || ~div_check(e.flatTime, system.gradRasterTime) || ~div_check(e.fallTime, system.gradRasterTime)
                ok=false;
            end
        end
        if ~ok
            is_ok=false;
            if ~isempty(text_error)
                text_error = [text_error ' '];
            end
            text_error = [text_error '[ '];
            if isfield(e, 'type')
                text_error = [text_error 'type:' e.type ' ' ];
            end
            if isfield(e, 'delay')
                text_error = [text_error 'delay:' num2str(e.delay*1e6) 'us ' ];
            end
            if isfield(e, 'type') && strcmp(e.type,'trap')
                text_error = [text_error 'riseTime:' num2str(e.riseTime*1e6) 'us flatTime:' num2str(e.flatTime*1e6) 'us fallTime:' num2str(e.fallTime*1e6) 'us '];
            end
            text_error = [text_error ']'];
        end
    end
end

function out = div_check(a, b)
%  checks wheher a can be divided by b to an accuracy of 1e-9
    c = a / b;
    out = (abs( c - round(c) ) < 1e-9);
end

