function [ is_ok, text_error, total_dur ] = checkTiming( system, varargin )
%checkTiming(sys, objects, ...) 
%   Function checks whether timing of the specified evets is aligned to
%   the corresponding raster

    if (isempty(varargin))
        text_error=['empty or damaged block detected' ];
        is_ok=false;
        return;
    end
    total_dur=mr.calcDuration(varargin{:});
    is_ok=div_check(total_dur,system.blockDurationRaster);
    if (is_ok)
        text_error=[];
    else
        text_error=['total duration:' num2str(total_dur*1e6) 'us' ];
    end
    for i=1:length(varargin)
        e=varargin{i};
        if isnumeric(e) % special handling for blockDuration
            continue;
        end
        assert(isstruct(e), 'wrong format of the variable aguments, list of structures is expected');
        ok=true;
        if length(e)>1
            % for now this is only the case for arrays of extensions, but
            % we actually cannot check extensons anyway...
            continue;
        end
        if isfield(e, 'type') && (strcmp(e.type,'adc') || strcmp(e.type,'rf') || strcmp(e.type,'output'))
            raster=system.rfRasterTime;
        else
            raster=system.gradRasterTime;
        end
        if isfield(e, 'delay')
            if e.delay<-eps
                ok=false;
            end
            if ~div_check(e.delay,raster)
                ok=false;
            end
        end
        if isfield(e, 'duration')
            if ~div_check(e.duration,raster)
                ok=false;
            end
        end
        if isfield(e, 'dwell') % special case ADC
            if e.dwell<system.adcRasterTime || abs(round(e.dwell/system.adcRasterTime)*system.adcRasterTime - e.dwell)>1e-10
                ok=false;
            end
        end
        if isfield(e, 'type') && strcmp(e.type,'rf')
            % check time vector
            if ~div_check(e.shape_dur,system.rfRasterTime)
                ok=false;
            end
            if length(e.t)>=4
                rt=e.t/system.rfRasterTime;
                drt=diff(rt);
                if all(abs(drt(2:end)-drt(1))<1e-9/system.rfRasterTime) % 1ns -- /system.rfRasterTime is necessary because 'rt' is in RF raster units
                    % equal stepping case -- constant dwell time
                    e.dwell=e.t(2)-e.t(1); % add dummy dwell so that it can be logged in case of error
                    if ~div_check(e.dwell,system.adcRasterTime) % should we check against rfRasterTime or adcRasterTime?
                        ok=false;
                    end
                else
                    % "extended"-shape -- all points should be on RF raster edges
                    if any(abs(rf-round(rt))>1e-6)
                        ok=false; % TODO: add a meaninfull error message, for now it will look very cryptic
                    end
                end
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
            if isfield(e, 'duration')
                text_error = [text_error 'duration:' num2str(e.duration*1e6) 'us ' ];
            end
            if isfield(e, 'shape_dur')
                text_error = [text_error 'shape_dur:' num2str(e.shape_dur*1e6) 'us ' ];
            end            
            if isfield(e, 'dwell')
                text_error = [text_error 'dwell:' num2str(e.dwell*1e9) 'ns ' ];
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

