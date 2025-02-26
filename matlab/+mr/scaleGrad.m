function [grad] = scaleGrad(grad, scale, system)
% scaleGrad just scales the gradint with the scalar
% if the optional parameter 'system' is provided the result is schecked
% against the system limits (amolityde and slew rate)
    if strcmp(grad.type,'trap')
        grad.amplitude=grad.amplitude*scale;
        grad.area=grad.area*scale;
        grad.flatArea=grad.flatArea*scale;
        if nargin>2
            if system.maxGrad<abs(grad.amplitude)
                error("mr.scaleGrad: maximum amplitude exceeded (%g %%)",100*abs(grad.amplitude)/system.maxGrad);
            end
            if abs(grad.amplitude)>eps && system.maxSlew<abs(grad.amplitude)/min(grad.riseTime,grad.fallTime)
                error("mr.scaleGrad: maximum slew rate exceeded (%g %%)",100*abs(grad.amplitude)/min(grad.riseTime,grad.fallTime)/system.maxSlew);
            end
        end
    else
        grad.waveform=grad.waveform*scale;
        grad.first=grad.first*scale;
        grad.last=grad.last*scale;
        if nargin>2
            if system.maxGrad<max(abs(grad.waveform))
                error("mr.scaleGrad: maximum amplitude exceeded (%g %%)",100*max(abs(grad.waveform))/system.maxGrad);
            end
            if max(abs(grad.waveform))>eps
                grad_max_abs_slew=max(abs(diff(grad.waveform)./diff(grad.tt)));
                if system.maxSlew<grad_max_abs_slew
                    error("mr.scaleGrad: maximum slew rate exceeded (%g %%)",100*grad_max_abs_slew/system.maxSlew);
                end
            end
        end
    end
    if isfield(grad,'id')
        grad = rmfield(grad,'id'); 
    end
end

