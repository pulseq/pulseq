function [grad] = scaleGrad(grad, scale)
% scaleGrad just scales the gradint with the scalar
    if strcmp(grad.type,'trap')
        grad.amplitude=grad.amplitude*scale;
        grad.area=grad.area*scale;
        grad.flatArea=grad.flatArea*scale;
    else
        grad.waveform=grad.waveform*scale;
        grad.first=grad.first*scale;
        grad.last=grad.last*scale;
    end
    if isfield(grad,'id')
        grad = rmfield(grad,'id'); 
    end
end

