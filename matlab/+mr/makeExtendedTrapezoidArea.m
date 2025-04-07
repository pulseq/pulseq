function [grad, times, amplitudes] = makeExtendedTrapezoidArea(channel, Gs, Ge, A, sys)
% make shortest possible extended trapezoid with a given area
% which starts and ends (optionally) as non-zero gradient values

% This is an upgraded version of the makeExtendedTrapezoidArea function
% modified by Mehmet Emin Öztürk in 2024,
% who is a Master student from Bilkent University, Turkey.

SR=sys.maxSlew*0.99; % otherwise we run into rounding errors during resampling

Tp=0;
obj1=@(x) (A-testGA(x,0,SR,sys.gradRasterTime,Gs,Ge)).^2;

% we seem to need a multistart...
[Gp(1),obj1val(1),exitf(1)] = fminsearch(obj1,-sys.maxGrad);
[Gp(2),obj1val(2),exitf(2)] = fminsearch(obj1,0);
[Gp(3),obj1val(3),exitf(3)] = fminsearch(obj1,sys.maxGrad);
[~,imin]=min(obj1val);
Gp=Gp(imin);
obj1val=obj1val(imin);
exitf=exitf(imin);

if  obj1val>1e-3 || ... % the search did not converge
    abs(Gp)> sys.maxGrad
    Gp=sys.maxGrad*sign(Gp);
    obj2=@(x) (A-testGA(Gp,x,SR,sys.gradRasterTime,Gs,Ge)).^2;
    [T,obj2val,exitf] = fminsearch(obj2,0);
    Tp=ceil(T/sys.gradRasterTime)*sys.gradRasterTime;
    
    % fix the ramps
    Tru=ceil(abs(Gp-Gs)/SR/sys.gradRasterTime)*sys.gradRasterTime;
    Trd=ceil(abs(Gp-Ge)/SR/sys.gradRasterTime)*sys.gradRasterTime;
    obj3=@(x) (A-testGA1(x,Tru,Tp,Trd, Gs,Ge)).^2;

    %obj3=@(x) (A-testGA(x,Tp,SR,sys.gradRasterTime,Gs,Ge)).^2;
    [Gp,obj3val,exitf] = fminsearch(obj3,Gp);
    assert(obj3val<1e-3); % did the final search converge?
end

if Tp<=0
    Tp=-sys.gradRasterTime;
    SR = sys.maxSlew*0.99;
    area = 0;
    while abs(area-A)>1e-3
        Tp = Tp+sys.gradRasterTime;
        obj3=@(x) (A-testGA1(x,ceil(abs(x-Gs)/SR/sys.gradRasterTime)*sys.gradRasterTime, ...
            Tp,ceil(abs(x-Ge)/SR/sys.gradRasterTime)*sys.gradRasterTime, ...
            Gs,Ge)).^2;
        [Gp,~,~] = fminsearch(obj3,Gp);
        if abs(Gp)>sys.maxGrad
            Gp = sys.maxGrad*sign(Gp); 

        end
        Tru=ceil(abs(Gp-Gs)/SR/sys.gradRasterTime)*sys.gradRasterTime;
        Trd=ceil(abs(Gp-Ge)/SR/sys.gradRasterTime)*sys.gradRasterTime;

        area=testGA1(Gp, Tru, Tp, Trd, Gs, Ge);
    end
end
assert(Tp>=0); % this was a nasty error condition when some of the above did not converge

if Tp>0
    times=cumsum([0 Tru Tp Trd]);
    amplitudes=[Gs Gp Gp Ge];
else
    Tru=ceil(abs(Gp-Gs)/SR/sys.gradRasterTime)*sys.gradRasterTime;
    Trd=ceil(abs(Gp-Ge)/SR/sys.gradRasterTime)*sys.gradRasterTime;
    
    if Trd>0
        if Tru>0
            times=cumsum([0 Tru Trd]);
            amplitudes=[Gs Gp Ge];
        else
            times=cumsum([0 Trd]);
            amplitudes=[Gs Ge];
        end
    else
        times=cumsum([0 Tru]);
        amplitudes=[Gs Ge];
    end
end


grad=mr.makeExtendedTrapezoid(channel,'system',sys,'times',times, 'amplitudes', amplitudes);
grad.area=testGA1(Gp, Tru, Tp, Trd, Gs, Ge);
assert(abs(grad.area-A)<1e-3);
assert(SR <= sys.maxSlew)
assert(all(abs(amplitudes)<=sys.maxGrad))
end

function ga = testGA(Gp, Tp, SR, dT, Gs, Ge)
Tru=ceil(abs(Gp-Gs)/SR/dT)*dT;
Trd=ceil(abs(Gp-Ge)/SR/dT)*dT;
%ga=0.5*Tru*(Gp+Gs)+Gp*Tp+0.5*(Gp+Ge)*Trd;
ga=testGA1(Gp, Tru, Tp, Trd, Gs, Ge);
end

function ga = testGA1(Gp, Tru, Tp, Trd, Gs, Ge)
ga=0.5*Tru.*(Gp+Gs)+Gp.*Tp+0.5*(Gp+Ge).*Trd;
end


function ga = testGA2(Gp, SR, Gs, Ge, dT)
Tru=ceil(abs(Gp-Gs)/SR/dT)*dT;
Trd=ceil(abs(Gp-Ge)/SR/dT)*dT;
ga=0.5*Tru*(Gp+Gs)+0.5*(Gp+Ge)*Trd;
end