function [tc ic]=calcRfCenter(rf)
%calcRfCenter Calculate the 'center' of the RF pulse
%   Returns the time point of the effective rotation calculated as the peak
%   of the RF amplitude for the shaped pulses and the center of the pulse
%   for the block pulses. Zeropadding in the RF pulse is considered as a
%   part of the shape. Delay field of the RF object is not taken into
%   account.
%
%   The main return value is the time point of the center, the optional
%   return value is the corresponding position in the rf envelope array.
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

    % we detect the excitation peak and if it is a plato we take its center
    rfmax=max(abs(rf.signal));
    ipeak=find(abs(rf.signal)>=rfmax*0.99999);
    tc=(rf.t(ipeak(1))+rf.t(ipeak(end)))/2;
    ic=ipeak(round(end/2));
    
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