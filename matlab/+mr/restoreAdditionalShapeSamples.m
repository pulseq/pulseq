function [tt_chg, waveform_chg] = restoreAdditionalShapeSamples(tt,waveform,first,last,gradRasterTime,iBlock)
% restore shape: if we had a
% trapezoid converted to shape we have to find
% the "corners" and we can eliminate internal
% samples on the straight segments
% but first we have to restore samples on the
% edges of the gradient raster intervals
% for that we need the first sample
    max_abs=max(abs(waveform));
    odd_step1=[first 2*waveform'];
    odd_step2=odd_step1.*(mod(1:length(odd_step1),2)*2-1);
    waveform_odd_rest=(cumsum(odd_step2).*(mod(1:length(odd_step2),2)*2-1))';
    waveform_odd_interp=[first; 0.5*(waveform(1:end-1)+waveform(2:end)); last];
    if abs(waveform_odd_rest(end)-last)>2e-5*max_abs % what's the reasonable threshold? 
        blInfo='';
        if exist('iBlock')
            blInfo=['[block ' num2str(iBlock) '] '];
        end
        warning('mr:restoreShape',[blInfo 'Last restored point ' ...
                 'differs too much from the recorded last, skipping the shape restoration step; ' ... 
                 'deviation: ' num2str(abs(waveform_odd_rest(end)-last)) 'Hz/m (' num2str(abs(waveform_odd_rest(end)-last)/max_abs*100) '%%); ' ...
                 'No reason to panic, this is typical for spirals. ' ... 
                 'To disable this warning (assuming you know what you are doing) add \n' ...
                 '`warning(''OFF'', ''mr:restoreShape'')´ to your sequence-generating script.']); 
        tt_chg=[0 tt' tt(end)+gradRasterTime/2];
        waveform_chg = [first waveform' last];
        return;
    end
    %figure; plot([0,10e-6+grad.t'],waveform_odd_rest-waveform_odd_interp);
    waveform_odd_mask=abs(waveform_odd_rest-waveform_odd_interp)<=eps+2e-5*max_abs; % threshold ???
    waveform_odd=waveform_odd_interp.*waveform_odd_mask+waveform_odd_rest.*(1-waveform_odd_mask);

    % combine odd & even
    comb=[ 0 waveform' ; waveform_odd' ];
    waveform_os=comb(2:end)';

    tt_odd=(0:(length(waveform_odd_rest)-1))*gradRasterTime;
    tt_os=(0:(length(waveform_os)-1))*gradRasterTime*0.5;

    waveform_even_reint=0.5*(waveform_odd_rest(1:end-1)+waveform_odd_rest(2:end));

    maskChanges = abs([1; diff(waveform_os,2); 1])>1e-8;   % TRUE if values change
    waveform_chg = waveform_os(maskChanges)';                     % Elements without repetitions
    tt_chg=tt_os(maskChanges);
    %figure;plot(grad.tt,grad.waveform);hold on; plot(tt_chg,waveform_chg); plot(tt_chg,waveform_chg,'o');
end                                