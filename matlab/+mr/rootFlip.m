function [rfw,bw] = rootFlip(b,d1,flip,tbw,peakConstraint,rfInd)

% root flip the envelope to minimize peak B1/SAR
b = b./max(abs(freqz(b)));
b = b*sin(flip/2 + atan(d1*2)/2); % scale to target flip angle
rfOrigMax = max(abs(b2rf(b))); % for checking constraint
r = roots(b); % calculate roots of single-band beta polynomial
% sort the roots by phase angle
[~,idx] = sort(angle(r));
r = r(idx);r = r(:);
r = leja_fast(r);

candidates = abs(1-abs(r)) > 0.001 & abs(angle(r)) < tbw/length(b)*pi;
maxRF = Inf;
maxPhsSlope = Inf;
if exist('rfInd','var')
    iiMin = min(rfInd);iiMax = max(rfInd);
else
    iiMin = 1;iiMax = 2^sum(candidates);
end
if exist('peakConstraint','var')
    if peakConstraint > 0
        % get a filter that smooths out the beta phase profile with a target
        % FWHM
        pbMask = abs(fftshift(fft(fftshift(b)))) > 0.1;% threshold the beta amp profile to get a passband mask        
    end
else
    peakConstraint = -1; % switch for no smoothing
end

for ii = iiMin:iiMax
    
    % convert this ii to binary pattern
    doFlipStr = fliplr(dec2bin(ii-1));
    doFlip = zeros(1,length(doFlipStr));
    for jj = 1:length(doFlipStr)
       doFlip(jj) = str2num(doFlipStr(jj)); 
    end
    % add zeros for the rest of the passbands
    doFlip = [doFlip(:); zeros(sum(candidates)-length(doFlip),1)];
    
    % embed in all-roots vector
    tmp = zeros(length(b)-1,1);
    tmp(candidates) = doFlip;
    doFlip = tmp;
    % flip those indices
    rt = r(:);
    rt(doFlip == 1) = conj(1./rt(doFlip == 1));
    % get polynomial coefficients back from flipped roots
    tmp = poly(rt);
    % important: normalize before modulating! doesn't work
    % for negative shifted slices otherwise
    tmp = tmp./max(abs(freqz(tmp)));
    % store the rf
    bt = tmp(:)*sin(flip/2 + atan(d1*2)/2);
    rft = b2rf(bt);
    if peakConstraint < 0
        % typical case: take RF with min peak amplitude
        if max(abs(rft)) < maxRF
            maxRF = max(abs(rft));
            rfw = rft;
            bw = bt;
        end
    else
        % OR, among pulses that meet peak constraint, take the one with 
        % lowest max phase slope in the passband
        if max(abs(rft)) < peakConstraint*rfOrigMax
            phs = unwrap(angle(fftshift(fft(fftshift(bt)))));
            phs = phs(pbMask > 0);
            if max(abs(diff(phs))) < maxPhsSlope 
                maxRF = max(abs(rft));
                maxPhsSlope = max(abs(diff(phs))) 
                rfw = rft;
                bw = bt;
            end
        end
    end
    
end

