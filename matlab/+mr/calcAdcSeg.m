function [adcSegments,adcSamplesPerSegment] = calcAdcSeg(numSamples,dwell,system,mode)
%mr.calcAdcSeg : Calculate splitting of the ADC in segments
%    On some scanners, notably Siemens, ADC objects that exceed a certain 
%    sample length (8192 samples on Siemens) should be splittable to N 
%    equal parts, each of which aligned to the gradient raster. Each 
%    segment, however, needs to have the number of samples smaller than 
%    system.adcSamplesLimit and divisible by system.adcSamplesDivisor to be
%    executable on the scanner. The optional parameter mode can be either
%    'shorten' or 'lengthen'.

if system.adcSamplesLimit<=0
    adcSamplesPerSegment=numSamples;
    adcSegments=1;
    return;
end

if ~exist('mode', 'var')
    mode='shorten';
end

switch mode
    case 'shorten'
        iMode=-1;
    case 'lenghen'
        iMode=1;
    otherwise
        error('In mr.calcAdcSeg(...,mode) mode should be either ''shorten'' or ''lenghen''');
end

t_eps=1e-9; % TODO: shift it to the system parameters???

iGR=round(system.gradRasterTime/system.adcRasterTime);
assert(abs(system.gradRasterTime/system.adcRasterTime-iGR)<t_eps);

iDwell=round(dwell/system.adcRasterTime);
assert(abs(dwell/system.adcRasterTime-iDwell)<t_eps);

iCommon=lcm(iGR,iDwell); % least common multiplier
samplesStep=iCommon/iDwell;

% Siemens-specific: number of samples should be divisible by system.adcSamplesDivisor
gcd_adcDiv=gcd(samplesStep,system.adcSamplesDivisor);
if gcd_adcDiv~=system.adcSamplesDivisor
    samplesStep=samplesStep*system.adcSamplesDivisor/gcd_adcDiv;
end

if mode=='shorten'
    numSamplesStepped=floor(numSamples/samplesStep); 
else
    numSamplesStepped=ceil(numSamples/samplesStep); 
end

while numSamplesStepped>0 && numSamplesStepped<2*numSamples/samplesStep
    adcSegmentFactors=factor(numSamplesStepped);
    adcSegments=1;        
    if(length(adcSegmentFactors)>1) 
         % we try all permutations and pick the smallest number of segments
        adcSegmentFactorsPerm=perms(adcSegmentFactors);
        adcSegmentFactorsPermProd=cumprod(adcSegmentFactorsPerm');
        adcSegmentCandidates=unique(adcSegmentFactorsPermProd(:)); % this sorts the sequence
        for i=1:(length(adcSegmentCandidates)-1)
            adcSegments=adcSegmentCandidates(i);
            adcSamplesPerSegment=numSamplesStepped*samplesStep/adcSegments;
            if (adcSamplesPerSegment<=system.adcSamplesLimit && adcSegments<=128)
                break
            end
        end
    else
        adcSamplesPerSegment=numSamplesStepped*samplesStep;
    end
    if (adcSamplesPerSegment<=system.adcSamplesLimit && adcSegments<=128)
        break
    end
    if mode=='shorten'
        numSamplesStepped=numSamplesStepped-1; % try again with a smaller number of samples
    else
        numSamplesStepped=numSamplesStepped+1; % try again with a greater number of samples
    end
        
end
assert(numSamplesStepped>0); % we could not find a suitable segmentation...
assert(adcSamplesPerSegment>0); 
assert(adcSegments<=128);
end

