function s=compressShape(w, forceCompression)
%compressShape Compress a gradient or pulse shape.
%   s=compressShape(w) Compress the waveform using a run-length compression
%   scheme on the derivative. This strategy encodes constant and linear
%   waveforms with very few samples. A structure is returned with the
%   fields: 
%     num_samples - the number of samples in the uncompressed waveform
%     data - containing the compressed waveform
%
%   See also decompressShape

if nargin<2
    forceCompression=false;
end

if any(~isfinite(w))
    error('compressShape() received infinite samples');
end

if ~forceCompression && length(w) <= 4 % avoid compressing very short shapes
    s.num_samples=length(w);
    s.data = w(:)';
    return;
end


% %MZ: old code with implicit quantization
% data = [w(1); diff(w(:))];
% maskChanges = [true; abs(diff(data))>1e-8];   % TRUE if values change
% vals = data(maskChanges);                     % Elements without repetitions

% MZ: explicit quantization with error correction
quant_fac=1e-7; % single precision floating point has ~7.25 decimal places 
ws=w./quant_fac;
datq=round([ws(1); diff(ws(:))]);
qerr=ws(:)-cumsum(datq);
qcor=[0; diff(round(qerr))];
datd=datq+qcor;
maskChanges=[true; diff(datd)~=0];
vals=datd(maskChanges).*quant_fac;            % Elements without repetitions

k = find([maskChanges', true]);               % Indices of changes
n = diff(k)';                                 % Number of repetitions

% Encode in Pulseq format
nExtra=n-2;
vals2=vals; 
vals2(nExtra<0)=nan;
nExtra(nExtra<0)=nan;
v=[vals vals2 nExtra]';
v=v(isfinite(v));
v(abs(v)<1e-10)=0;
s.num_samples = length(w);
% decide whether compression makes sense, otherwise store the original
if forceCompression || s.num_samples > length(v)
    s.data = v';
else
    s.data = w(:)';
end


