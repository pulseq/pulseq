function s=compressShape_mat(w)
%compressShape Compress a gradient or pulse shape.
%   s=compressShape(w) Compress the waveform using a run-length compression
%   scheme on the derivative. This strategy encodes constant and linear
%   waveforms with very few samples. A structure is returned with the
%   fields: 
%     num_samples - the number of samples in the uncompressed waveform
%     data - containing the compressed waveform
%
%   See also decompressShape


data = [w(1); diff(w(:))];

maskChanges = [true; abs(diff(data))>1e-8];   % TRUE if values change
vals = data(maskChanges);                     % Elements without repetitions
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
s.data = v.';

