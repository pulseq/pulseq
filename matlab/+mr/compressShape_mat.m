function s=compressShape_mat(w, forceCompression)
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

s=mr.compressShape(w(:), forceCompression);

end