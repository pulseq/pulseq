function s=compressShape(w)
%compressShape Compress a gradient or pulse shape.
%   s=compressShape(w) Compress the waveform using a run-length compression
%   scheme on the derivative. This strategy encodes constant and linear
%   waveforms with very few samples. A structure is returned with the
%   fields: 
%     num_samples - the number of samples in the uncompressed waveform
%     data - containing the compressed waveform
%
%   See also decompressShape

diff = w(2:end)-w(1:end-1);
diff =[w(1); diff(:)];

idx=1;
v=[];
while idx<=length(diff)
    curr=diff(idx);
    count=0;
    while idx<=length(diff)-1 && abs(curr-diff(idx+1))<1e-6
        count=count+1;
        idx=idx+1;
    end
    
    if count==0

        v=[v; curr];
    else
        v=[v; curr; curr; count-1]; % Encode extra repetitions (not including two values already read)
        
    end
    idx=idx+1;
end

s.num_samples = length(w);
s.data = v.';