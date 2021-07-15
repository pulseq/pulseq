function w = decompressShape(shape, forceDecompression)
%decompressShape Decompress a gradient or pulse shape.
%   w=decompressShape(shape) Decompress the shape compressed with a run-length
%   compression scheme on the derivative. The given shape is structure with
%   the following fields:
%     num_samples - the number of samples in the uncompressed waveform
%     data - containing the compressed waveform
%
%   See also compressShape


dataPack = shape.data;
dataPackLen = length(dataPack);
numSamples=shape.num_samples;

if nargin<2
    forceDecompression=false;
end

if ~forceDecompression && numSamples==dataPackLen
    % uncompressed shape
    w=dataPack';
    return;
end

w= zeros(1, numSamples) ;                 % pre-allocate the result matrix
                                          % dimensons: (1,length of the data set)
                                       
% decompression starts here

dataPackDiff = dataPack(2:end) - dataPack(1:end-1);

% when dataPackDiff == 0 the subsequent samples are equal ==> marker for
% repeats (run-length encoding)
dataPackMarkers=find(dataPackDiff==0.0);

countPack= 1;                                               % counter 1: points to the current compressed sample
countUnpack= 1;                                             % counter 2: points to the current uncompressed sample

for i=1:length(dataPackMarkers)
    nextPack=dataPackMarkers(i); % careful, this index may have "false positives" , e.g. if the value 3 repeats 3 times, then we will have 3 3 3
    currUnpackSamples=nextPack-countPack;
    if currUnpackSamples < 0 % this rejects false positives
        continue;        
    elseif currUnpackSamples > 0 % do we have an unpacked block to copy?
        w(countUnpack:(countUnpack+currUnpackSamples-1)) = dataPack(countPack:(nextPack-1));
        countPack = countPack + currUnpackSamples;
        countUnpack = countUnpack + currUnpackSamples;
    end
    % now comes the packed/repeated section
    rep=dataPack(countPack+2)+2;
    w(countUnpack:(countUnpack+rep-1))=dataPack(countPack);
    countPack= countPack + 3;
    countUnpack= countUnpack + rep;
end

% samples left?
if (countPack<=dataPackLen)
    assert(dataPackLen-countPack==numSamples-countUnpack);
    % copy the rest of the shape, it is unpacked
    w(countUnpack:end)= dataPack(countPack:end);
end

w = cumsum(w);
w=w(:);

% % test the new function against the old slow version
% w1=mr.decompressShape0(shape);
% assert(all(w==w1));
