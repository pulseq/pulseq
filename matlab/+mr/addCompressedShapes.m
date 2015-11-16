function s=addCompressedShapes(s1,s2)
%addCompressedShapes Join two compressed gradient or pulse shapes.
%   s=addCompressedShapes(s1,s2) Joint the compressed shape structures s1
%   and s2 to create new shape s.
%
%   *Important* This function assumes the two shapes start and end with a
%   values of 0.
%
%   The result is *logically* the same if the two shapes were uncompressed,
%   concatenated and then recompressed; however, the implementation is
%   memory efficient so very long shapes can be combined in compressed
%   form.
%
%   A structure is returned with the
%   fields: 
%     num_samples - the number of samples in the uncompressed waveform
%     data - containing the compressed waveform
%
%   See also compressShape, decompressShape

x=s1.data;
y=s2.data;
s.num_samples = s1.num_samples + s2.num_samples;

% Test if s1 ends with a run of zeros
isZeroRun = length(x)>=3 && x(end-1)==0 && x(end-2)==0;
if isZeroRun
    numEndingReps=x(end);
    % Test if s2 starts with a run of zeros
    isZeroRun2 = length(y)>=3 && y(1)==0 && y(2)==0;
    if isZeroRun2
        numStartingReps = y(3);
        data = [x(1:end-3), 0, 0, numEndingReps+numStartingReps+2, y(4:end)];
        s.data = data;
        return
    end
end
% Other cases a simple concat is most compressed
s.data = [x,y];
