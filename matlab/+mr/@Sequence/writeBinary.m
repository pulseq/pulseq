function writeBinary(obj,filename)
%WRITEBINARY Write sequence to file in binary format.
%   WRITEBINARY(seqObj, filename) Write the sequence data to the given
%   filename using the binary version of the Pulseq open file format for MR
%   sequences. The file specification is available at 
%   http://pulseq.github.io
%
%   Examples:
%   Write the sequence file to the sequences directory
%
%       writeBinary(seqObj,'sequences/gre.bin')
%
% See also  readBinary

binaryCodes = obj.getBinaryCodes();
fid=fopen(filename,'w');
fwrite(fid,binaryCodes.fileHeader);
fwrite(fid,binaryCodes.version_major,'int64');
fwrite(fid,binaryCodes.version_minor,'int64');
fwrite(fid,binaryCodes.version_revision,'int64');

if ~isempty(obj.definitions)
    fwrite(fid,binaryCodes.section.definitions,'int64');
    keys = obj.definitions.keys;
    values = obj.definitions.values;
    fwrite(fid,length(keys),'int64');
    for i=1:length(keys)
        fwrite(fid,[keys{i} 0]);
        fwrite(fid,2,'int8');
        val = values{i};
        fwrite(fid,length(val),'int8');
        fwrite(fid,val,'float64');
    end
end

fwrite(fid,binaryCodes.section.blocks,'int64');
fwrite(fid,size(obj.blockEvents,2),'int64');
for i = 1:length(obj.blockEvents)
    fwrite(fid,obj.blockEvents{i}(:),'int32');
end

if ~isempty(obj.rfLibrary.keys)
    keys=obj.rfLibrary.keys;
    fwrite(fid,binaryCodes.section.rf,'int64');
    fwrite(fid,length(keys),'int64');
    for k=keys
        data=obj.rfLibrary.data(k).array;
        fwrite(fid,k,'int32');
        fwrite(fid,data(1),'float64');  % amp
        fwrite(fid,data(2:3),'int32');  % mag, phase shape ids
        fwrite(fid,data(4:5),'float64');  % freq, phase offsets
    end
end

arbGradMask=obj.gradLibrary.type=='g';
trapGradMask=obj.gradLibrary.type=='t';

if any(arbGradMask)
    keys=obj.gradLibrary.keys;
    fwrite(fid,binaryCodes.section.gradients,'int64');
    fwrite(fid,length(keys(arbGradMask)),'int64');
    for k=keys(arbGradMask)
        data=obj.gradLibrary.data(k).array;
        fwrite(fid,k,'int32');
        fwrite(fid,data(1),'float64');  % amp
        fwrite(fid,data(2),'int32');    % mag shape id
    end
end

if any(trapGradMask)
    keys=obj.gradLibrary.keys;
    fwrite(fid,binaryCodes.section.trapezoids,'int64');
    fwrite(fid,length(keys(trapGradMask)),'int64');
    for k=keys(trapGradMask)
        data=obj.gradLibrary.data(k).array;
        data(2:end)=round(1e6*data(2:end));
        fwrite(fid,k,'int32');
        fwrite(fid,data(1),'float64');  % amp
        fwrite(fid,data(2:4),'int64');  % rise, flat, fall
    end
end

if ~isempty(obj.adcLibrary.keys)
    keys=obj.adcLibrary.keys;
    fwrite(fid,binaryCodes.section.adc,'int64');
    fwrite(fid,length(keys),'int64');
    for k=keys
        data=obj.adcLibrary.data(k).array(1:5).*[1 1e9 1e6 1 1];
        fwrite(fid,k,'int32');
        fwrite(fid,data(1:3),'int64');      % number, dwell, delay
        fwrite(fid,data(4:5),'float64');    % freq, phase offset
    end
end

if ~isempty(obj.delayLibrary.keys)
    keys=obj.delayLibrary.keys;
    fwrite(fid,binaryCodes.section.delays,'int64');
    fwrite(fid,length(keys),'int64');
    for k=keys
        data = round(1e6*obj.delayLibrary.data(k).array);
        fwrite(fid,k,'int32');
        fwrite(fid,data,'int64');
    end
end


if ~isempty(obj.shapeLibrary.keys)
    keys=obj.shapeLibrary.keys;
    fwrite(fid,binaryCodes.section.shapes,'int64');
    fwrite(fid,length(keys),'int64');
    for k=keys
        shape = obj.shapeLibrary.data(k).array;
        num_samples=shape(1);
        data=shape(2:end);
        fwrite(fid,k,'int32');
        fwrite(fid,num_samples,'int64');  % num uncompressed
        fwrite(fid,length(data),'int64'); % num compressed
        fwrite(fid,data,'float64');
    end
end

fclose(fid);
end
