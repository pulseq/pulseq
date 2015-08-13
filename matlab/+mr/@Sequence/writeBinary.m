function writeBinary(obj,filename)
%WRITEBINARY Write sequence to file in binary format.
%   WRITE(seqObj, filename) Write the sequence data to the given
%   filename using the Pulseq open file format for MR sequences.
%
%   Examples:
%   Write the sequence file to the sequences directory
%
%       write(seqObj,'sequences/gre.seq')
%
% See also  read

binaryCodes = mr.Sequence.getBinaryCodes();
fid=fopen(filename,'w');
fwrite(fid,binaryCodes.fileHeader);
fwrite(fid,binaryCodes.version,'int64');

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
fwrite(fid,size(obj.blockEvents,1),'int64');
for i=1:size(obj.blockEvents,1)
    fwrite(fid,obj.blockEvents(i,:),'int32');
end

if ~isempty(obj.rfLibrary)
    keys=cell2mat(obj.rfLibrary.keys);
    fwrite(fid,binaryCodes.section.rf,'int64');
    fwrite(fid,length(keys),'int64');
    for k=keys
        data=obj.rfLibrary(k).data;
        fwrite(fid,k,'int32');
        fwrite(fid,data(1),'float64');  % amp
        fwrite(fid,data(2:3),'int32');  % mag, phase shape ids
        fwrite(fid,data(4:5),'float64');  % freq, phase offsets
    end
end

arbGradMask=cellfun(@(x)strcmp(x.type,'grad'), obj.gradLibrary.values);
trapGradMask=cellfun(@(x)strcmp(x.type,'trap'), obj.gradLibrary.values);

if any(arbGradMask)
    keys=cell2mat(obj.gradLibrary.keys);
    fwrite(fid,binaryCodes.section.gradients,'int64');
    fwrite(fid,length(keys(arbGradMask)),'int64');
    for k=keys(arbGradMask)
        data=obj.gradLibrary(k).data;
        fwrite(fid,k,'int32');
        fwrite(fid,data(1),'float64');  % amp
        fwrite(fid,data(2),'int32');    % mag shape id
    end
end

if any(trapGradMask)
    keys=cell2mat(obj.gradLibrary.keys);
    fwrite(fid,binaryCodes.section.trapezoids,'int64');
    fwrite(fid,length(keys(trapGradMask)),'int64');
    for k=keys(trapGradMask)
        data=obj.gradLibrary(k).data;
        data(2:end)=round(1e6*data(2:end));
        fwrite(fid,k,'int32');
        fwrite(fid,data(1),'float64');  % amp
        fwrite(fid,data(2:4),'int64');  % rise, flat, fall
    end
end

if ~isempty(obj.adcLibrary)
    keys=cell2mat(obj.adcLibrary.keys);
    fwrite(fid,binaryCodes.section.adc,'int64');
    fwrite(fid,length(keys),'int64');
    for k=keys
        data=obj.adcLibrary(k).data.*[1 1e9 1e6 1 1];
        fwrite(fid,k,'int32');
        fwrite(fid,data(1:3),'int64');      % number, dwell, delay
        fwrite(fid,data(4:5),'float64');    % freq, phase offset
    end
end

if ~isempty(obj.delayLibrary)
    keys=cell2mat(obj.delayLibrary.keys);
    fwrite(fid,binaryCodes.section.delays,'int64');
    fwrite(fid,length(keys),'int64');
    for k=keys
        data = round(1e6*obj.delayLibrary(k).data);
        fwrite(fid,k,'int32');
        fwrite(fid,data,'int64');
    end
end


if ~isempty(obj.shapeLibrary)
    keys=cell2mat(obj.shapeLibrary.keys);
    fwrite(fid,binaryCodes.section.shapes,'int64');
    fwrite(fid,length(keys),'int64');
    for k=keys
        shape = obj.shapeLibrary(k);
        fwrite(fid,k,'int32');
        fwrite(fid,shape.num_samples,'int64');  % num uncompressed
        fwrite(fid,length(shape.data),'int64'); % num compressed
        fwrite(fid,shape.data,'float64');
    end
end

fclose(fid);
end
