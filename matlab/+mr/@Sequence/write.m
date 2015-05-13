function write(obj,filename)
%WRITE Write sequence to file.
%   WRITE(seqObj, filename) Write the sequence data to the given
%   filename using the open file format for MR sequences.
%
%   Examples:
%   Write the sequence file to the my_sequences directory
%
%       write(seqObj,'my_sequences/se.txt')
%
% See also  read

fid=fopen(filename,'w');
fprintf(fid,'# Sequence Blocks\n');
fprintf(fid,'# Created by MATLAB open sequence\n\n');

if ~isempty(obj.definitions)
    fprintf(fid,'[DEFINITIONS]\n');
    keys = obj.definitions.keys;
    values = obj.definitions.values;
    for i=1:length(keys)
        fprintf(fid,'%s ',keys{i});
        fprintf(fid,'%g ',values{i});
        fprintf(fid,'\n');
    end
end

fprintf(fid,'# Format of blocks:\n');
fprintf(fid,'##  D RF  GX  GY  GZ ADC\n');
fprintf(fid,'[BLOCKS]\n');
idFormatWidth=length(num2str(size(obj.blockEvents,1)));
idFormatStr=['%' num2str(idFormatWidth) 'd'];
for i=1:size(obj.blockEvents,1)
    fprintf(fid,[idFormatStr ' %2d %2d %3d %3d %3d %2d\n'],[i obj.blockEvents(i,:)]);
end
fprintf(fid,'\n');

if ~isempty(obj.rfLibrary)
    fprintf(fid,'# Format of RF events:\n');
    fprintf(fid,'# id amplitude mag_id phase_id freq phase\n');
    fprintf(fid,'# ..        Hz   ....     ....   Hz   rad\n');
    fprintf(fid,'[RF]\n');
    keys=cell2mat(obj.rfLibrary.keys);
    for k=keys
        fprintf(fid,'%d %12g %d %d %g %g\n',[k obj.rfLibrary(k).data]);
    end
    fprintf(fid,'\n');
end

arbGradMask=cellfun(@(x)strcmp(x.type,'arb'), obj.gradLibrary.values);
trapGradMask=cellfun(@(x)strcmp(x.type,'trap'), obj.gradLibrary.values);

if any(arbGradMask)
    fprintf(fid,'# Format of arbitrary gradients:\n');
    fprintf(fid,'# id amplitude shape_id\n');
    fprintf(fid,'# ..      Hz/m     ....\n');
    fprintf(fid,'[GRADIENTS]\n');
    keys=cell2mat(obj.gradLibrary.keys);
    for k=keys(arbGradMask)
        fprintf(fid,'%d %12g %d \n',[k obj.gradLibrary(k).data]);
    end
    fprintf(fid,'\n');
end

if any(trapGradMask)
    fprintf(fid,'# Format of trapezoid gradients:\n');
    fprintf(fid,'# id amplitude rise flat fall\n');
    fprintf(fid,'# ..      Hz/m   us   us   us\n');
    fprintf(fid,'[TRAP]\n');
    keys=cell2mat(obj.gradLibrary.keys);
    for k=keys(trapGradMask)
        data=obj.gradLibrary(k).data;
        data(2:end)=round(1e6*data(2:end));
        fprintf(fid,'%2d %12g %3d %4d %3d\n',[k data]);
    end
    fprintf(fid,'\n');
end

if ~isempty(obj.adcLibrary)
    fprintf(fid,'# Format of ADC events:\n');
    fprintf(fid,'# id num dwell delay freq phase\n');
    fprintf(fid,'# ..  ..    ns    us   Hz   rad\n');
    fprintf(fid,'[ADC]\n');
    keys=cell2mat(obj.adcLibrary.keys);
    for k=keys
        data=obj.adcLibrary(k).data.*[1 1e9 1e6 1 1];
        fprintf(fid,'%2d %3d %6d %3d %g %g\n',[k data]);
    end
    fprintf(fid,'\n');
end

if ~isempty(obj.delayLibrary)
    fprintf(fid,'# Format of delays:\n');
    fprintf(fid,'# id delay (us)\n');
    fprintf(fid,'[DELAYS]\n');
    keys=cell2mat(obj.delayLibrary.keys);
    for k=keys
        fprintf(fid,'%d %d\n',[k round(1e6*obj.delayLibrary(k).data)]);
    end
    fprintf(fid,'\n');
end


if ~isempty(obj.shapeLibrary)
    fprintf(fid,'# Sequence Shapes\n');
    fprintf(fid,'[SHAPES]\n\n');
    keys=cell2mat(obj.shapeLibrary.keys);
    for k=keys
        shape = obj.shapeLibrary(k);
        fprintf(fid,'shape_id %d\n',k);
        fprintf(fid,'num_samples %d\n',shape.num_samples);
        fprintf(fid,'%g\n',shape.data);
        fprintf(fid,'\n');
    end
end

fclose(fid);
end
